#include <cub/cub.cuh>
#include <cuda_runtime.h>

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

static void check_cuda(cudaError_t status, const char* what) {
    if (status != cudaSuccess) {
        std::cerr << what << ": " << cudaGetErrorString(status) << "\n";
        std::exit(2);
    }
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: cuda_kmer_reduce_runs <kmers.u64> <runs.bin>\n";
        return 2;
    }

    std::ifstream input(argv[1], std::ios::binary | std::ios::ate);
    if (!input) {
        std::cerr << "failed to open " << argv[1] << "\n";
        return 2;
    }
    const auto bytes = input.tellg();
    if (bytes < 0 || (static_cast<uint64_t>(bytes) % sizeof(uint64_t)) != 0) {
        std::cerr << "input length is not a u64 multiple\n";
        return 2;
    }
    const size_t n = static_cast<size_t>(bytes) / sizeof(uint64_t);
    if (n > static_cast<size_t>(INT32_MAX)) {
        std::cerr << "too many k-mers for one experimental CUB batch: " << n << "\n";
        return 2;
    }
    input.seekg(0);
    std::vector<uint64_t> host_keys(n);
    input.read(reinterpret_cast<char*>(host_keys.data()), bytes);
    if (!input) {
        std::cerr << "failed to read input\n";
        return 2;
    }

    uint64_t* device_keys_in = nullptr;
    uint64_t* device_keys_sorted = nullptr;
    uint64_t* device_unique = nullptr;
    int* device_counts = nullptr;
    int* device_unique_count = nullptr;
    void* sort_temp = nullptr;
    void* reduce_temp = nullptr;
    size_t sort_temp_bytes = 0;
    size_t reduce_temp_bytes = 0;

    check_cuda(cudaMalloc(&device_keys_in, n * sizeof(uint64_t)), "cudaMalloc device_keys_in");
    check_cuda(cudaMalloc(&device_keys_sorted, n * sizeof(uint64_t)), "cudaMalloc device_keys_sorted");
    check_cuda(cudaMalloc(&device_unique, n * sizeof(uint64_t)), "cudaMalloc device_unique");
    check_cuda(cudaMalloc(&device_counts, n * sizeof(int)), "cudaMalloc device_counts");
    check_cuda(cudaMalloc(&device_unique_count, sizeof(int)), "cudaMalloc device_unique_count");

    check_cuda(cub::DeviceRadixSort::SortKeys(
                   nullptr, sort_temp_bytes, device_keys_in, device_keys_sorted, static_cast<int>(n)),
               "cub DeviceRadixSort temp sizing");
    check_cuda(cudaMalloc(&sort_temp, sort_temp_bytes), "cudaMalloc sort_temp");
    check_cuda(cub::DeviceRunLengthEncode::Encode(
                   nullptr, reduce_temp_bytes, device_keys_sorted, device_unique, device_counts,
                   device_unique_count, static_cast<int>(n)),
               "cub DeviceRunLengthEncode temp sizing");
    check_cuda(cudaMalloc(&reduce_temp, reduce_temp_bytes), "cudaMalloc reduce_temp");

    check_cuda(cudaMemcpy(device_keys_in, host_keys.data(), n * sizeof(uint64_t),
                          cudaMemcpyHostToDevice),
               "cudaMemcpy host keys to device");
    check_cuda(cub::DeviceRadixSort::SortKeys(
                   sort_temp, sort_temp_bytes, device_keys_in, device_keys_sorted, static_cast<int>(n)),
               "cub DeviceRadixSort");
    check_cuda(cub::DeviceRunLengthEncode::Encode(
                   reduce_temp, reduce_temp_bytes, device_keys_sorted, device_unique, device_counts,
                   device_unique_count, static_cast<int>(n)),
               "cub DeviceRunLengthEncode");

    int unique_count_i32 = 0;
    check_cuda(cudaMemcpy(&unique_count_i32, device_unique_count, sizeof(int),
                          cudaMemcpyDeviceToHost),
               "cudaMemcpy unique_count to host");
    const size_t unique_count = static_cast<size_t>(unique_count_i32);
    std::vector<uint64_t> unique(unique_count);
    std::vector<int> counts(unique_count);
    check_cuda(cudaMemcpy(unique.data(), device_unique, unique_count * sizeof(uint64_t),
                          cudaMemcpyDeviceToHost),
               "cudaMemcpy unique keys to host");
    check_cuda(cudaMemcpy(counts.data(), device_counts, unique_count * sizeof(int),
                          cudaMemcpyDeviceToHost),
               "cudaMemcpy counts to host");

    std::ofstream output(argv[2], std::ios::binary);
    if (!output) {
        std::cerr << "failed to create " << argv[2] << "\n";
        return 2;
    }
    for (size_t i = 0; i < unique_count; ++i) {
        const uint32_t count = static_cast<uint32_t>(counts[i]);
        output.write(reinterpret_cast<const char*>(&unique[i]), sizeof(uint64_t));
        output.write(reinterpret_cast<const char*>(&count), sizeof(uint32_t));
    }
    if (!output) {
        std::cerr << "failed to write " << argv[2] << "\n";
        return 2;
    }

    std::cerr << "kmers\t" << n << "\n";
    std::cerr << "unique\t" << unique_count << "\n";

    cudaFree(reduce_temp);
    cudaFree(sort_temp);
    cudaFree(device_unique_count);
    cudaFree(device_counts);
    cudaFree(device_unique);
    cudaFree(device_keys_sorted);
    cudaFree(device_keys_in);
    return 0;
}
