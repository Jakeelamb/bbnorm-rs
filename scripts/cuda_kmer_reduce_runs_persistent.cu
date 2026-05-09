#include <cub/cub.cuh>
#include <cuda_runtime.h>

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

static constexpr uint64_t STOP_SENTINEL = UINT64_MAX;

static void check_cuda(cudaError_t status, const char* what) {
    if (status != cudaSuccess) {
        std::cerr << what << ": " << cudaGetErrorString(status) << "\n";
        std::exit(2);
    }
}

static void reduce_chunk(const std::vector<uint64_t>& host_keys) {
    const size_t n = host_keys.size();
    if (n > static_cast<size_t>(INT32_MAX)) {
        std::cerr << "too many k-mers for one experimental CUB batch: " << n << "\n";
        std::exit(2);
    }
    if (n == 0) {
        const uint64_t unique_count = 0;
        std::cout.write(reinterpret_cast<const char*>(&unique_count), sizeof(unique_count));
        std::cout.flush();
        return;
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

    const uint64_t unique_count_u64 = static_cast<uint64_t>(unique_count);
    std::cout.write(reinterpret_cast<const char*>(&unique_count_u64), sizeof(unique_count_u64));
    for (size_t i = 0; i < unique_count; ++i) {
        const uint32_t count = static_cast<uint32_t>(counts[i]);
        std::cout.write(reinterpret_cast<const char*>(&unique[i]), sizeof(uint64_t));
        std::cout.write(reinterpret_cast<const char*>(&count), sizeof(uint32_t));
    }
    std::cout.flush();

    cudaFree(reduce_temp);
    cudaFree(sort_temp);
    cudaFree(device_unique_count);
    cudaFree(device_counts);
    cudaFree(device_unique);
    cudaFree(device_keys_sorted);
    cudaFree(device_keys_in);
}

int main() {
    std::ios::sync_with_stdio(false);
    while (true) {
        uint64_t n = 0;
        std::cin.read(reinterpret_cast<char*>(&n), sizeof(n));
        if (!std::cin) {
            return 0;
        }
        if (n == STOP_SENTINEL) {
            return 0;
        }
        std::vector<uint64_t> keys(static_cast<size_t>(n));
        std::cin.read(reinterpret_cast<char*>(keys.data()), static_cast<std::streamsize>(n * sizeof(uint64_t)));
        if (!std::cin) {
            std::cerr << "truncated input chunk\n";
            return 2;
        }
        reduce_chunk(keys);
    }
}
