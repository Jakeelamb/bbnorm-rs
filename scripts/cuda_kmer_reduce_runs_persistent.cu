#include <cub/cub.cuh>
#include <cuda_runtime.h>

#include <climits>
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

class ReusableReducer {
  public:
    ~ReusableReducer() {
        cudaFree(reduce_temp_);
        cudaFree(sort_temp_);
        cudaFree(device_unique_count_);
        cudaFree(device_counts_);
        cudaFree(device_unique_);
        cudaFree(device_keys_sorted_);
        cudaFree(device_keys_in_);
    }

    void reduce_chunk(const std::vector<uint64_t>& host_keys) {
        const size_t n = host_keys.size();
        if (n > static_cast<size_t>(INT_MAX)) {
            std::cerr << "too many k-mers for one experimental CUB batch: " << n << "\n";
            std::exit(2);
        }
        if (n == 0) {
            const uint64_t unique_count = 0;
            std::cout.write(reinterpret_cast<const char*>(&unique_count), sizeof(unique_count));
            std::cout.flush();
            return;
        }

        reserve(n);

        check_cuda(cudaMemcpy(device_keys_in_, host_keys.data(), n * sizeof(uint64_t),
                              cudaMemcpyHostToDevice),
                   "cudaMemcpy host keys to device");
        check_cuda(cub::DeviceRadixSort::SortKeys(
                       sort_temp_, sort_temp_bytes_, device_keys_in_, device_keys_sorted_,
                       static_cast<int>(n)),
                   "cub DeviceRadixSort");
        check_cuda(cub::DeviceRunLengthEncode::Encode(
                       reduce_temp_, reduce_temp_bytes_, device_keys_sorted_, device_unique_,
                       device_counts_, device_unique_count_, static_cast<int>(n)),
                   "cub DeviceRunLengthEncode");

        int unique_count_i32 = 0;
        check_cuda(cudaMemcpy(&unique_count_i32, device_unique_count_, sizeof(int),
                              cudaMemcpyDeviceToHost),
                   "cudaMemcpy unique_count to host");
        const size_t unique_count = static_cast<size_t>(unique_count_i32);
        host_unique_.resize(unique_count);
        host_counts_.resize(unique_count);
        check_cuda(cudaMemcpy(host_unique_.data(), device_unique_, unique_count * sizeof(uint64_t),
                              cudaMemcpyDeviceToHost),
                   "cudaMemcpy unique keys to host");
        check_cuda(cudaMemcpy(host_counts_.data(), device_counts_, unique_count * sizeof(int),
                              cudaMemcpyDeviceToHost),
                   "cudaMemcpy counts to host");

        const uint64_t unique_count_u64 = static_cast<uint64_t>(unique_count);
        std::cout.write(reinterpret_cast<const char*>(&unique_count_u64), sizeof(unique_count_u64));
        for (size_t i = 0; i < unique_count; ++i) {
            const uint32_t count = static_cast<uint32_t>(host_counts_[i]);
            std::cout.write(reinterpret_cast<const char*>(&host_unique_[i]), sizeof(uint64_t));
            std::cout.write(reinterpret_cast<const char*>(&count), sizeof(uint32_t));
        }
        std::cout.flush();
    }

  private:
    void reserve(size_t n) {
        if (n > capacity_) {
            cudaFree(device_keys_in_);
            cudaFree(device_keys_sorted_);
            cudaFree(device_unique_);
            cudaFree(device_counts_);
            device_keys_in_ = nullptr;
            device_keys_sorted_ = nullptr;
            device_unique_ = nullptr;
            device_counts_ = nullptr;

            check_cuda(cudaMalloc(&device_keys_in_, n * sizeof(uint64_t)),
                       "cudaMalloc device_keys_in");
            check_cuda(cudaMalloc(&device_keys_sorted_, n * sizeof(uint64_t)),
                       "cudaMalloc device_keys_sorted");
            check_cuda(cudaMalloc(&device_unique_, n * sizeof(uint64_t)),
                       "cudaMalloc device_unique");
            check_cuda(cudaMalloc(&device_counts_, n * sizeof(int)), "cudaMalloc device_counts");
            capacity_ = n;
        }
        if (device_unique_count_ == nullptr) {
            check_cuda(cudaMalloc(&device_unique_count_, sizeof(int)),
                       "cudaMalloc device_unique_count");
        }

        size_t required_sort_temp_bytes = 0;
        check_cuda(cub::DeviceRadixSort::SortKeys(
                       nullptr, required_sort_temp_bytes, device_keys_in_, device_keys_sorted_,
                       static_cast<int>(n)),
                   "cub DeviceRadixSort temp sizing");
        if (required_sort_temp_bytes > sort_temp_bytes_) {
            cudaFree(sort_temp_);
            sort_temp_ = nullptr;
            check_cuda(cudaMalloc(&sort_temp_, required_sort_temp_bytes), "cudaMalloc sort_temp");
            sort_temp_bytes_ = required_sort_temp_bytes;
        }

        size_t required_reduce_temp_bytes = 0;
        check_cuda(cub::DeviceRunLengthEncode::Encode(
                       nullptr, required_reduce_temp_bytes, device_keys_sorted_, device_unique_,
                       device_counts_, device_unique_count_, static_cast<int>(n)),
                   "cub DeviceRunLengthEncode temp sizing");
        if (required_reduce_temp_bytes > reduce_temp_bytes_) {
            cudaFree(reduce_temp_);
            reduce_temp_ = nullptr;
            check_cuda(cudaMalloc(&reduce_temp_, required_reduce_temp_bytes),
                       "cudaMalloc reduce_temp");
            reduce_temp_bytes_ = required_reduce_temp_bytes;
        }
    }

    size_t capacity_ = 0;
    uint64_t* device_keys_in_ = nullptr;
    uint64_t* device_keys_sorted_ = nullptr;
    uint64_t* device_unique_ = nullptr;
    int* device_counts_ = nullptr;
    int* device_unique_count_ = nullptr;
    void* sort_temp_ = nullptr;
    void* reduce_temp_ = nullptr;
    size_t sort_temp_bytes_ = 0;
    size_t reduce_temp_bytes_ = 0;
    std::vector<uint64_t> host_unique_;
    std::vector<int> host_counts_;
};

int main() {
    std::ios::sync_with_stdio(false);
    ReusableReducer reducer;
    std::vector<uint64_t> keys;
    while (true) {
        uint64_t n = 0;
        std::cin.read(reinterpret_cast<char*>(&n), sizeof(n));
        if (!std::cin) {
            return 0;
        }
        if (n == STOP_SENTINEL) {
            return 0;
        }
        keys.resize(static_cast<size_t>(n));
        std::cin.read(reinterpret_cast<char*>(keys.data()), static_cast<std::streamsize>(n * sizeof(uint64_t)));
        if (!std::cin) {
            std::cerr << "truncated input chunk\n";
            return 2;
        }
        reducer.reduce_chunk(keys);
    }
}
