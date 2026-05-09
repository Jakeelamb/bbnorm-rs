#!/usr/bin/env python3
"""Probe whether CUDA sort/reduce is promising for BBNorm k-mer counting.

This is intentionally a standalone feasibility harness. It extracts real
canonical short k-mers from FASTQ/FASTQ.GZ inputs into a binary u64 stream, then
compiles and runs a tiny CUDA/CUB benchmark that compares CPU sort/reduce
against GPU host-to-device copy, radix sort, run-length encode, and
device-to-host copy.
"""

from __future__ import annotations

import argparse
import gzip
import os
import struct
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_R1 = ROOT / "tmp/human_benchmark_8threads/human_GRCh38_500k_R1.fq.gz"
DEFAULT_R2 = ROOT / "tmp/human_benchmark_8threads/human_GRCh38_500k_R2.fq.gz"

CUDA_SOURCE = r"""
#include <algorithm>
#include <cub/cub.cuh>
#include <chrono>
#include <cstring>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>

static void check_cuda(cudaError_t status, const char* what) {
    if (status != cudaSuccess) {
        std::cerr << what << ": " << cudaGetErrorString(status) << "\n";
        std::exit(2);
    }
}

static double elapsed_ms(std::chrono::steady_clock::time_point start,
                         std::chrono::steady_clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

static double event_ms(cudaEvent_t start, cudaEvent_t stop) {
    float ms = 0.0f;
    check_cuda(cudaEventElapsedTime(&ms, start, stop), "cudaEventElapsedTime");
    return static_cast<double>(ms);
}

static uint64_t checksum_pairs(const std::vector<uint64_t>& keys,
                               const std::vector<int>& counts) {
    uint64_t hash = 1469598103934665603ull;
    for (size_t i = 0; i < keys.size(); ++i) {
        hash ^= keys[i];
        hash *= 1099511628211ull;
        hash ^= static_cast<uint32_t>(counts[i]);
        hash *= 1099511628211ull;
    }
    return hash;
}

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "usage: cuda_kmer_sort_reduce <kmers.u64|-> [--skip-cpu-reference]\n";
        return 2;
    }
    const bool skip_cpu_reference = argc == 3 && std::string(argv[2]) == "--skip-cpu-reference";
    if (argc == 3 && !skip_cpu_reference) {
        std::cerr << "unknown option: " << argv[2] << "\n";
        return 2;
    }

    std::vector<uint64_t> host_keys;
    const std::string input_path(argv[1]);
    if (input_path == "-") {
        std::cin.sync_with_stdio(false);
        std::vector<char> raw(
            std::istreambuf_iterator<char>(std::cin),
            std::istreambuf_iterator<char>());
        if ((raw.size() % sizeof(uint64_t)) != 0) {
            std::cerr << "stdin length is not a u64 multiple\n";
            return 2;
        }
        host_keys.resize(raw.size() / sizeof(uint64_t));
        std::memcpy(host_keys.data(), raw.data(), raw.size());
    } else {
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
        input.seekg(0);
        host_keys.resize(n);
        input.read(reinterpret_cast<char*>(host_keys.data()), bytes);
        if (!input) {
            std::cerr << "failed to read input\n";
            return 2;
        }
    }
    const size_t n = host_keys.size();

    std::vector<uint64_t> cpu_unique;
    std::vector<int> cpu_counts;
    uint64_t cpu_checksum = 0;
    double cpu_sort_reduce_ms = 0.0;
    if (!skip_cpu_reference) {
    auto cpu_keys = host_keys;
    auto cpu_start = std::chrono::steady_clock::now();
    std::sort(cpu_keys.begin(), cpu_keys.end());
    cpu_unique.reserve(cpu_keys.size());
    cpu_counts.reserve(cpu_keys.size());
    for (size_t i = 0; i < cpu_keys.size();) {
        size_t j = i + 1;
        while (j < cpu_keys.size() && cpu_keys[j] == cpu_keys[i]) {
            ++j;
        }
        cpu_unique.push_back(cpu_keys[i]);
        cpu_counts.push_back(static_cast<int>(j - i));
        i = j;
    }
    auto cpu_end = std::chrono::steady_clock::now();
    cpu_checksum = checksum_pairs(cpu_unique, cpu_counts);
    cpu_sort_reduce_ms = elapsed_ms(cpu_start, cpu_end);
    }

    cudaEvent_t h2d_start, h2d_stop, sort_start, sort_stop, reduce_start, reduce_stop, d2h_start, d2h_stop;
    check_cuda(cudaEventCreate(&h2d_start), "cudaEventCreate");
    check_cuda(cudaEventCreate(&h2d_stop), "cudaEventCreate");
    check_cuda(cudaEventCreate(&sort_start), "cudaEventCreate");
    check_cuda(cudaEventCreate(&sort_stop), "cudaEventCreate");
    check_cuda(cudaEventCreate(&reduce_start), "cudaEventCreate");
    check_cuda(cudaEventCreate(&reduce_stop), "cudaEventCreate");
    check_cuda(cudaEventCreate(&d2h_start), "cudaEventCreate");
    check_cuda(cudaEventCreate(&d2h_stop), "cudaEventCreate");

    uint64_t* device_keys_in = nullptr;
    uint64_t* device_keys_sorted = nullptr;
    uint64_t* device_unique = nullptr;
    int* device_counts = nullptr;
    int* device_unique_count = nullptr;
    void* sort_temp = nullptr;
    void* reduce_temp = nullptr;
    size_t sort_temp_bytes = 0;
    size_t reduce_temp_bytes = 0;
    int unique_count_i32 = 0;

    check_cuda(cudaMalloc(&device_keys_in, n * sizeof(uint64_t)), "cudaMalloc device_keys_in");
    check_cuda(cudaMalloc(&device_keys_sorted, n * sizeof(uint64_t)), "cudaMalloc device_keys_sorted");
    check_cuda(cudaMalloc(&device_unique, n * sizeof(uint64_t)), "cudaMalloc device_unique");
    check_cuda(cudaMalloc(&device_counts, n * sizeof(int)), "cudaMalloc device_counts");
    check_cuda(cudaMalloc(&device_unique_count, sizeof(int)), "cudaMalloc device_unique_count");

    check_cuda(cub::DeviceRadixSort::SortKeys(
        nullptr,
        sort_temp_bytes,
        device_keys_in,
        device_keys_sorted,
        static_cast<int>(n)),
        "cub DeviceRadixSort temp sizing");
    check_cuda(cudaMalloc(&sort_temp, sort_temp_bytes), "cudaMalloc sort_temp");

    check_cuda(cub::DeviceRunLengthEncode::Encode(
        nullptr,
        reduce_temp_bytes,
        device_keys_sorted,
        device_unique,
        device_counts,
        device_unique_count,
        static_cast<int>(n)),
        "cub DeviceRunLengthEncode temp sizing");
    check_cuda(cudaMalloc(&reduce_temp, reduce_temp_bytes), "cudaMalloc reduce_temp");

    check_cuda(cudaEventRecord(h2d_start), "cudaEventRecord h2d_start");
    check_cuda(cudaMemcpy(
        device_keys_in,
        host_keys.data(),
        n * sizeof(uint64_t),
        cudaMemcpyHostToDevice),
        "cudaMemcpy host keys to device");
    check_cuda(cudaEventRecord(h2d_stop), "cudaEventRecord h2d_stop");
    check_cuda(cudaEventSynchronize(h2d_stop), "cudaEventSynchronize h2d_stop");

    check_cuda(cudaEventRecord(sort_start), "cudaEventRecord sort_start");
    check_cuda(cub::DeviceRadixSort::SortKeys(
        sort_temp,
        sort_temp_bytes,
        device_keys_in,
        device_keys_sorted,
        static_cast<int>(n)),
        "cub DeviceRadixSort");
    check_cuda(cudaEventRecord(sort_stop), "cudaEventRecord sort_stop");
    check_cuda(cudaEventSynchronize(sort_stop), "cudaEventSynchronize sort_stop");

    check_cuda(cudaEventRecord(reduce_start), "cudaEventRecord reduce_start");
    check_cuda(cub::DeviceRunLengthEncode::Encode(
        reduce_temp,
        reduce_temp_bytes,
        device_keys_sorted,
        device_unique,
        device_counts,
        device_unique_count,
        static_cast<int>(n)),
        "cub DeviceRunLengthEncode");
    check_cuda(cudaEventRecord(reduce_stop), "cudaEventRecord reduce_stop");
    check_cuda(cudaEventSynchronize(reduce_stop), "cudaEventSynchronize reduce_stop");
    check_cuda(cudaMemcpy(
        &unique_count_i32,
        device_unique_count,
        sizeof(int),
        cudaMemcpyDeviceToHost),
        "cudaMemcpy unique_count to host");
    const size_t unique_count = static_cast<size_t>(unique_count_i32);

    std::vector<uint64_t> gpu_unique(unique_count);
    std::vector<int> gpu_counts(unique_count);
    check_cuda(cudaEventRecord(d2h_start), "cudaEventRecord d2h_start");
    check_cuda(cudaMemcpy(
        gpu_unique.data(),
        device_unique,
        unique_count * sizeof(uint64_t),
        cudaMemcpyDeviceToHost),
        "cudaMemcpy unique keys to host");
    check_cuda(cudaMemcpy(
        gpu_counts.data(),
        device_counts,
        unique_count * sizeof(int),
        cudaMemcpyDeviceToHost),
        "cudaMemcpy counts to host");
    check_cuda(cudaEventRecord(d2h_stop), "cudaEventRecord d2h_stop");
    check_cuda(cudaEventSynchronize(d2h_stop), "cudaEventSynchronize d2h_stop");

    const uint64_t gpu_checksum = checksum_pairs(gpu_unique, gpu_counts);
    const bool ok = skip_cpu_reference || (unique_count == cpu_unique.size() && gpu_checksum == cpu_checksum);
    const double h2d_ms = event_ms(h2d_start, h2d_stop);
    const double sort_ms = event_ms(sort_start, sort_stop);
    const double reduce_ms = event_ms(reduce_start, reduce_stop);
    const double d2h_ms = event_ms(d2h_start, d2h_stop);
    const double gpu_total_ms = h2d_ms + sort_ms + reduce_ms + d2h_ms;

    std::cout << "kmers\t" << n << "\n";
    std::cout << "unique\t" << unique_count << "\n";
    std::cout << "gpu_backend\tcub_radix_sort_run_length_encode\n";
    std::cout << "cpu_reference\t" << (skip_cpu_reference ? "false" : "true") << "\n";
    if (!skip_cpu_reference) {
        std::cout << "cpu_sort_reduce_ms\t" << cpu_sort_reduce_ms << "\n";
    }
    std::cout << "gpu_h2d_ms\t" << h2d_ms << "\n";
    std::cout << "gpu_sort_ms\t" << sort_ms << "\n";
    std::cout << "gpu_reduce_ms\t" << reduce_ms << "\n";
    std::cout << "gpu_d2h_ms\t" << d2h_ms << "\n";
    std::cout << "gpu_total_timed_ms\t" << gpu_total_ms << "\n";
    if (!skip_cpu_reference) {
        std::cout << "cpu_checksum\t" << cpu_checksum << "\n";
    }
    std::cout << "gpu_checksum\t" << gpu_checksum << "\n";
    std::cout << "match\t" << (skip_cpu_reference ? "not_checked" : (ok ? "true" : "false")) << "\n";
    cudaFree(reduce_temp);
    cudaFree(sort_temp);
    cudaFree(device_unique_count);
    cudaFree(device_counts);
    cudaFree(device_unique);
    cudaFree(device_keys_sorted);
    cudaFree(device_keys_in);
    return ok ? 0 : 1;
}
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", type=Path, default=DEFAULT_R1)
    parser.add_argument("--r2", type=Path, default=DEFAULT_R2)
    parser.add_argument("--reads", type=int, default=50_000, help="paired records to read")
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=ROOT / "tmp" / "cuda_kmer_sort_reduce_probe",
    )
    parser.add_argument("--keep-kmers", action="store_true")
    parser.add_argument("--nvcc", default=os.environ.get("NVCC", "nvcc"))
    parser.add_argument("--force-rebuild", action="store_true")
    parser.add_argument(
        "--extractor",
        choices=("rust", "python"),
        default="rust",
        help="k-mer extraction path; rust uses the crate parser/k-mer semantics",
    )
    parser.add_argument(
        "--gzip-threads",
        type=int,
        default=4,
        help="gzip decoder threads for the Rust extractor; 0 disables threaded gzip",
    )
    parser.add_argument(
        "--stream",
        action="store_true",
        help="pipe Rust extractor stdout directly into the CUDA helper instead of writing kmers.u64",
    )
    parser.add_argument(
        "--skip-cpu-reference",
        action="store_true",
        help="skip helper-side CPU sort/reduce checksum after verified correctness runs",
    )
    return parser.parse_args()


def opener(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii", errors="strict")
    return path.open("r", encoding="ascii", errors="strict")


def fastq_records(path: Path):
    with opener(path) as handle:
        while True:
            name = handle.readline()
            if not name:
                return
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                raise ValueError(f"truncated FASTQ: {path}")
            yield seq.strip().encode("ascii")


BASE = {
    ord("A"): 0,
    ord("a"): 0,
    ord("C"): 1,
    ord("c"): 1,
    ord("G"): 2,
    ord("g"): 2,
    ord("T"): 3,
    ord("t"): 3,
}


def emit_kmers_for_seq(seq: bytes, k: int, out) -> int:
    mask = (1 << (2 * k)) - 1
    shift = 2 * (k - 1)
    forward = 0
    reverse = 0
    valid = 0
    emitted = 0
    pack = struct.pack
    for base in seq:
        bits = BASE.get(base)
        if bits is None:
            forward = 0
            reverse = 0
            valid = 0
            continue
        forward = ((forward << 2) | bits) & mask
        reverse = (reverse >> 2) | ((3 - bits) << shift)
        valid += 1
        if valid >= k:
            out.write(pack("<Q", min(forward, reverse)))
            emitted += 1
    return emitted


def extract_kmers_python(args: argparse.Namespace, kmers_path: Path) -> tuple[int, int, float]:
    started = time.perf_counter()
    reads = 0
    kmers = 0
    with kmers_path.open("wb") as out:
        r1_iter = fastq_records(args.r1)
        r2_iter = fastq_records(args.r2) if args.r2 else None
        for seq1 in r1_iter:
            if reads >= args.reads:
                break
            kmers += emit_kmers_for_seq(seq1, args.k, out)
            if r2_iter is not None:
                try:
                    seq2 = next(r2_iter)
                except StopIteration as exc:
                    raise ValueError(f"{args.r2} has fewer records than {args.r1}") from exc
                kmers += emit_kmers_for_seq(seq2, args.k, out)
            reads += 1
    elapsed = time.perf_counter() - started
    return reads, kmers, elapsed


def parse_tsv(stdout: str) -> dict[str, str]:
    parsed = {}
    for line in stdout.splitlines():
        if not line.strip() or "\t" not in line:
            continue
        key, value = line.split("\t", 1)
        parsed[key] = value
    return parsed


def rust_extractor_binary(args: argparse.Namespace) -> Path:
    binary = ROOT / "target" / "release" / "examples" / "cuda_kmer_extract"
    if args.force_rebuild or not binary.exists():
        subprocess.run(
            ["cargo", "build", "--release", "--example", "cuda_kmer_extract"],
            cwd=ROOT,
            check=True,
        )
    return binary


def extract_kmers_rust(args: argparse.Namespace, kmers_path: Path) -> tuple[int, int, float]:
    binary = rust_extractor_binary(args)
    cmd = [
        str(binary),
        "--r1",
        str(args.r1),
        "--out",
        str(kmers_path),
        "--reads",
        str(args.reads),
        "--k",
        str(args.k),
        "--gzip-threads",
        str(args.gzip_threads),
    ]
    if args.r2:
        cmd.extend(["--r2", str(args.r2)])
    completed = subprocess.run(
        cmd,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        sys.stderr.write(completed.stdout)
        sys.stderr.write(completed.stderr)
        raise subprocess.CalledProcessError(completed.returncode, cmd)
    parsed = parse_tsv(completed.stdout)
    return (
        int(parsed["reads"]),
        int(parsed["extracted_kmers"]),
        float(parsed["extract_seconds"]),
    )


def extract_kmers(args: argparse.Namespace, kmers_path: Path) -> tuple[int, int, float]:
    if args.extractor == "rust":
        return extract_kmers_rust(args, kmers_path)
    return extract_kmers_python(args, kmers_path)


def rust_extractor_command(args: argparse.Namespace, out: str) -> list[str]:
    binary = rust_extractor_binary(args)
    cmd = [
        str(binary),
        "--r1",
        str(args.r1),
        "--out",
        out,
        "--reads",
        str(args.reads),
        "--k",
        str(args.k),
        "--gzip-threads",
        str(args.gzip_threads),
    ]
    if args.r2:
        cmd.extend(["--r2", str(args.r2)])
    return cmd


def build_cuda_binary(args: argparse.Namespace, outdir: Path) -> Path:
    source = outdir / "cuda_kmer_sort_reduce.cu"
    binary = outdir / "cuda_kmer_sort_reduce"
    if args.force_rebuild or not binary.exists() or not source.exists():
        source.write_text(CUDA_SOURCE, encoding="utf-8")
        cmd = [
            args.nvcc,
            "-O3",
            "-std=c++17",
            str(source),
            "-o",
            str(binary),
        ]
        subprocess.run(cmd, cwd=ROOT, check=True)
    return binary


def cuda_helper_command(args: argparse.Namespace, binary: Path, input_name: str) -> list[str]:
    cmd = [str(binary), input_name]
    if args.skip_cpu_reference:
        cmd.append("--skip-cpu-reference")
    return cmd


def run_cuda_helper(
    args: argparse.Namespace, binary: Path, kmers_path: Path
) -> tuple[int, str, str, float]:
    started = time.perf_counter()
    completed = subprocess.run(
        cuda_helper_command(args, binary, str(kmers_path)),
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    return (
        completed.returncode,
        completed.stdout,
        completed.stderr,
        time.perf_counter() - started,
    )


def stream_extract_to_cuda(args: argparse.Namespace, binary: Path) -> tuple[int, int, float, int, str, str, float, str]:
    if args.extractor != "rust":
        raise SystemExit("--stream currently requires --extractor rust")
    extractor_cmd = rust_extractor_command(args, "-")
    started = time.perf_counter()
    helper = subprocess.Popen(
        cuda_helper_command(args, binary, "-"),
        cwd=ROOT,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert helper.stdin is not None
    extractor = subprocess.Popen(
        extractor_cmd,
        cwd=ROOT,
        stdout=helper.stdin,
        stderr=subprocess.PIPE,
    )
    helper.stdin.close()
    _, extractor_stderr = extractor.communicate()
    helper_stdout, helper_stderr = helper.communicate()
    helper_wall = time.perf_counter() - started

    extractor_report = extractor_stderr.decode("utf-8", errors="replace")
    if extractor.returncode != 0:
        sys.stderr.write(extractor_report)
        raise subprocess.CalledProcessError(extractor.returncode, extractor_cmd)
    parsed = parse_tsv(extractor_report)
    return (
        int(parsed["reads"]),
        int(parsed["extracted_kmers"]),
        float(parsed["extract_seconds"]),
        helper.returncode,
        helper_stdout.decode("utf-8", errors="replace"),
        helper_stderr.decode("utf-8", errors="replace"),
        helper_wall,
        extractor_report,
    )


def main() -> int:
    args = parse_args()
    if args.k <= 0 or args.k > 31:
        raise SystemExit("--k must be in 1..31 for this short-kmer CUDA probe")
    if not args.r1.exists():
        raise SystemExit(f"missing R1 input: {args.r1}")
    if args.r2 and not args.r2.exists():
        raise SystemExit(f"missing R2 input: {args.r2}")
    if args.stream and args.keep_kmers:
        raise SystemExit("--stream and --keep-kmers are mutually exclusive")
    args.outdir.mkdir(parents=True, exist_ok=True)
    kmers_path = args.outdir / f"kmers_k{args.k}_reads{args.reads}.u64"

    binary = build_cuda_binary(args, args.outdir)
    if args.stream:
        (
            reads,
            kmers,
            extract_seconds,
            helper_returncode,
            helper_stdout,
            helper_stderr,
            helper_wall,
            _extractor_report,
        ) = stream_extract_to_cuda(args, binary)
    else:
        reads, kmers, extract_seconds = extract_kmers(args, kmers_path)
        helper_returncode, helper_stdout, helper_stderr, helper_wall = run_cuda_helper(
            args, binary, kmers_path
        )
    pipeline_wall = helper_wall if args.stream else extract_seconds + helper_wall
    report = [
        f"r1\t{args.r1}",
        f"r2\t{args.r2 if args.r2 else ''}",
        f"extractor\t{args.extractor}",
        f"stream\t{str(args.stream).lower()}",
        f"skip_cpu_reference\t{str(args.skip_cpu_reference).lower()}",
        f"k\t{args.k}",
        f"reads\t{reads}",
        f"extracted_kmers\t{kmers}",
        f"extract_seconds\t{extract_seconds:.6f}",
        f"pipeline_wall_seconds\t{pipeline_wall:.6f}",
        f"cuda_helper_wall_seconds\t{helper_wall:.6f}",
        "",
        helper_stdout.strip(),
    ]
    if helper_stderr.strip():
        report.extend(["", "stderr:", helper_stderr.strip()])
    (args.outdir / "report.tsv").write_text("\n".join(report) + "\n", encoding="utf-8")
    sys.stdout.write("\n".join(report) + "\n")

    if not args.keep_kmers and not args.stream:
        kmers_path.unlink(missing_ok=True)
    return helper_returncode


if __name__ == "__main__":
    raise SystemExit(main())
