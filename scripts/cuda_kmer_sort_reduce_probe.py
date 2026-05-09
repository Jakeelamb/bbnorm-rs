#!/usr/bin/env python3
"""Probe whether CUDA sort/reduce is promising for BBNorm k-mer counting.

This is intentionally a standalone feasibility harness. It extracts real
canonical short k-mers from FASTQ/FASTQ.GZ inputs into a binary u64 stream, then
compiles and runs a tiny CUDA/Thrust benchmark that compares CPU sort/reduce
against GPU host-to-device copy, sort, reduce_by_key, and device-to-host copy.
"""

from __future__ import annotations

import argparse
import gzip
import os
import struct
import subprocess
import sys
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_R1 = ROOT / "tmp/human_benchmark_8threads/human_GRCh38_500k_R1.fq.gz"
DEFAULT_R2 = ROOT / "tmp/human_benchmark_8threads/human_GRCh38_500k_R2.fq.gz"

CUDA_SOURCE = r"""
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/reduce.h>
#include <thrust/sort.h>

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
                               const std::vector<uint32_t>& counts) {
    uint64_t hash = 1469598103934665603ull;
    for (size_t i = 0; i < keys.size(); ++i) {
        hash ^= keys[i];
        hash *= 1099511628211ull;
        hash ^= counts[i];
        hash *= 1099511628211ull;
    }
    return hash;
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: cuda_kmer_sort_reduce <kmers.u64>\n";
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
    input.seekg(0);
    std::vector<uint64_t> host_keys(n);
    input.read(reinterpret_cast<char*>(host_keys.data()), bytes);
    if (!input) {
        std::cerr << "failed to read input\n";
        return 2;
    }

    auto cpu_keys = host_keys;
    auto cpu_start = std::chrono::steady_clock::now();
    std::sort(cpu_keys.begin(), cpu_keys.end());
    std::vector<uint64_t> cpu_unique;
    std::vector<uint32_t> cpu_counts;
    cpu_unique.reserve(cpu_keys.size());
    cpu_counts.reserve(cpu_keys.size());
    for (size_t i = 0; i < cpu_keys.size();) {
        size_t j = i + 1;
        while (j < cpu_keys.size() && cpu_keys[j] == cpu_keys[i]) {
            ++j;
        }
        cpu_unique.push_back(cpu_keys[i]);
        cpu_counts.push_back(static_cast<uint32_t>(j - i));
        i = j;
    }
    auto cpu_end = std::chrono::steady_clock::now();
    const uint64_t cpu_checksum = checksum_pairs(cpu_unique, cpu_counts);

    cudaEvent_t h2d_start, h2d_stop, sort_start, sort_stop, reduce_start, reduce_stop, d2h_start, d2h_stop;
    check_cuda(cudaEventCreate(&h2d_start), "cudaEventCreate");
    check_cuda(cudaEventCreate(&h2d_stop), "cudaEventCreate");
    check_cuda(cudaEventCreate(&sort_start), "cudaEventCreate");
    check_cuda(cudaEventCreate(&sort_stop), "cudaEventCreate");
    check_cuda(cudaEventCreate(&reduce_start), "cudaEventCreate");
    check_cuda(cudaEventCreate(&reduce_stop), "cudaEventCreate");
    check_cuda(cudaEventCreate(&d2h_start), "cudaEventCreate");
    check_cuda(cudaEventCreate(&d2h_stop), "cudaEventCreate");

    check_cuda(cudaEventRecord(h2d_start), "cudaEventRecord h2d_start");
    thrust::device_vector<uint64_t> device_keys(host_keys.begin(), host_keys.end());
    check_cuda(cudaEventRecord(h2d_stop), "cudaEventRecord h2d_stop");
    check_cuda(cudaEventSynchronize(h2d_stop), "cudaEventSynchronize h2d_stop");

    check_cuda(cudaEventRecord(sort_start), "cudaEventRecord sort_start");
    thrust::sort(device_keys.begin(), device_keys.end());
    check_cuda(cudaEventRecord(sort_stop), "cudaEventRecord sort_stop");
    check_cuda(cudaEventSynchronize(sort_stop), "cudaEventSynchronize sort_stop");

    thrust::device_vector<uint64_t> reduced_keys(n);
    thrust::device_vector<uint32_t> reduced_counts(n);
    thrust::constant_iterator<uint32_t> ones(1);
    check_cuda(cudaEventRecord(reduce_start), "cudaEventRecord reduce_start");
    auto ends = thrust::reduce_by_key(
        device_keys.begin(), device_keys.end(),
        ones,
        reduced_keys.begin(),
        reduced_counts.begin());
    check_cuda(cudaEventRecord(reduce_stop), "cudaEventRecord reduce_stop");
    check_cuda(cudaEventSynchronize(reduce_stop), "cudaEventSynchronize reduce_stop");
    const size_t unique_count = static_cast<size_t>(ends.first - reduced_keys.begin());

    std::vector<uint64_t> gpu_unique(unique_count);
    std::vector<uint32_t> gpu_counts(unique_count);
    check_cuda(cudaEventRecord(d2h_start), "cudaEventRecord d2h_start");
    thrust::copy(reduced_keys.begin(), reduced_keys.begin() + unique_count, gpu_unique.begin());
    thrust::copy(reduced_counts.begin(), reduced_counts.begin() + unique_count, gpu_counts.begin());
    check_cuda(cudaEventRecord(d2h_stop), "cudaEventRecord d2h_stop");
    check_cuda(cudaEventSynchronize(d2h_stop), "cudaEventSynchronize d2h_stop");

    const uint64_t gpu_checksum = checksum_pairs(gpu_unique, gpu_counts);
    const bool ok = unique_count == cpu_unique.size() && gpu_checksum == cpu_checksum;
    const double h2d_ms = event_ms(h2d_start, h2d_stop);
    const double sort_ms = event_ms(sort_start, sort_stop);
    const double reduce_ms = event_ms(reduce_start, reduce_stop);
    const double d2h_ms = event_ms(d2h_start, d2h_stop);
    const double gpu_total_ms = h2d_ms + sort_ms + reduce_ms + d2h_ms;

    std::cout << "kmers\t" << n << "\n";
    std::cout << "unique\t" << unique_count << "\n";
    std::cout << "cpu_sort_reduce_ms\t" << elapsed_ms(cpu_start, cpu_end) << "\n";
    std::cout << "gpu_h2d_ms\t" << h2d_ms << "\n";
    std::cout << "gpu_sort_ms\t" << sort_ms << "\n";
    std::cout << "gpu_reduce_ms\t" << reduce_ms << "\n";
    std::cout << "gpu_d2h_ms\t" << d2h_ms << "\n";
    std::cout << "gpu_total_timed_ms\t" << gpu_total_ms << "\n";
    std::cout << "cpu_checksum\t" << cpu_checksum << "\n";
    std::cout << "gpu_checksum\t" << gpu_checksum << "\n";
    std::cout << "match\t" << (ok ? "true" : "false") << "\n";
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


def extract_kmers(args: argparse.Namespace, kmers_path: Path) -> tuple[int, int]:
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


def main() -> int:
    args = parse_args()
    if args.k <= 0 or args.k > 31:
        raise SystemExit("--k must be in 1..31 for this short-kmer CUDA probe")
    if not args.r1.exists():
        raise SystemExit(f"missing R1 input: {args.r1}")
    if args.r2 and not args.r2.exists():
        raise SystemExit(f"missing R2 input: {args.r2}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    kmers_path = args.outdir / f"kmers_k{args.k}_reads{args.reads}.u64"

    reads, kmers, extract_seconds = extract_kmers(args, kmers_path)
    binary = build_cuda_binary(args, args.outdir)
    completed = subprocess.run(
        [str(binary), str(kmers_path)],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    report = [
        f"r1\t{args.r1}",
        f"r2\t{args.r2 if args.r2 else ''}",
        f"k\t{args.k}",
        f"reads\t{reads}",
        f"extracted_kmers\t{kmers}",
        f"extract_seconds\t{extract_seconds:.6f}",
        "",
        completed.stdout.strip(),
    ]
    if completed.stderr.strip():
        report.extend(["", "stderr:", completed.stderr.strip()])
    (args.outdir / "report.tsv").write_text("\n".join(report) + "\n", encoding="utf-8")
    sys.stdout.write("\n".join(report) + "\n")

    if not args.keep_kmers:
        kmers_path.unlink(missing_ok=True)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
