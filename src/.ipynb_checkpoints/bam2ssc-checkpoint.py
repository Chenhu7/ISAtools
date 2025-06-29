#!/usr/bin/env python3
import argparse
import os
import sys
import gzip
import multiprocessing as mp
import tempfile
import subprocess
import pysam
from Bio.Seq import Seq
from collections import defaultdict
from functools import partial

# 设置日志
logging.basicConfig(level=logging.DEBUG, format='[DEBUG] %(message)s')
logger = logging.getLogger()

def parse_args():
    parser = argparse.ArgumentParser(description='Process BAM files in parallel.')
    parser.add_argument('fasta', help='Input FASTA file (can be gzipped)')
    parser.add_argument('bam_list', help='Comma-separated list of BAM files')
    parser.add_argument('out_dir', help='Output directory')
    parser.add_argument('threads', type=int, default=4, nargs='?', help='Number of threads (default: 4)')
    return parser.parse_args()

def count_single_bam(bam, threads_per_bam):
    bai_file = bam + '.bai'
    if not os.path.exists(bai_file):
        logger.info(f'No index found for {bam}, generating index...')
        try:
            subprocess.run(['samtools', 'index', '-@', str(threads_per_bam), bam], check=True)
            logger.info(f'Index created for {bam}')
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning(f'Failed to create index for {bam}, proceeding without index')

    try:
        result = subprocess.run(['samtools', 'view', '-c', '-@', str(threads_per_bam), bam],
                               capture_output=True, text=True, check=True)
        count = int(result.stdout.strip())
        # logger.info(f'BAM {bam} has {count} reads')
        return bam, count
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.warning(f'samtools view -c failed for {bam}: {e}, falling back to pysam')
        with pysam.AlignmentFile(bam, 'rb', threads=threads_per_bam) as bf:
            count = bf.count()
        logger.info(f'BAM {bam} has {count} reads (via pysam)')
        return bam, count

def get_bam_read_counts(bam_files, threads):
    num_bams = len(bam_files)
    threads_per_bam = max(1, threads // num_bams) if num_bams > 0 else 1
    # logger.info(f'Using {threads_per_bam} threads per BAM file for {num_bams} files')

    with mp.Pool(processes=num_bams) as pool:
        results = pool.starmap(count_single_bam, [(bam, threads_per_bam) for bam in bam_files])

    bam_lines = dict(results)
    total_lines = sum(bam_lines.values())
    if total_lines == 0:
        sys.exit('No BAM lines to process')
    return bam_lines, total_lines

def allocate_chunks(bam_files, bam_lines, total_lines, threads):
    chunk_allocations = []
    num_bams = len(bam_files)
    
    # 确保每个样本至少分配一个线程
    base_threads_per_bam = max(1, threads // num_bams)
    remaining_threads = threads - base_threads_per_bam * num_bams
    
    # logger.debug(f'Base threads per BAM: {base_threads_per_bam}, Remaining threads: {remaining_threads}')
    
    for bam in bam_files:
        bam_reads = bam_lines[bam]
        # 根据读取数比例分配线程，但至少分配 base_threads_per_bam
        proportional_threads = int((bam_reads / total_lines) * threads) if total_lines > 0 else 1
        chunks = max(base_threads_per_bam, proportional_threads)
        chunk_allocations.append((bam, chunks))
        logger.debug(f'Allocating {chunks} chunks for {bam} ({bam_reads} reads)')
    
    # 计算当前分配的总线程数
    total_allocated = sum(chunks for _, chunks in chunk_allocations)
    
    # 如果分配的线程总数超过可用线程，调整分配
    if total_allocated > threads:
        # logger.debug(f'Total allocated threads ({total_allocated}) exceeds available threads ({threads}), scaling down')
        # 按比例缩减
        scale_factor = threads / total_allocated
        chunk_allocations = [(bam, max(1, int(chunks * scale_factor))) for bam, chunks in chunk_allocations]
        total_allocated = sum(chunks for _, chunks in chunk_allocations)
    
    # 如果还有剩余线程，均匀分配到样本
    remaining_threads = threads - total_allocated
    if remaining_threads > 0:
        # logger.debug(f'Distributing {remaining_threads} remaining threads')
        for i in range(len(chunk_allocations)):
            if remaining_threads <= 0:
                break
            chunk_allocations[i] = (chunk_allocations[i][0], chunk_allocations[i][1] + 1)
            remaining_threads -= 1
            # logger.debug(f'Assigned extra chunk to {chunk_allocations[i][0]}')
    
    total_chunks = sum(chunks for _, chunks in chunk_allocations)
    # logger.debug(f'Total chunks across all BAMs: {total_chunks}')
    return chunk_allocations

def process_bam_chunk(bam, fasta_file, temp_dir, out_dir, threads, chunk_idx, start_read, end_read):
    bam_basename = os.path.splitext(os.path.basename(bam))[0]
    out1_tmp = os.path.join(temp_dir, f'out1_{bam_basename}_chunk_{chunk_idx}.txt')
    out2_tmp = os.path.join(temp_dir, f'out2_{bam_basename}_chunk_{chunk_idx}.txt')
    # logger.debug(f'Processing chunk {chunk_idx} for {bam} (reads {start_read} to {end_read})')

    id_count = defaultdict(int)
    ec = defaultdict(int)
    eci = defaultdict(float)
    seq_cache = {}  # 缓存序列
    processed_lines = 0
    written_out1 = 0
    written_out2 = 0

    with open(out1_tmp, 'w') as out1_fh, pysam.AlignmentFile(bam, 'rb', threads=threads) as bf, pysam.FastaFile(fasta_file) as fa:
        for i, read in enumerate(bf):
            if i < start_read:
                continue
            if i >= end_read:
                break
            if read.is_unmapped:
                continue
            if read.query_sequence is None:
                continue
            processed_lines += 1

            astrand = '-' if read.is_reverse else '+'
            xs = ts = erro = None
            for tag, value in read.tags:
                if tag == 'XS':
                    xs = value
                elif tag == 'ts':
                    ts = value
                elif tag == 'NM':
                    erro = value
            if ts and not xs and ts in ('+', '-'):
                xs = ('+' if ts == '-' else '-') if read.is_reverse else ts
            strand = xs or astrand

            pos = read.reference_start + 1
            cov = clip = gap = 0
            positions = []
            for op, length in read.cigartuples:
                if op in (4, 5):
                    clip += length
                elif op == 1:
                    cov += length
                elif op == 3:
                    end = pos + gap - 1
                    positions.extend([pos, end])
                    pos = end + length + 1
                    gap = 0
                elif op == 0:
                    cov += length
                    gap += length
                else:
                    gap += length
            end = pos + gap - 1
            positions.extend([pos, end])

            seqlen = len(read.query_sequence)
            identity = 1 - (erro / cov) if cov else 0
            coverage = (seqlen - clip) / seqlen if seqlen else 0

            id_count[read.query_name] += 1
            s1 = positions.pop(0)
            e1 = positions.pop(-1)
            str_pos = '-'.join(map(str, positions)) if positions else 'NA'
            out1_fh.write(f'{read.query_name}.m{id_count[read.query_name]}\t'
                         f'{read.reference_name}\t{strand}\t{s1}\t{e1}\t{str_pos}\t'
                         f'{identity}\t{coverage}\n')
            written_out1 += 1

            if str_pos != 'NA':
                key = f'{read.reference_name}\t{strand}\t{str_pos}'
                if key not in seq_cache:
                    b = str_pos.split('-')
                    ds = ''
                    if strand == '+':
                        for i, k1 in enumerate(b):
                            k1 = int(k1)
                            seq = fa.fetch(reference=read.reference_name, start=k1, end=k1+2) if i % 2 == 0 else fa.fetch(reference=read.reference_name, start=k1-3, end=k1-1)
                            ds += f'{seq}-' if i % 2 == 0 else f'{seq},'
                    else:
                        for i, k1 in enumerate(reversed(b)):
                            k1 = int(k1)
                            seq = fa.fetch(reference=read.reference_name, start=k1-3, end=k1-1) if i % 2 == 0 else fa.fetch(reference=read.reference_name, start=k1, end=k1+2)
                            seq = str(Seq(seq).reverse_complement())
                            ds += f'{seq}-' if i % 2 == 0 else f'{seq},'
                    ds = ds.rstrip(',')
                    seq_cache[key] = ds
                ec[key] += 1
                eci[key] += identity

    with open(out2_tmp, 'w') as out2_fh:
        for k in ec:
            # idmean = eci[k] / ec[k]
            # out2_fh.write(f'{k}\t{ec[k]}\t{idmean}\t{seq_cache[k]}\n')
            out2_fh.write(f'{k}\t{ec[k]}\t{seq_cache[k]}\n')
            written_out2 += 1

    # logger.debug(f'Chunk {chunk_idx} for {bam}: processed {processed_lines} lines, '
                 # f'wrote {written_out1} to {out1_tmp}, {written_out2} to {out2_tmp}')
    return out1_tmp, out2_tmp, bam

def merge_single_bam(bam, files, out_dir):
    bam_basename = os.path.splitext(os.path.basename(bam))[0]
    out1 = os.path.join(out_dir, f'{bam_basename}_flnc.ssc')
    out2 = os.path.join(out_dir, f'{bam_basename}_ssc.count')

    # 合并 out1
    # logger.debug(f'Merging out1 files for {bam}')
    out1_count = 0
    with open(out1, 'w') as out1_fh:
        for out1_tmp, _ in files:
            with open(out1_tmp, 'r') as in_fh:
                for line in in_fh:
                    out1_fh.write(line)
                    out1_count += 1
    # logger.debug(f'Wrote {out1_count} lines to {out1}')

    # 合并 out2
    # logger.debug(f'Merging out2 files for {bam}')
#     global_ec = defaultdict(int)
#     global_eci = defaultdict(float)
#     global_seq = {}
#     for _, out2_tmp in files:
#         with open(out2_tmp, 'r') as in_fh:
#             for line in in_fh:
#                 ref, strand, pos_str, count, idmean, ds = line.strip().split('\t', 5)
#                 key = f'{ref}\t{strand}\t{pos_str}'
#                 global_ec[key] += int(count)
#                 global_eci[key] += float(count) * float(idmean)
#                 global_seq[key] = ds

#     with open(out2, 'w') as out2_fh:
#         out2_count = 0
#         for k in sorted(global_ec):
#             idmean = global_eci[k] / global_ec[k]
#             out2_fh.write(f'{global_ec[k]}\t{k}\t{global_seq[k]}\t{idmean}\n')
#             out2_count += 1
#         # logger.debug(f'Wrote {out2_count} lines to {out2}')
#     return bam, out1_count, out2_count
    for _, out2_tmp in files:
        with open(out2_tmp, 'r') as in_fh:
            for line in in_fh:
                ref, strand, pos_str, count, ds = line.strip().split('\t', 4)
                key = f'{ref}\t{strand}\t{pos_str}'
                global_ec[key] += int(count)
                global_seq[key] = ds

    # 写入合并后的 out2 文件
    with open(out2, 'w') as out2_fh:
        out2_count = 0
        for k in sorted(global_ec):
            out2_fh.write(f'{global_ec[k]}\t{k}\t{global_seq[k]}\n')
            out2_count += 1
    

def merge_results(chunk_results, fasta_file, out_dir, threads):
    bam_groups = defaultdict(list)
    for out1_tmp, out2_tmp, bam in chunk_results:
        bam_groups[bam].append((out1_tmp, out2_tmp))

    with mp.Pool(processes=threads) as pool:
        results = pool.starmap(partial(merge_single_bam, out_dir=out_dir), bam_groups.items())

    # for bam, out1_count, out2_count in results:
    #     logger.info(f'Merged {bam}: {out1_count} lines in out1, {out2_count} lines in out2')

def main():
    args = parse_args()
    bam_files = args.bam_list.split(',')
    os.makedirs(args.out_dir, exist_ok=True)
    bam_lines, total_lines = get_bam_read_counts(bam_files, args.threads)
    chunk_allocations = allocate_chunks(bam_files, bam_lines, total_lines, args.threads)

    tasks = []
    for bam, chunk_count in chunk_allocations:
        total_lines = bam_lines[bam]
        lines_per_chunk = total_lines // chunk_count + 1 if chunk_count > 1 else total_lines
        for i in range(chunk_count):
            start_read = i * lines_per_chunk
            end_read = min((i + 1) * lines_per_chunk, total_lines)
            tasks.append((bam, args.fasta, tempfile.gettempdir(), args.out_dir, args.threads, i, start_read, end_read))

    with mp.Pool(processes=args.threads) as pool:
        chunk_results = pool.starmap(process_bam_chunk, tasks)

    merge_results(chunk_results, args.fasta, args.out_dir, args.threads)
    # print('Processing complete.')

if __name__ == '__main__':
    main()
