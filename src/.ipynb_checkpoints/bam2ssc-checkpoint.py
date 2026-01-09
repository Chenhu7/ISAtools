#!/usr/bin/env python3
import argparse
import os
import sys
import logging
import multiprocessing as mp
import tempfile
import subprocess
import pysam
from collections import defaultdict
from functools import partial

# 配置日志
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
logger = logging.getLogger()

# 高效反向互补映射
RC_MAP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def get_rc(seq):
    return seq.translate(RC_MAP)[::-1]

def parse_args():
    parser = argparse.ArgumentParser(description='Highly Optimized BAM processor v2.')
    parser.add_argument("--reference", "-r", type=str, required=True, help="Reference FASTA")
    parser.add_argument("--bam", "-b", type=str, nargs='+', required=True, help="Sorted/indexed BAMs")
    parser.add_argument("--output", "-o", type=str, default="isatools_output", help="Output folder")
    parser.add_argument("--threads", "-t", type=int, default=4, help="Threads")
    return parser.parse_args()

def count_single_bam(bam):
    """利用 idxstats 极速获取已比对序列数"""
    try:
        result = subprocess.run(['samtools', 'idxstats', bam], capture_output=True, text=True, check=True)
        count = sum(int(line.split('\t')[2]) for line in result.stdout.strip().split('\n') if line.strip())
        return bam, count
    except:
        with pysam.AlignmentFile(bam, 'rb') as bf:
            return bam, bf.mapped

def process_bam_chunk(bam, fasta_file, temp_dir, chunk_idx, start_read, end_read):
    bam_basename = os.path.splitext(os.path.basename(bam))[0]
    out1_tmp = os.path.join(temp_dir, f'o1_{bam_basename}_{chunk_idx}.tmp')
    out2_tmp = os.path.join(temp_dir, f'o2_{bam_basename}_{chunk_idx}.tmp')

    # 使用局部变量加速访问
    id_count = defaultdict(int)
    ec = defaultdict(int)
    seq_cache = {}
    
    # 减少函数查找开销
    get_id_count = id_count.__getitem__
    
    with pysam.AlignmentFile(bam, 'rb') as bf, \
         pysam.FastaFile(fasta_file) as fa, \
         open(out1_tmp, 'w', buffering=1024*1024) as o1_fh:
        
        # 优化读取循环
        for i, read in enumerate(bf):
            if i < start_read: continue
            if i >= end_read: break
            if read.is_unmapped or read.query_sequence is None: continue

            # 快速获取标签
            tags = dict(read.tags)
            xs = tags.get('XS')
            ts = tags.get('ts')
            if ts and not xs and ts in ('+', '-'):
                xs = ('+' if ts == '-' else '-') if read.is_reverse else ts
            strand = xs or ('-' if read.is_reverse else '+')

            # CIGAR 解析优化：只关注 M/I/N/S/H/D
            pos = read.reference_start + 1
            cov = clip = gap = 0
            positions = []
            
            for op, length in read.cigartuples:
                if op == 0: # MATCH
                    cov += length
                    gap += length
                elif op == 3: # SKIP (N)
                    e_p = pos + gap - 1
                    positions.append(pos)
                    positions.append(e_p)
                    pos = e_p + length + 1
                    gap = 0
                elif op == 1: # INS
                    cov += length
                elif op == 4 or op == 5: # CLIP
                    clip += length
                else: # DEL etc
                    gap += length
            
            e_p = pos + gap - 1
            positions.append(pos)
            positions.append(e_p)

            # 计算指标
            seqlen = len(read.query_sequence)
            identity = 1 - (tags.get('NM', 0) / cov) if cov else 0
            coverage = (seqlen - clip) / seqlen if seqlen else 0
            
            qname = read.query_name
            id_count[qname] += 1
            
            # 位置信息处理
            if len(positions) > 2:
                s1, e1 = positions[0], positions[-1]
                str_pos = '-'.join(map(str, positions[1:-1]))
            else:
                s1, e1 = (positions[0], positions[1]) if positions else ("NA", "NA")
                str_pos = 'NA'

            # 写入 o1
            o1_fh.write(f'{qname}.m{id_count[qname]}\t{read.reference_name}\t{strand}\t{s1}\t{e1}\t{str_pos}\t{identity:.4f}\t{coverage:.4f}\n')

            # 序列缓存逻辑优化
            if str_pos != 'NA':
                ref_n = read.reference_name
                key = (ref_n, strand, str_pos)
                ec[key] += 1
                if key not in seq_cache:
                    b = str_pos.split('-')
                    ds_list = []
                    if strand == '+':
                        for idx, k_val in enumerate(b):
                            k = int(k_val)
                            seq = fa.fetch(ref_n, k, k+2) if idx % 2 == 0 else fa.fetch(ref_n, k-3, k-1)
                            ds_list.append(seq)
                    else:
                        for idx, k_val in enumerate(reversed(b)):
                            k = int(k_val)
                            seq = fa.fetch(ref_n, k-3, k-1) if idx % 2 == 0 else fa.fetch(ref_n, k, k+2)
                            ds_list.append(get_rc(seq))
                    
                    # 组装 ds 字符串
                    res = ""
                    for j in range(0, len(ds_list), 2):
                        if j+1 < len(ds_list):
                            res += f"{ds_list[j]}-{ds_list[j+1]},"
                        else:
                            res += f"{ds_list[j]},"
                    seq_cache[key] = res.rstrip(',')

    # 写入 o2 临时文件
    with open(out2_tmp, 'w') as o2_fh:
        for k, count in ec.items():
            o2_fh.write(f'{count}\t{k[0]}\t{k[1]}\t{k[2]}\t{seq_cache[k]}\n')

    return out1_tmp, out2_tmp, bam

def merge_single_bam(bam, files, out_dir):
    bam_basename = os.path.splitext(os.path.basename(bam))[0]
    out1 = os.path.join(out_dir, f'{bam_basename}_flnc.ssc')
    out2 = os.path.join(out_dir, f'{bam_basename}_ssc.count')

    # 合并 o1：流式合并
    with open(out1, 'w', buffering=2*1024*1024) as out1_fh:
        for o1_tmp, _ in files:
            with open(o1_tmp, 'r') as f:
                while True:
                    chunk = f.read(1024*1024)
                    if not chunk: break
                    out1_fh.write(chunk)
            os.remove(o1_tmp)

    # 合并 o2：聚合计数
    global_ec = defaultdict(int)
    global_seq = {}
    for _, o2_tmp in files:
        with open(o2_tmp, 'r') as f:
            for line in f:
                c, r, s, p, ds = line.strip().split('\t')
                key = (r, s, p)
                global_ec[key] += int(c)
                global_seq[key] = ds
        os.remove(o2_tmp)

    with open(out2, 'w') as out2_fh:
        # 排序输出保证结果一致性
        for k in sorted(global_ec.keys()):
            out2_fh.write(f'{global_ec[k]}\t{k[0]}\t{k[1]}\t{k[2]}\t{global_seq[k]}\n')

def main():
    args = parse_args()
    if not os.path.exists(args.output): os.makedirs(args.output)
    
    # 1. 快速获取 BAM 计数
    with mp.Pool(processes=min(len(args.bam), args.threads)) as p:
        counts_res = p.map(count_single_bam, args.bam)
    bam_lines = dict(counts_res)
    total_reads = sum(bam_lines.values())
    
    # 2. 任务分配：改进分配算法，确保各核心负载均衡
    tasks = []
    temp_dir = os.path.join(args.output, "bam2ssc_tmp")
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    for bam in args.bam:
        b_count = bam_lines[bam]
        if b_count == 0: continue
        # 即使线程很多，每个 chunk 也不建议太小（至少 1000 行）
        n_chunks = max(1, min(b_count // 1000, int((b_count / total_reads) * args.threads) if total_reads > 0 else 1))
        r_per_chunk = (b_count // n_chunks) + 1
        for i in range(n_chunks):
            s = i * r_per_chunk
            e = min((i + 1) * r_per_chunk, b_count)
            if s < e:
                tasks.append((bam, args.reference, temp_dir, i, s, e))

    # 3. 多进程执行
    logger.info(f"Processing {len(tasks)} chunks on {args.threads} threads...")
    with mp.Pool(processes=args.threads) as pool:
        results = pool.starmap(process_bam_chunk, tasks)

    # 4. 按 BAM 分组并并行合并
    logger.info("Merging and cleaning up...")
    groups = defaultdict(list)
    for o1, o2, bam_path in results:
        groups[bam_path].append((o1, o2))
    
    with mp.Pool(processes=min(len(groups), args.threads)) as pool:
        pool.starmap(partial(merge_single_bam, out_dir=args.output), groups.items())
    
    logger.info("Process completed successfully.")

if __name__ == '__main__':
    main()
