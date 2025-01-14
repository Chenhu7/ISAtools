import pandas as pd
import numpy as np
import multiprocessing
import re
import subprocess
import os


# import pandas as pd
# import multiprocessing
# import re
# import subprocess
# from Bio import SeqIO

# def load_fasta(file):
#     """
#     使用 BioPython 高效加载 FASTA 文件，返回字典 {read_id: sequence}
#     """
#     return {record.id: str(record.seq) for record in SeqIO.parse(file, "fasta")}


# def process_bam_line(line):
#     """
#     解析BAM文件中的单行记录，返回处理后的结果
#     """
#     cols = line.strip().split("\t")
#     read_id = cols[0]
#     chr_name = cols[2]
#     strand = "+" if (int(cols[1]) & 0x10) == 0 else "-"
#     cigar = cols[5]
#     sequence = cols[9]

#     positions = []
#     current_pos = int(cols[3])
#     coverage = 0
#     mismatches = 0

#     for token in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
#         length, op = int(token[0]), token[1]
#         if op in "M=X":
#             coverage += length
#             current_pos += length
#         elif op == "N":
#             positions.append((current_pos, current_pos + length - 1))
#             current_pos += length
#         elif op == "I":
#             coverage += length

#     positions.append((int(cols[3]), current_pos - 1))
#     identity = 1 - (mismatches / coverage if coverage > 0 else 0)

#     exon_chain = "-".join(f"{start}-{end}" for start, end in positions)
#     return {
#         "read_id": read_id,
#         "Chr": chr_name,
#         "Strand": strand,
#         "exonChain": exon_chain,
#         "identity": identity,
#         "coverage": coverage / len(sequence) if len(sequence) > 0 else 0,
#     }


# def parse_bam_multiprocess(bam_file, num_processes=4):
#     """
#     并行解析BAM文件，返回解析后的DataFrame和exonChain统计信息
#     """
#     ec = {}
#     eci = {}

#     # 运行samtools view命令获取BAM内容
#     process = subprocess.Popen(
#         ["samtools", "view", bam_file],
#         stdout=subprocess.PIPE,
#         text=True,
#     )

#     # 多进程处理
#     pool = multiprocessing.Pool(num_processes)
#     results = pool.imap(process_bam_line, process.stdout)

#     parsed_data = []
#     for result in results:
#         parsed_data.append(result)
#         exon_key = f"{result['Chr']}\t{result['Strand']}\t{result['exonChain']}"
#         ec[exon_key] = ec.get(exon_key, 0) + 1
#         eci[exon_key] = eci.get(exon_key, 0) + result["identity"]

#     pool.close()
#     pool.join()

#     return pd.DataFrame(parsed_data), ec, eci


# def parse_and_merge(fasta_file, bam_file, num_processes=4):
#     """
#     解析FASTA和BAM文件，并合并解析结果，返回最终的DataFrame
#     """
#     fa = load_fasta(fasta_file)
#     bam_parsed, ec, eci = parse_bam_multiprocess(bam_file, num_processes)

#     # 处理junction统计
#     junction_data = []
#     for key, count in ec.items():
#         chr_name, strand, exon_chain = key.split("\t")
#         mean_identity = eci[key] / count
#         junction_data.append([chr_name, strand, exon_chain, count, mean_identity])

#     df_junction = pd.DataFrame(
#         junction_data, columns=["Chr", "Strand", "exonChain", "junction", "identity"]
#     )

#     # 合并数据
#     merged_df = bam_parsed.merge(
#         df_junction, on=["Chr", "Strand", "exonChain"], how="left"
#     )
#     return merged_df.dropna()


def run_bam2exonchain(reference, bam, output_flnc, output_count):
    """
    调用 BAM2exonChain.pl 脚本生成 exonChain 数据文件
    """
    cmd = [
        "perl",
        "src/BAM2exonChain.pl",
        reference,
        bam,
        output_flnc,
        output_count,
    ]
    try:
        subprocess.run(cmd, check=True)
        print(f"BAM2exonChain.pl 运行成功！生成文件：\n{output_flnc}\n{output_count}")
    except subprocess.CalledProcessError as e:
        print(f"运行 BAM2exonChain.pl 时出错: {e}")
        raise


def process_data(flnc_path, count_path):
    """
    读取和处理 exonChain 数据
    """
    # 读取 flnc 文件
    df_read = pd.read_csv(flnc_path, sep="\t", header=None)
    df_read = df_read.dropna()  # 去除缺失值
    df_read = df_read.drop(columns=[0, 6, 7])  # 删除无关列
    df_read.columns = ["Chr", "Strand", "TrStart", "TrEnd", "exonChain"]

    # 按照 Chr, Strand, exonChain 分组，聚合 TrStart 和 TrEnd
    df_read = (
        df_read.groupby(["Chr", "Strand", "exonChain"], as_index=False)
        .agg({"TrStart": list, "TrEnd": list})
        .reset_index(drop=True)
    )

    # 读取 junction count 文件
    df_junction = pd.read_csv(count_path, sep="\t", header=None)
    df_junction = df_junction.drop(columns=[0])  # 删除无关列
    df_junction.columns = ["Chr", "Strand", "exonChain", "junction", "identity"]

    # 合并数据
    df = df_read.merge(df_junction, on=["Chr", "Strand", "exonChain"], how="left")
    df = df.dropna()  # 删除合并后可能出现的缺失值

    return df


def load_data(reference, bam, output):
    """
    加载 BAM 文件并生成处理后的 exonChain 数据

    参数:
        reference: 参考基因组文件路径（FASTA 格式）
        bam: BAM 文件路径
        output: 输出目录，用于存放生成的文件

    返回:
        df: 处理后的 Pandas DataFrame
    """
    
    # 创建输出目录
    process_dir = os.path.join(output, "process")
    os.makedirs(process_dir, exist_ok=True)

    # 输出文件路径
    output_flnc = os.path.join(process_dir, "flnc.exonChain")
    output_count = os.path.join(process_dir, "exonChain.count")

    # 调用 BAM2exonChain.pl 脚本
    run_bam2exonchain(reference, bam, output_flnc, output_count)

    # 加载和处理数据
    df = process_data(output_flnc, output_count)
    df['frequency'] = df['TrStart'].apply(len)
    
    return df



def junction_screening(df, junction_freq_ratio, conservative_base=None):
    """
    过滤掉包含非保守碱基的isoform。
    参数:
        df: 输入的DataFrame。
        junction_freq_ratio: 非保守碱基的频率比率阈值。
        conservative_base: 可选，保守碱基集合，默认为 {'GT-AG', 'AT-AC', 'GC-AG'}。
    返回:
        过滤后的DataFrame。
    """
    if conservative_base is None:
        conservative_base = {'GT-AG', 'AT-AC', 'GC-AG'}
    
    df = df.copy()

    # 判断junction是否包含非保守碱基
    def contains_non_conservative(junction_str):
        return not set(junction_str.split(',')).issubset(conservative_base)
    
    # 标记是否包含非保守碱基
    df['contains_non_conservative'] = df['junction'].apply(contains_non_conservative)
    
    # 计算Group频率总和和频率比率
    df['Group_freq'] = df.groupby(['Chr', 'Strand', 'Group'])['frequency'].transform('sum')
    df['freq_ratio'] = df['frequency'] / df['Group_freq']
    
    # 过滤包含非保守碱基且频率比率小于阈值的行
    df = df[~((df['contains_non_conservative']) & (df['freq_ratio'] <= junction_freq_ratio))]
    
    # 删除临时列
    df.drop(columns=['contains_non_conservative', 'Group_freq', 'freq_ratio'], inplace=True)
    return df


def identity_filtering(df, threshold_identity):
    """
    根据identity值过滤数据。
    参数:
        df: 输入的DataFrame。
        threshold_identity: identity的最小阈值。
    返回:
        过滤后的DataFrame。
    """
    return df[df['identity'] >= threshold_identity].copy()


def filter_two_exon(df, threshol_twoExon_bp = 50):
    """
    过滤具有两个外显子并满足特定条件的记录
    :param df: 输入数据框
    :param threshol_twoExon_bp: 最小外显子长度阈值（单位：bp）
    """
    conservative_base = {'GT-AG', 'AT-AC', 'GC-AG'}

    # 判断junction是否包含非保守碱基
    def contains_non_conservative(junction_str):
        return not set(junction_str.split(',')).issubset(conservative_base)
    
    # 标记是否包含非保守碱基
    df['contains_non_conservative'] = df['junction'].apply(contains_non_conservative)
    
    
    # 计算总的平均频率
    total_meanfreq = df['frequency'].sum() / len(df)

    # 计算 TrStart 和 TrEnd 的均值
    df['TrStart_mean'] = df['TrStart'].apply(np.mean)
    df['TrEnd_mean'] = df['TrEnd'].apply(np.mean)
    
    # 解析 exonChain 为整型列表
    df['exonChain2'] = df['exonChain'].apply(lambda x: list(map(int, x.split('-'))))
    
    # 计算最小转录长度 (Tr_length_min)
    df['Tr_length_min'] = df.apply(
        lambda row: np.inf if len(row['exonChain2']) > 2 else min([
            (row['exonChain2'][0] - row['TrStart_mean']),
            (row['TrEnd_mean'] - row['exonChain2'][1])
        ]),
        axis=1
    )
    
    # 过滤：外显子最小长度小于阈值，或包含非保守位点
    df = df[~(((df['Tr_length_min'] < threshol_twoExon_bp) & (df['frequency'] < 0.01 * total_meanfreq)) |
              (df['contains_non_conservative']) & (df['frequency'] < 0.01 * total_meanfreq))]
    
    return df
