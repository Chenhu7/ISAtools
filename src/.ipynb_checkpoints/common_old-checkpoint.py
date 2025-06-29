import pandas as pd
import numpy as np
import multiprocessing
import re
import subprocess
import os
from src.gene_grouping import GeneClustering


def run_bam2exonchain(reference, bam, output_flnc, output_count,num_threads):
    """
    调用 BAM2exonChain.pl 脚本生成 exonChain 数据文件
    """
    current_dir = os.path.dirname(os.path.realpath(__file__))  # 获取当前脚本所在路径
    bam2exonchain_script = os.path.join(current_dir, 'BAM2exonChain.pl')
    cmd = [
        "perl",
        bam2exonchain_script,
        reference,
        bam,
        output_flnc,
        output_count,
        str(num_threads)
    ]
    # print(cmd)
    try:
        subprocess.run(cmd, check=True)
        # print(f"BAM2exonChain.pl 运行成功！生成文件：\n{output_flnc}\n{output_count}")
    except subprocess.CalledProcessError as e:
        print(f"运行 BAM2exonChain.pl 时出错: {e}")
        raise


def run_Ref2exonChain(gtf_anno, output):
    """
    调用 processRef2exonChain.pl 脚本生成 exonChain 数据文件
    """
    # 创建输出目录
    process_dir = os.path.join(output, "process")
    os.makedirs(process_dir, exist_ok=True)

    # 输出文件路径
    output_exonchain = os.path.join(process_dir, "anno.exonChain")

    # 构造命令
    current_dir = os.path.dirname(os.path.realpath(__file__))
    process_script = os.path.join(current_dir,'processRef2exonChain.pl')
    cmd = f"perl {process_script} {gtf_anno} | awk '$8!=\"NA\"' > {output_exonchain}"
    
    try:
        # 执行命令
        subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
        # print(f"processRef2exonChain.pl 运行成功！生成文件：{output_exonchain}")
    except subprocess.CalledProcessError as e:
        print(f"运行 processRef2exonChain.pl 时出错: {e}")
        raise


def process_data_filter(flnc_path, count_path, min_aln_coverage, min_aln_identity):
    """
    读取和处理、过滤 exonChain 数据
    """
    # 读取 flnc 文件
    df_read = pd.read_csv(flnc_path, sep="\t", header=None)
    df_read = df_read[(df_read[6] >= min_aln_identity) & (df_read[7] >= min_aln_coverage)]
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
    df_junction.columns = ["Chr", "Strand", "exonChain", "junction"]

    # 合并数据
    df = df_read.merge(df_junction, on=["Chr", "Strand", "exonChain"], how="left")
    df = df.dropna()  # 删除合并后可能出现的缺失值

    return df


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
    df_junction.columns = ["Chr", "Strand", "exonChain", "junction"]

    # 合并数据
    df = df_read.merge(df_junction, on=["Chr", "Strand", "exonChain"], how="left")
    df = df.dropna()  # 删除合并后可能出现的缺失值

    return df


def load_data(reference, bam, output,num_threads, min_aln_coverage, min_aln_identity):
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
    run_bam2exonchain(reference, bam, output_flnc, output_count,num_threads)

    # 加载和处理数据
    df = process_data_filter(output_flnc, output_count,min_aln_coverage, min_aln_identity)
    df['frequency'] = df['TrStart'].apply(len)

    df_raw = process_data(output_flnc, output_count)
    df_raw['frequency'] = df_raw['TrStart'].apply(len)
    
    return df,df_raw



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
    else:
        conservative_base = set(conservative_base.split(','))
    
    df = df.copy()

    # 判断junction是否包含非保守碱基 忽略大小写
    def contains_non_conservative(junction_str):
        junctions = {junc.upper() for junc in junction_str.split(',')}
        return not junctions.issubset(conservative_base)
    
    # # 判断junction是否包含非保守碱基
    # def contains_non_conservative(junction_str):
    #     return not set(junction_str.split(',')).issubset(conservative_base)
    
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


def correct_NNC(df,ref_anno,num_processes = 20):
    """
    矫正单个位点错配10bp内的NNC--> FSM 或(NNC截断--> ISM)
    """
    df_anno = ref_anno.copy()
    df = df.copy()
    df_anno['source'] = 'ref'
    df['source'] = 'data'
    df_anno['exonChain_sites'] = df_anno['exonChain'].apply(lambda x: list(map(int, x.split('-'))))
    df['exonChain_sites'] = df['exonChain'].apply(lambda x: list(map(int, x.split('-'))))
    
    df_merged = pd.concat([df, df_anno], axis=0, ignore_index=True)
    
    # 将data和ref聚类
    gene_clustering = GeneClustering(num_processes=num_processes)
    df_merged = gene_clustering.cluster(df_merged)
    
    # 去除全是ref的group
    df_merged = df_merged.groupby(['Chr', 'Strand', 'Group']).filter(lambda group: not (group['source'] == 'ref').all())
    
    # 只留下NNC
    df_merged2 = df_merged[(df_merged['category'] == "NNC") | (df_merged['category'].isna())]
    df_merged2 = df_merged2.groupby(['Chr', 'Strand', 'Group']).filter(lambda group: not (group['source'] == 'ref').all())
    
    # 进行矫正
    for _,df_group in df_merged2.groupby(['Chr', 'Strand', 'Group']):
        df_data = df_group[df_group['source'] == 'data'].copy()
        df_ref = df_group[df_group['source'] == 'ref'].copy()
        
        ref_siets = list(set(df_ref['exonChain_sites'].explode().tolist()))
        
        for index,row in df_data.iterrows():
            # print(len(row['exonChain_sites']) - len(set(row['exonChain']) & set(ref_siets)))
            if len(row['exonChain_sites']) - len(set(row['exonChain_sites']) & set(ref_siets)) == 1:
                data_site = list(set(row['exonChain_sites']) - set(ref_siets))[0]
                
                # 计算差值，提取差值最小的ref位点 并替换
                differences = [abs(data_site - ref_site) for ref_site in ref_siets]
                
                # 若含有差值小于10的ref sites则替换
                if min(differences) <= 10 and data_site not in ref_siets:
                    # 找到最小差值的索引
                    min_diff_index = differences.index(min(differences))
                    # 获取最小差值的数
                    ref_site = ref_siets[min_diff_index]
                    
                    # df_merged.loc[index, 'exonChain_sites'] = [ref_site if x == data_site else x for x in df_merged.loc[index, 'exonChain_sites']]
                    df_merged.loc[index, 'exonChain'] = df_merged.loc[index, 'exonChain'].replace(str(data_site),str(ref_site))

    df_return = df_merged[df_merged['source'] == 'data'].copy()
    
    df_return = df_return.explode(['TrStart', 'TrEnd'], ignore_index=True)
    df_return = df_return.groupby(['Chr', 'Strand', 'Group', 'exonChain'],as_index=False).agg({'TrStart': list,'TrEnd': list}).reset_index(drop=True)
    df_return['frequency'] = df_return['TrStart'].apply(len)
    
    return df_return
