import pandas as pd
import numpy as np
import multiprocessing
import re
import subprocess
import os
import gc
import glob
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
    process_dir = os.path.join(output, "temp")
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


def process_data(flnc_path, count_path, df_raw_path, min_aln_coverage=None, min_aln_identity=None):
    """
    使用 Pandas 处理 exonChain 数据，优化内存
    参数:
        flnc_path: flnc.exonChain 文件路径
        count_path: exonChain.count 文件路径
        min_aln_coverage: 最小覆盖率（可选）
        min_aln_identity: 最小 identity（可选）
    返回:
        处理后的 Pandas DataFrame
    """

    # 定义 flnc.exonChain 数据类型
    dtypes_flnc = {
        1: "category",  # Chr
        2: "category",  # Strand
        3: "int32",     # TrStart_reads
        4: "int32",     # TrEnd_reads
        5: str,         # exonChain
        6: "float32",   # identity
        7: "float32"    # coverage
    }

    # 读取 flnc.exonChain
    df_flnc = pd.read_csv(
        flnc_path,
        sep="\t",
        header=None,
        dtype=dtypes_flnc,
        usecols=[1, 2, 3, 4, 5, 6, 7],low_memory=True
    )
    df_flnc.columns = ["Chr", "Strand", "TrStart_reads", "TrEnd_reads", "exonChain", "identity", "coverage"]

    # 应用过滤
    if min_aln_coverage is not None and min_aln_identity is not None:
        df_flnc = df_flnc[(df_flnc["identity"] >= min_aln_identity) & (df_flnc["coverage"] >= min_aln_coverage)]
    
    # 删除无用列和 NaN
    df_flnc = df_flnc.drop(columns=["identity", "coverage"]).dropna()

    # 分组和聚合
    df_grouped = (
        df_flnc.groupby(["Chr", "Strand", "exonChain"], observed=True)
        .agg({"TrStart_reads": list, "TrEnd_reads": list})
        .reset_index()
    )
    df_grouped["TrStart_reads"] = df_grouped["TrStart_reads"].apply(lambda x: np.array(x, dtype=np.int32))
    df_grouped["TrEnd_reads"] = df_grouped["TrEnd_reads"].apply(lambda x: np.array(x, dtype=np.int32))

    df_grouped["frequency"] = df_grouped["TrStart_reads"].apply(len)

    # 将df_group写入文件
    df_grouped.to_parquet(df_raw_path)

    # 释放内存
    del df_flnc
    gc.collect()

    # 读取 exonChain.count
    df_junction = pd.read_csv(
        count_path,
        sep="\t",
        header=None,
        dtype={1: "category", 2: "category", 3: str, 4: str},
        usecols=[1, 2, 3, 4]
    )
    df_junction.columns = ["Chr", "Strand", "exonChain", "junction"]

    # 合并
    df = df_grouped.merge(df_junction, on=["Chr", "Strand", "exonChain"], how="inner").dropna()

    del df_grouped
    gc.collect()
    # return df,df_grouped
    return df


def load_data(reference, bam, output, num_threads, min_aln_coverage, min_aln_identity):
    """
    加载BAM文件并生成处理后的exonChain数据
    参数:
        reference: 参考基因组文件路径
        bam: BAM文件路径
        output: 输出目录
        num_threads: 线程数
        min_aln_coverage: 最小覆盖率
        min_aln_identity: 最小identity
    返回:
        df: 过滤后的DataFrame
        df_raw: 未过滤的DataFrame（延迟创建）
    """
    process_dir = os.path.join(output, "temp")
    os.makedirs(process_dir, exist_ok=True)
    sample = os.path.splitext(os.path.basename(bam))[0] 
    output_flnc = os.path.join(process_dir, f"{sample}.flnc.exonChain")
    output_count = os.path.join(process_dir, f"{sample}.exonChain.count")
    df_raw_path = os.path.join(process_dir, f"{sample}.exonChain_flnc.parquet")

    # 调用Perl脚本
    run_bam2exonchain(reference, bam, output_flnc, output_count, num_threads)

    # 处理过滤后的数据
    df = process_data(
        flnc_path=output_flnc,
        count_path=output_count,
        df_raw_path=df_raw_path,
        min_aln_coverage=min_aln_coverage,
        min_aln_identity=min_aln_identity
    )

    # return df,df_raw
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
    df['Group_freq'] = df.groupby(['Chr', 'Strand', 'Group'], observed=True)['frequency'].transform('sum')
    df['freq_ratio'] = df['frequency'] / df['Group_freq']
    
    # 过滤包含非保守碱基且频率比率小于阈值的行
    df = df[~((df['contains_non_conservative']) & (df['freq_ratio'] <= junction_freq_ratio))]
    
    # 删除临时列
    df.drop(columns=['contains_non_conservative', 'Group_freq', 'freq_ratio'], inplace=True)
    return df


def filter_fragmentary_transcript(df, threshold_fragmentary_transcript_bp = 50):
    """
    过滤具有两个外显子并满足特定条件的记录
    :param df: 输入数据框
    :param threshold_fragmentary_transcript_bp: 最小外显子长度阈值（单位：bp）
    """
    conservative_base = {'GT-AG', 'AT-AC', 'GC-AG'}

    # 判断junction是否包含非保守碱基
    def contains_non_conservative(junction_str):
        return not set(junction_str.split(',')).issubset(conservative_base)
    
    # 标记是否包含非保守碱基
    df['contains_non_conservative'] = df['junction'].apply(contains_non_conservative)
    
    
    # 计算总的平均频率
    total_meanfreq = df['frequency'].sum() / len(df)

    # 计算 TrStart_reads 和 TrEnd_reads 的均值
    df['TrStart_mean'] = df['TrStart_reads'].apply(np.mean)
    df['TrEnd_mean'] = df['TrEnd_reads'].apply(np.mean)
    
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
    df = df[~(((df['Tr_length_min'] < threshold_fragmentary_transcript_bp) & (df['frequency'] < 0.01 * total_meanfreq)) |
              (df['contains_non_conservative']) & (df['frequency'] < 0.01 * total_meanfreq))]
    
    return df


def annotate_transcripts_novel_only(df_result):
    df = df_result.copy()
    df['Group'] = df['Chr'].astype(str) + '_' + df['Strand'].astype(str) + '_' + df['Group'].astype(str)

    Gi = 1  # Gene 编号
    gene_id_map = {}
    tr_id_map = {}

    # 全部 group 分组处理
    for _, group_df in df.groupby('Group'):
        Ti = 1  # Transcript 编号
        first_row = group_df.iloc[0]
        gene_id = f"{first_row['Chr']}.{first_row['Strand']}.novel_gene{Gi}"
        Gi += 1

        for idx in group_df.index:
            gene_id_map[idx] = gene_id
            tr_id_map[idx] = f"{gene_id}.transcript{Ti}"
            Ti += 1

    # 应用结果
    df['GeneID'] = pd.Series(gene_id_map)
    df['TrID'] = pd.Series(tr_id_map)

    return df.drop(columns = 'Group')

def annotate_transcripts(df_result, anno_exonChain):
    """
    给 df_result 添加 TrID 和 GeneID 列
    1. 若 exonChain 与 anno_exonChain 完全匹配，选起止点误差最小的注释
    2. 对未注释的按 Group 补全：
       - 有已知 GeneID：共享 GeneID，TrID 为 .nnic_transcript{Ti}
       - 全新 group：生成 GeneID 为 chr.strand.gene{Gi}，TrID 为 novel_gene_transcript{Ti}
    """
    df = df_result.copy()
    anno = anno_exonChain.copy()

    # === 创建匹配索引 ===
    anno_key = anno[['Chr', 'Strand', 'exonChain', 'TrStart', 'TrEnd', 'TrID', 'GeneID']].copy()
    anno_key['match_key'] = anno_key['Chr'].astype(str) + '|' + anno_key['Strand'] + '|' + anno_key['exonChain']
    df['match_key'] = df['Chr'].astype(str) + '|' + df['Strand'] + '|' + df['exonChain']

    # === 合并匹配注释 ===
    merged = df.merge(anno_key, on='match_key', suffixes=('', '_ref'))
    merged['diff'] = (merged['TrStart'] - merged['TrStart_ref']).abs() + (merged['TrEnd'] - merged['TrEnd_ref']).abs()
    merged_sorted = merged.sort_values(['match_key', 'diff'])
    best_match = merged_sorted.drop_duplicates(subset=df.columns.tolist(), keep='first')

    # === 合并注释结果 ===
    df = df.merge(
        best_match[['Chr', 'Strand', 'exonChain', 'TrStart', 'TrEnd', 'TrID', 'GeneID']],
        on=['Chr', 'Strand', 'exonChain', 'TrStart', 'TrEnd'],
        how='left'
    )

    # === 处理未注释 ===
    Gi = 1
    mask_none = df['GeneID'].isna()
    df['Group'] = df['Chr'].astype(str) + df['Strand'].astype(str) + df['Group'].astype(str)
    grouped = df[mask_none].groupby('Group')

    new_gene_ids = {}
    new_tr_ids = {}

    for group, sub_df in grouped:
        Ti = 1
        indices = sub_df.index.tolist()
        annotated_gene_ids = df.loc[(df['Group'] == group) & (~df['GeneID'].isna()), 'GeneID'].unique()

        if len(annotated_gene_ids) > 0:
            gene_id = annotated_gene_ids[0]
            for idx in indices:
                new_gene_ids[idx] = gene_id
                new_tr_ids[idx] = f"{gene_id}.novel_transcript{Ti}"
                Ti += 1
        else:
            row = sub_df.iloc[0]
            gene_id = f"{row['Chr']}.{row['Strand']}.novel_gene{Gi}"
            Gi += 1
            for idx in indices:
                new_gene_ids[idx] = gene_id
                new_tr_ids[idx] = f"{gene_id}.transcript{Ti}"
                Ti += 1

    # === 应用注释 ===
    df.loc[list(new_gene_ids.keys()), 'GeneID'] = pd.Series(new_gene_ids)
    df.loc[list(new_tr_ids.keys()), 'TrID'] = pd.Series(new_tr_ids)
    df = df.drop(columns=['match_key'])

    return df.drop(columns = 'Group')


def convert_to_gtf(df_result, output_file):
    # 获取转录本所有位点
    df_result['sites'] = df_result.apply(lambda row: sorted(
        list(map(int, row['exonChain'].split('-'))) + [int(row['TrStart']),int(row['TrStart'])]
    ), axis=1)

    # 只保留所需列
    df_result = df_result[['Chr', 'Strand', 'sites','TrID','GeneID']].copy()

    # 生成 GTF 文件内容
    gtf_lines = []
    for idx, row in df_result.iterrows():
        chr_name = row["Chr"]
        strand = row["Strand"]
        sites = row["sites"]
        transcript_id = row["TrID"]
        gene_id = row["GeneID"]
        
        gtf_lines.append(
            f'{chr_name}\tISAtools\ttranscript\t{sites[0]}\t{sites[-1]}\t.\t{strand}\t.\tgene_id "{gene_id}"; transcript_id "{transcript_id}";'
        )

        # 生成 exon 记录
        for i in range(0, len(sites), 2):
            exon_start = sites[i]
            exon_end = sites[i+1]
            exon_number = i // 2 + 1
            gtf_lines.append(
                f'{chr_name}\tISAtools\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\tgene_id "{gene_id}"; transcript_id "{transcript_id}"; exon_number "{exon_number}";'
            )

    # 保存 GTF 文件
    with open(output_file, "w") as f:
        f.write("\n".join(gtf_lines) + "\n")
        

def compute_tr_length(sites):
    return sum(sites[i+1] - sites[i] for i in range(0, len(sites), 2))

def prepare_quantification(df):
    df = df[['GeneID', 'TrID', 'quantification', 'sites']].copy()
    df['length'] = df['sites'].apply(compute_tr_length)
    df['rpk'] = df['quantification'] / (df['length'] / 1000)
    sum_rpk = df['rpk'].sum()
    df['TPM'] = (df['rpk'] / sum_rpk * 1e6).round(2)
    return df

def save_quantification(df_result, output_dir, sample=None):
    df_quant = prepare_quantification(df_result)

    # 动态文件名前缀
    prefix = f'temp/isatools.{sample}_' if sample else 'isatools_'

    # transcript counts
    df_quant[['GeneID', 'TrID', 'quantification']] \
        .rename(columns={'quantification': 'count'}) \
        .to_csv(f'{output_dir}/{prefix}transcript_counts.tsv', sep='\t', index=False)

    # transcript tpm
    df_quant[['GeneID', 'TrID', 'TPM']] \
        .to_csv(f'{output_dir}/{prefix}transcript_tpm.tsv', sep='\t', index=False)

    # gene counts
    df_gene_count = df_quant.groupby('GeneID', as_index=False)['quantification'] \
        .sum().rename(columns={'quantification': 'count'})
    df_gene_count.to_csv(f'{output_dir}/{prefix}gene_counts.tsv', sep='\t', index=False)

    # gene tpm
    df_gene_tpm = df_quant.groupby('GeneID', as_index=False)['TPM'].sum().round(2)
    df_gene_tpm.to_csv(f'{output_dir}/{prefix}gene_tpm.tsv', sep='\t', index=False)
    

def save_anno_quantification_single(df_result,output_dir,gtf = None):
    if gtf:
        anno_exonChain = pd.read_csv(f"{output_dir}/temp/anno.exonChain",sep='\t')
        df_result = annotate_transcripts(df_result, anno_exonChain)
    else:
        df_result = annotate_transcripts_novel_only(df_result)
    
    # 保存SSC
    df_result[['Chr','Strand','GeneID','TrID','quantification','TrStart','TrEnd','exonChain']].to_csv(f'{output_dir}/isatools_exonChain.tsv', sep='\t',index = False)
    
    # 保存gtf
    convert_to_gtf(df_result, f'{output_dir}/isatools_transcript.gtf')
    
    # 保存定量
    save_quantification(df_result,output_dir)

    
def extract_sample_name(filename, type_suffix):
    basename = os.path.basename(filename)
    pattern = rf'^isatools\.(.+?)_{type_suffix}\.tsv$'
    match = re.match(pattern, basename)
    return match.group(1) if match else basename


def merge_expression_files(pattern, id_cols, type_suffix, output_file):
    files = glob.glob(pattern)
    combined_df = None

    for file in files:
        sample_name = extract_sample_name(file, type_suffix)
        df = pd.read_csv(file, sep="\t")

        value_col = [col for col in df.columns if col not in id_cols][0]  # 自动找 count 或 TPM 列
        df = df.rename(columns={value_col: sample_name})

        if combined_df is None:
            combined_df = df
        else:
            combined_df = pd.merge(combined_df, df, on=id_cols, how="outer")

    combined_df = combined_df.fillna(0)

    combined_df.to_csv(output_file, sep="\t", index=False)


def save_anno_quantification_multi(df_result,output_dir,sample,gtf = None):
    if gtf:
        anno_exonChain = pd.read_csv(f"{output_dir}/temp/anno.exonChain",sep='\t')
        df_result = annotate_transcripts(df_result, anno_exonChain)
    else:
        df_result = annotate_transcripts_novel_only(df_result)
    
    # 保存SSC
    df_result[['Chr','Strand','GeneID','TrID','quantification','TrStart','TrEnd','exonChain']].to_csv(f'{output_dir}/isatools.{sample}_exonChain.tsv', sep='\t',index = False)
    
    # 保存gtf
    convert_to_gtf(df_result, f'{output_dir}/isatools.{sample}_transcript.gtf')
    
    # 保存多样本定量
    save_quantification(df_result, output_dir, sample)
    
    # 合并多样本定量
    merge_expression_files(f"{output_dir}/temp/*_transcript_counts.tsv", ["GeneID", "TrID"], "transcript_counts", f"{output_dir}/combined_transcript_counts.tsv")
    merge_expression_files(f"{output_dir}/temp/*_transcript_tpm.tsv",    ["GeneID", "TrID"], "transcript_tpm",    f"{output_dir}/combined_transcript_tpm.tsv")
    merge_expression_files(f"{output_dir}/temp/*_gene_counts.tsv",       ["GeneID"],        "gene_counts",       f"{output_dir}/combined_gene_counts.tsv")
    merge_expression_files(f"{output_dir}/temp/*_gene_tpm.tsv",          ["GeneID"],        "gene_tpm",         f"{output_dir}/combined_gene_tpm.tsv")
    
    



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
    
    # df_return = df_return.explode(['TrStart_reads', 'TrEnd_reads'], ignore_index=True)
    # df_return = df_return.groupby(['Chr', 'Strand', 'Group', 'exonChain'],as_index=False).agg({'TrStart_reads': list,'TrEnd_reads': list}).reset_index(drop=True)
    def merge_TssTes(series):
        return sum(series, [])
    df_return.groupby(['Chr', 'Strand', 'Group', 'exonChain'],as_index=False, observed=True).agg({'TrStart_reads': merge_TssTes,'TrEnd_reads': merge_TssTes}).reset_index(drop=True)
    
    df_return['frequency'] = df_return['TrStart_reads'].apply(len)
    
    return df_return
