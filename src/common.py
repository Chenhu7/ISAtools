import pandas as pd
import numpy as np
import multiprocessing
import re
import subprocess
import os
import gc
import glob
from src.gene_grouping import GeneClustering

def run_bam2ssc(reference, bam, output_ssc,num_threads):
    """
    bam2ssc
    """
    current_dir = os.path.dirname(os.path.realpath(__file__))
    bam2SSC_script = os.path.join(current_dir, 'bam2ssc.py')
    cmd = ["python", bam2SSC_script,
        "-r", reference,
        "-b", *bam,
        "-o", output_ssc,
        "-t", str(num_threads)]
    
    subprocess.run(cmd, check=True)

def run_Ref2SSC(gtf_anno, output,num_threads):
    """
    anno2ssc
    """
    process_dir = os.path.join(output, "temp")
    os.makedirs(process_dir, exist_ok=True)

    output_SSC = os.path.join(process_dir, "anno.ssc")
    current_dir = os.path.dirname(os.path.realpath(__file__))
    gtf2SSC_script = os.path.join(current_dir,'gtf2ssc.py')
    cmd = ["python", gtf2SSC_script,
        "-i", gtf_anno,
        "-o", output_SSC,
        "-w", str(num_threads)]
    subprocess.run(cmd, check=True)

def process_data(flnc_path, count_path, df_raw_path, min_aln_coverage=None, min_aln_identity=None):
    """
    read preprocessing
    """
    dtypes_flnc = {
        1: "category",  # Chr
        2: "category",  # Strand
        3: "int32",     # TrStart_reads
        4: "int32",     # TrEnd_reads
        5: str,         # SSC
        6: "float32",   # identity
        7: "float32"    # coverage
    }

    df_flnc = pd.read_csv(
        flnc_path,
        sep="\t",
        header=None,
        dtype=dtypes_flnc,
        usecols=[1, 2, 3, 4, 5, 6, 7],low_memory=True
    )
    df_flnc.columns = ["Chr", "Strand", "TrStart_reads", "TrEnd_reads", "SSC", "identity", "coverage"]

    if min_aln_coverage is not None and min_aln_identity is not None:
        df_flnc = df_flnc[(df_flnc["identity"] >= min_aln_identity) & (df_flnc["coverage"] >= min_aln_coverage)]
    
    df_flnc = df_flnc.drop(columns=["identity", "coverage"]).dropna()

    df_grouped = (
        df_flnc.groupby(["Chr", "Strand", "SSC"], observed=True)
        .agg({"TrStart_reads": list, "TrEnd_reads": list})
        .reset_index()
    )

    df_grouped["TrStart_reads"] = df_grouped["TrStart_reads"].apply(lambda x: np.array(x, dtype=np.int32))
    df_grouped["TrEnd_reads"] = df_grouped["TrEnd_reads"].apply(lambda x: np.array(x, dtype=np.int32))
    df_grouped["frequency"] = df_grouped["TrStart_reads"].apply(len)
    df_grouped.to_parquet(df_raw_path)

    del df_flnc
    gc.collect()

    df_junction = pd.read_csv(
        count_path,
        sep="\t",
        header=None,
        dtype={1: "category", 2: "category", 3: str, 4: str},
        usecols=[1, 2, 3, 4]
    )

    df_junction.columns = ["Chr", "Strand", "SSC", "junction"]
    df = df_grouped.merge(df_junction, on=["Chr", "Strand", "SSC"], how="inner").dropna()

    del df_grouped
    gc.collect()

    return df

def load_data(reference, bam, output, num_threads, min_aln_coverage, min_aln_identity):
    """
    load ssc data
    """
    process_dir = os.path.join(output, "temp")
    os.makedirs(process_dir, exist_ok=True)
    sample = os.path.splitext(os.path.basename(bam))[0] 
    output_flnc = os.path.join(process_dir, f"{sample}_flnc.ssc")
    output_count = os.path.join(process_dir, f"{sample}_ssc.count")
    df_raw_path = os.path.join(process_dir, f"{sample}.ssc_flnc.parquet")

    df = process_data(
        flnc_path=output_flnc,
        count_path=output_count,
        df_raw_path=df_raw_path,
        min_aln_coverage=min_aln_coverage,
        min_aln_identity=min_aln_identity
    )

    return df

def junction_screening(df, junction_freq_ratio, conservative_base=None):
    """
    fiter non-canonical splice motifs
    """
    if conservative_base is None:
        conservative_base = {'GT-AG', 'AT-AC', 'GC-AG'}
    else:
        conservative_base = set(conservative_base.split(','))
    
    df = df.copy()
    def contains_non_conservative(junction_str):
        junctions = {junc.upper() for junc in junction_str.split(',')}
        return not junctions.issubset(conservative_base)

    df['contains_non_conservative'] = df['junction'].apply(contains_non_conservative)
    df['Group_freq'] = df.groupby(['Chr', 'Strand', 'Group'], observed=True)['frequency'].transform('sum')
    df['freq_ratio'] = df['frequency'] / df['Group_freq']

    df = df[~((df['contains_non_conservative']) & (df['freq_ratio'] <= junction_freq_ratio))]
    df.drop(columns=['contains_non_conservative', 'Group_freq', 'freq_ratio'], inplace=True)

    return df


def filter_fragmentary_transcript(df, threshold_fragmentary_transcript_bp = 50):
    """
    filter fragmentary transcript
    """
    conservative_base = {'GT-AG', 'AT-AC', 'GC-AG'}

    def contains_non_conservative(junction_str):
        return not set(junction_str.split(',')).issubset(conservative_base)
    
    df['contains_non_conservative'] = df['junction'].apply(contains_non_conservative)
    total_meanfreq = df['frequency'].sum() / len(df)

    df['TrStart_mean'] = df['TrStart_reads'].apply(np.mean)
    df['TrEnd_mean'] = df['TrEnd_reads'].apply(np.mean)
    df['SSC2'] = df['SSC'].apply(lambda x: list(map(int, x.split('-'))))
    
    df['Tr_length_min'] = df.apply(
        lambda row: np.inf if len(row['SSC2']) > 2 else min([
            (row['SSC2'][0] - row['TrStart_mean']),
            (row['TrEnd_mean'] - row['SSC2'][1])
        ]),
        axis=1
    )

    df = df[~(((df['Tr_length_min'] < threshold_fragmentary_transcript_bp) & (df['frequency'] < 0.01 * total_meanfreq)) |
              (df['contains_non_conservative']) & (df['frequency'] < 0.01 * total_meanfreq))]
    
    return df
