import sys
import logging
from traceback import print_exc
from io import StringIO
import argparse
import os
from src.common import *
from src.consensus import ConsensusFilter
from src.gene_grouping import GeneClustering
from src.isoform_classify import IsoformClassifier
from src.remove_lowWeight_edges import GraphEdgeFilter
from src.NNC_NIC_filter import NncNicGraphProcessor
from src.filter_mutually_exclusive_exons import MutuallyExclusiveExonFilter
from src.ISM_filter import TruncationProcessor
from src.get_terminal_sites import TerminalSitesProcessor
from src.gtf_rescue_TSS import GTFRescueTSS
from src.exonChain_collapse import ExonChainCollapse
from src.filter_lowFreq_junction import FilterLowfreqJunction
from src.exonChain_assign import ExonChainAssign
from src.isoform_quantification import Quantification
import psutil

# 初始化日志
logger = logging.getLogger("ISAtools")
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def run_pipeline(args):
    """
    运行主管道流程,包括数据预处理、分组、isoform 分析等。
    """
    # def print_memory_usage(note=""):
    #     process = psutil.Process(os.getpid())
    #     mem = process.memory_info().rss / 1024 ** 2  # rss = Resident Set Size
    #     print(f"[{note}] 内存使用: {mem:.2f} MB")

    # 示例用法
    #print_memory_usage("初始 ")
    
    logger.info(" === ISAtools pipeline started === ")

    # Stage 1: 数据预处理
    logger.info("Stage 1: Preprocessing input data...")
    df = load_data(reference = args.reference, bam = args.bam, output = args.output, num_threads = args.threads, min_aln_coverage = args.min_aln_coverage, min_aln_identity = args.min_aln_identity)
    logger.info(f"Loaded input data with {len(df)} records.")
    
    #print_memory_usage("数据预处理 ")
    #print(df.columns)
    
    len_freq1 = len(df)
    
    # frequency过滤
    logger.info(f"Stage 2: Filtering frequency {args.filter_freq}...")
    df = df[df['frequency'] >= args.filter_freq]
    logger.info(f"Filtered data with {len(df)} records.")
    len_freq2 = len(df)
    
    #print_memory_usage("frequency过滤 ")
    #print(df.columns)
    
    # 基因分组
    logger.info("Stage 3: Performing gene grouping...")
    gene_clustering = GeneClustering(num_processes=args.threads)  # 实例化 GeneClustering 类
    df = gene_clustering.cluster(df)
    logger.info(f"Gene grouping completed. Result: {len(df)} records.")
    
    #print_memory_usage("grouping ")
    #print(df.columns)
 
    # Junction 筛选
    logger.info("Stage 4: Screening junctions...")
    df = junction_screening(df, junction_freq_ratio = args.junction_freq_ratio)
    logger.info(f"Junction screening completed. Result: {len(df)} records.")
    
    #print_memory_usage("Junction ")
    #print(df.columns)
    
    # Consensus 策略
    logger.info("Stage 5: Applying consensus strategy...")
    consensusfilter = ConsensusFilter(consensus_bp = args.consensus_bp,
                                      consensus_multiple = args.consensus_multiple,
                                      num_processes=args.threads)
    df = consensusfilter.consensus(df)
    logger.info(f"Consensus strategy applied. Result: {len(df)} records.")
    
    #print_memory_usage("Consensus ")
    #print(df.columns)
    
    # remove lowWeight edges
    logger.info("Stage 6: Prune low-confidence splice edges...")
    if len_freq2/len_freq1 < 0.55: # 测序噪音过多，进行低权重边的裁剪
        graphedgefilter = GraphEdgeFilter(threshold_lowWeight_edges = args.threshold_lowWeight_edges,
                                          num_processes=args.threads)
        df = graphedgefilter.remove_edges(df)
    logger.info(f"Prune low-confidence splice edges. Result: {len(df)} records.")
    
    #print_memory_usage("remove lowWeight edges ")
    #print(df.columns)

    # Fragmentary transcript
    logger.info("Stage 7:Fragmentary transcript filter...")
    df = filter_fragmentary_transcript(df, threshold_fragmentary_transcript_bp = args.threshold_fragmentary_transcript_bp)
    logger.info(f"Fragmentary transcript filter applied. Result: {len(df)} records.")

    # NNC/NIC
    logger.info("Stage 8: NNC/NIC filter...")
    nncnicgraphprocessor = NncNicGraphProcessor(exon_excursion_diff_bp=args.exon_excursion_diff_bp, error_sites_diff_bp=args.error_sites_diff_bp,
                                                error_sites_multiple=args.error_sites_multiple, little_exon_bp=args.little_exon_bp,
                                                little_exon_mismatch_diff_bp=args.little_exon_mismatch_diff_bp, 
                                                Nolittle_exon_mismatch_diff_bp=args.Nolittle_exon_mismatch_diff_bp,
                                                little_exon_jump_multiple=args.little_exon_jump_multiple, 
                                                num_processes = args.threads)
    df = nncnicgraphprocessor.nnc_nic_graph(df)
    logger.info(f"NNC/NIC filter applied. Result: {len(df)} records.")
    
    #print_memory_usage("NNC/NIC ")
    #print(df.columns)
    
    # Stage9: ISM filter
    logger.info("Stage 9: filter truncation...")
    truncationprocessor = TruncationProcessor(delta_Notrun_bp = args.delta_Notrun_bp, 
                                              threshold_logFC_truncation_source_freq = args.threshold_logFC_truncation_source_freq, 
                                              threshold_logFC_truncation_group_freq = args.threshold_logFC_truncation_group_freq,
                                              num_processes = args.threads)
    df = truncationprocessor.get_truncation(df)
    df = truncationprocessor.filter_truncation_logFC(df)
    logger.info(f"ISM filter applied. Result: {len(df)} records.")
    
    #print_memory_usage("ISM filter ")
    #print(df.columns)
    
    # Stage10: TS prediction
    if args.gtf_anno is None:
        logger.info("Stage 10: TS prediction...")
        terminalsitesprocessor = TerminalSitesProcessor(cluster_group_size = args.cluster_group_size,
                                                        eps = args.eps,
                                                        min_samples = args.min_samples,
                                                        num_processes = args.threads)
        df = terminalsitesprocessor.get_terminal_sites(df)
        logger.info(f"TS prediction applied. Result: {len(df)} records.")
        
        df = df.drop(columns = ['TrStart','TrEnd'])
        
        #print_memory_usage("TS prediction ")
        #print(df.columns)
        
    else:
        logger.info("Stage 10: Annotation-based filtering and transcript start prediction...")
        run_Ref2exonChain(args.gtf_anno, args.output)
        dtypes = {
            'Chr': 'category',
            'Strand': 'category',
            'TrStart': 'int32',
            'TrEnd': 'int32',
            'exonChain': str
        }
        ref_anno = pd.read_csv(
            f"{args.output}/process/anno.exonChain",
            sep='\t',
            usecols=['Chr', 'Strand', 'TrStart', 'TrEnd', 'exonChain'],
            dtype=dtypes,
            low_memory=True
        )
        ref_anno = ref_anno.dropna().drop_duplicates()

        df_raw = pd.read_parquet(os.path.join(args.output, "process/exonChain_flnc.parquet"))
        gtfrescuetss = GTFRescueTSS(little_exon_bp=args.little_exon_bp,
                                mismatch_error_sites_bp=args.mismatch_error_sites_bp, 
                                mismatch_error_sites_groupfreq_multiple=args.mismatch_error_sites_groupfreq_multiple,
                                exon_excursion_diff_bp=args.exon_excursion_diff_bp,
                                fake_exon_group_freq_multiple=args.fake_exon_group_freq_multiple,
                                fake_exon_bp=args.fake_exon_group_freq_multiple,

                                # fragmentary_transcript_FreqRatio = args.fragmentary_transcript_FreqRatio,
                                ism_freqRatio_notrun = args.ism_freqRatio_notrun,
                                nic_freqratio_mean = args.nic_freqratio_mean, nic_freqratio_group = args.nic_freqratio_group,
                                nnc_freqratio_mean = args.nnc_freqratio_mean, nnc_freqratio_group = args.nnc_freqratio_group,
                                cluster_group_size = args.cluster_group_size, eps = args.eps, min_samples = args.min_samples,
                                num_processes = args.threads)
        df = gtfrescuetss.anno_process(df_raw,df,ref_anno)
        logger.info(f"Anno filter and TS prediction applied. Result: {len(df)} records.")
        
        # df['frequency'] = df['TrStart'].apply(len)
        df = df.drop(columns = ['TrStart','TrEnd'])
        
        #print_memory_usage("TS prediction ")
        #print(df.columns)
        
    # quantification
    logger.info("Stage 11: Quantification...")
    quantifier = Quantification(assign_delta = args.assign_delta,assign_exonChain_weight = args.assign_exonChain_weight,
                        truncation_weight = args.truncation_weight,num_processes = args.threads)
    df = quantifier.assign_quantify_exonChain(df,os.path.join(args.output, "process/exonChain.count"))
    logger.info(f"Quantification applied. Result: {len(df)} records.")
    
    #print_memory_usage("Quantification ")
    #print(df.columns)
    
    df.to_csv(f'{args.output}/isatools.exonChain.tsv',index = False)
    df.to_parquet(f'{args.output}/isatools.exonChain.parquet')
    logger.info(" === ISAtools pipeline finished === ")


def main(cmd_args):
    """
    主函数，解析命令行参数并调用管道运行函数。
    """
    parser = argparse.ArgumentParser(description="ISAtools: A tool for isoform analysis and quantification.")

    parser.add_argument("--reference", "-r",type=str, required=True, help="reference genome in FASTA format (can be gzipped)")
    parser.add_argument("--bam", "-b",type=str, required=True, help="sorted and indexed BAM file")
    parser.add_argument("--output", "-o", help="output folder, will be created automatically [default=isatools_output]",
                        type=str, default="isatools_output")
    parser.add_argument("--threads", "-t", help="number of threads to use [default=16]", type=int,default=16)
    
    # 数据预处理
    parser.add_argument("--min_aln_identity", type=float, default=0.97, help="min_aln_identity")
    parser.add_argument("--min_aln_coverage", type=float, default=0.99, help="min_aln_coverage")
    
    parser.add_argument("--collapse_max_diff", type=float, default=10, help="collapse_max_diff")
    parser.add_argument("--filter_freq", type=float, default=2, help="Filtering exonChain frequency")
    parser.add_argument("--min_junction_freq", type=float, default=2, help="min_junction_freq")
    parser.add_argument("--junction_freq_ratio", type=float, default=0.25, help="Frequency ratio for junction screening")
    parser.add_argument("--consensus_bp", type=int, default=10, help="Consensus correction threshold in bp")
    parser.add_argument("--consensus_multiple", type=float, default=0.1, help="Consensus multiple for correction")
    
    # NNC/NIC
    parser.add_argument("--threshold_lowWeight_edges", type=float, default=0.05, help="Frequency ratio for removing low weight edges")
    
    parser.add_argument("--exon_excursion_diff_bp", type=int, default=20, help="exon_excursion_diff_bp")
    parser.add_argument("--error_sites_diff_bp", type=int, default=10, help="error_sites_diff_bp")
    parser.add_argument("--error_sites_multiple", type=float, default=0.01, help="error_sites_multiple")
    parser.add_argument("--little_exon_bp", type=int, default=30, help="little_exon_bp")
    parser.add_argument("--little_exon_mismatch_diff_bp", type=int, default=10, help="little_exon_mismatch_diff_bp")
    parser.add_argument("--Nolittle_exon_mismatch_diff_bp", type=int, default=20, help="Nolittle_exon_mismatch_diff_bp")
    parser.add_argument("--little_exon_jump_multiple", type=float, default=0.1, help="little_exon_jump_multiple")

    # ISM
    parser.add_argument("--threshold_logFC_truncation_source_freq", type=float, default = 0, help="threshold_logFC_truncation_source_freq")
    parser.add_argument("--threshold_logFC_truncation_group_freq", type=float, default = -100, help="threshold_logFC_truncation_group_freq")
    parser.add_argument("--delta_Notrun_bp", type=int, default = 50, help="delta_Notrun_bp") # remove
    
    # two exon
    parser.add_argument("--threshold_fragmentary_transcript_bp", type=int, default=100, help="threshold_fragmentary_transcript_bp")
    
    # anno
    parser.add_argument("--gtf_anno", '-g', type=str, default = None, help = "gene database in gffutils DB format or GTF/GFF [optional]")
    parser.add_argument("--mismatch_error_sites_bp", type=int, default=20, help="mismatch_error_sites_bp")
    parser.add_argument("--mismatch_error_sites_groupfreq_multiple", type=float, default=0.25, help="mismatch_error_sites_groupfreq_multiple")
    # parser.add_argument("--exon_excursion_diff_bp", type=int, default=20, help="exon_excursion_diff_bp")
    parser.add_argument("--fake_exon_group_freq_multiple", type=float, default=0.1, help="fake_exon_group_freq_multiple")
    parser.add_argument("--fake_exon_bp", type=int, default=50, help="fake_exon_bp")
    
    # parser.add_argument("--fragmentary_transcript_FreqRatio", type=float, default=0.00001, help="fragmentary_transcript_FreqRatio")
    parser.add_argument("--ism_freqRatio_notrun", type=float, default=1, help="ism_freqRatio_notrun")
    parser.add_argument("--nic_freqratio_mean", type=float, default=0, help="nic_freqratio_mean")
    parser.add_argument("--nic_freqratio_group", type=float, default=0, help="nic_freqratio_group")
    parser.add_argument("--nnc_freqratio_mean", type=float, default=0, help="nnc_freqratio_mean")
    parser.add_argument("--nnc_freqratio_group", type=float, default=0, help="nnc_freqratio_group")

    # TS prediction
    parser.add_argument("--cluster_group_size", type=int, default = 1500, help="cluster_group_size")
    parser.add_argument("--eps", type=int, default = 10, help="eps")
    parser.add_argument("--min_samples", type=int, default = 20, help="min_samples")
    
    # quantification
    parser.add_argument("--assign_delta", type=int, default = 5, help="assign_delta")
    parser.add_argument("--assign_exonChain_weight", type=float, default=0, help="assign_exonChain_weight")
    parser.add_argument("--truncation_weight", type=float, default=1, help="truncation_weight")

    args = parser.parse_args(cmd_args)

    try:
        run_pipeline(args)
    except KeyboardInterrupt:
        logger.error("Pipeline interrupted by user.")
    except Exception as e:
        strout = StringIO()
        print_exc(file=strout)
        error_message = strout.getvalue()
        logger.critical(f"Pipeline failed with error: {error_message}")
        sys.exit(1)


if __name__ == "__main__":
    main(sys.argv[1:])
