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

# 初始化日志
logger = logging.getLogger("ISAtools")
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def run_pipeline(args):
    """
    运行主管道流程，包括数据预处理、分组、isoform 分析等。
    """
    logger.info(" === ISAtools pipeline started === ")

    # # Step 1: 数据预处理
    # logger.info("Step 1: Preprocessing input data...")
    # df = load_data(reference = args.reference, bam = args.bam, output = args.output)
    # logger.info(f"Loaded input data with {len(df)} records.")
    
    # df.to_parquet('SIRV_PB_0.parquet')
    
    df = pd.read_parquet('SIRV_PB_0.parquet')
    
    # Step 2: 基因分组
    logger.info("Step 2: Performing gene grouping...")
    gene_clustering = GeneClustering(num_processes=args.threads)  # 实例化 GeneClustering 类
    df = gene_clustering.cluster(df)
    logger.info(f"Gene grouping completed. Result: {len(df)} records.")
    
    # Step 3: Junction 筛选
    logger.info("Step 3: Screening junctions...")
    df = junction_screening(df, junction_freq_ratio = args.junction_freq_ratio)
    logger.info(f"Junction screening completed. Remaining records: {len(df)}.")
    
    # Step 4: Consensus 策略
    logger.info("Step 4: Applying consensus strategy...")
    consensusfilter = ConsensusFilter(consensus_bp = args.consensus_bp,
                                      consensus_multiple = args.consensus_multiple,
                                      num_processes=args.threads)
    df = consensusfilter.consensus(df)
    logger.info(f"Consensus strategy applied. Result: {len(df)} records.")
    
    # Step 5: remove lowWeight edges
    logger.info("Step 5: remove lowWeight edges...")
    graphedgefilter = GraphEdgeFilter(threshold_lowWeight_edges = args.threshold_lowWeight_edges,
                                      num_processes=args.threads)
    df = graphedgefilter.remove_edges(df)
    logger.info(f"Consensus strategy applied. Result: {len(df)} records.")
    
    # Step 6: NNC/NIC
    logger.info("Step 6: NNC/NIC filter...")
    nncnicgraphprocessor = NncNicGraphProcessor(exon_excursion_diff_bp = args.exon_excursion_diff_bp, 
                                                little_exon_bp = args.little_exon_bp, 
                                                little_exon_jump_multiple = args.little_exon_jump_multiple, 
                                                num_processes = args.threads)
    df = nncnicgraphprocessor.nnc_nic_graph(df)
    logger.info(f"NNC/NIC filter applied. Result: {len(df)} records.")
    
    # Step 7: filter_mutually_exclusive_exons
    logger.info("Step 7: filter_mutually_exclusive_exons...")
    mutuallyexclusiveexonfilter = MutuallyExclusiveExonFilter(mutually_exclusive_ref_ratio = args.mutually_exclusive_ref_ratio,
                                                              num_processes = args.threads)
    df = mutuallyexclusiveexonfilter.filter_mxes(df)
    logger.info(f"MXEs filter applied. Result: {len(df)} records.")
    
    # step8: two exon
    logger.info("Step 8: two exon filter...")
    df = filter_two_exon(df, threshol_twoExon_bp = args.threshol_twoExon_bp)
    logger.info(f"two exon filter applied. Result: {len(df)} records.")
    
    
    # step9: ISM filter
    logger.info("Step 9: filter truncation...")
    truncationprocessor = TruncationProcessor(delta_Notrun_bp = args.delta_Notrun_bp, 
                                              threshold_logFC_truncation_source_freq = args.threshold_logFC_truncation_source_freq, 
                                              threshold_logFC_truncation_group_freq = args.threshold_logFC_truncation_group_freq,
                                              num_processes = args.threads)
    df = truncationprocessor.get_truncation(df)
    df = truncationprocessor.filter_truncation_logFC(df)
    df = truncationprocessor.filter_truncation_delta(df)
    logger.info(f"ISM filter applied. Result: {len(df)} records.")

    
    # # Step 5: Isoform 分类
    # logger.info("Step 5: Classifying isoforms...")
    # isoformclassifier = IsoformClassifier(num_processes=args.threads)
    # df = isoformclassifier.add_category(df)  # 使用 isoform_classify 模块
    # logger.info(f"Isoform classification completed. Total isoforms: {len(df)}.")
    # print(df.value_count('category'))

    # # 保存结果
    # output_file = f"{args.output}/final_result.csv"
    # df.to_csv(output_file, index=False)
    # logger.info(f"Results saved to {output_file}.")

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
    parser.add_argument("--junction_freq_ratio", type=float, default=0.25, help="Frequency ratio for junction screening")
    parser.add_argument("--consensus_bp", type=int, default=10, help="Consensus correction threshold in bp")
    parser.add_argument("--consensus_multiple", type=float, default=0.1, help="Consensus multiple for correction")
    
    # NNC/NIC
    parser.add_argument("--threshold_lowWeight_edges", type=float, default=0.001, help="Frequency ratio for removing low weight edges")
    parser.add_argument("--exon_excursion_diff_bp", type=int, default=5, help="exon_excursion_diff_bp")
    parser.add_argument("--little_exon_bp", type=int, default=20, help="little_exon_bp")
    parser.add_argument("--little_exon_jump_multiple", type=float, default=0.01, help="little_exon_jump_multiple")
    parser.add_argument("--mutually_exclusive_ref_ratio", type=float, default=0.001, help="mutually_exclusive_ref_ratio")
    parser.add_argument("--threshol_twoExon_bp", type=int, default=50, help="threshol_twoExon_bp")
    
    # ISM
    parser.add_argument("--threshold_logFC_truncation_source_freq", type=float, default = -0.602, help="threshold_logFC_truncation_source_freq")
    parser.add_argument("--threshold_logFC_truncation_group_freq", type=float, default = -100, help="threshold_logFC_truncation_group_freq")
    parser.add_argument("--delta_Notrun_bp", type=int, default = 50, help="delta_Notrun_bp")
    
    
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