import sys
import logging
from traceback import print_exc
from io import StringIO
import argparse
import os
import shutil
import psutil
from src.common import *
from src.consensus import ConsensusFilter
from src.gene_grouping import GeneClustering
from src.isoform_classify import IsoformClassifier
from src.remove_lowConfidence_junction import SpliceConsensusFilter
from src.NNC_NIC_filter import NncNicGraphProcessor
from src.ISM_filter import TruncationProcessor
from src.get_terminal_sites import TerminalSitesProcessor
from src.gtf_rescue_filtering import GTFRescueFiltering
from src.SSC_assign import SSCAssign
from src.isoform_quantification import IsoformQuantifier
from src.generate_reports import IsoformAnnotator

def setup_logger(output_dir):
    logger = logging.getLogger("ISAtools")
    logger.setLevel(logging.INFO)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))

    log_file = os.path.join(output_dir, "isatools.log")
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))

    if logger.hasHandlers():
        logger.handlers.clear()

    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    return logger


def run_pipeline(bam,args,ref_anno = None):
    """
    Run the main ISAtools pipeline, including data preprocessing, gene grouping,
    isoform classification, and filtering steps.
    """
    sample = os.path.splitext(os.path.basename(bam))[0]
    logger = logging.getLogger("ISAtools")
    logger.info(f"\tProcessing sample: {bam}")
    
    # Stage 1: Preprocessing
    logger.info("Stage 1: Preprocessing input data...")
    df = load_data(reference = args.reference, bam = bam, output = args.output, num_threads = args.threads, min_aln_coverage = args.min_aln_coverage, min_aln_identity = args.min_aln_identity)
    len_freq1 = len(df)
    logger.info(f"Loaded input data with {len(df)} records.")

    # Stage 2: Filtering low-support SSCs
    logger.info("Stage 2: Filtering SSC with low read support...")
    df = df[df['frequency'] >= args.filter_freq]
    len_freq2 = len(df)
    logger.info(f"Low-support SSCs filtered out. Result: {len(df)} records.")
    
    # Stage 3: Gene grouping
    logger.info("Stage 3: Performing gene grouping...")
    gene_clustering = GeneClustering(num_processes=args.threads)
    df = gene_clustering.cluster(df)
    logger.info(f"Gene grouping completed. Result: {len(df)} records.")
    
    # Stage 4: Analyze conserved splice site dinucleotides (canonical and non-canonical)
    logger.info("Stage 4: Analyzing conserved splice site dinucleotides...")
    df = junction_screening(df, junction_freq_ratio=args.junction_freq_ratio)
    logger.info(f"Splice site analysis completed. Result: {len(df)} records.")

    # Stage 5: Consensus
    logger.info("Stage 5: Pruning low-confidence splice junctions based on supporting evidence...")
    consensusfilter = ConsensusFilter(consensus_bp = args.consensus_bp,
                                      consensus_multiple = args.consensus_multiple,
                                      num_processes=args.threads)
    df = consensusfilter.consensus(df)
    logger.info(f"Consensus strategy applied. Result: {len(df)} records.")
    
    # Stage7: Prune low-confidence splice junction
    logger.info("Stage 6: Prune low-confidence splice junction...")
    if len_freq2/len_freq1 < 0.55:
        splice_consensus_filter = SpliceConsensusFilter(threshold_lowWeight_edges = args.threshold_lowWeight_edges,
                                          num_processes=args.threads)
        df = splice_consensus_filter.remove_SJ(df)
    logger.info(f"Low-confidence junction pruning completed. Result: {len(df)} records.")
    
    # Stage7: Fragmentary transcript
    logger.info("Stage 7: Fragmentary transcript filter...")
    df = filter_fragmentary_transcript(df, threshold_fragmentary_transcript_bp = args.threshold_fragmentary_transcript_bp)
    logger.info(f"Fragmentary transcript filter applied. Result: {len(df)} records.")

    # Stage8: NNC/NIC filter
    logger.info("Stage 8: NNC/NIC filter...")
    nncnicgraphprocessor = NncNicGraphProcessor(exon_excursion_diff_bp=args.exon_excursion_diff_bp, error_sites_diff_bp=args.error_sites_diff_bp,
                                                error_sites_multiple=args.error_sites_multiple, little_exon_bp=args.little_exon_bp,
                                                little_exon_mismatch_diff_bp=args.little_exon_mismatch_diff_bp, 
                                                Nolittle_exon_mismatch_diff_bp=args.Nolittle_exon_mismatch_diff_bp,
                                                little_exon_jump_multiple=args.little_exon_jump_multiple, 
                                                num_processes = args.threads)
    df = nncnicgraphprocessor.nnc_nic_graph(df)
    logger.info(f"NNC/NIC filter applied. Result: {len(df)} records.")
    
    # Stage9: ISM filter
    logger.info("Stage 9: filter truncation...")
    truncationprocessor = TruncationProcessor(threshold_logFC_truncation_source_freq = args.threshold_logFC_truncation_source_freq, 
                                              threshold_logFC_truncation_group_freq = args.threshold_logFC_truncation_group_freq,
                                              num_processes = args.threads)
    df = truncationprocessor.get_truncation(df)
    df = truncationprocessor.filter_truncation_logFC(df)
    logger.info(f"ISM filter applied. Result: {len(df)} records.")

    # Stage10: TS prediction
    if ref_anno is None:
        logger.info("Stage 10: TS prediction...")
        terminalsitesprocessor = TerminalSitesProcessor(cluster_group_size = args.cluster_group_size,
                                                        eps = args.eps,
                                                        min_samples = args.min_samples,
                                                        num_processes = args.threads)
        df = terminalsitesprocessor.get_terminal_sites(df)
        logger.info(f"TS prediction applied. Result: {len(df)} records.")
        
    else:
        logger.info("Stage 10: Annotation-based filtering and transcript start prediction...")

        df_raw = pd.read_parquet(os.path.join(args.output, f"temp/{sample}.ssc_flnc.parquet"))
        gtf_rescue_filtering = GTFRescueFiltering(little_exon_bp=args.little_exon_bp,
                                mismatch_error_sites_bp=args.mismatch_error_sites_bp, 
                                mismatch_error_sites_groupfreq_multiple=args.mismatch_error_sites_groupfreq_multiple,
                                exon_excursion_diff_bp=args.exon_excursion_diff_bp,
                                fake_exon_group_freq_multiple=args.fake_exon_group_freq_multiple,
                                fake_exon_bp=args.fake_exon_group_freq_multiple,
                                ism_freqRatio_notrun = args.ism_freqRatio_notrun,
                                cluster_group_size = args.cluster_group_size, eps = args.eps, min_samples = args.min_samples,
                                num_processes = args.threads)
        df = gtf_rescue_filtering.anno_process(df_raw,df,ref_anno)
        logger.info(f"Anno filter and TS prediction applied. Result: {len(df)} records.")
    
    # Stage 11: quantification
    logger.info("Stage 11: Quantification...")
    quantifier = IsoformQuantifier(num_processes=args.threads)
    df = quantifier.assign_quantification(df, os.path.join(args.output, f"temp/{sample}.ssc_flnc.parquet"))
    logger.info(f"Quantification applied. Result: {len(df)} records.")

    df['sample'] = sample
    return df.drop(columns = 'Group')


def run_pipeline_multi(bam,args,ref_anno = None):
    result_list = []
    for sample_bam in bam:
        df_R = run_pipeline(sample_bam,args,ref_anno)
        result_list.append(df_R)

    return pd.concat(result_list)


def parse_args(cmd_args):
    parser = argparse.ArgumentParser(description="ISAtools: A pipeline for full-length RNA isoform reconstruction from long-read RNA-seq data.")

    # Required arguments
    parser.add_argument("--reference", "-r", type=str, required=True, help="Reference genome in FASTA format (can be gzipped). [Required]")
    parser.add_argument("--bam", "-b", type=str, nargs='+', required=True, help="Input sorted and indexed BAM file(s). [Required]")

    parser.add_argument("--output", "-o", type=str, default="isatools_output", help="Output directory. Default: isatools_output")
    parser.add_argument("--threads", "-t", type=int, default=16, help="Number of threads to use. Default: 16")
    parser.add_argument("--keep_temp", action="store_true", help="Keep intermediate files in the temp directory.")

    # Alignment filtering
    parser.add_argument("--min_aln_identity", type=float, default=0.97, help="Minimum alignment identity. Default: 0.97")
    parser.add_argument("--min_aln_coverage", type=float, default=0.99, help="Minimum alignment coverage. Default: 0.99")

    # SSC filtering
    parser.add_argument("--filter_freq", type=float, default=2, help="Minimum read support to retain an SSC. Default: 2")
    parser.add_argument("--min_junction_freq", type=float, default=2, help="Minimum read support for splice junctions. Default: 2")
    parser.add_argument("--junction_freq_ratio", type=float, default=0.25, help="Minimum frequency ratio for junction screening. Default: 0.25")

    # Consensus filtering
    parser.add_argument("--consensus_bp", type=int, default=10, help="Allowed deviation (bp) in consensus correction. Default: 10")
    parser.add_argument("--consensus_multiple", type=float, default=0.1, help="Supporting read ratio for consensus correction. Default: 0.1")

    # Graph edge filtering
    parser.add_argument("--threshold_lowWeight_edges", type=float, default=0.05, help="Threshold ratio for filtering weak edges. Default: 0.05")

    # NNC/NIC filtering
    parser.add_argument("--exon_excursion_diff_bp", type=int, default=20, help="Maximum exon position deviation allowed. Default: 20")
    parser.add_argument("--error_sites_diff_bp", type=int, default=10, help="Max position deviation for suspected error sites. Default: 10")
    parser.add_argument("--error_sites_multiple", type=float, default=0.01, help="Read ratio threshold for error site detection. Default: 0.01")
    parser.add_argument("--little_exon_bp", type=int, default=30, help="Maximum size for a small exon. Default: 30")
    parser.add_argument("--little_exon_mismatch_diff_bp", type=int, default=10, help="Position difference for mismatch in small exons. Default: 10")
    parser.add_argument("--Nolittle_exon_mismatch_diff_bp", type=int, default=20, help="Mismatch difference threshold for non-small exons. Default: 20")
    parser.add_argument("--little_exon_jump_multiple", type=float, default=0.1, help="Ratio threshold for exon skipping detection. Default: 0.1")

    # ISM truncation filtering
    parser.add_argument("--threshold_logFC_truncation_source_freq", type=float, default=0, help="logFC threshold of source SSCs for truncation. Default: 0")
    parser.add_argument("--threshold_logFC_truncation_group_freq", type=float, default=-100, help="logFC threshold for group truncation filtering. Default: -100")

    # Fragmentary transcript filtering
    parser.add_argument("--threshold_fragmentary_transcript_bp", type=int, default=100, help="Minimum length required to retain transcript. Default: 100")

    # Annotation-based filtering
    parser.add_argument("--gtf_anno", "-g", type=str, default=None, help="Optional GTF/GFF file for annotation-based transcript rescue and filtering.")
    parser.add_argument("--mismatch_error_sites_bp", type=int, default=20, help="Max deviation to define mismatched error sites. Default: 20")
    parser.add_argument("--mismatch_error_sites_groupfreq_multiple", type=float, default=0.25, help="Read ratio threshold for mismatch errors. Default: 0.25")
    parser.add_argument("--fake_exon_group_freq_multiple", type=float, default=0.1, help="Group ratio threshold for fake exon detection. Default: 0.1")
    parser.add_argument("--fake_exon_bp", type=int, default=50, help="Max length of a potential fake exon. Default: 50")
    parser.add_argument("--ism_freqRatio_notrun", type=float, default=0.5, help="Minimum ratio to retain non-truncated ISMs. Default: 0.5")

    # Transcription start/end prediction
    parser.add_argument("--cluster_group_size", type=int, default=1500, help="Max group size for TS clustering. Default: 1500")
    parser.add_argument("--eps", type=int, default=50, help="DBSCAN epsilon (distance threshold). Default: 50")
    parser.add_argument("--min_samples", type=int, default=20, help="Minimum samples for TS cluster. Default: 20")

    args = parser.parse_args(cmd_args)
    return args


def main(cmd_args):
    args = parse_args(cmd_args)
    os.makedirs(args.output, exist_ok=True)
    os.makedirs(os.path.join(args.output, "temp"), exist_ok=True)
    logger = setup_logger(args.output)
    logger.info(" === ISAtools pipeline started === ")

    try:
        logger.info(f"Processing BAM files...")
        output_ssc = os.path.join(args.output, "temp")
        run_bam2ssc(args.reference, args.bam, output_ssc,args.threads)
        logger.info(f"BAM files were converted to SSC format.")

        if args.gtf_anno:
            run_Ref2SSC(args.gtf_anno, args.output, args.threads)
            ref_anno = pd.read_csv(f"{args.output}/temp/anno.ssc",sep='\t')
            ref_anno = ref_anno[ref_anno['SSC'].notna()]
        else:
            ref_anno = None

        if len(args.bam) == 1:
            if ref_anno is None:
                df_result = run_pipeline(args.bam[0],args)
            else:
                df_result =run_pipeline(args.bam[0],args,ref_anno)

        else:
            if ref_anno is None:
                df_result = run_pipeline_multi(args.bam,args)
                  
            else:
                df_result = run_pipeline_multi(args.bam,args,ref_anno)

        annotator = IsoformAnnotator(num_processes=args.threads)
        annotator.save_results(df_result,args.output,ref_anno)

        if not args.keep_temp:
            temp_dir = os.path.join(args.output, "temp")
            shutil.rmtree(temp_dir, ignore_errors=True)

        logger.info(" === ISAtools pipeline finished === ")

    except Exception as e:
        logger.error("An error occurred:")
        logger.error(str(e))
        print_exc()
        sys.exit(1)


if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except KeyboardInterrupt:
        print("\nInterrupted by user. Exiting.")
        sys.exit(130)
