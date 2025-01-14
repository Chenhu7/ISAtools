import pandas as pd
import multiprocessing
import re
import subprocess
from Bio import SeqIO

def load_fasta(file):
    """
    使用 BioPython 高效加载 FASTA 文件，返回字典 {read_id: sequence}
    """
    return {record.id: str(record.seq) for record in SeqIO.parse(file, "fasta")}


def process_bam_line(line):
    """
    解析BAM文件中的单行记录，返回处理后的结果
    """
    cols = line.strip().split("\t")
    read_id = cols[0]
    chr_name = cols[2]
    strand = "+" if (int(cols[1]) & 0x10) == 0 else "-"
    cigar = cols[5]
    sequence = cols[9]

    positions = []
    current_pos = int(cols[3])
    coverage = 0
    mismatches = 0

    for token in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length, op = int(token[0]), token[1]
        if op in "M=X":
            coverage += length
            current_pos += length
        elif op == "N":
            positions.append((current_pos, current_pos + length - 1))
            current_pos += length
        elif op == "I":
            coverage += length

    positions.append((int(cols[3]), current_pos - 1))
    identity = 1 - (mismatches / coverage if coverage > 0 else 0)

    exon_chain = "-".join(f"{start}-{end}" for start, end in positions)
    return {
        "read_id": read_id,
        "Chr": chr_name,
        "Strand": strand,
        "exonChain": exon_chain,
        "identity": identity,
        "coverage": coverage / len(sequence) if len(sequence) > 0 else 0,
    }


def parse_bam_multiprocess(bam_file, num_processes=4):
    """
    并行解析BAM文件，返回解析后的DataFrame和exonChain统计信息
    """
    ec = {}
    eci = {}

    # 运行samtools view命令获取BAM内容
    process = subprocess.Popen(
        ["samtools", "view", bam_file],
        stdout=subprocess.PIPE,
        text=True,
    )

    # 多进程处理
    pool = multiprocessing.Pool(num_processes)
    results = pool.imap(process_bam_line, process.stdout)

    parsed_data = []
    for result in results:
        parsed_data.append(result)
        exon_key = f"{result['Chr']}\t{result['Strand']}\t{result['exonChain']}"
        ec[exon_key] = ec.get(exon_key, 0) + 1
        eci[exon_key] = eci.get(exon_key, 0) + result["identity"]

    pool.close()
    pool.join()

    return pd.DataFrame(parsed_data), ec, eci


def parse_and_merge(fasta_file, bam_file, num_processes=4):
    """
    解析FASTA和BAM文件，并合并解析结果，返回最终的DataFrame
    """
    fa = load_fasta(fasta_file)
    bam_parsed, ec, eci = parse_bam_multiprocess(bam_file, num_processes)

    # 处理junction统计
    junction_data = []
    for key, count in ec.items():
        chr_name, strand, exon_chain = key.split("\t")
        mean_identity = eci[key] / count
        junction_data.append([chr_name, strand, exon_chain, count, mean_identity])

    df_junction = pd.DataFrame(
        junction_data, columns=["Chr", "Strand", "exonChain", "junction", "identity"]
    )

    # 合并数据
    merged_df = bam_parsed.merge(
        df_junction, on=["Chr", "Strand", "exonChain"], how="left"
    )
    return merged_df.dropna()


def main(fasta_file, bam_file, merged_output, num_processes=10):
    # 解析和合并数据
    merged_df = parse_and_merge(fasta_file, bam_file, num_processes)
    # 保存结果到输出文件
    merged_df.to_csv(merged_output, sep="\t", index=False)
    # return merged_df


if __name__ == "__main__":
    main('../ISA-tools_debug/0data/Ref/GRCh38_ref_genome.fa', '../ISA-tools_debug/data/Mouse/Mouse.PB.simulated.sorted.bam', 't', num_processes=10)