import pandas as pd
import argparse

def convert_to_gtf(input_file, output_file):
    # 读取 parquet 文件
    exonChain = pd.read_parquet(input_file)

    # 解析 exonChain 并合并 pred_terminal_sites
    exonChain['sites'] = exonChain.apply(lambda row: sorted(
        list(map(int, row['exonChain'].split('-'))) + list(row['pred_terminal_sites'])
    ), axis=1)

    # 只保留所需列
    exonChain = exonChain[['Chr', 'Strand', 'sites']].copy()

    # 生成 GTF 文件内容
    gtf_lines = []
    for idx, row in exonChain.iterrows():
        chr_name = row["Chr"]
        strand = row["Strand"]
        sites = row["sites"]
        
        # 生成 transcript 记录
        transcript_id = f"ISAtoolsTx{idx+1}"
        gene_id = f"ISAtoolsGene{idx+1}"
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

    print(f"GTF 文件已保存至: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert exonChain parquet to GTF format")
    parser.add_argument("input_file", type=str, help="Path to the input parquet file")
    # parser.add_argument("output_file", type=str, help="Path to the output GTF file")

    args = parser.parse_args()
    
    convert_to_gtf(args.input_file,output_file="ISAtools_annotations.gtf")

