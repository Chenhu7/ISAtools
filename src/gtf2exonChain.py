import pandas as pd
import re
import os
import sys

def extract_ids(description):
    gene_id_match = re.search(r'gene_id "([^"]+)"', description)
    transcript_id_match = re.search(r'transcript_id "([^"]+)"', description)
    return gene_id_match.group(1) if gene_id_match else None, transcript_id_match.group(1) if transcript_id_match else None

def aggregate_sites(group):
    # 仅处理 TrStart 和 TrEnd 列
    sites = sorted(group['TrStart'].tolist() + group['TrEnd'].tolist())
    return pd.Series({"sites": sites})

def main(input_gtf):
    output = os.path.splitext(input_gtf)[0] + ".exonChain"

    # 读取 GTF 文件
    gtf = pd.read_csv(input_gtf, sep='\t', header=None)
    gtf = gtf.drop(columns=[1, 5, 7])
    gtf.columns = ['Chr', 'type', 'TrStart', 'TrEnd', 'Strand', 'description']

    # 提取 gene_id 和 transcript_id
    gtf[['gene_id', 'transcript_id']] = gtf['description'].apply(lambda x: pd.Series(extract_ids(x)))

    # 仅保留 exon 类型的数据
    gtf = gtf[gtf['type'] == 'exon']

    # Groupby 操作，显式排除分组列
    grouped_gtf = (
        gtf.groupby(["Chr", "Strand", "gene_id", "transcript_id"], group_keys=False)[['TrStart', 'TrEnd']]
        .apply(aggregate_sites)
        .reset_index()
    )

    # 处理位点信息
    grouped_gtf['TrStart'] = grouped_gtf['sites'].apply(lambda x: x[0])
    grouped_gtf['TrEnd'] = grouped_gtf['sites'].apply(lambda x: x[-1])
    grouped_gtf['exonChain'] = grouped_gtf['sites'].apply(lambda x: '-'.join(map(str, x[1:-1])))
    grouped_gtf = grouped_gtf.drop(columns='sites')

    # 保存为新的 exonChain 文件
    grouped_gtf.to_csv(output, index=False)
    print(f"exonChain file saved to: {output}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python gtf2exonChain.py <input_gtf>")
        sys.exit(1)
    input_gtf = sys.argv[1]
    main(input_gtf)

