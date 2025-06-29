import sys
import gzip
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
import argparse

def extract_ids(description):
    gene_match = re.search(r'gene_id "([^"]+)"', description)
    tx_match = re.search(r'transcript_id "([^"]+)"', description)
    gene_name_match = re.search(r'gene_name "([^"]+)"', description)
    gene_id = gene_match.group(1) if gene_match else "NA"
    transcript_id = tx_match.group(1) if tx_match else "NA"
    gene_name = gene_name_match.group(1) if gene_name_match else "NA"
    return gene_id, transcript_id, gene_name

def parse_exon_line(line):
    cols = line.strip().split('\t')
    if len(cols) < 9 or cols[2] != 'exon':
        return None
    gene_id, transcript_id, gene_name = extract_ids(cols[8])
    if not transcript_id or not gene_id:
        return None
    return {
        'TrID': transcript_id,
        'GeneID': gene_id,
        'GeneName': gene_name,
        'Chr': cols[0],
        'Strand': cols[6],
        'TrStart': int(cols[3]),
        'TrEnd': int(cols[4])
    }

def process_chunk(lines):
    transcripts = dict()
    for line in lines:
        exon = parse_exon_line(line)
        if exon is None:
            continue
        tr_id = exon['TrID']
        if tr_id not in transcripts:
            transcripts[tr_id] = {
                'Chr': exon['Chr'],
                'Strand': exon['Strand'],
                'GeneID': exon['GeneID'],
                'GeneName': exon['GeneName'],
                'starts': [],
                'ends': []
            }
        transcripts[tr_id]['starts'].append(exon['TrStart'])
        transcripts[tr_id]['ends'].append(exon['TrEnd'])
    return transcripts

def merge_dicts(dicts):
    merged = dict()
    for d in dicts:
        for tr_id, info in d.items():
            if tr_id not in merged:
                merged[tr_id] = info
            else:
                merged[tr_id]['starts'].extend(info['starts'])
                merged[tr_id]['ends'].extend(info['ends'])
    return merged

def read_gtf_chunks(input_file, chunk_size=10000):
    open_func = gzip.open if input_file.endswith('.gz') else open
    chunk = []
    with open_func(input_file, 'rt') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            chunk.append(line)
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
        if chunk:
            yield chunk

def main(input_gtf, workers, chunk_size, output_file=None):
    all_results = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = []
        for chunk in read_gtf_chunks(input_gtf, chunk_size=chunk_size):
            futures.append(executor.submit(process_chunk, chunk))
        for future in as_completed(futures):
            all_results.append(future.result())
    merged = merge_dicts(all_results)

    if output_file is None:
        output_file = input_gtf.rsplit('.', 1)[0] + ".SSC"
    with open(output_file, 'w') as out:
        out.write("TrID\tGeneID\tGeneName\tChr\tStrand\tTrStart\tTrEnd\tSSC\n")
        for tr_id, info in merged.items():
            sites = sorted(info['starts'] + info['ends'])
            start = sites[0]
            end = sites[-1]
            inner = '-'.join(str(x) for x in sites[1:-1]) if len(sites) > 2 else 'NA'
            out.write(f"{tr_id}\t{info['GeneID']}\t{info.get('GeneName','NA')}\t{info['Chr']}\t{info['Strand']}\t{start}\t{end}\t{inner}\n")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Convert GTF annotation to SSC format with multiprocessing, including gene_name")
    parser.add_argument("-i", "--input_gtf", help="Input GTF or GTF.gz file")
    parser.add_argument("-w", "--workers", type=int, default=4, help="Number of workers (default: 4)")
    parser.add_argument("-c", "--chunk_size", type=int, default=10000, help="Chunk size for processing lines (default: 10000)")
    parser.add_argument("-o", "--output", default=None, help="Output SSC file path (default: input filename + .SSC)")

    args = parser.parse_args()

    main(args.input_gtf, args.workers, args.chunk_size, args.output)
