
# ISAtools

**ISAtools** (**I**soform **S**equencing **A**nalysis tools) is a **sequencing data–driven framework** for full-length RNA isoform reconstruction and quantification from **PacBio circular consensus sequencing data**. Designed with annotation flexibility and biological fidelity in mind, ISAtools supports isoform identification with high precision and recall, and accurately resolves splice junctions and transcript boundaries directly from read evidence.

If reference annotations are available, ISAtools incorporates conserved, low-abundance isoforms through guided filtering and rescue steps, further enhancing transcriptome completeness.

---

## Quick Start

To run ISAtools with aligned reads, the minimum required parameters are:

```bash
python isatools.py -r /PATH/TO/reference_genome.fasta \
	-b /PATH/TO/sample1.bam /PATH/TO/sample2.bam /PATH/TO/sample3.bam
```

For example, using the toy data provided in this repository:

```bash
python isatools.py -r test/toy_data/test_genome.fasta \
	-b test/toy_data/test_aligned.bam
```

To run ISAtools with a reference gene annotation:

```bash
python isatools.py -r /PATH/TO/reference_genome.fasta \
	-b /PATH/TO/sample1.bam /PATH/TO/sample2.bam /PATH/TO/sample3.bam \
	-g /PATH/TO/gene_annotation.gtf \
	-o OUTPUT_FOLDER
```

Example with toy data:

```bash
python isatools.py -r test/toy_data/test_genome.fasta \
	-b test/toy_data/test_aligned.bam \
	-g test/toy_data/test_annotation.gtf \
	-o OUTPUT_FOLDER
```

> **Note:** It is strongly recommended to pre-trim adapter and polyA sequences. If using untrimmed reads, add `--min_aln_coverage 0` to avoid misfiltering due to alignment artifacts.

---

## Installation

ISAtools requires **Python 3.9 or higher** and `samtools` in your `$PATH`.

We recommend using [Conda](https://docs.conda.io/) to manage dependencies and environments.

### Create a Conda environment

```bash
# Create and activate the Conda environment
conda create -n isatools python=3.9 -y
conda activate isatools

# Install samtools
conda install -c bioconda samtools -y

# Clone ISAtools repository and install Python dependencies
git clone https://github.com/Chenhu7/ISAtools.git
cd ISAtools
pip install -r requirements.txt
```

### Use environment.yml (alternative setup)
```bash
# Clone ISAtools repository
git clone https://github.com/Chenhu7/ISAtools.git
cd ISAtools

conda env create -f environment.yml
conda activate isatools
```

### Verify installation by running the toy example:

```bash
python isatools.py -r test/toy_data/test_genome.fasta -b test/toy_data/test_aligned.bam
```

Upon successful execution, output files will appear in the default `isatools_output` folder.

---

## Parameters

| Category                                | Parameter                                   | Description                                                 | Default           |
| --------------------------------------- | ------------------------------------------- | ----------------------------------------------------------- | ----------------- |
| **Basic**                               | `-r, --reference`                           | Reference genome in FASTA format (supports gzipped files).  | Required          |
|                                         | `-b, --bam`                                 | Input sorted and indexed BAM file(s).                       | Required          |
|                                         | `-o, --output`                              | Output directory.                                           | `isatools_output` |
|                                         | `-t, --threads`                             | Number of threads to use.                                   | `4`               |
|                                         | `--keep_temp`                               | Retain intermediate files.                                  | Off               |
| **Alignment Filtering**                 | `--min_aln_identity`                        | Minimum alignment identity to retain reads.                 | `0.97`            |
|                                         | `--min_aln_coverage`                        | Minimum alignment coverage.                                 | `0.99`            |
| **SSC Filtering**                       | `--filter_freq`                             | Minimum read support to retain an SSC.                      | `2`               |
|                                         | `--min_junction_freq`                       | Minimum read support for splice junctions.                  | `2`               |
|                                         | `--junction_freq_ratio`                     | Minimum frequency ratio for junction filtering.             | `0.25`            |
| **Consensus Refinement**                | `--consensus_bp`                            | Allowed positional deviation (bp) for consensus correction. | `10`              |
|                                         | `--consensus_multiple`                      | Minimum read support ratio for consensus correction.        | `0.1`             |
| **Splice Structure Filtering**          | `--threshold_lowWeight_edges`               | Threshold ratio for filtering weak graph edges.             | `0.05`            |
| **NNC/NIC Identification & Correction** | `--exon_excursion_diff_bp`                  | Maximum allowed exon position deviation.                    | `20`              |
|                                         | `--error_sites_diff_bp`                     | Max deviation for suspected error sites.                    | `10`              |
|                                         | `--error_sites_multiple`                    | Read ratio threshold for error site detection.              | `0.01`            |
|                                         | `--little_exon_bp`                          | Maximum size defining a small exon.                         | `30`              |
|                                         | `--little_exon_mismatch_diff_bp`            | Allowed mismatch deviation for small exons.                 | `10`              |
|                                         | `--Nolittle_exon_mismatch_diff_bp`          | Allowed mismatch deviation for regular exons.               | `20`              |
|                                         | `--little_exon_jump_multiple`               | Ratio threshold for exon skipping detection.                | `0.1`             |
| **ISM Identification & Filtering**      | `--threshold_logFC_truncation_source_freq`  | logFC threshold for source SSC filtering.                   | `0`               |
|                                         | `--threshold_logFC_truncation_group_freq`   | logFC threshold for group-level truncation filtering.       | `-100`            |
|                                         | `--threshold_fragmentary_transcript_bp`     | Minimum transcript length to retain.                        | `100`             |
| **Optional Annotation Refinement**      | `-g, --gtf_anno`                            | Reference annotation (GTF/GFF) for rescue and filtering.    | None              |
|                                         | `--mismatch_error_sites_bp`                 | Max deviation for mismatch error detection.                 | `20`              |
|                                         | `--mismatch_error_sites_groupfreq_multiple` | Ratio threshold for mismatch grouping.                      | `0.25`            |
|                                         | `--fake_exon_group_freq_multiple`           | Ratio threshold for fake exon detection.                    | `0.1`             |
|                                         | `--fake_exon_bp`                            | Max length of a candidate fake exon.                        | `50`              |
|                                         | `--ism_freqRatio_notrun`                    | Minimum ratio to retain non-truncated ISMs.                 | `0.5`             |
| **TSS/TES Detection**                   | `--cluster_group_size`                      | Maximum group size for clustering.                          | `1500`            |
|                                         | `--eps`                                     | DBSCAN epsilon (distance threshold, bp).                    | `50`              |
|                                         | `--min_samples`                             | Minimum reads required to form a cluster.                   | `20`              |

---

## Output Files

### Main Output

- `isatools.filtered.transcript.gtf`: Final filtered transcript models (including known and novel isoforms).
- `isatools_gene_counts.tsv`: Gene-level raw read counts.
- `isatools_gene_tpm.tsv`: Gene-level TPM expression values.
- `isatools_transcript_counts.tsv`: Transcript-level read counts.
- `isatools_transcript_tpm.tsv`: Transcript-level TPMs.

### Optional Output (with `--keep_temp`)

- `*_flnc.ssc`: Read-level SSC file with alignment details.
- `*_ssc.count`: Unique SSCs with aggregated read counts.
- `anno.ssc`: Reference annotation converted into SSC format (if `-g` is provided).

---

## File Format Descriptions

### Expression Tables (`*_counts.tsv`, `*_tpm.tsv`)

| Column        | Description                                                  |
| ------------- | ------------------------------------------------------------ |
| `TrID`        | Transcript ID, uniquely defined by the combination of `TrStart`, `SSC`, and `TrEnd`. Isoforms sharing the same SSC but differing in transcript boundaries (either `TrStart` or `TrEnd`) are treated as distinct transcripts. If multiple TSS/TES pairs are detected for the same SSC, a suffix `_TssTes*` is appended to distinguish each variant. |
| `GeneID`      | Gene ID; if reference annotation is unavailable, ISAtools assigns a gene cluster ID. |
| `GeneName`    | Same as above; derived from annotation if available.         |
| `count`/`TPM` | Expression value for transcript or gene.                     |

### SSC File – Read Level (`*_flnc.ssc`)

| Column | Description        |
| ------ | ------------------ |
| 0      | Read ID            |
| 1      | Chromosome ID      |
| 2      | Strand             |
| 3      | 5′ alignment start |
| 4      | 3′ alignment end   |
| 5      | Splice site chain  |
| 6      | Alignment identity |
| 7      | Alignment coverage |

### SSC File – Unique SSC Level (`*_ssc.count`)

| Column | Description                   |
| ------ | ----------------------------- |
| 0      | Frequency (read count of SSC) |
| 1      | Chromosome                    |
| 2      | Strand                        |
| 3      | Splice site chain             |
| 4      | Donor-acceptor motif sequence |

### Annotation SSC File (`anno.ssc`)

| Column    | Description                              |
|-----------|------------------------------------------|
| `TrID`    | Transcript ID                             |
| `GeneID`  | Gene ID                                   |
| `GeneName`| Gene name                                 |
| `Chr`     | Chromosome                                |
| `Strand`  | Strand                                     |
| `TrStart` | Transcript start site                     |
| `TrEnd`   | Transcript end site                       |
| `SSC`     | Splice Site Chain                         |
