
# ISAtools

**ISAtools** (**I**soform **S**equencing **A**nalysis tools) is a **sequencing data–driven framework** for full-length RNA isoform reconstruction and quantification from long-read RNA sequencing data. Designed with annotation flexibility and biological fidelity in mind, ISAtools supports isoform identification with high precision and recall, and accurately resolves splice junctions and transcript boundaries directly from read evidence.

If reference annotations are available, ISAtools incorporates conserved, low-abundance isoforms through guided filtering and rescue steps, further enhancing transcriptome completeness.

---

## 🔧 Quick Start

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

## 📦 Installation

ISAtools requires **Python 3.9 or higher** and `samtools` in your `$PATH`.

We recommend using [Conda](https://docs.conda.io/) to manage dependencies and environments.

### ✅ Create a Conda environment

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

### ✅ Use environment.yml (alternative setup)
```bash
# Clone ISAtools repository
git clone https://github.com/Chenhu7/ISAtools.git
cd ISAtools

conda env create -f environment.yml
conda activate isatools
```

### 🧪 Verify installation by running the toy example:

```bash
conda activate isatools

python isatools.py -r test/toy_data/test_genome.fasta -b test/toy_data/test_aligned.bam
```

Upon successful execution, output files will appear in the default `isatools_output` folder.

---

## 📁 Output Files

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

## 📄 File Format Descriptions

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
