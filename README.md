# ONTAP

**ONT** **A**mplicon sequencing **P**hylotyper

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15181157.svg)](https://doi.org/10.5281/zenodo.15181157)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[[_TOC_]]

## Pipeline overview

ONTAP is a Nextflow DSL2 pipeline for basecalling, read mapping, QC, variant calling, and phylogenetic analysis of nanopore multiplex amplicon data. It is designed for STI and epidemic pathogen genomics workflows using Oxford Nanopore long reads.

The pipeline carries out the following steps:

1. **Basecalling** — converts raw FAST5/POD5 signal data to reads using Dorado, followed by demultiplexing by barcode
2. **Pre-mapping QC** — FastQC on raw and adapter-trimmed reads; PycoQC on sequencing summary
3. **Adapter and primer trimming** — removes adapters and primers with Cutadapt, filtering reads by length
4. **Mapping** — aligns trimmed reads to the reference using Minimap2; sorts and indexes with SAMtools
5. **Post-mapping filtering** — separates on-target from off-target reads using the target regions BED file
6. **Post-mapping QC** — coverage reporting over defined amplicon regions (BEDtools + SAMtools depth); read-length distribution
7. **Variant calling** — per-sample SNP calling with Clair3 (haploid mode, gVCF output); gVCF merging with BCFtools
8. **Consensus curation** — constructs per-sample consensus FASTA sequences from gVCF output
9. **Phylogenetics** — builds a maximum-likelihood phylogenetic tree with RAxML-NG (optional recombination removal)
10. **Reporting** — MultiQC report aggregating QC metrics across all samples

![flowchart](images/ampliseq_pipeline.jpg)

## Usage

### Quickstart

#### Installation

Before running ONTAP, install the following dependencies:

1. [Install Nextflow](https://www.nextflow.io/docs/latest/install.html) (version >= 24.04.2)

2. [Install Docker](https://docs.docker.com/engine/install/) (required for the `docker` and `laptop` profiles)

3. Download the appropriate Dorado installer from the [Dorado repository](https://github.com/nanoporetech/dorado#installation). The path to the executable will be `<path to downloaded folder>/bin/dorado`.

4. (Optional) Download the appropriate Dorado basecalling model from the [Dorado repository](https://github.com/nanoporetech/dorado/#available-basecalling-models):

   ```bash
   # Download all models
   dorado download --model all
   # Download a particular model
   dorado download --model <model>
   ```

   If a pre-downloaded model path is not provided, the model specified by `--basecall_model` will be downloaded automatically during the run.

5. Download the appropriate Clair3 model from the [Rerio repository](https://github.com/nanoporetech/rerio?tab=readme-ov-file#clair3-models) (requires Python 3):

   ```bash
   # Clone the Rerio repo
   git clone https://github.com/nanoporetech/rerio

   # Download all Clair3 models
   python3 download_model.py --clair3
   # Or download a particular model
   python3 download_model.py --clair3 clair3_models/<config>_model
   ```

   The downloaded model will be found at `clair3_models/<config>` within the Rerio directory. The recommended model for the default basecalling configuration is `r1041_e82_400bps_hac_v430`.

6. Clone the ONTAP repository with its required submodules:
   ```bash
   git clone --recurse-submodules <repo-url>
   ```

#### From source code

After completing the installation steps above, run the pipeline from the cloned repository. Replace the placeholder values with absolute paths appropriate for your system.

**With Docker** (`-profile docker`):

```bash
nextflow run ONTAP/main.nf \
    --raw_read_dir <directory containing FAST5/POD5 files> \
    --reference <reference FASTA> \
    --primers <FASTA containing primers> \
    --target_regions_bed <BED file of target regions> \
    --additional_metadata <CSV mapping sample IDs to barcodes> \
    --dorado_local_path <absolute path to Dorado executable> \
    --clair3_model <path to Clair3 model> \
    -profile docker
```

The [examples](examples) folder contains example input files.

**With the laptop profile** (`-profile laptop`):

The `laptop` profile enables Docker and supports offline operation. It sets default local paths for Dorado, the basecalling model, and the Clair3 model (under `/Users/Shared/ampliseq/`). You can override any of these defaults on the command line or via a custom config file (`-c my_custom.config`).

```bash
nextflow run ONTAP/main.nf \
    --raw_read_dir <directory containing FAST5/POD5 files> \
    --reference <reference FASTA> \
    --primers <FASTA containing primers> \
    --target_regions_bed <BED file of target regions> \
    --additional_metadata <CSV mapping sample IDs to barcodes> \
    -profile laptop
```

#### Demo

To run a short demo using a Zenodo dataset, follow the instructions in [demo/Demo.md](demo/Demo.md).

#### Using on the Sanger farm

Load the required modules:

```bash
module load nextflow ISG/singularity
```

Follow installation steps 5 and 6 above to download a Clair3 model and clone the repository.

Submit the Nextflow master process as an LSF job in the oversubscribed queue:

```bash
bsub -o output.o -e error.e -q oversubscribed -R "select[mem>4000] rusage[mem=4000]" -M4000 \
    nextflow run ONTAP/main.nf \
        --raw_read_dir <directory containing FAST5/POD5 files> \
        --reference <reference FASTA> \
        --primers <FASTA containing primers> \
        --target_regions_bed <BED file of target regions> \
        --additional_metadata <CSV mapping sample IDs to barcodes> \
        --clair3_model <path to Clair3 model> \
        -profile standard
```

Once the run has finished successfully and you have inspected the output, clean up intermediate files. The `work/` directory and `.nextflow.log` are useful for troubleshooting — do not delete them until you are satisfied the outputs are correct:

```bash
rm -rf work .nextflow*
```

Alternatively, use `nextflow clean` for more fine-grained control over which runs and intermediate files are removed.

### Input

The following inputs are required for every run:

| Parameter               | Description                                                                                                                 |
| ----------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| `--raw_read_dir`        | Directory containing raw FAST5 or POD5 files from the sequencer.                                                            |
| `--reference`           | Reference genome in FASTA format to align reads against.                                                                    |
| `--primers`             | FASTA file containing primer sequences used to generate the amplicons. Used by Cutadapt for primer trimming.                |
| `--target_regions_bed`  | BED file defining the amplicon target regions. Used for on-target filtering and coverage reporting.                         |
| `--additional_metadata` | CSV file mapping sample IDs to barcodes. Must contain at minimum `barcode_kit` and `barcode` columns.                       |
| `--clair3_model`        | Absolute path to a locally downloaded Clair3 model directory (see [Installation](#installation) step 5).                    |
| `--dorado_local_path`   | Absolute path to the Dorado executable. Required when using `docker` or `laptop` profiles. Not required on the Sanger farm. |

Optional inputs:

| Parameter               | Description                                                                                                |
| ----------------------- | ---------------------------------------------------------------------------------------------------------- |
| `--basecall_model_path` | Path to a pre-downloaded Dorado basecalling model. If not provided, the model is downloaded automatically. |
| `--multiqc_config`      | Path to a custom MultiQC configuration file.                                                               |

### Output

The pipeline writes all results to `--outdir` (default: `results`):

```
results/
  sequencing_summary/
    summary.tsv                             # Dorado sequencing summary
  fastqs/                                   # Basecalled FASTQ files per sample (if --save_fastqs)
    <sample_ID>.fastq.gz
  cutadapt/                                 # Adapter/primer-trimmed reads (if --save_trimmed / --save_too_short / --save_too_long)
    <sample_ID>_trimmed.fastq.gz
    <sample_ID>_too_short.fastq.gz
    <sample_ID>_too_long.fastq.gz
  sorted_ref/                               # Reference FASTA index
    <reference>.fai
  mapped_reads/                             # Sorted BAM per sample (if --keep_sorted_bam); BAI index (if --keep_bam_files)
    <sample_ID>_sorted.bam
    <sample_ID>.bai
  qc/
    fastqc_pre_trim/
      <sample_ID>_fastqc.html               # FastQC report on raw reads
      <sample_ID>_fastqc.zip
    fastqc_post_trim/
      <sample_ID>_trimmed_fastqc.html       # FastQC report on trimmed reads
      <sample_ID>_trimmed_fastqc.zip
    pycoqc/
      *                                     # PycoQC basecalling QC report
    <qc_stage>/
      readlengths/
        <sample_ID>.read-lengths.tsv        # Read-length distribution
      samtools_stats/
        <sample_ID>.stats                   # SAMtools stats
        <sample_ID>.flagstats               # SAMtools flagstats
      coverage/
        samtools_depth/
          <sample_ID>_samtools_depth.tsv    # Per-position depth (SAMtools)
        coverage_summary/
          *coverage_summary.tsv             # Per-amplicon coverage summary
          *.html                            # Coverage plot
        bedtools_genome_coverage/
          *.bedGraph                        # Genome-wide coverage (BEDtools)
        bedtools_coverage/
          *coverage.bed                     # Per-region coverage (BEDtools)
  variants/
    <sample_ID>_clair3.gvcf.gz              # Per-sample Clair3 gVCF
    <sample_ID>_clair3.vcf.gz              # Per-sample Clair3 VCF
    logs/
      <sample_ID>_clair3.log               # Clair3 run log
    merged_gvcf/
      <run>_<date>_merged.vcf.gz           # Merged gVCF across all samples
      <run>_<date>_merged.tsv              # Variant table from merged gVCF
  curated_consensus/
    <sample_ID>.fasta                       # Per-sample consensus FASTA (multi-locus)
    <sample_ID>_multi_locus.fasta           # Full multi-locus consensus
    <sample_ID>_wg.fasta                    # Whole-genome consensus
  snp_aln/
    merged.fasta.snp.aln                    # SNP-only alignment (snp-sites)
  tree/
    *.support                               # RAxML-NG phylogenetic tree with bootstrap support
  gubbins/                                  # Gubbins recombination removal outputs (if --remove_recombination)
    gubbins_out.*
  multiqc/
    multiqc_report.html                     # Aggregated MultiQC report
    multiqc_data/                           # MultiQC data directory
```

### Parameters

**Reference files (mandatory)**

| Option                  | Default | Description                                            |
| ----------------------- | ------- | ------------------------------------------------------ |
| `--raw_read_dir`        | `""`    | Directory containing raw FAST5/POD5 files.             |
| `--reference`           | `""`    | Path to the reference genome FASTA.                    |
| `--primers`             | `""`    | Path to the primer sequences FASTA.                    |
| `--target_regions_bed`  | `""`    | Path to the BED file defining target amplicon regions. |
| `--additional_metadata` | `""`    | Path to CSV mapping sample IDs to barcodes.            |

---

**Basecalling**

| Option                  | Default                              | Description                                                                             |
| ----------------------- | ------------------------------------ | --------------------------------------------------------------------------------------- |
| `--basecall`            | `true`                               | Enable basecalling.                                                                     |
| `--basecall_model`      | `dna_r10.4.1_e8.2_400bps_hac@v4.3.0` | Dorado basecalling model. Must match the flow cell and chemistry used.                  |
| `--basecall_model_path` | `""`                                 | Path to a pre-downloaded Dorado model. If empty, the model is downloaded automatically. |
| `--dorado_local_path`   | `""`                                 | Absolute path to a locally installed Dorado executable.                                 |
| `--trim_adapters`       | `all`                                | Adapter/primer trimming mode passed to Dorado.                                          |
| `--min_qscore`          | `9`                                  | Minimum Phred quality score for read filtering during basecalling.                      |
| `--read_format`         | `fastq`                              | Output format for basecalled reads.                                                     |

---

**QC**

| Option                            | Default                              | Description                                                           |
| --------------------------------- | ------------------------------------ | --------------------------------------------------------------------- |
| `--cutadapt_args`                 | `"-e 0.15 --no-indels --overlap 18"` | Additional arguments passed to Cutadapt for primer trimming.          |
| `--lower_read_length_cutoff`      | `450`                                | Minimum read length (bp) after primer trimming.                       |
| `--upper_read_length_cutoff`      | `800`                                | Maximum read length (bp) after primer trimming.                       |
| `--coverage_reporting_thresholds` | `"1,2,8,10,25,30,40,50,100"`         | Comma-separated depth thresholds for per-amplicon coverage reporting. |
| `--coverage_filtering_threshold`  | `"25"`                               | Minimum mean coverage depth for a sample to pass filtering.           |

---

**Variant calling**

| Option                  | Default | Description                                                                                 |
| ----------------------- | ------- | ------------------------------------------------------------------------------------------- |
| `--clair3_model`        | `""`    | Path to the locally downloaded Clair3 model directory.                                      |
| `--clair3_min_coverage` | `"8"`   | Minimum read depth required to call a variant with Clair3.                                  |
| `--masking_quality`     | `15`    | Phred quality score threshold for base masking. Bases below this score are replaced with N. |

---

**Tree building**

| Option                   | Default  | Description                                                        |
| ------------------------ | -------- | ------------------------------------------------------------------ |
| `--remove_recombination` | `false`  | Remove recombination events before building the phylogenetic tree. |
| `--raxml_base_model`     | `GTR+G4` | Substitution model used by RAxML-NG.                               |
| `--raxml_threads`        | `2`      | Number of threads allocated to RAxML-NG.                           |

### Dependencies

The following dependencies must be installed separately before running the pipeline:

- **Dorado** — the Oxford Nanopore basecaller. Download from the [Dorado GitHub releases](https://github.com/nanoporetech/dorado#installation) and supply the path via `--dorado_local_path`. On the Sanger farm, the containerised version is used automatically.

- **Clair3 model** — the statistical model used for variant calling. Download via the [Rerio repository](https://github.com/nanoporetech/rerio) and supply the path via `--clair3_model`. The recommended model for the default basecalling configuration is `r1041_e82_400bps_hac_v430`.

All other pipeline dependencies are containerised and pulled automatically.

## Software versions

| Tool        | Version | Container                                             |
| ----------- | ------- | ----------------------------------------------------- |
| bcftools    | 1.20    | `quay.io/biocontainers/bcftools:1.20--h8b25389_0`     |
| bedtools    | 2.31.1  | `quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1`   |
| clair3      | v1.0.9  | `hkubal/clair3:v1.0.9`                                |
| cutadapt    | 4.7     | `quay.io/biocontainers/cutadapt:4.7--py310h4b81fae_1` |
| cuda_dorado | 0.7.1   | `quay.io/sangerpathogens/cuda_dorado:0.7.1`           |
| fastqc      | 0.12.1  | `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0`     |
| minimap2    | 2.26    | `quay.io/biocontainers/minimap2:2.26--he4a0461_2`     |
| multiqc     | 1.22.2  | `quay.io/biocontainers/multiqc:1.22.2--pyhdfd78af_0`  |
| pod5        | 0.3.6   | `quay.io/sangerpathogens/pod5:0.3.6`                  |
| pycoqc      | 2.5.2   | `quay.io/biocontainers/pycoqc:2.5.2--py_0`            |
| samtools    | 1.19.2  | `quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1`   |
| seqtk       | 1.4     | `quay.io/biocontainers/seqtk:1.4--he4a0461_2`         |

## Troubleshooting

- **Runtime performance**: a full analysis run including basecalling typically takes 10–12 hours in laptop/Docker mode, and 40 minutes to 1.5 hours on the Sanger HPC with GPU access. Choosing a faster Dorado model (e.g. `fast`) will reduce basecalling time.
- **GPU support**: the pipeline runs without a GPU, but Dorado basecalling is substantially faster with GPU hardware. See the [Dorado documentation](https://github.com/nanoporetech/dorado?tab=readme-ov-file#platforms) for supported GPU platforms.
- **Offline operation**: use `-profile laptop` for offline runs. Ensure all models are pre-downloaded and their paths are supplied.
- **Resuming a failed run**: add `-resume` to your command to restart from cached intermediate results.
- For further help, check `.nextflow.log` and the per-process `.command.log` logs in the `work/` directory.

Sanger users may find [this page](https://ssg-confluence.internal.sanger.ac.uk/spaces/PaMI/pages/181078206/General+pipeline+info#Generalpipelineinfo-Troubleshootingafailedpipelinerunandsendingabugreport) useful for troubleshooting Nextflow pipeline runs.

## Issues and Contributions

**GitHub users:** if you find an issue with this pipeline, or would like to suggest an improvement, please log an issue or open a pull request on this repository.

**Sanger users:** if you need internal support, you can raise an issue on the PAM Freshservice portal: https://sanger.freshservice.com/support/catalog/items/426
