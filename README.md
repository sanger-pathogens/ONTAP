# Long-read Ampliseq

A Nextflow pipeline for basecalling, read mapping, QC, variant calling and analysis of nanopore multiplex amplicon data.

## Installation
1. [Install Nextflow](https://www.nextflow.io/docs/latest/install.html)

2. [Install Docker](https://docs.docker.com/engine/install/)

3. Download the appropriate Dorado installer from the [repo](https://github.com/nanoporetech/dorado#installation). The path to the executable will be ```<path to downloaded folder>/bin/dorado```

4. (Optional) Download the appropriate Dorado model from the [repo](https://github.com/nanoporetech/dorado/#available-basecalling-models)
    ```
    # Download all models
    dorado download --model all
    # Download particular model
    dorado download --model <model>
    ```
    If a pre-downloaded model path is not provided to the pipeline, the model specified by the `--basecall_model` parameter will be downloaded on the fly.

5. Download the appropriate Clair3 model from the [Rerio repo](https://github.com/nanoporetech/rerio?tab=readme-ov-file#clair3-models) (you will need Python3)
    
    First clone the repo:
    ```
    git clone https://github.com/nanoporetech/rerio
    ```
    This contains scripts to download the model(s) to ```clair3_models/<config>```
    ```
    #  Download all models
    python3 download_model.py --clair3
    #  Download particular model
    python3 download_model.py --clair3 clair3_models/<config>_model
    ```
    Each downloaded model can be found in the repo directory under ```clair3_models/<config>```

6. Clone the repository

    Please note you will also need access to the [assorted-sub-workflows](https://github.com/sanger-pathogens/assorted-sub-workflows/tree/b065b17b0ee663483fa14c09fc9b1dede9afa8ba) and [nextflowtool](https://github.com/sanger-pathogens/nextflowtool/tree/74b25a9346d243db662caccb777296061400b65a) submodule repos

    With SSH (will need an SSH key- instructions [here](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)):
    ```
    git clone --recurse-submodules git@github.com:sanger-pathogens/long-read-ampliseq.git
    ```

    With HTTPS:
    
    ```
    git clone --recurse-submodules https://github.com/sanger-pathogens/long-read-ampliseq.git
    ```
    You may need to enter a personal access token for authentication in place of a password- please see this [guide](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens). 

## Usage
```
nextflow run long-read-ampliseq/main.nf \
--raw_read_dir <directory containing FAST5/POD5 files> \
--reference <reference fasta> \
--primers <fasta containing primers> \
--target_regions_bed <BED file containing target regions> \
--additional_metadata <CSV mapping sample IDs to barcodes> \
--dorado_local_path <absolute path to Dorado executable> \
--clair3_model <path to Clair3 model> \
-profile docker
```
The [examples](examples) folder contains some example files

### Other parameters:

#### Basecalling
- --basecall = "true"
- --basecall_model = "dna_r10.4.1_e8.2_400bps_hac@v4.3.0"
- --basecall_model_path = ""
- --trim_adapters = "all"
- --barcode_kit_name = ["SQK-NBD114-24"] (currently this can only be edited via the config file)
- --read_format = "fastq"

#### Saving output files
- --keep_sorted_bam = true
- --save_fastqs = true
- --save_trimmed = true
- --save_too_short = true
- --save_too_long = true

#### QC
- --qc_reads = true
- --min_qscore = 9
- --cutadapt_args = "-e 0.15 --no-indels --overlap 18"
- --lower_read_length_cutoff = 450
- --upper_read_length_cutoff = 800
- --coverage_reporting_thresholds = "1,2,8,10,25,30,40,50,100"
- --coverage_filtering_threshold = "25"
- --multiqc_config = ""

#### Variant calling
- --clair3_min_coverage = "5"

###### Consensus curation
- --min_ref_gt_qual = 1
- --min_alt_gt_qual = 1

###### Tree building
- --remove_recombination = false
- --raxml_base_model = 'GTR+G4'
- --raxml_threads = 2


## Profiles

### -profile standard

This is the default profile and is intended to allow the pipeline to run (with internet access) on the Sanger HPC (farm). It ensures the pipeline can run with the LSF job scheduler and uses singularity images for dependencies management, as well as the latest versions of the pipeline base configuration (from [PaM Info common config file](https://github.com/sanger-pathogens/nextflow-commons/blob/master/configs/nextflow.config)) and Dorado models.

### -profile laptop

This profile has been provided specifically for use with an macOS laptop, on which a docker daemon has been installed via Docker Desktop. The profile has several features that allow the pipeline to be used offline:
- A local copy of a configuration file that is otherwise downloaded.
- Default local paths for the following parameters: `--dorado_local_path`, `--clair3_model`, `--basecall_model_path`

Should you need to run the pipeline offline, it is best to make use of pre-populated dependency caches. These can be created with any of the supported profiles (e.g. -profile docker) and involves running the pipeline once to completion.

You can override the default paths using the command line parameters directly when invoking nextflow or supplying an additional config file in which these parameters are set, using the `-c my_custom.config` nextflow option.

## Support
Please contact PaM Informatics for support through our [helpdesk portal](https://jira.sanger.ac.uk/servicedesk/customer/portal/16) or for external users please reach out by email: pam-informatics@sanger.ac.uk
