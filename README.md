# remap-gff3 : Python programs for updating gff3 coordinates to new assembly versions

## Background

It is reasonable to expect that genome assemblies will be updated, sometimes from entirely new data. A continuing problem is mapping existing datasets to the new assemblies. For this purpose, we generate a workflow to update gff3 coordinates to new assembly versions, tailored towards the manually curated annotations that the i5k Workspace@NAL facilitates. This workflow relies on existing services and programs - specifically, NCBI's whole-genome alignment service, and Crossmap - and packages them in a way that generates error-free gff3 models on the new assembly, for models that have 100% alignment between old and new assemblies.

## General Workflow

1. Generate whole-genome alignment file.
    * For now, the i5k Workspace is using NCBI's whole-genome alignments, available here: ftp://ftp.ncbi.nlm.nih.gov/pub/remap/

2. Filter whole-genome alignment results based on quality criteria.
    * For now, we will focus on perfect alignments.
    * Use NCBI's alignment results in gff3 format to filter out alignments with reciprocity=3 and pct_identity_gap=100
3. Generate a chain file for coordinate remap program based on filtered alignment results.
    * UCSC [chain format](https://genome.ucsc.edu/goldenpath/help/chain.html)
    * Generate a chain file using gff_to_chain.py
4. Use CrossMap to update gff3 coordinates
    * [CrossMap](http://crossmap.sourceforge.net/)
5. Post-process the output from CrossMap
    * remove all the not exact match features from the CrossMap output
    * re-construct the parent features for the models where all child features are perfectly re-mapped, but the parents aren't
    * run gff3_QC to generate QC report for re-constructed gff3 file
        * For now, we will remove all the incorrectly merged gene parents (Ema0009) and incorrectly split parents (Emr0002) from QC report, and then run gff3_fix to correct format errors.
    * run gff3_fix to correct GFF3 format errors
    * get updated and removed GFF3 files

## Dependencies

* python2.7
    * Packages/Modules:
        * CrossMap
        * gff3
        * gff3tool
* Perl
* zlib1g-dev
* liblzo2-dev


## Installation

There are two options for installing: directly from Github or Docker

### Install from Github

1. `git clone https://github.com/NAL-i5K/remap-gff3.git`
2. `cd remap-gff3`
3. `pip install .`

### Docker

#### Install Docker

Follow [instructions](https://docs.docker.com/install/) to install Docker for your environment.

#### Docker image

For container deployment, you can build a image from a Dockerfile or get a pre-built image from DockerHub.

##### Build a image from a Dockerfile

1. `git clone https://github.com/NAL-i5K/remap-gff3.git`
2. `cd remap-gff3`
3. `docker build -t remap-gff3-image .`
4. `docker run -itp 8000:8000 remap-gff3-image`

##### Get a pre-built image from DockerHub

1. `docker pull dytk2134/remap-gff3-image`
2. `docker run -itp 8000:8000 dytk2134/remap-gff3-image`

## Quick Start

* remap-gff3.py

``` shell
usage: remap-gff3.py [-h] -a ALIGNMENT_FILE -t_fa TARGET_FASTA -q_fa
                     QUERY_FASTA -g INPUT_GFF [INPUT_GFF ...] [-tmp_ID]
                     [-chain CHAIN_FILE] [-tmp] [-u UPDATED_POSTFIX]
                     [-r REMOVED_POSTFIX] [-v]

Quick start:
remap-gff3.py -a example_file/alignment.gff3 -t_fa example_file/target.fa -q_fa example_file/query.fa -dir output -tmp_ID -g example_file/example1.gff3 example_file/example2.gff3

optional arguments:
  -h, --help            show this help message and exit
  -a ALIGNMENT_FILE, --alignment_file ALIGNMENT_FILE
                        NCBI's whole-genome alignments(gff3 format).
  -t_fa TARGET_FASTA, --target_fasta TARGET_FASTA
                        Target genome assembly
  -q_fa QUERY_FASTA, --query_fasta QUERY_FASTA
                        Query genome assembly
  -dir OUT_DIR, --out_dir OUT_DIR
                        Output directory
  -g INPUT_GFF [INPUT_GFF ...], --input_gff INPUT_GFF [INPUT_GFF ...]
                        List one or more GFF3 files to be updated.
  -tmp_ID, --tmp_identifier
                        Generate a unique temporary identifier for all the
                        feature in the input gff3 files. (Default: False)
  -chain CHAIN_FILE, --chain_file CHAIN_FILE
                        Input a ready-made chain file.
  -tmp, --temp          Store all the intermediate files/temporary files into
                        [alignment_filename]_tmp/ directory. (Default: False)
  -u UPDATED_POSTFIX, --updated_postfix UPDATED_POSTFIX
                        The filename postfix for updated features (default:
                        "_updated")
  -r REMOVED_POSTFIX, --removed_postfix REMOVED_POSTFIX
                        The filename postfix for removed features (default:
                        "_removed")
  -v, --version         show program's version number and exit
```

## Additional Functionality

### generate file in UCSC chain format

* gff_to_chain.py
  * Quick start: `gff_to_chain.py -t_fa example_file/target.fa -q_fa example_file/query.fa -g alignment.gff3 -o chain.txt`

### re-construct missing features in gff3

* re_construct_gff3_deature.py
  * Quick start: `re_construct_gff3_deature.py -old_g old.gff3 -new_g new.gff3 -og re_construct.gff3 -re report.txt`

### get removed features

* get_remove_feature.py
  * Quick start: `get_remove_feature.py -old_g old.gff3 -new_g new.gff3 -og remove_features.gff3`
