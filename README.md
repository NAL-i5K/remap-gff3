# remap-gff3 : Python programs for updating gff3 coordinates to new assembly versions
## Background
It is reasonable to expect that these assemblies will be updated, sometimes from entirely new data. A continuing problem is mapping existing datasets (in particular the manually curated annotations that we facilitate) to the new assemblies. For this purpose, we generate a workflow to update gff3 coordinates to new assembly versions.

## General Workflow
1. Generate whole-genome alignment file.
    * For now, we will use NCBI's whole-genome alignments, available here: ftp://ftp.ncbi.nlm.nih.gov/pub/remap/
2. Filter whole-genome alignment results based on quality criteria.
    * For now, we will focus on perfect alignments. 
    * Use NCBI's alignment results in gff3 format to filter out alignments with reciprocity=3 and pct_identity_gap=100
3. Generate input file for coordinate remap program based on filtered alignment results.
    * generate file in UCSC chain format. https://genome.ucsc.edu/goldenpath/help/chain.html
    * Use crossmap to update gff3 coordinates. http://crossmap.sourceforge.net/
    * remove all the not exact match features from the CrossMap output
    * re-construct the parent features for the models where all child features are perfectly re-mapped, but the parents aren't
    * run gff3_QC to generate QC report for re-constructed gff3 file
    * remove all the incorrectly merged gene parents (Ema0009) and incorrectly split parents (Emr0002) from QC report
    * run gff3_fix to correct GFF3 format errors
    * get updated and removed GFF3 files

## Current functions
### generate file in UCSC chain format
* bin/gff_to_chain.py
    * Quick start: `python gff_to_chain.py -t_fa target.fasta -q_fa query.fasta -g example.gff3 -o chain.txt`
### re-construct missing features in gff3
* bin/re_construct_gff3_deature.py
    * Quick start: `python re_construct_gff3_deature.py -old_g old.gff3 -new_g new.gff3 -og re_construct.gff3 -re report.txt`
### get removed features
* bin/get_remove_feature.py
    * Quick start: `python get_remove_feature.py -old_g old.gff3 -new_g new.gff3 -og remove_features_.gff3`
### re-map gff3 files
* bin/remap-gff3.py
    * Quick start: `python remap-gff3.py -a alignment.gff3 -t_fa target.fa -q_fa query.fa -tmp_ID -g A.gff3 B.gff3 C.gff3`

## Installation
### Prerequisite
* python2.7
* Perl
### Install packages/modules via the setup.py
```
python setup.py install
```
### Install packages/modules in a Docker container
#### Install Docker
Follow [instructions](https://docs.docker.com/install/) to install Docker for your environment.
#### Docker image
Build a image from a Dockerfile  
```
git clone https://github.com/NAL-i5K/remap-gff3.git
cd remap-gff3
docker build -t remap-gff3-image .
docker run -itp 8000:8000 remap-gff3-image
```
Get a pre-built image from DockerHub
```
docker pull dytk2134/remap-gff3-image
docker run -itp 8000:8000 dytk2134/remap-gff3-image
```
