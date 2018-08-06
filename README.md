# remap-gff3 : Python programs for updating gff3 coordinates to new assembly versions

[![Build Status](https://travis-ci.org/NAL-i5K/remap-gff3.svg?branch=master)](https://travis-ci.org/NAL-i5K/remap-gff3)

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
        * To clean up errors that happen after model reconstruction, we will only fix the following errors in the QC report:
            * Ema0001: Parent feature start and end coordinates exceed those of child features
            * Ema0003: This feature is not contained within the parent feature coordinates
            * Ema0006: Wrong phase
            * Ema0007: CDS and parent feature on different strands
            * Emr0001: Duplicate transcript found
            * Esf0014: ##gff-version" missing from the first line
    * run gff3_fix to correct GFF3 format errors
    * get updated and removed GFF3 files
    * generate a summary report and a list of removed features
    * run gff3_QC to generate QC report for updated GFF3 files

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

`pip install git+https://github.com/NAL-i5K/remap-gff3.git`

#### Troubleshooting

If you install `remap-gff3` locally(e.g. `pip install . --user`), please make sure the **~/.local/bin** is in your **$PATH**.

### Docker

#### Install Docker

Follow [instructions](https://docs.docker.com/install/) to install Docker for your environment.

#### Docker image

For container deployment, you can build a image from a Dockerfile or get a pre-built image from DockerHub.

##### Build a image from a Dockerfile

1. `git clone https://github.com/NAL-i5K/remap-gff3.git`
2. `cd remap-gff3`
3. `docker build -t remap-gff3 .`
4. `docker run -it remap-gff3`

##### Get a pre-built image from DockerHub

1. `docker pull nali5k/remap-gff3`
2. `docker run -it nali5k/remap-gff3`

## Usage

* remap-gff3.py

``` shell
usage: remap-gff3.py [-h] -a ALIGNMENT_FILE -t_fa TARGET_FASTA -q_fa
                     QUERY_FASTA -g INPUT_GFF [INPUT_GFF ...] -s
                     SOURCE -b BUILDNAME [-tmp_ID]
                     [-chain CHAIN_FILE] [-tmp] [-u UPDATED_POSTFIX]
                     [-r REMOVED_POSTFIX] [-v]

Quick start:
remap-gff3.py -a example_file/alignment.gff3 -t_fa example_file/target.fa -q_fa example_file/query.fa -dir output -tmp_ID -g example_file/example1.gff3 example_file/example2.gff3 -s NCBI -b Hazt_2.0

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
  -s SOURCE, --source SOURCE
                        Source of the assembly (e.g. NCBI). This generates 
			a pragma line recommended by the gff3 specification
			(https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
  -b BUILDNAME, --buildName BUILDNAME
                        The genome assembly (build) name used for the
                        coordinates. This generates a pragma line 
			recommended by the gff3 specification 
			(https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
  -tmp_ID, --tmp_identifier
                        Generate a unique temporary identifier for all the
                        feature in the input gff3 files. (Default: False)
  -summary, --summary_report
                        Generate a document that summarizes the change in
                        feature types after remapping and lists the removed
                        features
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

### Example

* Quick start: `remap-gff3.py -a example_file/alignment.gff3 -t_fa example_file/target.fa -q_fa example_file/query.fa -dir output -tmp_ID -g example_file/example1.gff3 example_file/example2.gff3 -summary -s NCBI -b Hazt_2.0`

* NCBI's whole-genome  alignments(gff3 format)

```shell
NW_017236740.1  Genbank match   54416   55533   .       +       .       ID=448e258e-d822-4ed5-89de-f1d1e5086dfd;Target=KN240439.1 1812 2929 -;best_on_query=1;best_on_query_same_unit=1;best_on_subject=1;best_on_subject_same_unit=1;gap_count=0;genomic_to_genomic=1;num_ident=1118;num_mismatch=0;pct_coverage=21.1382;pct_coverage_hiqual=21.1382;pct_ident_quantized=98;pct_identity_gap=100;pct_identity_gapopen_only=100;pct_identity_ungap=100;reciprocity=3;same_unit_reciprocity=3
NW_017240107.1  Genbank match   1       457     .       +       .       ID=184c62dc-3a03-49fb-92c3-0a2dcb71cc57;Target=KN239297.1 1 457 +;best_on_query=1;best_on_query_same_unit=1;best_on_subject=1;best_on_subject_same_unit=1;gap_count=0;genomic_to_genomic=1;num_ident=457;num_mismatch=0;pct_coverage=2.64483;pct_coverage_hiqual=2.64483;pct_ident_quantized=98;pct_identity_gap=100;pct_identity_gapopen_only=100;pct_identity_ungap=100;reciprocity=3;same_unit_reciprocity=3;sequence_matched_component_align=1
NW_017237251.1  Genbank match   4820    5961    .       +       .       ID=894eaa1e-b8a4-4777-8ae6-feec1c20c870;Target=KN239297.1 16138 17279 +;best_on_query=1;best_on_query_same_unit=1;best_on_subject=1;best_on_subject_same_unit=1;gap_count=0;genomic_to_genomic=1;num_ident=1142;num_mismatch=0;pct_coverage=6.60918;pct_coverage_hiqual=6.60918;pct_ident_quantized=98;pct_identity_gap=100;pct_identity_gapopen_only=100;pct_identity_ungap=100;reciprocity=3;same_unit_reciprocity=3
NW_017235329.1	Genbank	match	115007	119836	.	+	.	ID=799b2a11-39b7-499e-aedd-126f543ef23b;Target=KN234834.1 557151 561964 +;gap_count=4;genomic_to_genomic=1;num_ident=4782;num_mismatch=28;pct_coverage=0.0905731;pct_coverage_hiqual=0.0905731;pct_ident_quantized=98;pct_identity_gap=98.9243;pct_identity_gapopen_only=99.3353;pct_identity_ungap=99.4179;reciprocity=1;same_unit_reciprocity=1;Gap=M820 D4 M1519 D15 M50 I4 M13 D1 M2408
```

* original gff3 file

```shell
KN234853.1      OGSv1.0 gene    612439  613408  .       +       .       ID=HAZT001062;Dbxref=I5KNAL:HAZT001062;Name=HAZT001062;method=
Maker2
KN234853.1      OGSv1.0 mRNA    612439  613408  .       +       .       Name=HAZT001062-RA;Parent=HAZT001062;ID=HAZT001062-RA;Dbxref=G
O:0004591,GO:0030976,GO:0006099,EC:1.2.4.2;method=Maker2
KN234853.1      OGSv1.0 exon    612439  612576  .       +       .       method=Maker2;Parent=HAZT001062-RA;ID=HAZT001062-RA-EXON01
KN234853.1      OGSv1.0 CDS     612442  612576  .       +       0       method=Maker2;Dbxref=GO:0004591,GO:0030976,GO:0006099,EC:1.2.4
.2;Parent=HAZT001062-RA;ID=HAZT001062-RA-CDS
###
KN234853.1      OGSv1.0 gene    961300  961512  .       -       .       ID=HAZT001066;Name=HAZT001066;Dbxref=I5KNAL:HAZT001066;method=
Maker2
KN234853.1      OGSv1.0 mRNA    961300  961512  .       -       .       ID=HAZT001066-RA;Name=HAZT001066-RA;Parent=HAZT001066;method=Maker2
KN234853.1      OGSv1.0 exon    961300  961512  .       -       .       ID=HAZT001066-RA-EXON01;Parent=HAZT001066-RA;method=Maker2
KN234853.1      OGSv1.0 CDS     961300  961512  .       -       0       ID=HAZT001066-RA-CDS;Parent=HAZT001066-RA;method=Maker2
###
```

* updated gff3 file

```shell
##genome-build NCBI Hazt_2.0
NW_017250073.1  OGSv1.0 gene    157594  157806  .       -       .       method=Maker2;Dbxref=I5KNAL:HAZT001066;ID=HAZT001066;Name=HAZT001066
NW_017250073.1  OGSv1.0 mRNA    157594  157806  .       -       .       method=Maker2;ID=HAZT001066-RA;Parent=HAZT001066;Name=HAZT001066-RA
NW_017250073.1  OGSv1.0 exon    157594  157806  .       -       .       ID=HAZT001066-RA-EXON01;Parent=HAZT001066-RA;method=Maker2
NW_017250073.1  OGSv1.0 CDS     157594  157806  .       -       0       ID=HAZT001066-RA-CDS;Parent=HAZT001066-RA;method=Maker2
```

* removed gff3 file

```shell
KN234853.1      OGSv1.0 gene    612439  613408  .       +       .       ID=HAZT001062;Dbxref=I5KNAL:HAZT001062;Name=HAZT001062;method=
Maker2
KN234853.1      OGSv1.0 mRNA    612439  613408  .       +       .       Name=HAZT001062-RA;Parent=HAZT001062;ID=HAZT001062-RA;Dbxref=G
O:0004591,GO:0030976,GO:0006099,EC:1.2.4.2;method=Maker2
KN234853.1      OGSv1.0 exon    612439  612576  .       +       .       method=Maker2;Parent=HAZT001062-RA;ID=HAZT001062-RA-EXON01
KN234853.1      OGSv1.0 CDS     612442  612576  .       +       0       method=Maker2;Dbxref=GO:0004591,GO:0030976,GO:0006099,EC:1.2.4
.2;Parent=HAZT001062-RA;ID=HAZT001062-RA-CDS
```

If `-summary` is given, remap-gff3.py will generate two extra summary files: **[input gff3 filename]_summary.tsv** and **[input gff3 filename]_removed.tsv**.

* summary report ([input gff3 filename]_summary.tsv)

```shell
Feature_type    Original_count  New_count       Retained(%)
gene            30              26              86.67%
mRNA            29              25              86.21%
exon            111             88              79.28%
CDS             107             85              79.44%
polypeptide     29              25              86.21%
three_prime_utr 1               1               100.00%
miRNA           1               1               100.00%
```

* removed features list ([input gff3 filename]_removed.tsv)

```shell
type             ID                     Name            owner
gene            HAZT001062              HAZT001062      NA
mRNA            HAZT001062-RA           HAZT001062-RA   NA
exon            HAZT001062-RA-EXON01    NA              NA
exon            HAZT001062-RA-EXON02    NA              NA
exon            HAZT001062-RA-EXON03    NA              NA
CDS             HAZT001062-RA-CDS001    NA              NA
CDS             HAZT001062-RA-CDS002    NA              NA
CDS             HAZT001062-RA-CDS003    NA              NA
polypeptide     HAZT001062-PA           NA              NA
```

## Additional Functionality

### generate file in UCSC chain format

* gff_to_chain.py
  * Quick start: `gff_to_chain.py -t_fa example_file/target.fa -q_fa example_file/query.fa -a example_file/alignment.gff3 -o chain.txt`
  * [Full documentation](docs/gff_to_chain.md)

### re-construct missing features in gff3

* re_construct_gff3_features.py
  * Quick start: `re_construct_gff3_features.py -old_g example_file/example1.gff3 -new_g example_file/example1_CrossMap_filtered.gff3 -og example1_re_construct.gff3 -re report.txt`
  * [Full documentation](docs/re_construct_gff3_features.md)

### get removed features

* get_remove_feature.py
  * Quick start: `get_remove_feature.py -old_g example_file/example1.gff3 -new_g example_file/example1_updated.gff3 -og example1_removed.gff3`
  * [Full documentation](docs/get_remove_feature.md)

