# generate file in UCSC chain format

## Introduction

Generate file in UCSC [chain format](https://genome.ucsc.edu/goldenpath/help/chain.html) is part of the gff3 coordinates update workflow. In order to avoid probable mapping error, we focus on perfect alignments and generate a chain file for coordinate remap program from the perfect alignments. For now, the i5k Workspace is using NCBI's whole-genome alignments, available here: ftp://ftp.ncbi.nlm.nih.gov/pub/remap/.

## quality criteria

- reciprocity=3
- pct_identity_gap=100
- The size of the alignment between old and new assemblies should be the same

## gff_to_chain.py

### Dependencies

    - python

### Method
Extract the following information from a NCBI's alignment files to build a Chain file.
#### Chain file
- Header
    - score = `num_ident` attributes
    - tName = Sequence id in `Target` attribute
    - tSize = The size of target sequence in target fasta file
    - tStrand = Strand in `Target` attribute
    - tStart = Start in `Target` attribute - 1 (GFF file is 1-based; chain file is 0-based)
    - tEnd = End in `Target` attribute
    - qName = Column 1(Sequence id)
    - qSize = The size of query sequence in query fasta file
    - qStrand = Column 7(Strand)
    - qStart = Column 4 (Start) - 1 (GFF file is 1-based; chain file is 0-based)
    - qEnd = Column 5 (End)
    - id = auto increment number
- alignment data line (assume all alignment are single ungapped segments)
    - size = tEnd - tStart
##### Note
If **tStrand** is on the negative strand, this script will transform it to positive strand and then transform the query to negative strand.
- qStart_new = qSize -qEnd_old
- qEnd_new = qSize - qStart_old

### Usage
#### help and usage messages

```shell
usage: gff_to_chain.py [-h] -a ALIGNMENT_FILE -t_fa TARGET_FASTA -q_fa
                       QUERY_FASTA -o OUTPUT [-v]

Quick start:
python gff_to_chain.py -t_fa example_file/target.fa -q_fa example_file/query.fa -g alignment.gff3 -o chain.txt

optional arguments:
  -h, --help            show this help message and exit
  -a ALIGNMENT_FILE, --alignment_file ALIGNMENT_FILE
                        NCBI's whole-genome alignments(gff3 format).
  -t_fa TARGET_FASTA, --target_fasta TARGET_FASTA
                        Target genome assembly
  -q_fa QUERY_FASTA, --query_fasta QUERY_FASTA
                        Query genome assembly
  -o OUTPUT, --output OUTPUT
                        output chain format file
  -v, --version         show program's version number and exit
```

#### example
##### NCBI's whole-genome  alignments(gff3 format)
```
NW_017236740.1  Genbank match   54416   55533   .       +       .       ID=448e258e-d822-4ed5-89de-f1d1e5086dfd;Target=KN240439.1 1812 2929 -;best_on_query=1;best_on_query_same_unit=1;best_on_subject=1;best_on_subject_same_unit=1;gap_count=0;genomic_to_genomic=1;num_ident=1118;num_mismatch=0;pct_coverage=21.1382;pct_coverage_hiqual=21.1382;pct_ident_quantized=98;pct_identity_gap=100;pct_identity_gapopen_only=100;pct_identity_ungap=100;reciprocity=3;same_unit_reciprocity=3
NW_017240107.1  Genbank match   1       457     .       +       .       ID=184c62dc-3a03-49fb-92c3-0a2dcb71cc57;Target=KN239297.1 1 457 +;best_on_query=1;best_on_query_same_unit=1;best_on_subject=1;best_on_subject_same_unit=1;gap_count=0;genomic_to_genomic=1;num_ident=457;num_mismatch=0;pct_coverage=2.64483;pct_coverage_hiqual=2.64483;pct_ident_quantized=98;pct_identity_gap=100;pct_identity_gapopen_only=100;pct_identity_ungap=100;reciprocity=3;same_unit_reciprocity=3;sequence_matched_component_align=1
NW_017237251.1  Genbank match   4820    5961    .       +       .       ID=894eaa1e-b8a4-4777-8ae6-feec1c20c870;Target=KN239297.1 16138 17279 +;best_on_query=1;best_on_query_same_unit=1;best_on_subject=1;best_on_subject_same_unit=1;gap_count=0;genomic_to_genomic=1;num_ident=1142;num_mismatch=0;pct_coverage=6.60918;pct_coverage_hiqual=6.60918;pct_ident_quantized=98;pct_identity_gap=100;pct_identity_gapopen_only=100;pct_identity_ungap=100;reciprocity=3;same_unit_reciprocity=3
NW_017235329.1	Genbank	match	115007	119836	.	+	.	ID=799b2a11-39b7-499e-aedd-126f543ef23b;Target=KN234834.1 557151 561964 +;gap_count=4;genomic_to_genomic=1;num_ident=4782;num_mismatch=28;pct_coverage=0.0905731;pct_coverage_hiqual=0.0905731;pct_ident_quantized=98;pct_identity_gap=98.9243;pct_identity_gapopen_only=99.3353;pct_identity_ungap=99.4179;reciprocity=1;same_unit_reciprocity=1;Gap=M820 D4 M1519 D15 M50 I4 M13 D1 M2408
```

#### Chain file
```
chain 1118 KN240439.1 5289 + 1811 2929 NW_017236740.1 124178 - 68645 69763 4
1118

chain 457 KN239297.1 17279 + 0 457 NW_017240107.1 457 + 0 457 5
457

chain 1142 KN239297.1 17279 + 16137 17279 NW_017237251.1 120242 + 4819 5961 6
1142

```
- The **pct_identity_gap** and **reciprocity** of the alignment, ID=799b2a11-39b7-499e-aedd-126f543ef23b did not match the quality criteria(pct_identity_gap=100 and reciprocity=3). Therefore, this alignment won't be converted to a chain alignment.