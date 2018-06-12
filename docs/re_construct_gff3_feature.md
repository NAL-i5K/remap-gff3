# Re-construct missing features in gff3

## Background

After using CrossMap to update gff3 coordinates, we found that there are some parent features in a model might be erroneously removed. Since the parent feature (e.g. gene, mRNA) is usually longer than the child feature(e.g. exon, CDS) and the alignment information from NCBI's alignment files contain both long and short sequence alignments, a parent feature might map more than two alignments in the chain file we generate from NCBI's alignment files. Those multiple map features will be removed. To make sure, there are no erronously removed feature, we developed the re_construct_gff3_feature.py to re-construct the parent features for the models where all child features are perfectly re-mapped, but the parents aren't.

## re_construct_gff3_feature.py

### Dependencies

    - python2.7
        - Packages/Modules:
            - gff3

### Method

1. collect all the child features from original and updated Gff3 file
2. sort those child features by Start/End coordinate and the feature type
3. map the identifier between original and updated Gff3 file
4. check if the parent feature match the re-construct criteria
    - re-construct criteria
        - the number of child features of the model in these two Gff3 file should be the same
        - all the child features of the model should on the same scaffold and on the same strand
        - after remapping, the order of child features of the model in original and updated Gff3 file should still be the same
5. re-construct the parent feature
    - Sequence ID: same with the child features
    - Source: collect from the feature in original Gff3 file
    - Type: collect from the feature in original Gff3 file
    - Start: get the minimum coordinate from the child features
    - End: get the maximum coordinate from the child features
    - Score: collect from the feature in original Gff3 file
    - Strand: same with the child features
    - Phase: collect from the feature in original Gff3 file
    - Attributes: collect from the feature in original Gff3 file
6. continuously re-construct parent feature until there is no feature pass the re-construct criteria
#### Note
Assume that each feature in the GFF3 file have an unique identifier.

### Usage
#### help and usage messages

```shell
usage: re_construct_gff3_features.py [-h] -old_g OLD_GFF -new_g NEW_GFF -og
                                     OUTPUT_GFF [-re RE_CONSTRUCT_FEATURES]
                                     [-tmp_ID] [-v]

Quick start:
python2.7 re_construct_gff3_features.py -old_g old.gff3 -new_g new.gff3 -og re_construct.gff3 -re report.txt

optional arguments:
  -h, --help            show this help message and exit
  -old_g OLD_GFF, --old_gff OLD_GFF
                        The original gff3 file
  -new_g NEW_GFF, --new_gff NEW_GFF
                        The updated gff3 file
  -og OUTPUT_GFF, --output_gff OUTPUT_GFF
                        output re-construct gff3 file
  -re RE_CONSTRUCT_FEATURES, --re_construct_features RE_CONSTRUCT_FEATURES
                        output re-construct features
  -tmp_ID, --tmp_identifier
                        Use unique temporary identifier (tmp_identifier attribute) as identifier. (Default: use ID attribute).
  -v, --version         show program's version number and exit
```

#### example
##### original GFF3 file
```shell
KN234853.1      OGSv1.0 gene    961300  961512  .       -       .       ID=HAZT001066;Name=HAZT001066;Dbxref=I5KNAL:HAZT001066;method=Maker2
KN234853.1      OGSv1.0 mRNA    961300  961512  .       -       .       ID=HAZT001066-RA;Name=HAZT001066-RA;Parent=HAZT001066;method=Maker2
KN234853.1      OGSv1.0 exon    961300  961512  .       -       .       ID=HAZT001066-RA-EXON01;Parent=HAZT001066-RA;method=Maker2
KN234853.1      OGSv1.0 CDS     961300  961512  .       -       0       ID=HAZT001066-RA-CDS;Parent=HAZT001066-RA;method=Maker2
KN234853.1      OGSv1.0 polypeptide     961300  961512  .       -       .       ID=HAZT001066-PA;Parent=HAZT001066-RA;method=Maker2
###
```
##### updated GFF3 file (Output from CrossMap)
```shell
NW_017250073.1  OGSv1.0 exon    157594  157806  .       -       .       ID=HAZT001066-RA-EXON01;Parent=HAZT001066-RA;method=Maker2
NW_017250073.1  OGSv1.0 CDS     157594  157806  .       -       0       ID=HAZT001066-RA-CDS;Parent=HAZT001066-RA;method=Maker2
###
```

##### re-constructed GFF3 file
```shell
NW_017250073.1  OGSv1.0 gene    157594  157806  .       -       .       method=Maker2;Dbxref=I5KNAL:HAZT001066;ID=HAZT001066;Name=HAZT001066
NW_017250073.1  OGSv1.0 mRNA    157594  157806  .       -       .       method=Maker2;ID=HAZT001066-RA;Parent=HAZT001066;Name=HAZT001066-RA
NW_017250073.1  OGSv1.0 exon    157594  157806  .       -       .       ID=HAZT001066-RA-EXON01;Parent=HAZT001066-RA;method=Maker2
NW_017250073.1  OGSv1.0 CDS     157594  157806  .       -       0       ID=HAZT001066-RA-CDS;Parent=HAZT001066-RA;method=Maker2
NW_017250073.1  OGSv1.0 polypeptide     157594  157806  .       -       .       ID=HAZT001066-PA;Parent=HAZT001066-RA;method=Maker2
###
```