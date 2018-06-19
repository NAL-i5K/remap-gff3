# get removed features

Compare the identifier of the feature between old and new gff3 file and output the removed features.

## get_remove_feature.py

### Dependencies

    - python

### Usage
#### help and usage messages

```shell
usage: get_remove_feature.py [-h] -old_g OLD_GFF -new_g NEW_GFF -og OUTPUT_GFF
                             [-tmp_ID] [-v]

Quick start:
get_remove_feature.py -old_g example_file/example1.gff3 -new_g example_file/example1_updated.gff3 -og example1_removed.gff3

optional arguments:
  -h, --help            show this help message and exit
  -old_g OLD_GFF, --old_gff OLD_GFF
                        The original gff3 file
  -new_g NEW_GFF, --new_gff NEW_GFF
                        The updated gff3 file
  -og OUTPUT_GFF, --output_gff OUTPUT_GFF
                        output removed feature gff3 file
  -tmp_ID, --tmp_identifier
                        Use unique temporary identifier (tmp_identifier attribute) as identifier. (Default: use ID attribute)
  -v, --version         show program's version number and exit
```

#### example
##### original gff3 file
```
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

#### updated gff3 file
```
NW_017250073.1  OGSv1.0 gene    157594  157806  .       -       .       method=Maker2;Dbxref=I5KNAL:HAZT001066;ID=HAZT001066;Name=HAZT001066
NW_017250073.1  OGSv1.0 mRNA    157594  157806  .       -       .       method=Maker2;ID=HAZT001066-RA;Parent=HAZT001066;Name=HAZT001066-RA
NW_017250073.1  OGSv1.0 exon    157594  157806  .       -       .       ID=HAZT001066-RA-EXON01;Parent=HAZT001066-RA;method=Maker2
NW_017250073.1  OGSv1.0 CDS     157594  157806  .       -       0       ID=HAZT001066-RA-CDS;Parent=HAZT001066-RA;method=Maker2
```

#### removed gff3 file
```
KN234853.1      OGSv1.0 gene    612439  613408  .       +       .       ID=HAZT001062;Dbxref=I5KNAL:HAZT001062;Name=HAZT001062;method=
Maker2
KN234853.1      OGSv1.0 mRNA    612439  613408  .       +       .       Name=HAZT001062-RA;Parent=HAZT001062;ID=HAZT001062-RA;Dbxref=G
O:0004591,GO:0030976,GO:0006099,EC:1.2.4.2;method=Maker2
KN234853.1      OGSv1.0 exon    612439  612576  .       +       .       method=Maker2;Parent=HAZT001062-RA;ID=HAZT001062-RA-EXON01
KN234853.1      OGSv1.0 CDS     612442  612576  .       +       0       method=Maker2;Dbxref=GO:0004591,GO:0030976,GO:0006099,EC:1.2.4
.2;Parent=HAZT001062-RA;ID=HAZT001062-RA-CDS
```
