# agp_to_chain.py

## Introduction

agp_to_chain.py is aim to convert an [AGP](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/) file to UCSC [chain](https://genome.ucsc.edu/goldenpath/help/chain.html) file.

### Dependencies

    - python

### Method
Extract the following information from the AGP files to build a Chain file.

#### Chain file
- Header
    - score = 1000
    - tName = component_id
    - tSize = The size of target sequence in target fasta file
    - tStrand = +
    - tStart = component_beg - 1 (AGP file is 1-based; chain file is 0-based)
    - tEnd = component_end
    - qName = object
    - qSize = The size of query sequence in query fasta file
    - qStrand = orientation
    - qStart = object_beg - 1 (if qStrand = +) or qSize - object_end + 1 (if qStrand = -)
    - qEnd = object_end (if qStrand = +) or qSize - object_beg (if qStrand = -)
    - id = auto increment number
- alignment data line
    - size = tEnd - tStart

##### Note
- This script will only handle nine-column, tab-delimited lines in agp file.
- gap line (component_type = N or U) will be ignored
- components with unknown orientation (?, 0 or na) are treated as if they had + orientation

### Usage
#### help and usage messages

```shell
usage: agp_to_chain.py [-h] -a AGP_FILE -t_fa TARGET_FASTA -q_fa QUERY_FASTA
                       -o OUTPUT [-v]

Quick start:
agp_to_chain.py -t_fa example_file/target.fa -q_fa example_file/query.fa -a example_file/example.agp -o chain.txt

optional arguments:
  -h, --help            show this help message and exit
  -a AGP_FILE, --agp_file AGP_FILE
                        Input agp file
  -t_fa TARGET_FASTA, --target_fasta TARGET_FASTA
                        Target genome assembly
  -q_fa QUERY_FASTA, --query_fasta QUERY_FASTA
                        Query genome assembly
  -o OUTPUT, --output OUTPUT
                        output chain format file
  -v, --version         show program's version number and exit
```

#### example
##### AGP file
```
NW_017236740.1  54415   55534   1   W   KN240439.1  1811    2929    -
NW_017236740.1  1   457 1   W   KN239297.1  1   457 +
NW_017237251.1  4819    5961    1   W   KN239297.1  16137   17279   +
```

#### Chain file
```
chain 1000 KN240439.1 5289 + 1811 2929 NW_017236740.1 124178 - 68645 69763 4
1118

chain 1000 KN239297.1 17279 + 0 457 NW_017240107.1 457 + 0 457 5
457

chain 1000 KN239297.1 17279 + 16137 17279 NW_017237251.1 120242 + 4819 5961 6
1142

```

