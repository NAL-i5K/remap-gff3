#! /usr/bin/env python
# Contributed by Li-Mei Chiang <dytk2134 [at] gmail [dot] com> (2018)
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

import sys
import re
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
if not logger.handlers:
    lh = logging.StreamHandler()
    lh.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    logger.addHandler(lh)

__version__ = '1.0'

def fasta_file_sequence_length(fasta_file):
    # get the length of the sequence in the fasta_file
    # sequence_length = {'SequenceID': {'length': Sequence_length, 'state': True}
    # True: finished; False: processing
    sequence_length = dict()
    sequence_id = None
    with open(fasta_file, 'rb') as fasta_file_f:
        for line in fasta_file_f:
            line = line.strip()
            if len(line) != 0:
                if line[0] == '>':
                    lines = line.split(' ')
                    if sequence_id == None:
                        # the first sequence
                        sequence_id = lines[0][1:]
                    elif sequence_id != lines[0][1:]:
                        # next sequence
                        sequence_length[sequence_id]['state'] = True
                        sequence_id = lines[0][1:]
                    if sequence_id not in sequence_length:
                        sequence_length[sequence_id] = {
                            'length': 0,
                            'state': False
                        }
                    else:
                        if sequence_length[sequence_id]['state'] == True:
                            logger.warning('Duplicate ID found! %s' % (sequence_id))
                else:
                    sequence_length[sequence_id]['length'] += len(line)
    return sequence_length
def main(alignment_file, target, query, output):
    logger.info('Reading target genome assembly: (%s)...\n', target)
    target_length_dict = fasta_file_sequence_length(target)

    logger.info('Reading query genome assembly: (%s)...\n', query)
    query_length_dict = fasta_file_sequence_length(query)

    chain_ID = 0
    total_alignments = 0
    filtered_alignments = 0
    removed_alignments = 0
    Size_different = 0

    out_f = open(output, 'w')

    with open(alignment_file, 'rb') as in_f:
        for line in in_f:
            line = line.strip()
            if len(line) != 0:
                if not line.startswith('#'):
                    total_alignments += 1
                    tokens = line.split('\t')
                    attribute = dict(re.findall('([^=;]+)=([^=;\n]+)', tokens[8]))
                    # focus on prefect alignments
                    if attribute['reciprocity'] != '3' and attribute['pct_identity_gap'] != '100':
                        removed_alignments += 1
                        continue
                    Target = attribute['Target'].split(' ')
                    # chain header
                    score = attribute['num_ident']
                    tName = Target[0]
                    try:
                        tSize = target_length_dict[Target[0]]['length']
                    except KeyError:
                        logger.warning('ID: %s not found in the target fasta file' % Target[0])
                        removed_alignments += 1
                        continue
                    tStrand = Target[3]
                    tStart = int(Target[1]) - 1
                    tEnd = int(Target[2])
                    qName = tokens[0]
                    try:
                        qSize = query_length_dict[tokens[0]]['length']
                    except KeyError:
                        logger.warning('ID: %s not found in the query fasta file' % tokens[0])
                        removed_alignments += 1
                        continue
                    qStrand = tokens[6]
                    if qStrand == '-':
                        qStart = qSize - int(tokens[4])
                        qEnd = qSize - int(tokens[3]) + 1
                    else:
                        qStart = int(tokens[3]) - 1
                        qEnd = int(tokens[4])

                    # check if target strand is '+'. If not, transforms the strand coordinates
                    if tStrand != '+':
                        tStrand = '+'
                        tStart = tStart
                        tEnd = tEnd
                        if qStrand == '-':
                            qStrand = '+'
                        else:
                            qStrand = '-'
                        qStart = qSize - int(tokens[4])
                        qEnd = qSize - int(tokens[3]) + 1

                    if (tEnd - tStart) == (qEnd - qStart):
                        chain_ID += 1
                        ID = str(chain_ID)
                        filtered_alignments += 1
                        header = ['chain', score, tName, str(tSize), tStrand, str(tStart), str(tEnd), qName, str(qSize), qStrand, str(qStart), str(qEnd), ID]
                        out_f.write(' '.join(header) + '\n')
                        size = str(tEnd - tStart)
                        out_f.write(size + '\n')
                        out_f.write('\n')
                    else:
                        Size_different += 1
    out_f.close()
    logger.info('Total alignment: %d', total_alignments)
    logger.info('Filtered alignment: %d', filtered_alignments)
    logger.info('Removed alignment: %d', removed_alignments)
    logger.info('size different alignment: %d', Size_different)

if __name__ == '__main__':
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\

    Quick start:
    %(prog)s -t_fa example_file/target.fa -q_fa example_file/query.fa -g alignment.gff3 -o chain.txt
    """))

    parser.add_argument('-a', '--alignment_file', type=str, help='NCBI\'s whole-genome alignments(gff3 format).', required=True)
    parser.add_argument('-t_fa', '--target_fasta', type=str, help='Target genome assembly', required=True)
    parser.add_argument('-q_fa', '--query_fasta', type=str, help='Query genome assembly', required=True)
    parser.add_argument('-o', '--output', type=str, help='output chain format file', required=True)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s' + __version__)

    args = parser.parse_args()
    main(alignment_file=args.alignment_file, target=args.target_fasta, query=args.query_fasta, output=args.output)