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

__version__ = '1.0'


def fasta_file_sequence_length(fasta_file, logger):
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

if __name__ == '__main__':
    logger_stderr = logging.getLogger(__name__+'stderr')
    logger_stderr.setLevel(logging.INFO)
    stderr_handler = logging.StreamHandler()
    stderr_handler.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    logger_stderr.addHandler(stderr_handler)
    logger_null = logging.getLogger(__name__+'null')
    null_handler = logging.NullHandler()
    logger_null.addHandler(null_handler)
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\

    Quick start:
    python %(prog)s -t_fa target.fasta -q_fa query.fasta -g example.gff3 -o chain.txt
    """))

    parser.add_argument('-g', '--gff', type=str, help='NCBI remap GFF3 file, gff3 format', required=True)
    parser.add_argument('-t_fa', '--target_fasta', type=str, help='target alignment fasta file', required=True)
    parser.add_argument('-q_fa', '--query_fasta', type=str, help='query alignment fasta file', required=True)
    parser.add_argument('-o', '--output', type=str, help='output chain format file', required=True)
    parser.add_argument('-r', '--report', type=str, help='output strand convert report')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()

    logger_stderr.info('Reading target fasta file: (%s)...\n', args.target_fasta)
    target_length_dict = fasta_file_sequence_length(args.target_fasta, logger_stderr)
    query_length_dict = fasta_file_sequence_length(args.query_fasta, logger_stderr)

    logger_stderr.info('Reading query fasta file: (%s)...\n', args.query_fasta)
    #error_dict example: {'Emr0001': [[15,16],[13]],'Esf0005': [[17]]}
    error_dict = {}
    #line_num_dict example: {3: ['Emr0001','Esf0003'], 15: ['Emr0026']}
    line_num_dict = {}
    chain_ID = 0
    total_alignments = 0
    filtered_alignments = 0
    removed_alignments = 0
    Size_different = 0
    out_f = open(args.output, 'w')
    out_r = None
    if args.report:
        out_r = open(args.report, 'w')

    with open(args.gff, 'r') as in_f:
        for line in in_f:
            Convert_strand = False
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
                        logger_stderr.warning('ID: %s not found in the target fasta file' % Target[0])
                        removed_alignments += 1
                        continue
                    tStrand = Target[3]
                    tStart = int(Target[1]) - 1
                    tEnd = int(Target[2])
                    qName = tokens[0]
                    try:
                        qSize = query_length_dict[tokens[0]]['length']
                    except KeyError:
                        logger_stderr.warning('ID: %s not found in the query fasta file' % tokens[0])
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
                        Convert_strand = True


                    if (tEnd - tStart) == (qEnd - qStart):
                        chain_ID += 1
                        ID = str(chain_ID)
                        filtered_alignments += 1
                        header = ['chain', score, tName, str(tSize), tStrand, str(tStart), str(tEnd), qName, str(qSize), qStrand, str(qStart), str(qEnd), ID]
                        out_f.write(' '.join(header) + '\n')
                        size = str(tEnd - tStart)
                        out_f.write(size + '\n')
                        out_f.write('\n')
                        if out_r != None and Convert_strand == True:
                            out_r.write(' '.join(header) + '\n')
                    else:
                        Size_different += 1
    out_f.close()
    if out_r != None:
        out_r.close()
    logger_stderr.info('Total alignment: %d', total_alignments)
    logger_stderr.info('Filtered alignment: %d', filtered_alignments)
    logger_stderr.info('Removed alignment: %d', removed_alignments)
    logger_stderr.info('size different alignment: %d', Size_different)
