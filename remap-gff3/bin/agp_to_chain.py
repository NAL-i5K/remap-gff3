#! /usr/bin/env python
import logging
import gzip

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
    if fasta_file.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open
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



def main(agp_file, target, query, output):
    logger.info('Reading target genome assembly: (%s)...\n', target)
    target_length_dict = fasta_file_sequence_length(target)

    logger.info('Reading query genome assembly: (%s)...\n', query)
    query_length_dict = fasta_file_sequence_length(query)
    # gap line will be ignored
    gap_type = set(['N', 'U'])
    chain_ID = 0
    component_line = 0 # total number of component lines
    removed_line = 0
    with open(output, 'w') as out_f:
        with open(agp_file, 'r') as in_f:
            for line in in_f:
                line = line.strip()
                if line:
                    if not line.startswith('#'):
                        tokens = line.split('\t')
                        if len(tokens) == 9:
                            # only handle component line
                            if tokens[4] not in gap_type:
                                component_line += 1
                                score = 1000
                                tName = tokens[5] # component_id
                                try:
                                    tSize = target_length_dict[tokens[5]]['length']
                                except KeyError:
                                    logger.warning('ID: %s not found in the target fasta file' % tokens[5])
                                    removed_line += 1
                                    continue
                                tStrand = '+'
                                tStart = int(tokens[6]) - 1 # component_beg - 1
                                tEnd = int(tokens[7]) # component_end
                                qName = tokens[0] # object
                                try:
                                    qSize = query_length_dict[tokens[0]]['length']
                                except KeyError:
                                    logger.warning('ID: %s not found in the query fasta file' % tokens[0])
                                    removed_line += 1
                                    continue
                                if tokens[8] == '-':
                                    qStrand = '-'
                                    qStart = qSize - int(tokens[2]) + 1
                                    qEnd = qSize - int(tokens[1])
                                else:
                                    # components with unknown orientation (?, 0 or na) are treated as if they had + orientation
                                    qStrand = '+'
                                    qStart = int(tokens[1]) - 1
                                    qEnd = int(tokens[2])
                                chain_ID += 1
                                header = ['chain', str(score), tName, str(tSize), tStrand, str(tStart), str(tEnd), qName, str(qSize), qStrand, str(qStart), str(qEnd), str(chain_ID)]
                                out_f.write(' '.join(header) + '\n')
                                size = str(tEnd - tStart)
                                out_f.write(size + '\n')
                                out_f.write('\n')
    logger.info('Total component lines: %d', component_line)
    logger.info('Removed lines: %d', removed_line)
    logger.info('Converted lines: %d', component_line - removed_line)


if __name__ == '__main__':
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\

    Quick start:
    %(prog)s -t_fa example_file/target.fa -q_fa example_file/query.fa -a example_file/example.agp -o chain.txt
    """))

    parser.add_argument('-a', '--agp_file', type=str, help='Input agp file', required=True)
    parser.add_argument('-t_fa', '--target_fasta', type=str, help='Target genome assembly', required=True)
    parser.add_argument('-q_fa', '--query_fasta', type=str, help='Query genome assembly', required=True)
    parser.add_argument('-o', '--output', type=str, help='output chain format file', required=True)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s' + __version__)

    args = parser.parse_args()
    main(agp_file=args.agp_file, target=args.target_fasta, query=args.query_fasta, output=args.output)
