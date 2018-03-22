#! /usr/bin/env python
# Contributed by Li-Mei Chiang <dytk2134 [at] gmail [dot] com> (2018)
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
if not logger.handlers:
    lh = logging.StreamHandler()
    lh.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    logger.addHandler(lh)

__version__ = '1.0'

if __name__ == '__main__':
    import os
    import subprocess
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\

    Quick start:
    python %(progs)s -a alignment.gff3 -t_fa target.fa -q_fa query.fa -g A.gff3 B.gff3 C.gff3
    """))

    parser.add_argument('-a', '--alignment_file', type=str, help='NCBI\'s whole-genome alignments(gff3 format).', required=True)
    parser.add_argument('-t_fa', '--target_fasta', type=str, help='Target genome assembly', required=True)
    parser.add_argument('-q_fa', '--query_fasta', type=str, help='Query genome assembly', required=True)
    parser.add_argument('-g', '--input_gff', nargs='+', type=str, help='List one or more GFF3 files to be updated.', required=True)
    parser.add_argument('-tmp_ID', '--tmp_identifier', action="store_true", help='Generate a unique temporary identifier for all the feature in the input gff3 files. (Default: False)', default=False)
    parser.add_argument('-chain', '--chain_file', type=str, help='Input a ready-made chain file.')
    parser.add_argument('-tmp', '--temp', action="store_false", help='Store all the intermediate files/temporary files into [alignment_filename]_tmp/ directory. (Default: False)')
    parser.add_argument('-u', '--updated_postfix', default='_updated', help='The filename postfix for updated features (default: "_updated")')
    parser.add_argument('-r', '--removed_postfix', default='_removed', help='The filename postfix for removed features (default: "_removed")')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()

    temp_dir = '_'.join([os.path.splitext(args.alignment_file)[0], 'tmp'])
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    chain_file = None
    if not args.chain_file:
        import gff3_to_chain
        # generate a chain file and store in the [alignment_filename]_tmp/ directory
        chain_file = '.'.join([os.path.splitext(args.alignment_file)[0], 'chain'])
        logger.info('Generate a chain file: (%s)', chain_file)
        gff3_to_chain.main(alignment_file=args.alignment_file, target=args.target_fasta, query=args.query_fasta, output=chain_file)
    else:
        chain_file = args.chain_file
    
    for gff3 in args.input_gff:
        gff3_filename, gff3_extension = os.path.splitext(gff3)
        CrossMap_mapped = 
        if args.tmp_identifier:
            # add tmp_identifier attribute to all the features in the input gff3 files
        # run CrossMap
        subprocess.Popen(['CrossMap.py', 'gff', chain_file, gff3, ).wait()





    