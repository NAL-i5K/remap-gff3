# ! /usr/bin/env python
# Contributed by Li-Mei Chiang <dytk2134 [at] gmail [dot] com> (2018)
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

import sys
import re
import logging
import copy

__version__ = '1.0'

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
    python2.7 filter_not_exact_match.py -g CrossMap_output.gff3 -log CrossMap.log -og filtered.gff3
    """))

    parser.add_argument('-g', '--gff', type=str, help='The output gff3 file from CrossMap', required=True)
    parser.add_argument('-og', '--output_gff', type=str, help='output filtered gff3 file', required=True)
    parser.add_argument('-log', '--log_f', type=str, help='output log from CrossMap')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()

    not_exact_match = set()
    with open(args.log_f, 'r') as log_file:
        for line in log_file:
            line = line.strip()
            if len(line) != 0:
                if line[0] != '#':
                    token = line.split("\t")
                    if len(token) < 10:
                        continue
                    match_state = token[9]
                    if 'not exact match' in match_state:
                        attribute_dict = dict(re.findall('([^=;]+)=([^=;\n]+)', token[8]))
                        if 'ID' in attribute_dict:
                            not_exact_match.add(attribute_dict['ID'])
    out_gff = open(args.output_gff, 'w')
    with open(args.gff, 'r') as gff_f:
        for line in gff_f:
            line = line.strip()
            if len(line) != 0:
                if line[0] != '#':
                    token = line.split("\t")
                    if len(token) != 9:
                        continue
                    else:
                        attribute_dict = dict(re.findall('([^=;]+)=([^=;\n]+)', token[8]))
                        if 'ID' in attribute_dict:
                            if attribute_dict['ID'] not in not_exact_match:
                                out_gff.write(line + '\n')

    out_gff.close()





