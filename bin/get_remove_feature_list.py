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
import copy
from gff3 import Gff3

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
    python2.7 get_remove_feature_list.py -old_g old.gff3 -new_g new.gff3 -og remove_features_.gff3
    """))

    parser.add_argument('-old_g', '--old_gff', type=str, help='The original gff3 file', required=True)
    parser.add_argument('-new_g', '--new_gff', type=str, help='The updated  gff3 file', required=True)
    parser.add_argument('-og', '--output_gff', type=str, help='output removed feature gff3 file', required=True)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()

    logger_stderr.info('Reading original GFF3 file: (%s)...\n', args.old_gff)
    old_gff3 = Gff3(gff_file=args.old_gff, logger=None)

    logger_stderr.info('Reading updated GFF3 file: (%s)...\n', args.new_gff)
    new_gff3 = Gff3(gff_file=args.new_gff, logger=None)

    with open(args.output_gff, 'w') as out_f:
        for feature_id in old_gff3.features:
            if feature_id not in new_gff3.features:
                for feature in old_gff3.features[feature_id]:
                    out_f.write(feature['line_raw'])

