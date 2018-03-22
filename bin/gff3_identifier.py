#! /usr/local/bin/python2.7
# Contributed by Li-Mei Chiang <dytk2134 [at] gmail [dot] com> (2017)
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

def idgenerator(prefix, lastnumber, digitlen):
    lastnumber += 1
    idnum = str(lastnumber)
    if len(idnum) < digitlen:
        adddigit = digitlen-len(idnum)
        for i in range(adddigit):
            idnum = str(0) + idnum
    result={}
    result['ID'] = prefix + idnum
    result['maxnum'] = lastnumber
    return(result)


def missing_ID(gff3, logger):
    missing_ID_dict = dict()
    for line in gff3.lines:
        if line['line_type'] == 'feature':
            if 'attributes' in line and 'ID' in line['attributes']:
                continue
            else:
                if 'Parent' in line['attributes']:
                    for parent in line['attributes']['Parent']:
                        if parent not in missing_ID_dict:
                            missing_ID_dict[parent] = [line]
                        else:
                            missing_ID_dict[parent].append(line)
                # all the parent features should has a ID attributes

    return missing_ID_dict

def re_assign(gff3, logger):
    missing_ID_dict = missing_ID(gff3, logger)
    for ID in gff3.features:
        features = []
        for feature in gff3.features[ID]:
            if 'Parent' in feature['attributes']:
                for parent in feature['attributes']['Parent']:
                    if parent in missing_ID_dict:
                        features.extend(missing_ID_dict[parent])
        features.extend(gff3.features[ID])
        sort_features = []
        # not unique within the scope of the GFF file
        if len(features) > 1:
            if features[0]['strand'] == '-':
                sort_features = sorted(features, key=lambda x: x['end'], reverse=True)
            else:
                sort_features = sorted(features, key=lambda x: x['start'])
        digitlen = 3
        IDnumber_dict = dict()
        for feature in sort_features:
            if feature['type'] not in IDnumber_dict:
                IDnumber_dict[feature['type']] = 0
            if 'ID' not in feature['attributes']:
                prefix = '-'.join([feature['attributes']['Parent'][0], feature['type']])
 
            else:
                prefix = feature['attributes']['ID']
            newid = idgenerator(prefix, IDnumber_dict[feature['type']], digitlen)
            IDnumber_dict[feature['type']] = newid['maxnum']
            feature['attributes']['ID'] = newid['ID']

def write_features(line, out_f):
    field_keys = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase']
    field_list = [str(line[k]) for k in field_keys]
    attribute_list = []
    for k, v in line['attributes'].items():
        if isinstance(v, list):
            v = ','.join(v)
        attribute_list.append('%s=%s' % (str(k), str(v)))
    field_list.append(';'.join(attribute_list))
    out_f.write('\t'.join(field_list) + '\n')

def write_gff3(gff3, out_f):
    wrote_lines = set()
    root_lines = []
    with open(out_f, 'w') as out_gff:
        for line in gff3.lines:
            if line['line_type'] == 'feature':
                parent_list = [x for x in line['parents'] if x]
                # feature that don't have 'Parent' attribute
                # feature that have 'Parent' attribute, but its parent feature not in the gff3 file
                if not parent_list:
                    root_lines.append(line)
        for root in root_lines:
            if root['line_index'] in wrote_lines:
                continue
            else:
                try:
                    root_feature = gff3.features[root['attributes']['ID']]
                except KeyError:
                    root_feature = [root]
                for line_data in root_feature:
                    write_features(line_data, out_gff)
                descendants = gff3.descendants(root)
                for descendant in descendants:
                    if descendant['line_index'] in wrote_lines:
                        continue
                    write_features(descendant, out_gff)

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
    python2.7 re_construct_gff3_features.py -old_g old.gff3 -new_g new.gff3 -og re_construct.gff3 -r report.txt
    """))

    parser.add_argument('-gff', '--gff', type=str, help='input gff3 file', required=True)
    parser.add_argument('-og', '--output_gff', type=str, help='output gff3 file', required=True)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()

    logger_stderr.info('Reading input GFF3 file: (%s)...\n', args.gff)
    gff3 = Gff3(gff_file=args.gff, logger=None)

    re_assign(gff3=gff3, logger=logger_stderr)
    write_gff3(gff3, args.output_gff)



