#! /usr/bin/env python

import re
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
if not logger.handlers:
    lh = logging.StreamHandler()
    lh.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    logger.addHandler(lh)

__version__ = '1.0'

def get_features_ID(gff, tmp_identifier):
    features_ID_set = set()
    ID_key = 'ID'
    if tmp_identifier:
        ID_key = 'tmp_identifier'
    with open(gff, 'rb') as in_f:
        for line in in_f:
            line = line.strip()
            if len(line) != 0:
                if not line.startswith('#'):
                    tokens = line.split('\t')
                    attributes = dict(re.findall('([^=;]+)=([^=;\n]+)', tokens[8]))
                    try:
                        features_ID_set.add(attributes[ID_key])
                    except KeyError:
                        continue
    return features_ID_set


def output_remove_features(old_gff, new_gff, output_gff, tmp_identifier):
    # get features ID
    new_features_ID_set = get_features_ID(new_gff, tmp_identifier)
    out_f = open(output_gff, 'w')
    ID_key = 'ID'
    if tmp_identifier:
        ID_key = 'tmp_identifier'
    with open(old_gff, 'rb') as in_f:
        for line in in_f:
            line = line.strip()
            if len(line) != 0:
                if not line.startswith('#'):
                    tokens = line.split('\t')
                    attributes = dict(re.findall('([^=;]+)=([^=;\n]+)', tokens[8]))
                    try:
                        if attributes[ID_key] not in new_features_ID_set:
                            if tmp_identifier:
                                del attributes['tmp_identifier']
                            attributes_list = list()
                            for key in attributes:
                                attributes_list.append('%s=%s' % (str(key), str(attributes[key])))
                            # update cloumn 9
                            tokens[8] = ';'.join(attributes_list)
                            out_f.write('\t'.join(tokens) + '\n')
                    except KeyError:
                        continue

if __name__ == '__main__':
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\

    Quick start:
    %(prog)s -old_g example_file/example1.gff3 -new_g example_file/example1_updated.gff3 -og example1_removed.gff3
    """))

    parser.add_argument('-old_g', '--old_gff', type=str, help='The original gff3 file', required=True)
    parser.add_argument('-new_g', '--new_gff', type=str, help='The updated  gff3 file', required=True)
    parser.add_argument('-og', '--output_gff', type=str, help='output removed feature gff3 file', required=True)
    parser.add_argument('-tmp_ID', '--tmp_identifier', action="store_true", help='Use unique temporary identifier (tmp_identifier attribute) as identifier. (Default: use ID attribute).', default=False)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()
    output_remove_features(args.old_gff, args.new_gff, args.output_gff, args.tmp_identifier)

