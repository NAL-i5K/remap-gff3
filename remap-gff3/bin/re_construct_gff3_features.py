#! /usr/bin/env python
import copy
from gff3tool.lib.gff3 import Gff3
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
if not logger.handlers:
    lh = logging.StreamHandler()
    lh.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    logger.addHandler(lh)

__version__ = '1.0'

def build_parentID_dict(gff3):
    # {parentID:[aaa,bbb,ccc]}
    parentID_dict = dict()
    for line in gff3.lines:
        if line['line_type'] == 'feature':
            if 'attributes' in line and 'Parent' in line['attributes']:
                for parent in line['attributes']['Parent']:
                    if parent not in parentID_dict:
                        parentID_dict[parent] = [line]
                    else:
                        parentID_dict[parent].append(line)

    return parentID_dict
def polypeptide_re_construct(old_gff3, new_gff3, tmp_identifier=False, report=None):
    old_CDS = dict()
    new_CDS = dict()
    old_polypeptide = dict()
    new_polypeptide = set()
    for line in old_gff3.lines:
        if line['line_type'] == 'feature':
            if line['type'] == 'polypeptide':
                try:
                    for parent in line['attributes']['Parent']:
                        old_polypeptide[parent] = line
                except KeyError:
                    continue
            elif line['type'] == 'CDS':
                try:
                    for parent in line['attributes']['Parent']:
                        if parent not in old_CDS:
                            old_CDS[parent] = [line]
                        else:
                            old_CDS[parent].append(line)
                except KeyError:
                    continue
    if len(old_polypeptide) != 0:
        for line in new_gff3.lines:
            if line['line_type'] == 'feature':
                if line['type'] == 'polypeptide':
                    try:
                        for parent in line['attributes']['Parent']:
                            new_polypeptide.add(parent)
                    except KeyError:
                        continue
                elif line['type'] == 'CDS':
                    try:
                        for parent in line['attributes']['Parent']:
                            if parent not in new_CDS:
                                new_CDS[parent] = [line]
                            else:
                                new_CDS[parent].append(line)
                    except KeyError:
                        continue
        eofindex = len(new_gff3.lines) - 1
        for parent_ID in new_CDS:
            if parent_ID not in new_polypeptide:
                new_CDS_features = new_CDS[parent_ID]
                old_CDS_features = old_CDS[parent_ID]
                if len(new_CDS_features) == len(old_CDS_features):
                    newpolypeptide = copy.deepcopy(old_polypeptide[parent_ID])
                    cPos = []
                    strand_set = set()
                    seqid_set = set()
                    new_CDS_features_sort = []
                    old_CDS_features_sort = []
                    try:
                        if new_CDS_features[0]['strand'] == '-':
                            new_CDS_features_sort = sorted(new_CDS_features, key=lambda x: x['end'], reverse=True)
                        else:
                            new_CDS_features_sort = sorted(new_CDS_features, key=lambda x: x['start'])
                        if old_CDS_features[0]['strand'] == '-':
                            old_CDS_features_sort = sorted(old_CDS_features, key=lambda x: x['end'], reverse=True)
                        else:
                            old_CDS_features_sort = sorted(old_CDS_features, key=lambda x: x['start'])
                    except IndexError:
                        continue
                    flag = False

                    for idx, child in enumerate(new_CDS_features_sort):
                        try:
                            if tmp_identifier:
                                if child['attributes']['tmp_identifier'] != old_CDS_features_sort[idx]['attributes']['tmp_identifier']:
                                    flag = True
                                    break
                            else:
                                if child['attributes']['ID'] != old_CDS_features_sort[idx]['attributes']['ID']:
                                    flag = True
                                    break
                        except KeyError:
                            flag = True
                        cPos.append(child['start'])
                        cPos.append(child['end'])
                        strand_set.add(child['strand'])
                        seqid_set.add(child['seqid'])
                    if flag == True:
                        continue
                    # update strand/start/end
                    if len(strand_set) == 1 and len(seqid_set) == 1:
                        eofindex += 1
                        newpolypeptide['line_index'] = eofindex
                        newpolypeptide['seqid'] = list(seqid_set)[0]
                        newpolypeptide['strand'] = list(strand_set)[0]
                        newpolypeptide['start'] = min(cPos)
                        newpolypeptide['end'] = max(cPos)
                        newpolypeptide['parents'] = []
                        try:
                            for parent in newpolypeptide['attributes']['Parent']:
                                newpolypeptide['parents'].append(new_gff3.features[parent])
                        except KeyError:
                            newpolypeptide['parents'] = []
                        new_gff3.features[newpolypeptide['attributes']['ID']].append(newpolypeptide)
                        new_gff3.lines.append(newpolypeptide)
                        if report != None:
                            write_features(newpolypeptide, report)

def re_construct(old_gff3, new_gff3, tmp_identifier=False, report=None):
    old_parentID_dict = build_parentID_dict(old_gff3)
    new_parentID_dict = build_parentID_dict(new_gff3)
    eofindex = len(new_gff3.lines) - 1
    new_parent = set()
    for parent_ID in new_parentID_dict:
        try:
            features = new_gff3.features[parent_ID]
        except KeyError:
            features = []
        # parent feature not in the gff3 file
        # might need to re-construct
        if len(features) == 0:
            new_children = new_parentID_dict[parent_ID]
            old_children = old_parentID_dict[parent_ID]
            if len(new_children) == len(old_children):
                # all child feature are prefectly re-mapped, but the parents aren't
                # get parent information from original gff3 file
                try:
                    old_features = old_gff3.features[parent_ID]
                except KeyError:
                    old_features = []
                # the ID of a parent feature should be unique
                if len(old_features) == 1:
                    newparent = copy.deepcopy(old_features[0])
                    # clean old child feature and update the child's parent list
                    newparent['children'] = []
                    cPos = []
                    strand_set = set()
                    seqid_set = set()
                    new_children_sort = []
                    old_children_sort = []
                    try:
                        if new_children[0]['strand'] == '-':
                            new_children_sort = sorted(new_children, key=lambda x: (x['end'],x['type']), reverse=True)
                        else:
                            new_children_sort = sorted(new_children, key=lambda x: (x['start'],x['type']))
                        if old_children[0]['strand'] == '-':
                            old_children_sort = sorted(old_children, key=lambda x: (x['end'],x['type']), reverse=True)
                        else:
                            old_children_sort = sorted(old_children, key=lambda x: (x['start'],x['type']))
                    except IndexError:
                        continue
                    flag = False
                    for idx, child in enumerate(new_children_sort):
                        try:
                            if tmp_identifier:
                                if child['attributes']['tmp_identifier'] != old_children_sort[idx]['attributes']['tmp_identifier']:
                                    flag = True
                                    break
                            else:
                                if child['attributes']['ID'] != old_children_sort[idx]['attributes']['ID']:
                                    flag = True
                                    break
                        except KeyError:
                            flag = True
                        cPos.append(child['start'])
                        cPos.append(child['end'])
                        strand_set.add(child['strand'])
                        seqid_set.add(child['seqid'])
                        newparent['children'].append(child)
                    if flag == True:
                        continue
                    # update strand/start/end
                    if len(strand_set) == 1 and len(seqid_set) == 1:
                        eofindex += 1
                        newparent['line_index'] = eofindex
                        newparent['seqid'] = list(seqid_set)[0]
                        newparent['strand'] = list(strand_set)[0]
                        newparent['start'] = min(cPos)
                        newparent['end'] = max(cPos)
                        newparent['children'] = sorted(newparent['children'], key=lambda x: (x['line_index']))
                        newparent['parents'] = []
                        try:
                            for parent in newparent['attributes']['Parent']:
                                new_parent.add(parent)
                                newparent['parents'].append(new_gff3.features[parent])
                        except KeyError:
                            newparent['parents'] = []
                        new_gff3.features[parent_ID].append(newparent)
                        new_gff3.lines.append(newparent)
                        for child in new_children:
                            if len(child['parents']) == 0:
                                child['parents'] = []
                            else:
                                child['parents'] = [x for x in child['parents'] if x]
                            child['parents'].append(new_gff3.features[parent_ID])

                        if report != None:
                            write_features(newparent, report)
                    else:
                        logger.warning('child features of %s are not on the same scaffold/strand.' % parent_ID)
    if len(new_parent) != 0:
        re_construct(old_gff3, new_gff3,logger, report)

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
                    wrote_lines.add(line_data['line_index'])
                descendants = gff3.descendants(root)

                for descendant in descendants:
                    if descendant['line_index'] in wrote_lines:
                        continue
                    write_features(descendant, out_gff)
                    wrote_lines.add(descendant['line_index'])
def main(old_gff, new_gff, output_gff, re_construct_features, tmp_identifier):
    logger.info('Reading original GFF3 file: (%s)...\n', old_gff)
    old_gff3 = Gff3(gff_file=old_gff)

    logger.info('Reading updated GFF3 file: (%s)...\n', new_gff)
    new_gff3 = Gff3(gff_file=new_gff)

    if re_construct_features:
        out_f = open(re_construct_features, 'w')
    else:
        out_f = None
    polypeptide_re_construct(old_gff3=old_gff3, new_gff3=new_gff3, tmp_identifier=tmp_identifier, report=out_f)
    re_construct(old_gff3=old_gff3, new_gff3=new_gff3, tmp_identifier=tmp_identifier, report=out_f)
    logger.info('Generating the re-constructed gff3 file: (%s)...\n', output_gff)
    write_gff3(new_gff3, output_gff)

    if re_construct_features:
        out_f.close()

if __name__ == '__main__':
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\

    Quick start:
    %(prog)s -old_g example_file/example1.gff3 -new_g example_file/example1_CrossMap_filtered.gff3 -og example1_re_construct.gff3 -re report.txt
    """))

    parser.add_argument('-old_g', '--old_gff', type=str, help='The original gff3 file', required=True)
    parser.add_argument('-new_g', '--new_gff', type=str, help='The updated  gff3 file', required=True)
    parser.add_argument('-og', '--output_gff', type=str, help='output re-construct gff3 file', required=True)
    parser.add_argument('-re', '--re_construct_features', type=str, help='output re-construct features')
    parser.add_argument('-tmp_ID', '--tmp_identifier', action="store_true", help='Use unique temporary identifier (tmp_identifier attribute) as identifier. (Default: use ID attribute).', default=False)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()
    main(args.old_gff, args.new_gff, args.output_gff, args.re_construct_features, args.tmp_identifier)


