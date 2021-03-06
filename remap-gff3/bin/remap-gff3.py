#! /usr/bin/env python
import os
import re
import sys
import subprocess
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
if not logger.handlers:
    lh = logging.StreamHandler()
    lh.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    logger.addHandler(lh)

__version__ = '1.0'

def tmp_identifier(in_gff, out_gff):
    import uuid
    out_f = open(out_gff, 'w')
    with open(in_gff, 'r') as in_f:
        for line in in_f:
            line = line.strip()
            if line:
                if not line.startswith('#'):
                    tokens = line.split('\t')
                    if len(tokens) != 9:
                        logger.warning('Features should contain 9 fields, got %d: %s' % (len(tokens) - 1, repr(tokens[1:])))
                        continue
                    attributes = dict(re.findall('([^=;]+)=([^=;\n]+)', tokens[8]))
                    # add tmp_identifier attribute to feature in gff3 file
                    temp_id = str(uuid.uuid1())
                    attributes['tmp_identifier'] = temp_id
                    attributes_list = list()
                    for key in attributes:
                        attributes_list.append('%s=%s' % (str(key), str(attributes[key])))
                    # update cloumn 9
                    tokens[8] = ';'.join(attributes_list)
                    out_f.write('\t'.join(tokens) + '\n')
                elif line.startswith('##FASTA'):
                    # ignore embedded ##FASTA sequence
                    break
                else:
                    out_f.write(line + '\n')

def remove_tmpID(in_gff, out_gff):
    out_f = open(out_gff, 'w')
    with open(in_gff, 'r') as in_f:
        for line in in_f:
            line = line.strip()
            if line:
                if not line.startswith('#'):
                    tokens = line.split('\t')
                    if len(tokens) != 9:
                        logger.warning('Features should contain 9 fields, got %d: %s' % (len(tokens) - 1, repr(tokens[1:])))
                        continue
                    attributes = dict(re.findall('([^=;]+)=([^=;\n]+)', tokens[8]))
                    del attributes['tmp_identifier']
                    attributes_list = list()
                    for key in attributes:
                        attributes_list.append('%s=%s' % (str(key), str(attributes[key])))
                    # update cloumn 9
                    tokens[8] = ';'.join(attributes_list)
                    out_f.write('\t'.join(tokens) + '\n')
                elif line.startswith('##FASTA'):
                    # ignore embedded ##FASTA sequence
                    break
                else:
                    out_f.write(line + '\n')

def filter_not_exact_match(CrossMap_mapped_file, CrossMap_log_file, filtered_file,tmp_identifier):
    not_exact_match = set()
    ID_set = set()
    with open(CrossMap_log_file, 'r') as log:
        for line in log:
            line = line.strip()
            if line:
                if line[0] != '#':
                    token = line.split("\t")
                    if len(token) < 10:
                        continue
                    match_state = token[9]
                    attribute_dict = dict(re.findall('([^=;]+)=([^=;\n]+)', token[8]))
                    try:
                        if not tmp_identifier:
                            if attribute_dict['ID'] in ID_set:
                                logger.error('All features should have an unique ID. Please run the command with -tmp_ID.')
                                sys.exit(0)
                            elif 'not exact match' not in match_state:
                                ID_set.add(attribute_dict['ID'])
                    except KeyError:
                        logger.error('All features should have an unique ID. Please run the command with -tmp_ID.')
                        sys.exit(0)
                    if 'not exact match' in match_state:
                        if tmp_identifier:
                            if 'tmp_identifier' in attribute_dict:
                                not_exact_match.add(attribute_dict['tmp_identifier'])
                        else:
                            if 'ID' in attribute_dict:
                                not_exact_match.add(attribute_dict['ID'])

    out_f = open(filtered_file, 'w')
    with open(CrossMap_mapped_file, 'r') as in_f:
        for line in in_f:
            line = line.strip()
            if line:
                if line[0] != '#':
                    token = line.split("\t")
                    if len(token) != 9:
                        continue
                    else:
                        attribute_dict = dict(re.findall('([^=;]+)=([^=;\n]+)', token[8]))
                        if tmp_identifier:
                            if 'tmp_identifier' in attribute_dict and attribute_dict['tmp_identifier'] not in not_exact_match:
                                out_f.write(line + '\n')
                        else:
                            if 'ID' in attribute_dict and attribute_dict['ID'] not in not_exact_match:
                                out_f.write(line + '\n')
    out_f.close()

def get_summary(update_gff, remove_gff, summary_report, removed_list):
    # summary_dict = {'gene': {'Feature_type': 'gene','Original_count':100, 'New_count': 90, 'Retained(%)': 90}}
    summary_dict = dict()
    type_order = list()
    with open(update_gff, 'r') as in_gff:
        for line in in_gff:
            line = line.strip()
            if len(line) != 0:
                if line[0] != '#':
                    token = line.split("\t")
                    # type
                    if token[2] not in summary_dict:
                        type_order.append(token[2])
                        summary_dict[token[2]] = {
                            'Feature_type': token[2],
                            'Original_count': 0,
                            'New_count': 0,
                            'Retained(%)': 0
                        }
                    summary_dict[token[2]]['Original_count'] += 1
                    summary_dict[token[2]]['New_count'] += 1
    # removed = [{'type': 'gene', 'ID': 'MODEL0000001', 'Name': 'MODEL', 'owner': 'i5k'}]
    # list of dict
    removed = list()
    with open(remove_gff, 'r') as in_gff:
        for line in in_gff:
            line = line.strip()
            if len(line) != 0:
                if line[0] != '#':
                    token = line.split("\t")
                    attribute_dict = dict(re.findall('([^=;]+)=([^=;\n]+)', token[8]))
                    # type
                    if token[2] not in summary_dict:
                        type_order.append(token[2])
                        summary_dict[token[2]] = {
                            'Feature_type': token[2],
                            'Original_count': 0,
                            'New_count': 0,
                            'Retained(%)': 0
                        }
                    summary_dict[token[2]]['Original_count'] += 1
                    remove_dict = {
                        'type': token[2],
                        'ID': 'NA',
                        'Name': 'NA',
                        'owner': 'NA'
                    }
                    if 'ID' in attribute_dict:
                        remove_dict['ID'] = attribute_dict['ID']
                    if 'Name' in attribute_dict:
                        remove_dict['Name'] = attribute_dict['Name']
                    if 'owner' in attribute_dict:
                        remove_dict['owner'] = attribute_dict['owner']
                    removed.append(remove_dict)
    # calculate % retained
    for feature_type in summary_dict:
        Retained = 100 * float(summary_dict[feature_type]['New_count'])/float(summary_dict[feature_type]['Original_count'])
        summary_dict[feature_type]['Retained(%)'] = '%.2f%%' % Retained

    with open(summary_report, 'w') as out_f:
        header = ['Feature_type', 'Original_count', 'New_count', 'Retained(%)']
        # write out header
        outline = '\t'.join(header)
        out_f.write(outline + '\n')
        for feature_type in type_order:
            values = [str(summary_dict[feature_type][key]) for key in header]
            outline = '\t'.join(values)
            out_f.write(outline + '\n')
    with open(removed_list, 'w') as out_f:
        header = ['type', 'ID', 'Name', 'owner']
        outline = '\t'.join(header)
        out_f.write(outline + '\n')
        for rm in removed:
            values = [rm[key] for key in header]
            outline = '\t'.join(values)
            out_f.write(outline + '\n')
def add_string_to_file(in_f, add_line, line_num):
    current_line = 0
    lines = list()
    with open(in_f, 'r') as out_f:
        lines = out_f.readlines()

    with open(in_f, 'w') as out_f:
        for line in lines:
            current_line += 1
            if current_line == line_num:
                out_f.write(add_line + '\n')
                out_f.write(line)
            else:
                out_f.write(line)

if __name__ == '__main__':
    import argparse
    from textwrap import dedent
    import re_construct_gff3_features
    import get_remove_feature
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\

    Quick start:
    %(prog)s -a example_file/alignment.gff3 -t_fa example_file/target.fa -q_fa example_file/query.fa -dir output -tmp_ID -g example_file/example1.gff3 example_file/example2.gff3  -s NCBI -b Hazt_2.0
    """))

    parser.add_argument('-a', '--alignment_file', type=str, help='NCBI\'s whole-genome alignments(gff3 format).')
    parser.add_argument('-t_fa', '--target_fasta', type=str, help='Target genome assembly', required=True)
    parser.add_argument('-q_fa', '--query_fasta', type=str, help='Query genome assembly', required=True)
    parser.add_argument('-dir', '--out_dir', type=str, help='Output directory', required=True)
    parser.add_argument('-g', '--input_gff', nargs='+', type=str, help='List one or more GFF3 files to be updated.', required=True)
    parser.add_argument('-s', '--source', type=str, help='Source of the assembly.', required=True)
    parser.add_argument('-b', '--buildName', type=str, help='The genome assembly build name used for the coordinates.', required=True)
    parser.add_argument('-tmp_ID', '--tmp_identifier', action="store_true", help='Generate a unique temporary identifier for all the feature in the input gff3 files. (Default: False)', default=False)
    parser.add_argument('-summary', '--summary_report', action="store_true", help='Generate a document that summarizes the change in feature types after remapping and lists the removed features',default=False)
    parser.add_argument('-chain', '--chain_file', type=str, help='Input a ready-made chain file.')
    parser.add_argument('-tmp', '--temp', action="store_false", help='Store all the intermediate files/temporary files into [alignment_filename]_tmp/ directory. (Default: False)')
    parser.add_argument('-u', '--updated_postfix', default='_updated', help='The filename postfix for updated features (default: "_updated")')
    parser.add_argument('-r', '--removed_postfix', default='_removed', help='The filename postfix for removed features (default: "_removed")')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    if args.alignment_file:
        temp_dir_name = '%s_%s' % (os.path.splitext(os.path.basename(args.alignment_file))[0], 'tmp')
    elif args.chain_file:
        temp_dir_name = '%s_%s' % (os.path.splitext(os.path.basename(args.chain_file))[0], 'tmp')
    else:
        logger.error('please specify NCBI\'s whole-genome alignment file with -a or input a ready-made chain file with -chain.')
        sys.exit(0)

    temp_dir = os.path.join(args.out_dir, temp_dir_name)
    rm_tmp_list = []

    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    chain_file = None
    if not args.chain_file:
        import gff_to_chain
        # generate a chain file and store in the [alignment_filename]_tmp/ directory
        chain_file_name = '%s.%s' % (os.path.splitext(os.path.basename(args.alignment_file))[0], 'chain')
        chain_file = os.path.join(temp_dir, chain_file_name)
        logger.info('===== Generate a chain file: (%s) =====', chain_file)
        gff_to_chain.main(alignment_file=args.alignment_file, target=args.target_fasta, query=args.query_fasta, output=chain_file)
        rm_tmp_list.append(chain_file)
    else:
        chain_file = args.chain_file

    for gff3 in args.input_gff:
        gff3_filename, gff3_extension = os.path.splitext(os.path.basename(gff3))
        in_gff = gff3
        if args.tmp_identifier:
            # add tmp_identifier attribute to all the features in the input gff3 files
            out_gff = '%s/%s_%s%s' % (temp_dir, gff3_filename, 'mod', gff3_extension)
            tmp_identifier(in_gff, out_gff)
            in_gff = out_gff
            rm_tmp_list.append(out_gff)
        # run CrossMap
        logger.info('===== Use CrossMap to update gff3 coordinates: (%s) =====', gff3)
        CrossMap_mapped_file = '%s/%s_CrossMap%s' % (temp_dir, gff3_filename, gff3_extension)
        CrossMap_log_file = '%s/%s_CrossMap%s' % (temp_dir, gff3_filename, '.log')
        subprocess.Popen(['CrossMap.py', 'gff', chain_file, in_gff, CrossMap_mapped_file]).wait()
        log_file = open(CrossMap_log_file, 'w')
        subprocess.Popen(['CrossMap.py', 'gff', chain_file, in_gff], stdout=log_file).wait()
        log_file.close()

        logger.info('===== Remove all the not exact match features from the CrossMap output =====')
        # remove all the not exact match features from the CrossMap output
        filtered_file = '%s/%s_CrossMap_filtered%s' % (temp_dir, gff3_filename, gff3_extension)
        filter_not_exact_match(CrossMap_mapped_file, CrossMap_log_file, filtered_file, args.tmp_identifier)

        # re-construct the parent features
        logger.info('===== Re-construct the parent features for the models where all child features are perfectly re-mapped =====')
        re_construct_file = '%s/%s_re_construct%s' % (temp_dir, gff3_filename, gff3_extension)
        re_construct_report = '%s/%s_re_construct%s' % (temp_dir, gff3_filename, '.report')
        re_construct_gff3_features.main(in_gff, filtered_file, re_construct_file, re_construct_report, args.tmp_identifier)
        # run gff3_QC to generate QC report for re-constructed gff3 file
        logger.info('===== Run gff3_QC to generate QC report for re-constructed gff3 file =====')
        re_construct_QC = '%s/%s_re_construct_QC.report' % (temp_dir, gff3_filename)
        subprocess.Popen(['gff3_QC', '-g', re_construct_file, '-f', args.query_fasta, '-o', re_construct_QC]).wait()
        # list of errors
        # only fix the errors in the list
        fix_errors = ['Ema0001', 'Ema0003', 'Ema0006', 'Ema0007', 'Emr0001', 'Esf0014']
        re_construct_QC_filtered = '%s/%s_re_construct_QC_filtered.report' % (temp_dir, gff3_filename)
        log_file = open(re_construct_QC_filtered, 'w')
        pattern = 'NR==1 || /%s/ {print}'% ('/ || /'.join(fix_errors))
        subprocess.Popen(['awk', pattern, re_construct_QC], stdout=log_file).wait()
        log_file.close()
        rm_tmp_list.extend([CrossMap_mapped_file, CrossMap_mapped_file + '.unmap',CrossMap_log_file, filtered_file, re_construct_file, re_construct_report, re_construct_QC, re_construct_QC_filtered])
        # run gff3_fix
        update_gff = os.path.join(args.out_dir, '%s%s%s' % (gff3_filename, args.updated_postfix, gff3_extension))
        remove_gff = os.path.join(args.out_dir, '%s%s%s' % (gff3_filename, args.removed_postfix, gff3_extension))
        update_gff_QC = os.path.join(args.out_dir, '%s_QC%s' % (gff3_filename, '.txt'))
        logger.info('===== Run gff3_fix to correct GFF3 format errors =====')
        if args.tmp_identifier:
            tmp_update_gff = '%s/%s%s_tmp%s' % (temp_dir, os.path.basename(in_gff), args.updated_postfix, gff3_extension)
            subprocess.Popen(['gff3_fix', '-qc_r', re_construct_QC_filtered, '-g', re_construct_file, '-og', tmp_update_gff]).wait()
            get_remove_feature.output_remove_features(in_gff, tmp_update_gff, remove_gff, tmp_identifier)
            remove_tmpID(tmp_update_gff, update_gff)
            rm_tmp_list.append(tmp_update_gff)
        else:
            subprocess.Popen(['gff3_fix', '-qc_r', re_construct_QC_filtered, '-g', re_construct_file, '-og', update_gff]).wait()
            get_remove_feature.output_remove_features(in_gff, update_gff, remove_gff, tmp_identifier)
        logger.info('===== Get updated and removed GFF3 files =====')
        if args.summary_report:
            logger.info('===== Generate a summary report =====')
            summary_report = os.path.join(args.out_dir, '%s%s' % (gff3_filename, '_summary.tsv'))
            removed_list = os.path.join(args.out_dir, '%s%s' % (gff3_filename, '_removed.tsv'))
            get_summary(update_gff, remove_gff, summary_report, removed_list)
        # add ##genome-build source buildName to the second line
        genome_build = '##genome-build %s %s' % (args.source, args.buildName)
        add_string_to_file(update_gff, genome_build, 2)
        logger.info('===== Run gff3_QC to generate QC report for updated GFF3 files =====')
        subprocess.Popen(['gff3_QC', '-g', update_gff, '-f', args.query_fasta, '-o', update_gff_QC]).wait()
    if args.temp:
        for rmfile in rm_tmp_list:
            os.remove(rmfile)
        try:
            os.rmdir(temp_dir)
        except OSError:
            pass
