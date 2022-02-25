import sys
import os
import glob
import pandas as pd
sys.path.insert(0, "../Utils")
from utils import *
import re


def arg_parsing():
        parser = argparse.ArgumentParser(description='Pipeline for mQTLs')
        parser.add_argument('yaml_file', , type=str
                            help='File storing the parameters of the analysis')
        parser.add_argument('-s', '--steps' , type=str, default='all',
                            help='Step to perform analysis')
        args = parser.parse_args()
        print(str(time.asctime()))
        print('The arguments are: ', vars(args), '\n')
        sys.stderr.write('\n' + str(time.asctime())
                         + ' -- Arguments: ' + str(vars(args)) + '\n')
        return args


def correspondance_list_files(list1, list2):
    list1_base = set([os.path.basename(file).split('.')[0] for file in list1])
    list2_base = set([os.path.basename(file).split('.')[0] for file in list2])
    list1 = [file for file in list1 if os.path.basename(file).split('.')[0] in list1_base & list2_base]
    list2 = [file for file in list2 if os.path.basename(file).split('.')[0] in list1_base & list2_base]
    return list1, list2


def select_SNP_per_pvalue(file, pval_col, dist_bp=500000):
    snp_df = pd.read_table(file, usecols=['#CHR', 'POS', 'SNP', pval_col])
    best_snp = []
    for chr in snp_df['#CHR'].unique():
        chr_df = snp_df[snp_df['#CHR'] == chr]
        while chr_df.shape[0] > 0:
            best_pval = chr_df[pval_col].min()
            best_snp += chr_df[chr_df[pval_col]
                               == best_pval]['SNP'].unique().tolist()
            pos = chr_df[chr_df[pval_col]
                         == best_pval]['POS'].unique().tolist()[-1]
            chr_df = chr_df[(chr_df['POS'] < pos - dist_bp)
                            | (chr_df['POS'] > pos + dist_bp)]
    return best_snp


def find_mem_request(nano, base):
    from math import ceil
    nano_size = os.path.getsize(nano) / 1000000
    base_size = os.path.getsize(base)
    if (nano_size < 1) or (base_size < 1):
        mem = 5000
    else:
        mem = (ceil(nano_size / 1000) + 5) * 1000
    return mem


def main(yaml_file, steps='all'):

    os.system(f"sed -e 's/^#//' {yaml_file} > {yaml_file}_new")
    # Extracting the arguments in
    yml = yaml_parser(f'{yaml_file}_new')
    yml_ls = [f'{key}: {value}' for key,value in yml['other'].items()]
    print(f'# VARIABLES\n {"\n".join(yml_ls)}')

    # Creating the output drectory
    dir, n = yml['output_directory'], 0
    while os.path.exists(dir):
        n = n + 1
        dir = yml['output_directory'] + str(n)
    os.makedirs(dir)
    os.makedirs(os.path.join(dir, 'temp'))
    temp_dir = os.path.join(dir, 'temp')
    os.system(f"cp {yaml_file}_new {dir}")

    # Verification of the list of files bam and nano
    nano_files = directory_to_file_list(yml['nanopolish_directory'])
    basecal_files = directory_to_file_list(yml['basecalling_directory'])
    initial_len_ls = (len(nano_files), len(basecal_files))
    nano_files, basecal_files = correspondance_list_files(nano_files, basecal_files)
    print(f'\nNumber of corresponding nano-basecalling files: {len(nano_files)} (input: {initial_len_ls})')

    # Finemapping when required
    target_snp = yml['snp_type']
    if yml['snp_list']:
            assert os.path.exist(yml['snp_list']), 'Error with the file containing the snps'
    if yml['finemapping']:
        # TODO: Implement the selection of some SNP if necessary
        # TODO: verificationin cpg_analysis that we have only the SNPs we targetted
        snp_ls = select_SNP_per_pvalue(yml['finemapping']['file'],
                                       yml['finemapping']['pval_col'],
                                       dist_bp=yml['finemapping']['distance_bp'])

    # Defining the steps of the pipeline to run
    # TODO: Change the way the steps are managed
    if ('1' in steps) or ('all' in steps):
        # Bam-basecalling file genotyped/filtered/merged with corresponding nanopolish file
        os.chdir(temp_dir)
        for n, (base_call, nano) in enumerate(zip(nano_files, basecal_files)):
            print(n, basecall, nano)
        #     assert os.path.basename(base_call).split('.')[0] == os.path.basename(nano).split('.')[0], f'The files are not corresponding {base_call, nano}'
        #     mem = find_mem_request(nano, base_call)
        #     os.system(f'bsub -Jbamnano{n} -M{mem} -ebamnano{n}.out -obamnano{n}.out "python3 ../../bam_nano_filtering.py  {base_call} {nano} {target_snp}"')
        #
        #     # Open the output and rerun all the LSF memory errors
        #     rerun_more_mem = f'bsub -Jbamnano{n} -M{mem + 5000} -ebamnano{n}.out -obamnano{n}.out "python3 ../../bam_nano_filtering.py  {base_call} {nano} {target_snp}"'
        #     os.system(f'bsub -w"done(bamnano{n})" -Jbamnano{n} python3 quality_analysis.py bamnano{n}.out {rerun_more_mem}')
        #
        # # Merging all files together
        # os.system(f'bsub -w"done(bamnano0)" -Jhead -e /dev/null -o /dev/null "head -n1 Filtered_nano_bam_files_{os.path.basename(basecal_files[0]).split('.')[0]}.csv > ../Filtered_nano_bam_files.csv"')
        # os.system(f'bsub -w"done(bamnano*)" -Jmerge -emerge.out -omerge.out "tail -n+2 -q Filtered_nano_bam_files_* >> ../Filtered_nano_bam_files.csv"')
        # os.system(f' bsub -w"post_done(merge)" rm -f Filtered_nano_bam_files_*')

    if ('2' in steps) or ('all' in steps):
        print('STEP2')
        # Statistics
        # if '1' in steps:
        #     os.system('bsub -w"done(merge)" -Jstats -M {} "python3 cpg_snp_analysis.py   "')
        # else:
        #     os.system('bsub -Jstats -M {} "python3 cpg_snp_analysis.py   "')

    if ('3' in steps) or ('all' in steps):
        print('STEP3')
        # Plotting
        # if '3' in steps:
        #     os.system('bsub -w"done(stats)" -Jplots -M {} "python3 methylation_analysis.py   "')
        # else:
        #     os.system('bsub -Jplots -M {} "python3 methylation_analysis.py   "')

if __name__ == '__main__':
    main(**vars(arg_parsing()))
