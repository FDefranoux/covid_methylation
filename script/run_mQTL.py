import sys
import os
import socket
import argparse
import pandas as pd
host = socket.gethostname()
if 'Fanny' in host:
    PATH_UTILS = '/home/fanny/Work/EBI/Utils'
    ABS_PATH = '/home/fanny/Work/EBI/covid_nanopore'
else:
    PATH_UTILS = '/nfs/research/birney/users/fanny/Utils'
    ABS_PATH = '/nfs/research/birney/users/fanny/covid_nanopore'

sys.path.insert(0, PATH_UTILS)
from utils import *
import time
import re
from math import ceil
import subprocess



def arg_parsing():
    parser = argparse.ArgumentParser(description='Pipeline for mQTLs')
    parser.add_argument('yaml_file', type=str,
                        help='File storing the parameters of the analysis')
    parser.add_argument('-s', '--steps' , type=str, default='all',
                        help='Step to perform analysis')
    args = parser.parse_args()
    print(str(time.asctime()))
    print('The arguments are: ', vars(args), '\n', flush=True)
    return args


def correspondance_list_files(list1, list2, details=True):
    list1_base = {os.path.basename(file).split('.')[0] for file in list1}
    list2_base = {os.path.basename(file).split('.')[0] for file in list2}
    list1 = [file for file in list1 if os.path.basename(file).split('.')[0] in
        (list1_base & list2_base)]
    list2 = [file for file in list2 if os.path.basename(file).split('.')[0] in
        (list1_base & list2_base)]
    if details:
        print('Non corrresponding files: ', (list1_base - list2_base) | (list2_base - list1_base))
    return sorted(list1), sorted(list2)


def select_SNP_per_pvalue(file, pval_col, dist_bp=500000):
    snp_df = pd.read_table(file, usecols=['#CHR', 'POS', 'SNP', pval_col])
    best_snp = []
    for chrom in snp_df['#CHR'].unique():
        chr_df = snp_df[snp_df['#CHR'] == chrom]
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
    nano_size = os.path.getsize(nano) / 1000000
    base_size = os.path.getsize(base)
    print(nano, nano_size, base, base_size)
    if (nano_size < 1) or (base_size < 1):
        mem = 5000
    else:
        mem = (ceil(nano_size / 1000) + 5) * 1000
    return mem


def main(yaml_file, steps='all'):

    os.system(f"sed -e 's/^# //' {yaml_file} > {yaml_file[:-4]}_new.yml")
    # Extracting the arguments in
    yml = yaml_parser(f'{yaml_file[:-4]}_new.yml')
    yml = yml['other']
    yml_ls = [f'{key}: {value}' for key,value in yml.items()]
    print(f'# VARIABLES\n {"   ---    ".join(yml_ls)}', flush=True)

    # Creating the output drectory
    out_dir, n = yml['output_directory'], 0
    while os.path.exists(out_dir):
        n = '-'.join([str(n) for n in time.localtime(time.time())[0:5]])
        out_dir = yml['output_directory'] + n
    os.makedirs(out_dir)
    subdir_ls = ['per_file_filtering', 'per_chr_analysis', 'pval_plots', 'individual_plots']
    for n, dir in enumerate(subdir_ls):
        subdir_ls[n] = os.path.join(out_dir, dir)
        os.makedirs(subdir_ls[n])
    os.system(f"mv {yaml_file[:-4]}_new.yml {out_dir}")

    # Verification of the list of files bam and nano
    nano_files = directory_to_file_list(yml['nanopolish_directory'])
    basecal_files = directory_to_file_list(yml['basecalling_directory'])
    nano_files, basecal_files = correspondance_list_files(nano_files, basecal_files, details=True)
    print(f'\nNumber of corresponding nano-basecalling files: {len(nano_files)}')

    # Finemapping when required
    target_snp = yml['snp_type']
    # if yml['snp_list']:
    #     assert os.path.exist(yml['snp_list']), 'Error with the file containing the snps'
    # if yml['finemapping']:
    #     # TODO: Implement the selection of some SNP if necessary
    #     # TODO: verificationin cpg_analysis that we have only the SNPs we targetted
    #     snp_ls = select_SNP_per_pvalue(yml['finemapping']['file'],
    #                                    yml['finemapping']['pval_col'],
    #                                    dist_bp=yml['finemapping']['distance_bp'])

    # Defining the steps of the pipeline to run
    # TODO: Change the way the steps are managed
    if ('1' in steps) or ('all' in steps):
        print('# STEP1 -- BAM AND NANO FILES FILTERING')
        # Bam-basecalling file genotyped/filtered/merged with corresponding nanopolish file
        os.chdir(subdir_ls[0])
        for n, (nano, base_call) in enumerate(zip(nano_files, basecal_files)):
            assert os.path.basename(base_call).split('.')[0] == os.path.basename(nano).split('.')[0], f'The files are not corresponding {base_call, nano}'
            file_title = os.path.basename(base_call).split('.')[0]
            mem = find_mem_request(nano=nano, base=base_call)
            print(n, base_call, nano, mem, flush=True)
            os.system(f'bsub -Jbamnano{n} -M{mem} -obamnano_{file_title}.out "python3 {ABS_PATH}/bam_nano_filtering.py  {base_call} {nano} {target_snp}"')

        # Merging all files together
        first = os.path.basename(basecal_files[0]).split('.')[0]
        os.system(f'bsub -w"done(bamnano0)" -Jhead -o /dev/null "head -n1 {subdir_ls[0]}/Filtered_nano_bam_files_{first}.csv > {out_dir}/Filtered_nano_bam_files.csv"')
        os.system(f'bsub -w"ended(bamnano*)" -Jverif -o /dev/null "grep LSBATCH bamnano_*.out -A5 | grep -v LSBATCH | grep -v python3 > {out_dir}/SUMMARY_filtering.out"')
        os.system(f'bsub -w"done(verif)" -M1000 -Jmerge -o merge.out "tail -n+2 -q {subdir_ls[0]}/Filtered_nano_bam_files_* >> {out_dir}/Filtered_nano_bam_files.csv"')
        os.system(f'bsub -w"done(merge)" -M1000 -Jrm "rm -f {subdir_ls[0]}/Filtered_nano_bam_files_*"')

    # if ('2' in steps) or ('all' in steps):
    #     list_snp_file = yml['list_snp_file']
    #     print('# STEP2 -- STATISTICAL ANALYSIS')
    #     os.chdir(subdir_ls[1])
    #     # TODO: Reorganize the output files
    #     # First we need to split the file into chromosomes
    #     os.system(f'bsub -w"done(merge)" -M 1000 -Jsplit -osplit.out "bash {ABS_PATH}/split_filter.sh {out_dir}/Filtered_nano_bam_files.csv {subdir_ls[1]}"')
    #     for chr in range(1,25):
    #         os.system(f'bsub -w"done(split)" -M 8000 -Jcpg_{chr} -o{subdir_ls[1]}/cpg_{chr}.out "python3 {ABS_PATH}/cpg_snp_analysis.py {subdir_ls[1]}/Filtered_nano_bam_files.csv_chr_{chr}.csv -o {subdir_ls[1]} -u {target_snp} -f {list_snp_file}"')
    #     # TODO implement the merging of the analysis files
    #     # TODO manage if file is empty
    #     os.system(f'bsub -w"ended(cpg*)" -Jverifcpg -o/dev/null "grep LSBATCH {subdir_ls[1]}/cpg_*.out -A5 | grep -v LSBATCH | grep -v python3 > {out_dir}/SUMMARY_analysis.out"')
    #     for file in ['Counts_Diff_means.csv', 'Mann_Whitney.csv', 'Spearmann_corr.csv' ]:
    #         os.system(f'bsub -w"done(verifcpg)" -Jhead{file} -o/dev/null "head -n1 {subdir_ls[1]}/Filtered_nano_bam_files.csv_chr_1_{file} > {out_dir}/All_{file}"')
    #         os.system(f'bsub -w"done(head{file})" -M1000 -Jmerge{file} -o/dev/null "tail -n+2 -q {subdir_ls[1]}/Filtered_nano_bam_files.csv_chr_*_{file} >> {out_dir}/All_{file}"')
    #         # TODO: generalize following line
    #         # os.system("bsub -M 1000 -o /dev/null 'cut -d, -f3,14,15,17 Filtered_nano_bam_files.csv | sort -u > snp_associated_cpg.csv'")

    # if ('3' in steps) or ('all' in steps):
    #     print('# STEP3 -- PLOTTING')
    #     file_ls = ['All_Mann_Whitney.csv', 'All_Spearmann_corr.csv']
    #     file_ls = [os.path.join(out_dir, file) for file in file_ls]
    #     files = '--'.join(file_ls)
    #     ABS_PATH + '/methylation_analysis.py'
    #
    #     os.system(f'bsub -Jplots -M10000  "python3 /nfs/research/birney/users/fanny/covid_nanopore/methylation_analysis.py All_Mann_Whitney.csv--All_Spearmann_corr.csv -s control_snp -p 6.5"')
    #
    #
    #     os.system('bsub -w"done(merge*)" -Jplots -M {} "python3 methylation_analysis.py" ')


if __name__ == '__main__':
    main(**vars(arg_parsing()))
