import pandas as pd
import os
import glob

# pd.options.display.max_rows = 2
# TODO: For each analysis, print the parameters used ! Report the logs ? The filtering steps ?
# TODO: add arg parser + file to log the entries and analyses
dir = '/hps/nobackup/birney/projects/gel_methylation/control_snps/reads/gcc*'
# dir = 'example.txt'
nanopolish_input = '/hps/nobackup/birney/projects/gel_methylation/nanopolish'
title = '_grep_step_control'
file_snps = 'finemapped'
lsb = True
target_snp = 'control_snp'
# fine_mapping = 'significant_hits_COVID19_HGI_A2_ALL_leave_23andme_20210607.txt'


def lsf_arrray(file_list):
    lsb_index = int(os.environ['LSB_JOBINDEX'])-1
    file_array = file_list[lsb_index]
    return file_array


def grep_target_readnames(file, list_readnames, nanopolish_input, output='', control=True):
    with open(f'{os.path.join(output, file + "_readnames.temp")}', 'w') as f:
        f.writelines('\n'.join(list_readnames))
        f.writelines(['\n'])
    nano_file = os.path.join(nanopolish_input, file + '.tsv.gz')
    os.system(
        f'zcat {nano_file} | grep -f {os.path.join(output, file + "_readnames.temp")} >> {os.path.join(output, file + ".txt")}')
    if control:
        os.system(f'cut -f5 {os.path.join(output, file + ".txt")} | sort -u > {os.path.join(output, file + "_recognized_readnames.temp")}')
        os.system(
            f'grep -v -f {os.path.join(output, file + "_recognized_readnames.temp")} {os.path.join(output, file + "_readnames.temp")} > {os.path.join(output, file + "_notrecognized_readnames.txt")}')
    # os.remove(f'{os.path.join(output, file + "_readnames.temp")}')
    # os.remove(f'{os.path.join(output, file + "_recognized_readnames.temp")}')


def merge_basefile_and_nanofile(base_file, colnames_basefile, target_snp,
                                 nanopolish_input, list_snp=[], title='',
                                 fine_mapping=True):
    # Opening base calling file
    base_df = pd.read_table(base_file, header=None, names=colnames_basefile)
    base_df[['chromosome', 'pos', 'ref', 'alt']
            ] = base_df[target_snp].str.split(':', expand=True)
    list_readnames = list(base_df['read_name'].unique())
    base_file = os.path.basename(base_file).split('.')[0]
    print(base_file, flush=True)
    # Recuperation of the read names and grepping nanopolish file
    grep_target_readnames(base_file, list_readnames, nanopolish_input=nanopolish_input,
                          output=f'nanopolish_greped{title}', control=True)


def main(dir, nanopolish_input, target_snp, file_snps='', lsb=False, fine_mapping=False, title=''):
    # TODO: Check wheter if the basefile for covid and control snp could be the same or associated ...
    # if target_snp = 'covid_snp'
    # colnames_basefile = ['sample_id', 'covid_snp', 'read_name', 'base_called']
    if (os.path.isdir(dir)) or ('*' in dir):
        list_files = glob.glob(dir)
    else:
        try:
            list_files = pd.read_table(dir, header=None)[0].tolist()
        except Exception as err:
            print('Could not read input dir or file list\n', err)

    colnames_basefile = ['sample_id', 'control_snp', 'covid_snp',
                         'control_snp_rs', 'covid_snp_rs', 'read_name', 'base_called']

    if file_snps:
        list_snp = pd.read_table(file_snps, header=None)[0].to_list()
    else:
        list_snp = []
    os.system(f'mkdir nanopolish_greped{title}')
    if not lsb:
        for base_file in list_files:
            print(base_file, flush=True)
            try:
                merge_basefile_and_nanofile(base_file, colnames_basefile,
                                            target_snp, nanopolish_input,
                                            list_snp, title, fine_mapping=fine_mapping)

            except Exception as err:
                print(base_file, 'ERROR', err, flush=True)
    else:
        # LSF job array management
        base_file = lsf_arrray(list_files)
        print(base_file, flush=True)
        try:
            merge_basefile_and_nanofile(base_file, colnames_basefile,
                                        target_snp, nanopolish_input,
                                        list_snp, title, fine_mapping=fine_mapping)

        except Exception as err:
            print(base_file, 'ERROR', err, flush=True)


if __name__ == '__main__':
    main(dir=dir, nanopolish_input=nanopolish_input, target_snp = target_snp,
         file_snps=file_snps, lsb=lsb, fine_mapping=False)


def new_function(basefile, snp='', finemapped='', outdir='.'):
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    table = pd.read_table(basefile)
    print(table.head())
    if finemapped:
        fine_ls = pd.read_table(finemapped, header=None)[0].to_list()
    print(fine_ls)
    for file in table['sample_id'].unique():
        if finemapped:
            table[(table['sample_id'] == file) & (table[snp].isin(fine_ls))].to_csv(os.path.join(outdir, f'{file}.txt'), index=False, sep='\t')
        else:
            table[table['sample_id'] == file].to_csv(os.path.join(outdir, f'{file}.txt'), index=False, sep='\t')





import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
before = pd.read_table('finemapped_readchr_count.txt', header=None, names=['file', 'before'])
before[['file', 'chr']] = before['file'].str[2:-2].str.split("', '", expand=True)
after = pd.read_table('finemapped_readchr_count_aftergrep.txt', header=0, names=['chr', 'after']).reset_index()
after['file'] = after['index'].str[:-4].str.split('/', expand=True)[1]
after.drop('index', axis=1, inplace=True)
after = after[after['chr'].isin([str(n) for n in range(30)] + ['X', 'Y'])]

n_bef = before.groupby('file').nunique()['chr'].sort_values()
n_aft = after.groupby('file').nunique()['chr'].sort_values()

bef_valid = n_bef[n_bef >= 22].index
aft_valid = n_aft[n_aft >= 22].index

len(set(bef_valid) - set(aft_valid))
len(set(aft_valid) & set(bef_valid))

err = pd.read_table('error_files_filtering.txt', header=None)[0].tolist()
err = [ l[:-4] for l in err]
wcoln = pd.read_table('wrong_col_numb_filtering.txt', header=None)[0].tolist()
wcoln = [ l[:-4] for l in wcoln]

all = pd.merge(after, before, on=['chr', 'file'])
all['err'] = 'None'
all.loc[all['file'].isin(err), 'err'] = 'other_err'
all.loc[all['file'].isin(wcoln), 'err'] = 'wrong_col_numb'


# NB if we wanto to control the quality of the analysis we need to control the bam files, the nanopolish (before and after grep) and their filtered
