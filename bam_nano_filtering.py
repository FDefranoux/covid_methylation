import pandas as pd
import os
import glob

# pd.options.display.max_rows = 2
# TODO: For each analysis, print the parameters used ! Report the logs ? The filtering steps ?
# TODO: add arg parser + file to log the entries and analyses
# dir = '/hps/nobackup/birney/projects/gel_methylation/control_snps/reads/gcc*'
dir = 'redo_lsf.txt'
nanopolish_input = '/hps/nobackup/birney/projects/gel_methylation/nanopolish'
title = '_control_finemapped'
file_snps = 'finemapped'
lsb = True
target_snp = 'control_snp'
# fine_mapping = 'significant_hits_COVID19_HGI_A2_ALL_leave_23andme_20210607.txt'


def lsf_arrray(file_list):
    lsb_index = int(os.environ['LSB_JOBINDEX'])-1
    file_array = file_list[lsb_index]
    return file_array


def grep_target_readnames(file, list_readnames, nanopolish_input_dir, output='', control=True):
    new_file = os.path.join(output, file)
    with open(f'{new_file + "_readnames.temp")}', 'w') as f:
        f.writelines('\n'.join(list_readnames))
        f.writelines(['\n'])
    nano_file = os.path.join(nanopolish_input_dir, file) + '.tsv.gz'
    os.system(f'zcat {nano_file} | head -n 1 > {new_file)}_greped.txt')
    os.system(f'zcat {nano_file} | grep -f {new_file}_readnames.temp) >> {new_file}_greped.txt')
    # Removing the lines containg extra fields:
    os.system(f'cut -f12,13,14,15,16,17,18,19,20 {new_file}_greped.txt > {new_file} _extracols.txt')
    os.system(f'sed -i "/^[[:space:]]*$/d" {new_file}_extracols.txt')
    line_extracols = len(open(new_file + "_extracols.txt")).readlines())
        if line_extracols != 0:
            print('This file has too many fields !')
            os.system(f'grep -f {new_file}_extracols.txt {new_file}_greped.txt > {new_file}_greped_extracols)
            line_greped = len(open(new_file + "_greped_extracols")).readlines())
            if line_extracols != line_greped:
                print(f'ERROR removed too much line corresponding to extra fields lines {line_extracols - line_greped}')
            os.system(f'grep -vf {new_file}_extracols.txt {new_file}_greped.txt > {new_file}_new.txt')
            os.remove(f'{new_file}_greped.txt ')
            os.system(f'mv {new_file}_new.txt {new_file}_greped.txt')

    if control:
        # TODO: Call or redo Tom's script to count the readnames per chromosome instead of total?
        os.system(
            f'grep -v -f {new_file}_readnames.temp) {new_file}_greped.txt > {new_file}_notrecognized_readnames.txt')
    os.remove(f'{new_file}_readnames.temp')
    os.remove(f'{new_file}_extracols.txt')
    return f'{new_file}_greped.txt'


def genotype_frombasecalling(df, t=0.90, print_counts=False):
    df[['chromosome', 'pos', 'ref', 'alt']
       ] = df['control_snp'].str.split(':', expand=True)

    # Association of the haplotype types
    df['haplotype'] = 'other'
    df.loc[df['base_called'] == df['ref'], 'haplotype'] = 'ref'
    df.loc[df['base_called'] == df['alt'], 'haplotype'] = 'alt'
    genome = df.groupby(['sample_id', 'control_snp',
                        'haplotype']).size().unstack(fill_value=0)

    # If one haplotype is not represented, set the whole column to zero.
    for haplo in ['ref', 'alt', 'other']:
        if df[df['haplotype'] == haplo].shape[0] == 0:
            genome[haplo] = 0

    # Transforming the values in percentage of total
    sum = genome.sum(axis=1)
    genome['0'] = genome['ref'].astype(int) / sum
    genome['1'] = genome['alt'].astype(int) / sum
    genome['2'] = genome['other'].astype(int) / sum

    # Discard the unsignificant values
    genome = genome[genome > 1 - t]

    # Setting up the genome as the two predominant allele called
    genome['Allele1'] = genome[['0', '1', '2']].T.apply(lambda x: x.idxmax())
    for col in ['2', '1', '0']:
        genome.loc[(genome['Allele1'] == col), col] = -1
    genome['Allele2'] = genome[['0', '1', '2']].T.apply(lambda x: x.idxmax())
    alleles = genome[['Allele1', 'Allele2']].astype(int).copy()
    genome['Genotype'] = alleles.min(axis=1).astype(
                        str) + '/' + alleles.max(axis=1).astype(str)

    # Reset the percentage if printing
    if print_counts:
        genome['0'] = genome['ref'].astype(int) / sum
        genome['1'] = genome['alt'].astype(int) / sum
        genome['2'] = genome['other'].astype(int) / sum
        print(genome.to_markdown(), '\n')

    df = pd.merge(df, genome[['Allele1', 'Allele2', 'Genotype']].reset_index(),
                  on=['sample_id', 'control_snp'])
    return df


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


def filtering_datas(df, list_snp=[], fine_mapping=False):
    log = {}
    log['initial shape'] = str(df.shape)

    # Filtering out the deletions
    new_df = df[df['ref'].str.len() < 2].copy()
    log['deletions'], shape_diff = int(df.shape[0] - new_df.shape[0]), new_df.shape[0]

    # Filtering out the miss-called alleles
    # (ie haplotype not corresponding nor to ref nor to alt)
    new_df = new_df[new_df['Genotype'].str.contains('2') == False]
    log['miss-called'], shape_diff = int(shape_diff - new_df.shape[0]), new_df.shape[0]

    # Filtering out the reads_names associated with several CHR
    double_chr = new_df.groupby(
        ['read_name', 'chromosome']).size().unstack().dropna(thresh=2).index
    if len(double_chr) > 0: new_df.drop(double_chr, inplace=True)
    log['read_name in several chr'], shape_diff = int(shape_diff - new_df.shape[0]), new_df.shape[0]

    # Filtering for specific SNPs
    if fine_mapping:
        # TODO: Add kwargs to modify the finemapping or remove ??
        list_snp = select_SNP_per_pvalue(fine_mapping,
                                         pval_col='all_inv_var_meta_p',
                                         dist_bp=500000)

    if list_snp:
        new_df = new_df[new_df['covid_snp'].isin(list_snp)]
        log['snp_filtering'], shape_diff = int(shape_diff - new_df.shape[0]), new_df.shape[0]

    log['final shape'] = str(new_df.shape)
    return new_df, log


def nanopolish_formatting(nano):
    nano.loc[nano['sample_id'].str.contains('PROM1'), 'phenotype'] = 'Severe'
    nano.loc[nano['sample_id'].str.contains(
        'PROM1') == False, 'phenotype'] = 'Mild'
    nano['cpg'] = nano['chromosome'].astype(str) + ':' + nano['start'].astype(
        str) + ':' + nano['num_motifs'].astype(str)
    nano.loc[nano['num_motifs'] == 1, 'distance_cpg_snp'] = abs(
        nano['pos'].astype(int) - nano['start'].astype(int))


def merge_basefile_and_nanofile(file, colnames_basefile, target_snp,
                                nanopolish_input, list_snp=[], title='',
                                fine_mapping=True):
    # Opening base calling file
    base_df = pd.read_table(file, header=None, names=colnames_basefile, dtype='object')
    base_df[['chromosome', 'pos', 'ref', 'alt']
            ] = base_df[target_snp].str.split(':', expand=True)
    file = os.path.basename(file).split('.')[0]
    print(file, flush=True)

    base_df = genotype_frombasecalling(base_df, t=0.90, print_counts=False)
    # Filtering datas
    base_df, log = filtering_datas(base_df, list_snp=list_snp, fine_mapping=fine_mapping)
    print(log)

    # Recuperation of the read names and grepping nanopolish file
    greped_nano = grep_target_readnames(file, list(base_df['read_name'].unique()), nanopolish_input=nanopolish_input,
                          output=f'nanopolish_greped{title}', control=True)

    # Reading and formatting nanopolish file
    nano = pd.read_table(greped_nano, dtype='object')
    nano['sample_id'] = file
    nanopolish_formatting(nano)

    # Merging nanopolish file with basecalling file
    merge = pd.merge(nano, base_df, on=['chromosome', 'read_name', 'sample_id'],
                     copy=False)

    #TODO: verification of the merge results
    if base_df.shape[0] != merge.shape[0]:
        print('\n\n Shapes after merging',
              f'Base calling file: {base_df.shape}',
              f'nanopolish file: {nano.shape}',
              f'merged file: {merge.shape}', flush=True)
    del nano, base_df

    return merge


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
                merge = merge_basefile_and_nanofile(base_file, colnames_basefile,
                                                    target_snp, nanopolish_input, list_snp, title, fine_mapping=False)
                merge.drop_duplicates(inplace=True)
                if base_file in list_files[0]:
                    merge.to_csv(f'nanopolish_greped{title}/Filtered_nano_bam_files{title}.csv', mode='w',
                                 header=True, index=False)
                else:
                    merge.to_csv(f'nanopolish_greped{title}/Filtered_nano_bam_files{title}.csv', mode='a',
                                 header=False, index=False)
            except Exception as err:
                print(base_file, 'ERROR', err, flush=True)
    else:
        # LSF job array management
        base_file = lsf_arrray(list_files)
        print(base_file, flush=True)
        try:
            merge = merge_basefile_and_nanofile(base_file, colnames_basefile,
                                                target_snp, nanopolish_input, list_snp, title, fine_mapping=fine_mapping)
            merge.drop_duplicates(inplace=True)
            # Saving
            merge.to_csv(f'nanopolish_greped{title}/Filtered_nano_bam_files{title}_{os.path.basename(base_file)[:-4]}.csv', mode='w',
                         header=True, index=False)
        except Exception as err:
            print(base_file, 'ERROR', err, flush=True)

        # Ending message
        # TODO: Create small bash script to be called to perform those after completion of job array
        print('\n\n------------------------\nTo gather the files run:')
        print('------------------------\n')
        print(
            f'head -n1 nanopolish_greped{title}/Filtered_nano_bam_files{title}_{os.path.basename(base_file)[:-4]}.csv > nanopolish_greped{title}/Filtered_nano_bam_files{title}.csv')
        print(f'\ntail -n+2 -q Filtered_nano_bam_files{title}_* >> Filtered_nano_bam_files{title}.csv')
        print(
            f'\nrm -f nanopolish_greped{title}/Filtered_nano_bam_files{title}_*')

# TODO remove duplicated lines if any
if __name__ == '__main__':
    main(dir=dir, nanopolish_input=nanopolish_input, target_snp=target_snp,
         file_snps=file_snps, lsb=lsb, fine_mapping=False)
