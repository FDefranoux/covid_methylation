import pandas as pd
import os
import glob

# pd.options.display.max_rows = 2
# TODO: TEST LSF JOB indexing possibility !
# TODO: For each analysis, print the parameters used ! Report the logs ? The filtering steps ?

nanopolish_input = '/hps/nobackup/birney/projects/gel_methylation/nanopolish'
dir = '/hps/nobackup/birney/projects/gel_methylation/control_snps/reads/gcc0085*'
title = '_control_finemapped'
lsb = True


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
        f'zcat {nano_file} | head -n 1 > {os.path.join(output, file)}.txt')
    os.system(
        f'zcat {nano_file} | grep -f {os.path.join(output, file + "_readnames.temp")} >> {os.path.join(output, file + ".txt")}')
    if control:
        os.system(
            f'grep -f -v {os.path.join(output, file + "_readnames.temp")} {os.path.join(output, file + ".txt")} > {os.path.join(output, file + "_notrecognized_readnames.txt")}')
    # os.remove(f'{os.path.join(output, file + '_readnames.temp')}')


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
    log['deletions'] = int(df.shape[0] - new_df.shape[0])
    shape_diff = new_df.shape[0]

    # Filtering out the miss-called alleles
    # (ie haplotype not corresponding nor to ref nor to alt)
    new_df = new_df[new_df['Genotype'].str.contains('2') == False]
    log['miss-called'] = int(shape_diff - new_df.shape[0])
    shape_diff = new_df.shape[0]

    # Filtering out the reads_names associated with several CHR
    double_chr = new_df.groupby(
        ['read_name', 'chromosome']).size().unstack().dropna(thresh=2).index
    if len(double_chr) > 0:
        new_df.drop(double_chr, inplace=True)
    log['read_name in several chr'] = int(shape_diff - new_df.shape[0])
    shape_diff = new_df.shape[0]

    # Filtering for specific SNPs
    if fine_mapping:
        # TODO: Add kwargs to modify the finemapping
        list_snp = select_SNP_per_pvalue(fine_mapping,
                                         pval_col='all_inv_var_meta_p',
                                         dist_bp=500000)

    if list_snp:
        new_df = new_df[new_df['covid_snp'].isin(list_snp)]
        log['snp_filtering'] = int(shape_diff - new_df.shape[0])

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
    base_df = genotype_frombasecalling(base_df, t=0.90, print_counts=False)
    # Filtering datas
    base_df, log = filtering_datas(
        base_df, list_snp=list_snp, fine_mapping=fine_mapping)
    print(log)
    base_df = base_df.astype(object)

    # Reading and merging nanopolish file
    nano = pd.read_table(f'nanopolish_greped{title}/{base_file}.txt')
    nano['sample_id'] = base_file
    nano['chromosome'] = nano['chromosome'].astype(str)
    nano = nano.astype(object)
    merge = pd.merge(nano, base_df, on=['chromosome', 'read_name', 'sample_id'],
                     copy=False)

    #TODO: verification of the merge results
    if base_df.shape[0] != merge.shape[0]:
        print('\n\n Shapes after merging',
              f'Base calling file: {base_df.shape}',
              f'nanopolish file: {nano.shape}',
              f'merged file: {merge.shape}', flush=True)
    del nano, base_df

    nanopolish_formatting(merge)
    return merge


def main(dir, nanopolish_input, title='', file_snps='', lsb=False, fine_mapping=True):
    # colnames_basefile = ['sample_id', 'covid_snp', 'read_name', 'base_called']
    # snp = 'covid_snp'
    colnames_basefile = ['sample_id', 'control_snp', 'covid_snp',
                         'control_snp_rs', 'covid_snp_rs', 'read_name', 'base_called']
    target_snp = 'control_snp'
    if file_snps:
        list_snp = pd.read_table(file_snps, header=None)[0].to_list()
    else:
        list_snp = []
    os.system(f'mkdir nanopolish_greped{title}')
    if not lsb:
        for base_file in glob.glob(dir):
            print(base_file, flush=True)
            try:
                merge = merge_basefile_and_nanofile(base_file, colnames_basefile,
                                                    target_snp, nanopolish_input, list_snp, title)
                if base_file in glob.glob(dir)[0]:
                    merge.to_csv(f'Filtered_nano_bam_files{title}.csv', mode='w',
                                 header=True, index=False)
                else:
                    merge.to_csv(f'Filtered_nano_bam_files{title}.csv', mode='a',
                                 header=False, index=False)
            except Exception as err:
                print(base_file, 'ERROR', err, flush=True)
    else:
        # LSF job array management
        base_file = lsf_arrray(glob.glob(dir))
        print(base_file, flush=True)
        merge = merge_basefile_and_nanofile(base_file, colnames_basefile,
                                            target_snp, nanopolish_input, list_snp, title, fine_mapping=fine_mapping)
        # Saving
        merge.to_csv(f'Filtered_nano_bam_files{title}_{os.path.basename(base_file)[:-4]}.csv', mode='w',
                     header=True, index=False)
        print('\n\n------------------------\nTo gather the files run:')
        print('------------------------\n\n')
        print(
            f'head -n1 Filtered_nano_bam_files{title}_{os.path.basename(base_file)[:-4]}.csv > Filtered_nano_bam_files{title}.csv')
        print(
            f'cat Filtered_nano_bam_files{title}_* >> Filtered_nano_bam_files{title}.csv')


if __name__ == '__main__':
    main(dir=dir, title=title, file_snps='finemapped', lsb=lsb,
         nanopolish_input=nanopolish_input, fine_mapping=True)
