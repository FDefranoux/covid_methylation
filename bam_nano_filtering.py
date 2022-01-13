import pandas as pd
import os
import glob

# pd.options.display.max_rows = 2

dir = '/hps/nobackup/birney/projects/gel_methylation/control_snps/reads/*'
title = '_control_finemapped'

def grep_target_readnames(file, list_readnames, output='nanopolish_greped'):
    with open(f'{file}_readnames.temp', 'w') as f:
        f.writelines('\n'.join(list_readnames))
        f.writelines(['\n'])
    nano_file = os.path.join('nanopolish_indexed', file + '.tsv.gz')
    os.system(f'zcat {nano_file} | head -n 1 > {file}.txt')
    os.system(f'mkdir {output}')
    os.system(f'zcat {nano_file} | grep -f {file}_readnames.temp >> {output}/{file}.txt')
    os.remove(f'{file}_readnames.temp')


def genotype_frombasecalling(df, t=0.90, print_counts=False):
    df[['chromosome', 'pos', 'ref', 'alt']] = df['control_snp'].str.split(':', expand=True)

    # Association of the haplotype types
    df['haplotype'] = 'other'
    df.loc[df['base_called'] == df['ref'], 'haplotype'] = 'ref'
    df.loc[df['base_called'] == df['alt'], 'haplotype'] = 'alt'
    genome = df.groupby(['sample_id', 'control_snp', 'haplotype']).size().unstack(fill_value=0)

    # If one haplotype is not represented, set the whole column to zero.
    for haplo in ['ref', 'alt', 'other']:
        if df[df['haplotype'] == haplo].shape[0] == 0:
            percent[haplo] = 0
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
        genome.loc[(genome['Allele1'] == col) , col] = -1
    genome['Allele2' ] = genome[['0', '1', '2']].T.apply(lambda x: x.idxmax())
    alleles = genome[['Allele1', 'Allele2']].astype(int).copy()
    genome['Genotype'] = alleles.min(axis=1).astype(
                        str) + '/' + alleles.max(axis=1).astype(str)

    # Reset the percentage if printing
    if print_counts:
        genome['0'] = genome['ref'].astype(int) / sum
        genome['1'] = genome['alt'].astype(int) / sum
        genome['2'] = genome['other'].astype(int) / sum
        print(genome.to_markdown(), '\n')

    df = pd.merge(df, genome[['Allele1', 'Allele2', 'Genotype']].reset_index(), on=['sample_id', 'control_snp'])
    return df


def filtering_datas(df, list_snp=[]):
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
    if len(double_chr) > 0: new_df.drop(double_chr, inplace=True)
    log['read_name in several chr'] = int(shape_diff - new_df.shape[0])
    shape_diff = new_df.shape[0]

    # Filtering for specific SNPs
    if list_snp:
        new_df = new_df[new_df['covid_snp'].isin(list_snp)]
        log['snp_filtering'] = int(shape_diff - new_df.shape[0])

    log['final shape'] = str(new_df.shape)
    return new_df, log


def nanopolish_formatting(nano):
    nano.loc[nano['sample_id'].str.contains('PROM1'), 'phenotype'] = 'Severe'
    nano.loc[nano['sample_id'].str.contains('PROM1') == False, 'phenotype'] = 'Mild'
    nano['cpg'] = nano['chromosome'].astype(str) + ':' + nano['start'].astype(
        str) + ':' + nano['num_motifs'].astype(str)
    nano.loc[nano['num_motifs'] == 1, 'distance_cpg_snp'] = abs(
        nano['pos'].astype(int) - nano['start'].astype(int))


def main(dir, title=''):
    # colnames_basefile = ['sample_id', 'covid_snp', 'read_name', 'base_called']
    # snp = 'covid_snp'
    colnames_basefile = ['sample_id', 'control_snp', 'covid_snp', 'control_snp_rs', 'covid_snp_rs', 'read_name', 'base_called']
    target_snp = 'control_snp'
    nano_cols = ['chromosome', 'strand', 'start', 'end', 'read_name',
                     'log_lik_ratio', 'log_lik_methylated', 'log_lik_unmethylated',
                     'num_calling_strands', 'num_motifs', 'sequence']

    for file in glob.glob(dir):
        base_file = pd.read_table(file, header=None, names=colnames_basefile)
        base_file[['chromosome', 'pos', 'ref', 'alt']] = base_file[target_snp].str.split(':', expand=True)
        list_readnames = list(base_file['read_name'].unique())
        file = os.path.basename(file).split('.')[0]
        # grep_target_readnames(file, list_readnames)
        base_file = genotype_frombasecalling(base_file, t=0.90, print_counts=False)
        finemapped_snp = pd.read_table('finemapped', header=None)[0].to_list()
        base_file, log = filtering_datas(base_file, list_snp=finemapped_snp)
        # TODO: Print something with the logs ?
        base_file = base_file.astype(object)

        nano = pd.read_table(f'nanopolish_greped/{file}.txt', header=None, names=nano_cols)
        nano['sample_id'] = file
        nano['chromosome'] = nano['chromosome'].astype(str)
        nano = nano.astype(object)
        merge = pd.merge(nano, base_file, on=['chromosome', 'read_name', 'sample_id'], copy=False)
        if base_file.shape[0] != merge.shape[0]:
            print('\n\n Shapes after merging',
                  f'Base calling file: {base_file.shape}',
                  f'nanopolish file: {nano.shape}',
                  f'merged file: {merge.shape}', flush=True)
        del nano, base_file

        nanopolish_formatting(merge)

        if file in glob.glob(dir)[0]:
            merge.to_csv(f'Filtered_nano_bam_files{title}.csv', mode='w',
                         header=True, index=False)
        else:
            merge.to_csv(f'Filtered_nano_bam_files{title}.csv', mode='a',
                         header=False, index=False)


if __name__ == '__main__':
    main(dir, title=title)
