import pandas as pd
import os
import argparse
from math import ceil


def find_mem_request(nano, base):
    nano_size = os.path.getsize(nano) / 1000000
    base_size = os.path.getsize(base)
    print(nano, nano_size, base, base_size)
    if (nano_size < 1) or (base_size < 1):
        mem = 5000
    else:
        mem = (ceil(nano_size / 1000) + 5) * 1000
    return mem

################################################################################
############################  FUNCTIONS  #######################################
################################################################################

## For manipulation of the BAM_base_calling_files
def genotype_frombasecalling(file, target_snp, t=0.90, **kwargs):
    try:
        df = pd.read_table(file, dtype='object')
        df[['chromosome', 'pos', 'ref', 'alt']
           ] = df[target_snp].str.split(':', expand=True)

        # Association of the haplotype types
        df['haplotype'] = 'other'
        df.loc[df['base_called'] == df['ref'], 'haplotype'] = 'ref'
        df.loc[df['base_called'] == df['alt'], 'haplotype'] = 'alt'
        genome = df.groupby(['sample_id', target_snp,
                            'haplotype']).size().unstack(fill_value=0)

        # If one haplotype is not represented, set the whole column to zero.
        sum = genome.sum(axis=1)
        for h, haplo in enumerate(['ref', 'alt', 'other']):
            if df[df['haplotype'] == haplo].shape[0] == 0:
                genome[haplo] = 0
            # Transforming the values in percentage of total
            genome[str(h)] = genome[haplo].astype(int) / sum

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

        # Association of genotypes to the dataframe
        df = pd.merge(df, genome[['Allele1', 'Allele2', 'Genotype']].reset_index(),
                      on=['sample_id', target_snp])
        return df
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def filtering_datas(df, target_snp, min=5):
    log = {}
    log['initial nrows'] = df.shape[0]

    # Filtering out the deletions
    new_df = df[df['ref'].str.len() < 2].copy()
    log['deletions'], shape_diff = int(df.shape[0] - new_df.shape[0]), new_df.shape[0]

    # Removing cpgs with not enough reads
    n_reads = df.groupby(['sample_id', target_snp])['read_name'].nunique()
    low_reads_snps = n_reads[n_reads < min].dropna(how='all').reset_index()[target_snp]
    if len(low_reads_snps) > 0: new_df = new_df[new_df[target_snp].isin(low_reads_snps) == False]
    log['not enough reads'], shape_diff = int(shape_diff - new_df.shape[0]), new_df.shape[0]

    # Filtering out the miss-called alleless
    # (ie haplotype not corresponding nor to ref nor to alt)
    new_df = new_df[new_df['Genotype'].str.contains('2') == False]
    log['miss-called (genotype=2)'], shape_diff = int(shape_diff - new_df.shape[0]), new_df.shape[0]

    # Filtering out the reads_names associated with several CHR
    double_chr = new_df.groupby(
        ['read_name', 'chromosome']).size().unstack().dropna(thresh=2).reset_index()['read_name']
    if len(double_chr) > 0: new_df = new_df[new_df['read_name'].isin(double_chr) == False]
    log['read_name in several chr'], shape_diff = int(shape_diff - new_df.shape[0]), new_df.shape[0]

    log['final nrows'] = new_df.shape[0]
    print('\nFiltering log:\n', '\n'.join([key + ': ' + str(val) for key, val in log.items()]))
    return new_df



def nanopolish_formatting(nanopolish_file, name):
    # Reading and formatting nanopolish file
    nano = pd.read_table(nanopolish_file, dtype='object')
    if not nano.empty:
        nano['sample_id'] = name
        nano.loc[nano['sample_id'].str.contains('PROM1'), 'phenotype'] = 'Severe'
        nano.loc[nano['sample_id'].str.contains('PROM1') == False, 'phenotype'] = 'Mild'
        nano['cpg'] = nano['chromosome'].astype(str) + ':' + nano['start'].astype(
            str) + ':' + nano['num_motifs'].astype(str)
    return nano


def perc_methylation(df, thresh=0.5, n_min=5):
    df['methylation'] = None
    # df.loc[abs(df['log_lik_ratio'].astype(float)) <= thresh, 'methylation'] = 'unsure'
    df.loc[df['log_lik_ratio'].astype(float) > thresh, 'methylation'] = 'methylated'
    df.loc[df['log_lik_ratio'].astype(float) <= thresh, 'methylation'] = 'unmethylated'
    print(df['log_lik_ratio'].min(), df['log_lik_ratio'].median(), df['log_lik_ratio'].max() )
    print(df['methylation'].value_counts().to_markdown())
    print('blou')
    perc_meth = df.groupby(['cpg', 'sample_id', 'haplotype', 'methylation']).size().unstack()
    perc_meth['tot'] = perc_meth.sum(axis=1)
    perc_meth = perc_meth[perc_meth['tot'] >= n_min]
    perc_meth['%meth'] = perc_meth['methylated'] / perc_meth['tot']
    # perc_meth['%unmeth'] = perc_meth['unmethylated'] / perc_meth['tot']
    # perc_meth['%unsure'] = perc_meth['unsure'] / perc_meth['tot']
    return perc_meth


def merge_verifications(df, list_snp_file=None):
    # To see if there is any header left in the df
    print('Colnames in df: ', df[(df == df.columns)].dropna(how='all').shape[0], flush=True)
    df = df[(df == df.columns) == False].dropna(how='all')

    # Print genotype misscalling
    print('\nNumber of rows with wrong genotype: ', df[df['Genotype'].isin(['0/0', '0/1', '1/1']) == False].shape[0], flush=True)
    df = df[df['Genotype'].isin(['0/0', '0/1', '1/1'])]

    # Print SNPs outside our selection
    if list_snp_file:
        snp_ls = pd.read_table(list_snp_file, header=None)[0].tolist()
        print('\nNumber of rows with wrong SNPs: ', df[(df[unit].isin(snp_ls) == False)].shape[0], flush=True)
        df = df[(df[unit].isin(snp_ls))]
        del snp_ls

    # Print duplicated lines
    print('\nNumber of duplicated lines: ', df[df.duplicated(keep=False)].shape[0], flush=True)
    df = df[df.duplicated() == False]


################################################################################
###############################  MAIN  #########################################
################################################################################
def main(basecalling_file, nanopolish_file, target_snp, name='', ratio=0.90, min=5, thresh=0.5, list_snp_file=None, print_counts=False, filtering=True):

    # Opening base calling file and generating genotype
    called_base_df = genotype_frombasecalling(basecalling_file, target_snp, t=ratio)
    if not called_base_df.empty:

        # Filtering
        if filtering:
            called_base_df = filtering_datas(called_base_df, target_snp, min=min)

        # Count results
        if print_counts:
            print('\nGenome Counts\n', called_base_df['Genotype'].value_counts().to_markdown(), '\n')

        called_base_df[['chromosome', 'pos', 'ref', 'alt']] = called_base_df[target_snp].str.split(':', expand=True)

        # Selecting readnames of interest from nanopolish file+formatting the table
        nano_df = nanopolish_formatting(nanopolish_file, name)

        # Merging nanopolish file with basecalling file
        if not nano_df.empty:
            print(nano_df.head(2).to_markdown(), called_base_df.head(2).to_markdown())
            merge = pd.merge(nano_df, called_base_df, on=['chromosome', 'read_name', 'sample_id'], copy=False)
            called_base_df_shape = called_base_df.shape[0]
            nano_df_shape = nano_df.shape[0]
            del nano_df, called_base_df
            merge.drop_duplicates(inplace=True)
            assert merge.shape[0] > 0, 'Merging failed'

            # Calculation of the %methylation
            merge['log_lik_ratio'] = merge['log_lik_ratio'].astype(float)
            perc_meth = perc_methylation(merge, n_min=min, thresh=thresh)

            merge = merge.groupby(['phenotype', 'sample_id', 'chromosome', 'cpg',
                                    target_snp, 'Genotype', 'haplotype']).median().reset_index()
            merge = pd.merge(merge, perc_meth['%meth'].reset_index(), on=['cpg', 'sample_id', 'haplotype'], how='outer')
            del perc_meth
            merge_verifications(merge, list_snp_file=list_snp_file)
            print(merge.dtypes.to_markdown())

            # Saving individual files
            # QUESTION: As median or not ?
            merge.to_csv(f'Filtered_nano_bam_files_{name}.csv', mode='w', header=True, index=False)

            #TODO: verification of the merge results
            if called_base_df_shape != merge.shape[0]:
                print('\nShapes after merging\n',
                      f'Base calling file: {called_base_df_shape}\n',
                      f'nanopolish file: {nano_df_shape}\n',
                      f'merged file: {merge.shape}\n', flush=True)
        else:
            print('ERROR Nanopolish greped file is EMPTY!')
    else:
        print('ERROR Basecalling file is EMPTY!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='STEP1 - Pipeline for mQTLs'
        + ' - Association of the basecalling file from BAM with nanopolish files')
    parser.add_argument('basecalling_file', type=str, help='basecalling file')
    parser.add_argument('nanopolish_file', type=str, help='nanopolish file')
    parser.add_argument('-n', '--name', type=str, help='filename', default='')
    parser.add_argument('-s', '--target_snp', type=str, help='which SNP (covid or control)', default='covid_snp')
    parser.add_argument('-r', '--ratio', type=float, help='ratio for deciding genome', default=0.9)
    parser.add_argument('-m', '--min', type=float, help='minimal number of reads for validating a cpg', default=5)
    parser.add_argument('-t', '--thresh', type=float, help='threshold of log_methylation', default=0.5)
    parser.add_argument('-l', '--list_snp_file', type=str, help='File list of SNPs to find in the final df', default=None)
    parser.add_argument('-f', '--filtering', type=bool, help='Filter the datas ?', default=True)
    parser.add_argument('-c', '--print_counts', type=bool, help='Print counts of genotype', default=True)
    args = parser.parse_args()
    print('Arguments (bam_filtering.py): ', ' - '.join([str(x)+': '+str(y) for x,y in vars(args).items()]), '\n')
    main(**vars(args))
