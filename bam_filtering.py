import pandas as pd
import os
import argparse


################################################################################
############################  FUNCTIONS  #######################################
################################################################################

## For manipulation of the BAM_base_calling_files
def genotype_frombasecalling(file, target_snp, t=0.90, **kwargs):
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


# TODO: Include function to calculate the %methylation


def nanopolish_formatting(nanopolish_file):

    # Reading and formatting nanopolish file
    nano = pd.read_table(nanopolish_file, dtype='object')
    nano.loc[nano['sample_id'].str.contains('PROM1'), 'phenotype'] = 'Severe'
    nano.loc[nano['sample_id'].str.contains('PROM1') == False, 'phenotype'] = 'Mild'
    nano['cpg'] = nano['chromosome'].astype(str) + ':' + nano['start'].astype(
        str) + ':' + nano['num_motifs'].astype(str)
    return nano


################################################################################
###############################  MAIN  #########################################
################################################################################
def main(basecalling_file, nanopolish_file, target_snp, name='', ratio=0.90, min=5, print_counts=False, filtering=True):

    # Opening base calling file and generating genotype
    called_base_df = genotype_frombasecalling(basecalling_file, target_snp, t=ratio)

    # Filtering
    if filtering:
        called_base_df = filtering_datas(called_base_df, target_snp, min=min)

    # Count results
    if print_counts:
        print('\nGenome Counts\n', called_base_df['Genotype'].value_counts().to_markdown(), '\n')

    called_base_df[['chromosome', 'pos', 'ref', 'alt']] = called_base_df[target_snp].str.split(':', expand=True)

    # Selecting readnames of interest from nanopolish file+formatting the table
    nano_df = nanopolish_formatting(nanopolish_file)

    # Merging nanopolish file with basecalling file

    nano_df['sample_id'] = name
    merge = pd.merge(nano_df, called_base_df, on=['chromosome', 'read_name', 'sample_id'],
                     copy=False)
    merge.drop_duplicates(inplace=True)
    # Saving individual files
    # QUESTION: As median or not ?
    merge.to_csv(f'Filtered_nano_bam_files_{name}.csv', mode='w', header=True, index=False)

    #TODO: verification of the merge results
    if called_base_df.shape[0] != merge.shape[0]:
        print('\nShapes after merging\n',
              f'Base calling file: {called_base_df.shape}\n',
              f'nanopolish file: {nano_df.shape}\n',
              f'merged file: {merge.shape}\n', flush=True)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='STEP1 - Pipeline for mQTLs'
        + ' - Association of the basecalling file from BAM with nanopolish files')
    parser.add_argument('basecalling_file', type=str, help='basecalling file')
    parser.add_argument('nanopolish_file', type=str, help='nanopolish file')
    parser.add_argument('-n', '--name', type=str, help='filename', default='')
    parser.add_argument('-s', '--target_snp', type=str, help='which SNP (covid or control)', default='covid_snp')
    parser.add_argument('-t', '--ratio', type=float, help='ratio for deciding genome', default=0.9)
    parser.add_argument('-m', '--min', type=float, help='minimal number of reads for validating a cpg', default=5)
    parser.add_argument('-f', '--filtering', type=bool, help='Filter the datas ?', default=True)
    parser.add_argument('-c', '--print_counts', type=bool, help='Print counts of genotype', default=True)
    args = parser.parse_args()
    print('## Bam filtering. \nThe arguments are: ', vars(args), '\n')
    main(**vars(args))
