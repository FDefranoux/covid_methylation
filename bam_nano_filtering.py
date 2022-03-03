import pandas as pd
import os
import argparse


################################################################################
############################  FUNCTIONS  #######################################
################################################################################

## For manipulation of the BAM_base_calling_files
def genotype_frombasecalling(file, target_snp, t=0.90, print_counts=False, filtering=True, **kwargs):
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
    if print_counts:
        print('\nGenome Counts', genome.to_markdown(), '\n')

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

    # Filtering
    if filtering:
        df = filtering_datas(df, **kwargs)
    df[['chromosome', 'pos', 'ref', 'alt']] = df[target_snp].str.split(':', expand=True)

    return df


def filtering_datas(df):
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

    log['final shape'] = str(new_df.shape)
    print(log)
    return new_df


## For manipulation of the nanopolish file
def grep_target_readnames(nano_file, list_readnames, control=True):
    print('GREP READNAMES STEP!!')
    out_file = os.path.basename(nano_file).split('.')[0]
    # Saving the read_names in temp file
    with open(f'{out_file}_readnames.temp', 'w') as f:
        f.writelines('\n'.join(list_readnames))
        f.writelines(['\n'])

    # Extraction of the lines from raw nanopolish file with corresponding readnames
    os.system(f'zcat {nano_file} | head -n 1 > {out_file}_greped.txt')
    os.system(f'zcat {nano_file} | grep -f {out_file}_readnames.temp >> {out_file}_greped.txt')
    if control:
        # TODO: Call or redo Tom's script to count the readnames per chromosome instead of total?
        os.system(
            f'grep -v -f {out_file}_readnames.temp {out_file}_greped.txt > {out_file}_notrecognized_readnames.txt')
    os.remove(f'{out_file}_readnames.temp')

    # Recuperating the lines contaning extra fields:
    os.system(f'cut -f12,13,14,15,16,17,18,19,20 {out_file}_greped.txt > {out_file}_extracols.txt')
    os.system(f'sed -i "/^[[:space:]]*$/d" {out_file}_extracols.txt') # del blank lines
    line_extracols = len(open(out_file + "_extracols.txt").readlines())
    if line_extracols != 0:
        print(f'\nThis file has too many fields ! {out_file}')
        # Verification that we will grep out not more than the extra-field lines
        os.system(f'grep -f {out_file}_extracols.txt {out_file}_greped.txt > {out_file}_greped_extracols')
        line_greped = len(open(out_file + "_greped_extracols").readlines())
        if line_extracols != line_greped:
            print(f'ERROR removed too much line corresponding to extra fields lines {line_extracols - line_greped}')
        os.system(f'grep -vf {out_file}_extracols.txt {out_file}_greped.txt > {out_file}_new.txt')
        os.remove(f'{out_file}_greped_extracols')
        os.system(f'mv {out_file}_new.txt {out_file}_greped.txt')

    os.remove(f'{out_file}_extracols.txt')
    return f'{out_file}_greped.txt'


def nanopolish_formatting(nanopolish_file, list_readnames, control=True):
    # Recuperation of the read names and grepping nanopolish file
    greped_nano = grep_target_readnames(nanopolish_file, list_readnames, control=control)

    # Reading and formatting nanopolish file
    nano = pd.read_table(greped_nano, dtype='object')
    nano['sample_id'] = os.path.basename(nanopolish_file).split('.')[0]
    nano.loc[nano['sample_id'].str.contains('PROM1'), 'phenotype'] = 'Severe'
    nano.loc[nano['sample_id'].str.contains('PROM1') == False, 'phenotype'] = 'Mild'
    nano['cpg'] = nano['chromosome'].astype(str) + ':' + nano['start'].astype(
        str) + ':' + nano['num_motifs'].astype(str)
    return nano


def merge_basefile_and_nanofile(basecalling_file, nanopolish_file, target_snp):

    # Opening base calling file and generating genotype (+filtering)
    called_base_df = genotype_frombasecalling(basecalling_file, target_snp, t=0.90, print_counts=False, filtering=True)

    # Selecting readnames of interest from nanopolish file+formatting the table
    nano_df = nanopolish_formatting(nanopolish_file, list(called_base_df['read_name'].unique()), control=True)

    # Merging nanopolish file with basecalling file
    merge = pd.merge(nano_df, called_base_df, on=['chromosome', 'read_name', 'sample_id'],
                     copy=False)
    merge.drop_duplicates(inplace=True)

    #TODO: verification of the merge results
    if called_base_df.shape[0] != merge.shape[0]:
        print('\n\n Shapes after merging',
              f'Base calling file: {called_base_df.shape}',
              f'nanopolish file: {nano_df.shape}',
              f'merged file: {merge.shape}', flush=True)

    return merge


################################################################################
###############################  MAIN  #########################################
################################################################################
def main(basecalling_file, nanopolish_file, target_snp):
    try:
        merge = merge_basefile_and_nanofile(basecalling_file, nanopolish_file, target_snp)
        # Saving individual files
        merge.to_csv(f'Filtered_nano_bam_files_{os.path.basename(basecalling_file).split(".")[0]}.csv',
                    mode='w', header=True, index=False)
    except Exception as err:
        print(basecalling_file, 'ERROR', err, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='STEP1 - Pipeline for mQTLs'
        + ' - Association of the basecalling file from BAM with nanopolish files')
    parser.add_argument('basecalling_file', type=str, help='basecalling file')
    parser.add_argument('nanopolish_file', type=str, help='nanopolish file')
    parser.add_argument('target_snp', type=str, help='which SNP (covid or control)')
    args = parser.parse_args()
    main(**vars(args))


# if not lsb:
#     for base_file in list_files:
#         print(base_file, flush=True)
#         try:
#             merge = merge_basefile_and_nanofile(base_file,
#                                                 target_snp, nanopolish_input, title, fine_mapping=False)
#             if base_file in list_files[0]:
#                 merge.to_csv(f'nanopolish_greped{title}/Filtered_nano_bam_files{title}.csv', mode='w',
#                              header=True, index=False)
#             else:
#                 merge.to_csv(f'nanopolish_greped{title}/Filtered_nano_bam_files{title}.csv', mode='a',
#                              header=False, index=False)
#         except Exception as err:
#             print(base_file, 'ERROR', err, flush=True)
# else:
#         # LSF job array management
#     file = lsf_arrray(list_files)
#     print(base_file, flush=True)
#     try:
#         merge = merge_basefile_and_nanofile(basecalling_file, nanopolish_file, target_snp)
#         # Saving individual files
#         merge.to_csv(f'Filtered_nano_bam_files_{os.path.basename(base_file).split('.')[0]}.csv',
#                     mode='w', header=True, index=False)
#     except Exception as err:
#         print(base_file, 'ERROR', err, flush=True)
