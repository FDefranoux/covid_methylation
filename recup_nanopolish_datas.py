import pandas as pd
import pysam
import os
import sys

# TODO: duplication of rows SOLVED ??
# TODO: SNPs that are in the range but do not appear as row SOLVED ??


# Variable definition
file_list = 'all_files.txt'
hits_table = 'significant_hits_COVID19_HGI_A2_ALL_leave_23andme_20210607.txt'
output_name = 'Sign_all.csv'
# output_name = 'SNPs_across_nanopolish_files.csv'


def region_select_fromSNP(pos_df):
    pos_df = pos_df.astype(int).copy()
    pos_df['new_POS'] = (pos_df['POS'] - 11).astype(str).str[:-6]
    pos_df['POS_end'] = (pos_df['new_POS'].astype(int) + 1) * 1000000
    pos_df['new_POS'] = pos_df['new_POS'] + '000000'
    pos_df['#CHR'] = pos_df['#CHR'].astype(str)
    pos_df[['new_POS', 'POS', 'POS_end']] = pos_df[[
        'new_POS', 'POS', 'POS_end']].astype(int)
    pos_df.set_index('POS', inplace=True)
    pos_df['tuple'] = pos_df.set_index(
        ['#CHR', 'new_POS', 'POS_end']).index.tolist()
    return pos_df['tuple']


def tabix_listregion_listfiles(list_region, file, error_files, cols=[]):
    df_final = pd.DataFrame(columns=cols)
    error_files[file] = []
    try:
        tbx = pysam.TabixFile(file)
        for n, region_n in enumerate(list_region):
            try:
                rows_n = [x for x in tbx.fetch(*region_n)]
                row_array = [r.split('\t') for r in rows_n]
                df_reg = pd.DataFrame(row_array, columns=cols)
                df_reg['region'] = str(region_n)
                df_final = pd.concat([df_final, df_reg], axis=0)

            except Exception as err:
                error_files[file].append(f'Error with {region_n}\n{err}')

    except OSError:
        error_files[file] = f'{file}.tbi not found\n'
    return df_final[df_final.duplicated() == False]


def associate_snp_to_sequence(df_final, snps):
    df = df_final[df_final.duplicated() == False].copy()
    # If working with regions, how to associate SNP value to region asked
    df['start_sequence'] = df['start'].astype(int) - 5
    df['end_sequence'] = df['start_sequence'].astype(
        int) + df['sequence'].str.len()
    df_out = pd.DataFrame()
    for snp in set(snps):
        df_int = df[(df['start_sequence'] <= snp)
                    & (df['end_sequence'] >= snp)].copy()
        df_int['SNP'] = snp
        df_out = pd.concat([df_out, df_int], axis=0)
    return df_out


def associate_snp_alleles(df_final, snps_refs):
    df = df_final[df_final.duplicated() == False].copy()
    df['rel_snp'] = df['SNP'] - df['start_sequence']
    df['allele'] = '0'
    df['#CHR'] = df['#CHR'].astype(int)

    for n in df['rel_snp'].unique():
        df.loc[df['rel_snp'] == n, 'allele'] = df['sequence'].str[n]
    allele_df = pd.merge(snps_refs, df, right_on=[
                         'SNP', '#CHR'], left_on=['POS', '#CHR'])

    # Genotype association
    allele_df.loc[allele_df['REF'] == allele_df['allele'], 'Genotype'] = '0'
    allele_df.loc[allele_df['ALT'] == allele_df['allele'], 'Genotype'] = '1'
    allele_df.loc[(allele_df['ALT'] != allele_df['allele']) & (
        allele_df['REF'] != allele_df['allele']), 'Genotype'] = '2'
    return allele_df


def main(file_list, hits_table):
    # Remove existing file
    if output_name in os.listdir():
        print('Removing previous file')
        os.remove(output_name)

    # Reading the tables
    file_ls = pd.read_table(file_list, header=0).iloc[:, 0].tolist()
    hits_df = pd.read_table(hits_table)

    # Find a way to select columns directly from the file
    cols = ['#CHR', 'strand', 'start', 'end', 'read_name', 'log_lik_ratio',
            'log_lik_methylated', 'log_lik_unmethylated', 'num_calling_strands',
            'num_motifs', 'sequence']
    list_region = set(region_select_fromSNP(hits_df[['#CHR', 'POS']]).values)
    error_files = {}
    for file in file_ls:
        error_files[file] = []
        # print(pd.read_table(file, sep='\t', nrows=1).columns.tolist() == cols)
        df_final = tabix_listregion_listfiles(
            list_region, file, error_files, cols=cols)
        df_final['#CHR'] = df_final['#CHR'].astype(int)

        if not df_final.empty:
            df_final['File'] = file
            df_snps = associate_snp_to_sequence(df_final[['sequence', 'start',
                                                          '#CHR']],
                                                hits_df['POS'])
            df_final = pd.merge(df_final, df_snps, on=['sequence', 'start',
                                                       '#CHR'])

            snp_genotype = associate_snp_alleles(df_final[['start_sequence',
                                                           'SNP', 'sequence',
                                                           '#CHR']],
                                                 hits_df[['#CHR', 'POS',
                                                          'REF', 'ALT']])
            df_final = pd.merge(df_final, snp_genotype, on=['#CHR',
                                                            'start_sequence',
                                                            'sequence', 'SNP'])
            if 'PROM1' in file:
                df_final['Phenotype'] = 'Severe'
            else:
                df_final['Phenotype'] = 'Mild'
            pos_snp_diff = df_final[df_final['POS'] != df_final['SNP']]
            if not pos_snp_diff.empty:
                print('ERRRROR, SNP and POS columns not equals')
                print(pos_snp_diff)
                # df_final.drop(pos_snp_diff.index, axis=0)
            del df_final['POS']
            if file == file_ls[0]:
                df_final.to_csv(output_name, mode='a',
                                header=True, index=False)
            else:
                df_final.to_csv(output_name, mode='a',
                                header=False, index=False)
        elif error_files[file] == []:
            error_files = 'No DFs'
    print('DICT OF ERRORS', file=sys.stderr)
    print([(file, err) for (file, err) in zip(
        error_files.keys(), error_files.values()) if err != []],
         file=sys.stderr)


if __name__ == '__main__':
    main(file_list, hits_table)
