import pandas as pd
import pysam
import os
import sys

# Variable definition
file_list = ['bams/gcc13846.bam']
print(file_list)
hits_table = 'significant_hits_COVID19_HGI_A2_ALL_leave_23andme_20210607.txt'
output_name = 'Bam_test.csv'

# TODO: Ajust the generalization and structure of the code
#   - Function for adding the relevant information in the summary
#   - Gathering function to create list of region and reading of hits_df
#       --> relevant info gathered in analysis script
# TODO: Verification that all RNAMEs are in the header
# TODO: Adjust function to have larger region selection


def region_select_fromSNP_large(pos_df):
    # SEE how to implement this correctly and usefully
    # TODO: Assignement of SNP for each region
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
    return pos_df[['tuple', 'POS']].set_index('tuple')


def region_select_fromSNP(pos_df):
    pos_df = pos_df.astype(int).copy()
    pos_df['#CHR'] = pos_df['#CHR'].astype(str)
    # Creating temporary columns
    pos_df['new_POS'] = (pos_df['POS'] - 1)
    pos_df['POS_end'] = (pos_df['POS'] + 1)
    # Regions to look up
    pos_df['tuple'] = pos_df.set_index(
        ['#CHR', 'new_POS', 'POS_end']).astype(str).index.tolist()
    # Concat of CHR and POS as potentiel index
    pos_df['SNP_POS'] = pos_df['#CHR'] + '_' + pos_df['SNP_POS'].astype(str)
    # Removing temporary columns
    pos_df.drop(['new_POS', 'POS_end', 'POS'])
    return pos_df


def tabix_listregion_listfiles(list_region, file, error_files, cols=[]):
    df_final = pd.DataFrame(columns=cols)
    error_files[file] = []
    try:
        bam = pysam.AlignmentFile(file)
        SQs = pd.DataFrame([el.split('\t')
                            for el in str(bam.header).split('\n')])
        SQs = SQs[SQs[0] == '@SQ'][[1, 2]]
        ref_namelenght = SQs[1].str.cat(SQs[2], '_')
        for n, region_n in enumerate(list_region):
            try:
                rows_n = [x for x in bam.fetch(*region_n)]
                row_array = [str(r).split('\t') for r in rows_n]
                df_reg = pd.DataFrame(row_array, columns=cols)
                # Adding
                df_reg['region'] = str(region_n)
                df_reg['file'] = file
                if 'PROM1' in file:
                    df_reg['Phenotype'] = 'Severe'
                else:
                    df_reg['Phenotype'] = 'Moderate'
                df_final = pd.concat([df_final, df_reg], axis=0)
            except Exception as err:
                error_files[file].append(f'Error with {region_n}\n{err}')

    except OSError:
        error_files[file] = f'{file} or its index not found\n'
    return df_final.loc[df_final.duplicated() == False, [
        'region', 'QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'TLEN', 'SEQ', 'CIGAR']], ref_namelenght


def main(file_list, hits_table):
    # Remove existing file
    if output_name in os.listdir():
        print('Removing previous file')
        os.remove(output_name)

    # Reading the tables
    # file_ls = pd.read_table(file_list, header=0).iloc[:, 0].tolist()
    hits_df = pd.read_table(hits_table)
    assert hits_df[hits_df[['#CHR', 'POS']].duplicated(
        )].empty, 'Chromosome and SNP Position not enough to create unique indexes'
    table_region = region_select_fromSNP(
        hits_df[['#CHR', 'POS', 'REF', 'ALT']])
    list_region = set(table_region['tuple'].astype(str))
    # Do we need to add info there ?

    ref_table = pd.DataFrame()
    error_files = {}
    # TODO: recuperate column names automatically
    bam_cols = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ',
                'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', '?']
    for file in file_list:
        print(file)
        error_files[file] = []
        df_final, ref_table[file] = tabix_listregion_listfiles(
            list_region, file, error_files, cols=bam_cols)
        if file == file_list[0]:
            df_final.to_csv(output_name, mode='a',
                            header=True, index=False)
        else:
            df_final.to_csv(output_name, mode='a',
                            header=False, index=False)

    print('DICT OF ERRORS', file=sys.stderr)
    print([(file, err) for (file, err) in zip(
        error_files.keys(), error_files.values()) if err != []],
         file=sys.stderr)
    ref_table.to_csv('Ref_table.csv')


if __name__ == '__main__':
    main(file_list, hits_table)
