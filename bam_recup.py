import pandas as pd
import pysam
import os
import sys

# Variable definition
file_list = ['bams/gcc13846.bam']
print(file_list)
hits_table = 'significant_hits_COVID19_HGI_A2_ALL_leave_23andme_20210607.txt'
output_name = 'Bam_test.csv'


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
    pos_df['new_POS'] = (pos_df['POS'] - 1)
    pos_df['POS_end'] = (pos_df['POS'] + 1)
    pos_df['#CHR'] = pos_df['#CHR'].astype(str)
    pos_df['tuple'] = pos_df.set_index(
        ['#CHR', 'new_POS', 'POS_end']).index.tolist()
    return pos_df[['tuple', 'POS']].set_index('tuple')


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
                df_reg['region'] = str(region_n)
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

    table_region = region_select_fromSNP(hits_df[['#CHR', 'POS']])
    list_region = set(table_region.index).to_list()

    ref_table = pd.DataFrame()
    error_files = {}
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

# ## POST-ANALYSIS
# bam = pd.read_csv('Bam_test_SEQ.csv')
# bam.head(2)
# # TODO: Work on the CIGAR sequence to get the right allele
# # NEED: Reference file!
# # Function to change seq for reverse reads when FLAG == 16
# trans = "ATGC".maketrans("ATGC", "TACG")
# bam.loc[bam['FLAG'] == 16, ['SEQ']] = bam['SEQ'].str.translate(trans)
# bam.loc[bam['FLAG'] == 16, ['CIGAR']] = bam['CIGAR'][::-1]
