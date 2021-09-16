import pandas as pd
import pysam
import os
import sys


# Variable definition
file_list = os.listdir('bams')[:-1]
hits_table = 'significant_hits_COVID19_HGI_A2_ALL_leave_23andme_20210607.txt'
output_name = 'Bam_test.csv'


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


def main(file_list, hits_table):
    # Remove existing file
    if output_name in os.listdir():
        print('Removing previous file')
        os.remove(output_name)

    # Reading the tables
    file_ls = pd.read_table(file_list, header=0).iloc[:, 0].tolist()
    hits_df = pd.read_table(hits_table)

    list_region = set(region_select_fromSNP(hits_df[['#CHR', 'POS']]).values)
    error_files = {}
    for file in file_ls:
        print(file)
        error_files[file] = []
        # print(pd.read_table(file, sep='\t', nrows=1).columns.tolist() == cols)
        df_final = tabix_listregion_listfiles(
            list_region, file, error_files, cols=None)
        print(df_final.head(5))

        if file == file_ls[0]:
            df_final.to_csv(output_name, mode='a',
                            header=True, index=False)
        else:
            df_final.to_csv(output_name, mode='a',
                            header=False, index=False)

    print('DICT OF ERRORS', file=sys.stderr)
    print([(file, err) for (file, err) in zip(
        error_files.keys(), error_files.values()) if err != []],
         file=sys.stderr)
