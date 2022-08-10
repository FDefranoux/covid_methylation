import pandas as pd
import pysam
import os
import sys
import seaborn as sns
# TODO: duplication of rows SOLVED ??
# TODO: SNPs that are in the range but do not appear as row SOLVED ??

# Variable definition
file_list = 'nano_list.txt'
hits_table = 'significant_hits_COVID19_HGI_A2_ALL_leave_23andme_20210607.txt'
output = {'datas': 'Nanopolish_per_significant_region.csv'}


class SamFiles:
    # def __init__(self):
    #     self.error_files = {}
    def region(file):
        return file.fetch

    def reads(file):
        return file.find

    def sam_iterators(func, args, cols=None):
        try:
            args = tuple(args)
            rows = [x for x in func(*args)]
            array = [str(r).split('\t') for r in rows]
            df = pd.DataFrame(array, columns=cols)
        except Exception as err:
            print(f'Error with {args}\n{err}')
            # error_files[file].append(f'Error with {args}\n{err}')
            df = pd.DataFrame()
        return df

    def bam_ref_names(bam):
        # Bam file header saving
        try:
            bam_head = pd.DataFrame([el.split('\t')
                                     for el in str(bam.header).split('\n')])
            bam_head = bam_head[bam_head[0] == '@SQ'][[1, 2]]
            ref_namelenght = pd.DataFrame(
                bam_head[1].str.cat(bam_head[2], '_'))
            return ref_namelenght
        except Exception as err:
            print(err)
            return pd.DataFrame()

    def open(file):
        try:
            if '.bam' in file:
                sam_file = pysam.AlignmentFile(file, 'rb')
            else:
                sam_file = pysam.TabixFile(file)
            return sam_file
        except OSError as err:
            print(file)
            print(err)
            # error_files[file] = f'{type} file {file} or its index were not found'


def region_select_fromSNP(pos_df):
    pos_df = pos_df.copy()
    pos_df['#CHR'] = pos_df['#CHR'].astype(str)
    # Creating temporary columns
    pos_df['new_POS'] = (pos_df['POS'] - 1)
    pos_df['POS_end'] = (pos_df['POS'] + 1)
    # Regions to look up
    pos_df['tuple'] = pos_df.set_index(
        ['#CHR', 'new_POS', 'POS_end']).astype(str).index.tolist()
    # Concat of CHR and POS as potentiel index
    pos_df['SNP_POS'] = pos_df['#CHR'] + '_' + pos_df['POS'].astype(str)
    # Removing temporary columns
    pos_df.drop(['new_POS', 'POS_end', 'POS'], axis=1)
    return pos_df


def region_select_fromSNP_minmax(pos_df):
    # SEE how to implement this correctly and usefully
    # TODO: Assignement of SNP for each region
    pos_df = pos_df.copy()
    pos_df['new_POS'] = (pos_df['POS'] - 11).astype(str).str[:-6]
    pos_df['new_POS'] = pos_df['new_POS'] + '000000'
    minmax_table = pos_df.groupby(['#CHR', 'new_POS']).min(['POS']) - 500000
    minmax_table['max'] = pos_df.groupby(
        ['#CHR', 'new_POS']).max()['POS'] + 500000
    minmax_table.reset_index('#CHR', inplace=True)
    min_list = minmax_table.set_index(['#CHR', 'POS', 'max']).index.tolist()
    for n in range(len(min_list)-1):
        if min_list[n][0] == min_list[n+1][0]:
            if min_list[n][2] > min_list[n+1][1]:
                mini = min(min_list[n][1], min_list[n+1][1])
                maxi = max(min_list[n][2], min_list[n+1][2])
                min_list[n] = (min_list[n][0], mini, maxi)
                min_list[n+1] = (min_list[n][0], mini, maxi)
    minmax_table['tuple'] = min_list
    return minmax_table


def recup_nanopolish(file_ls, region_list, output_names={'datas': 'Bam_summary.csv',
                                                         'reference': 'Ref_table.csv'}):
    nano_cols = ['#CHR', 'strand', 'start', 'end', 'read_name', 'log_lik_ratio',
                 'log_lik_methylated', 'log_lik_unmethylated', 'num_calling_strands',
                 'num_motifs', 'sequence']
    for file in file_ls:
        nano_file = SamFiles.open(file)
        for region in region_list:
            nano_df = pd.DataFrame()
            try:
                nano_df = SamFiles.sam_iterators(
                    SamFiles.region(nano_file), region, cols=nano_cols)
                nano_df['SNP_hit'] = str(
                    region[0]) + '_' + str(region[1] + 1)  # WORKING?
                nano_df['file'] = file
            except:
                print(f'Error with iterating over file {file}')

            try:
                if file == file_ls[0]:
                    nano_df.to_csv(output_names['datas'], mode='a',
                                   header=True, index=False)
                else:
                    nano_df.to_csv(output_names['datas'], mode='a',
                                   header=False, index=False)
            except:
                print(f'unexpected error during the saving of file {file}')


def main(file_list, hits_table):
    # Remove existing file
    for name in output.values():
        if name in os.listdir():
            print(f'Removing previous file {name}')
            os.remove(name)

    # Reading Hits table
    file_ls = pd.read_table(file_list, header=None).iloc[:, 0].tolist()
    hits_df = pd.read_table(hits_table)
    assert hits_df[hits_df[['#CHR', 'POS']].duplicated(
        )].empty, 'CHR and SNP Position not enough to create unique indexes'
    # table_region = region_select_fromSNP(hits_df[['#CHR', 'POS']])
    table_region = region_select_fromSNP_minmax(hits_df[['#CHR', 'POS']])
    list_region = set(table_region['tuple'])
    print('list region len:', len(list_region))

    # Look up bam files
    recup_nanopolish(file_ls, list_region, output_names=output)


if __name__ == '__main__':
    main(file_list, hits_table)
