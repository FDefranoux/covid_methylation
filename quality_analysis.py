import pandas as pd
import seaborn as sns
import os

file = 'Filtered_nano_bam_files_control_finemapped.csv'
# file = 'Filtered_nano_bam_files_test.csv'


def main(file):

    # Opening the allele_table
    df = pd.read_csv(file)
    os.system('mkdir quality')
    for name in df['sample_id'].unique()[0]:
        print(file, flush=True)
        df[df['sample_id'] == name].groupby('chromosome').nunique().to_csv(f'quality_test/Nunique_{name}.csv')
        df[df['sample_id'] == name].groupby('chromosome').size().to_csv(f'quality_test/Size_{name}.csv')

    # df.columns = ['name', 'SNP', 'read_name', 'Allele']
    # df[['CHR', 'pos', 'ref', 'alt']] = df['SNP'].str.split(':', expand=True)
    # #Getting nanopolish file
    # nano_cols = ['CHR', 'strand', 'start', 'end', 'read_name',
    #              'log_lik_ratio', 'log_lik_methylated', 'log_lik_unmethylated',
    #              'num_calling_strands', 'num_motifs', 'sequence']
    # file_ls = set(df['name'].to_list())
    # for file in file_ls:
    #     try:
    #         nano_file = os.path.join('nanopolish_grep_reads', file + '.txt')
    #         nano_df = pd.read_table(nano_file, header=None, names=nano_cols)
    #         df['CHR'] = df['CHR'].astype(int)
    #         nano_all = pd.merge(
    #             nano_df, df[df['name'] == file], on=['read_name', 'CHR'])
    #         nano_all[['pos', 'start', 'end']] = nano_all[[
    #             'pos', 'start', 'end']].astype(int)
    #
    #         count_df = count_table(nano_all)
    #         merge = pd.merge(nano_all, count_df[['SNP', 'name', 'Genotype']],
    #                          on=['SNP', 'name'], copy=False)
    #         verify_counts(merge)
    #         if nano_all.shape[0] != merge.shape[0]:
    #             print(f'{file} SNP-ID is not unique !')
    #         merge.to_csv(os.path.join('nanopolish_grep_reads',
    #                                   file + '_genotyped.csv'), index=False)
    #     except Exception as err:
    #         try:
    #             print('nano_df\n', nano_df.head(2))
    #         except:
    #             print(file, err)


if __name__ == '__main__':
    main(file)
