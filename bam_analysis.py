import pandas as pd
import seaborn as sns
import os

file_allele = 'all_nanopore_base_at_significant_snps_oncovmeth_05102021.txt'


def count_table(df, t=0.90):
    df['Gen'] = 'other'
    df.loc[df['Allele'] == df['ref'], 'Gen'] = 'ref'
    df.loc[df['Allele'] == df['alt'], 'Gen'] = 'alt'
    table_counts = df.groupby(['name', 'SNP', 'Gen']
                              ).size().unstack(fill_value=0)
    for haplo in ['ref', 'alt', 'other']:
        if df[df['Gen'] == haplo].shape[0] == 0:
            table_counts[haplo] = 0
    sum = table_counts[['alt', 'other', 'ref']].sum(axis=1)
    percent = pd.DataFrame()
    table_counts['sum'] = sum
    percent['0'] = table_counts['ref'].astype(int) / sum
    percent['1'] = table_counts['alt'].astype(int) / sum
    percent['2'] = table_counts['other'].astype(int) / sum
    percent = percent[percent > 1 - t]
    for col in percent.columns:
        percent.loc[(percent[col] == percent.max(axis=1)), 'Allele1'] = col
        percent.loc[(percent[col] == percent.min(axis=1))
                    & (percent['Allele1'] != col), 'Allele2'] = col
        percent.loc[(percent.max(axis=1) == percent.min(axis=1))
                    & (percent[col] == percent.max(axis=1))
                    & (percent['Allele2'].isna()), 'Allele2'] = col
    alleles = percent[['Allele1', 'Allele2']].astype(int).copy()
    percent['Genotype'] = alleles.min(axis=1).astype(
        str) + '/' + alleles.max(axis=1).astype(str)
    return percent['Genotype'].reset_index()


def verify_counts(df):
    counts = df.groupby(['name', 'Genotype', 'Gen']).size().unstack()
    sum = counts.sum(axis=1)
    counts['alt'] = counts['alt'] / sum
    counts['ref'] = counts['ref'] / sum
    counts['other'] = counts['other'] / sum
    print(counts.to_markdown())


def main(file_allele):

    # Opening the allele_table
    df = pd.read_table(file_allele, sep='\t',
                       error_bad_lines=False, header=None)
    df.columns = ['name', 'SNP', 'read_name', 'Allele']
    df[['CHR', 'pos', 'ref', 'alt']] = df['SNP'].str.split(':', expand=True)
    #Getting nanopolish file
    nano_cols = ['CHR', 'strand', 'start', 'end', 'read_name',
                 'log_lik_ratio', 'log_lik_methylated', 'log_lik_unmethylated',
                 'num_calling_strands', 'num_motifs', 'sequence']
    file_ls = set(df['name'].to_list())
    file_ls = ['gcc00022_PROM1']
    for file in file_ls:
        try:
            nano_file = os.path.join('nanopolish_grep_reads', file + '.txt')
            nano_df = pd.read_table(nano_file, header=None, names=nano_cols)
            df['CHR'] = df['CHR'].astype(int)
            nano_all = pd.merge(
                nano_df, df[df['name'] == file], on=['read_name', 'CHR'])
            nano_all[['pos', 'start', 'end']] = nano_all[[
                'pos', 'start', 'end']].astype(int)
            nano_all = nano_all[(nano_all['pos'] < nano_all['start']-5)
                                & (nano_all['pos'] < nano_all['end']+5)]

            count_df = count_table(nano_all)
            merge = pd.merge(nano_all, count_df[['SNP', 'name', 'Genotype']],
                             on=['SNP', 'name'], copy=False)
            verify_counts(merge)
            if nano_all.shape[0] != merge.shape[0]:
                print(f'{file} SNP-ID is not unique !')
            merge.to_csv(os.path.join('nanopolish_grep_reads',
                                      file + '_genotyped.csv'), index=False)
        except Exception as err:
            # print(file, err)
            try:
                print('nano_df\n', nano_df.head(2))
                print('\nnano_all\n', nano_all.head(2))
                print('\ncount_df\n', count_df.head(2))
                print('\nmerge\n', merge.head(2))
            except:
                pass


if __name__ == '__main__':
    main(file_allele)

# # Graphs
# blou = df_global[df_global != 0]
# blou['cut'] = pd.qcut(blou['sum'], q=4)
# f = sns.displot(data=blou, x=blou['ref'], kde=True, col=blou['cut'], col_wrap=2, facet_kws={
#                 'sharex': False, 'sharey': False})
# f.savefig('histogram_ref_proportions_SNP_alleles_fromBAM.png')
# f.savefig('histogram_ref_proportions_SNP_alleles_fromBAM.png')
