import pandas as pd
import seaborn as sns
import os

file_allele = 'all_nanopore_base_at_significant_snps_oncovmeth_05102021.txt'


def count_table(df, t=0.95):
    df['Gen'] = 'other'
    df.loc[df['Allele'] == df['ref'], 'Gen'] = 'ref'
    df.loc[df['Allele'] == df['alt'], 'Gen'] = 'alt'
    table_counts = df.groupby(['name', 'SNP', 'Gen']
                              ).size().unstack(fill_value=0)
    for haplo in ['ref', 'alt', 'other']:
        if df[df['Gen'] == haplo].shape[0] == 0:
            table_counts[haplo] = 0
    sum = table_counts[['alt', 'other', 'ref']].sum(axis=1)
    table_counts['sum'] = sum
    table_counts['alt%'] = table_counts['alt'].astype(int) / sum
    table_counts['other%'] = table_counts['other'].astype(int) / sum
    table_counts['ref%'] = table_counts['ref'].astype(int) / sum

    table_counts.loc[(table_counts['alt'] > 0.5 * t)
                     & (table_counts['ref'] < 1 * t), 'Genotype'] = '0/1'
    table_counts.loc[(table_counts['ref'] > 0.5 * t)
                     & (table_counts['alt'] < 1 * t), 'Genotype'] = '0/1'
    table_counts.loc[table_counts['ref'] > 1.0 * t, 'Genotype'] = '0/0'
    table_counts.loc[table_counts['alt'] > 1.0 * t, 'Genotype'] = '1/1'
    table_counts.loc[table_counts['other'] > 1.0 * t, 'Genotype'] = '2/2'
    return table_counts


def main(file_allele):

    # Opening the allele_table
    df = pd.read_table(file_allele, sep='\t',
                       error_bad_lines=False, header=None)
    df.columns = ['name', 'SNP', 'read_name', 'Allele']
    df[['chr', 'pos', 'ref', 'alt']] = df['SNP'].str.split(
            ':', expand=True)

    #Getting nanopolish file
    nano_cols = ['CHR', 'strand', 'start', 'end', 'read_name',
                 'log_lik_ratio', 'log_lik_methylated', 'log_lik_unmethylated',
                 'num_calling_strands', 'num_motifs', 'sequence']
    file_ls = set(df['name'].to_list())
    for file in file_ls:
        try:
            print(file)
            nano_file = os.path.join('nanopolish_grep_reads', file + '.txt')
            nano_df = pd.read_table(nano_file, header=None, names=nano_cols)
            nano_all = pd.merge(
                nano_df, df[df['name'] == file], on='read_name')
            nano_all[['pos', 'start', 'end']] = nano_all[[
                'pos', 'start', 'end']].astype(int)
            nano_all = nano_all[(nano_all['pos'] < nano_all['start']-5)
                                & (nano_all['pos'] < nano_all['end']+5)]

            count_df = count_table(nano_all).reset_index()
            merge = pd.merge(nano_all, count_df[['SNP', 'name', 'Genotype']], on=[
                             'SNP', 'name'], copy=False)
            if nano_all.shape[0] != merge.shape[0]:
                print(f'{file} SNP-ID is not unique !')
            merge.to_csv(os.path.join('nanopolish_grep_reads',
                                      file + '_genotyped.csv'), index=False)
        except Exception as err:
            print(file, err)


if __name__ == '__main__':
    main(file_allele)

# # Graphs
# blou = df_global[df_global != 0]
# blou['cut'] = pd.qcut(blou['sum'], q=4)
# f = sns.displot(data=blou, x=blou['ref'], kde=True, col=blou['cut'], col_wrap=2, facet_kws={
#                 'sharex': False, 'sharey': False})
# f.savefig('histogram_ref_proportions_SNP_alleles_fromBAM.png')
# f.savefig('histogram_ref_proportions_SNP_alleles_fromBAM.png')
