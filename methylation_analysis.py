import pandas as pd


hits_table = 'significant_hits_COVID19_HGI_A2_ALL_leave_23andme_20210607.txt'
hits_df = pd.read_table(hits_table)
hits_df[['#CHR', 'POS']] = hits_df[['#CHR', 'POS']].astype(int)

file = 'Sign_all.csv'
df = pd.read_table(file, header=0, dtype=object, sep=',')
df.head(2)
df.columns
len(set(df.File))


def associate_snp_alleles(df_final, snps_refs):
    df = df_final[df_final.duplicated() == False].copy()
    df['rel_snp'] = df['SNP'].astype(int) - df['start'].astype(int)
    df['allele'] = '0'
    df[['SNP', '#CHR']] = df[['SNP', '#CHR']].astype(int)

    for n in df['rel_snp'].unique():
        df.loc[df['rel_snp'] == n, 'allele'] = df['sequence'].str[n]
    allele_df = pd.merge(snps_refs, df, right_on=[
                         'SNP', '#CHR'], left_on=['POS', '#CHR'])
    return allele_df


associate_snp_alleles(df, hits_df[['#CHR', 'POS', 'REF', 'ALT']])[
                      'allele'].value_counts()


# pd.options.display.html.table_schema = True
# pd.options.display.max_rows = None


# TABLE COUNTS
# TODO: genotype homozygote ref or homozygote alt --> table of counts
df.loc[df['strand'] == '+', ['File', '#CHR', 'SNP',
                             'sequence', 'start', 'end', 'REF', 'ALT', 'allele']]


df[['File', '#CHR', 'SNP', 'allele']].sort_values(
    ['#CHR', 'SNP']).value_counts().sort_index()
df.groupby(['File', 'SNP', 'Genotype']).count()
# Until you solve the problem of the columsn names
df = df[df['SNP'] == df['POS']]


#  GC MOTIFS
# TODO: Add gestion of the number of GC motifs to handle
# TODO: Measure the distance between GC motif and SNP
# # To find the other motifs
# df_final['distance_snp_motif'] = df_final['start'].astype(int) - df_final['SNP']
# df_final[df_final['num_motifs'].astype(int) > 1]['sequence'].str.find('GC')
##     df[df['num_motifs'] != '1']['sequence']
#     df['num_motifs'].astype(int).max()
#     df['sequence'].str.rfind('GC').max()
#     GC_motifs = pd.DataFrame()
#     for motif_n in range(1, df['num_motifs'].astype(int).max()):
#     df['sequence'].str.rfind('GC')
#     blou.value_counts()
#     df = df.copy()

# PLOTS
# TODO: Plot the right things
# import seaborn as sns
# sns.boxplot(x = df_final['Phenotype'], y=df_final['log_lik_methylated'].astype(float), hue=df_final['Genotype'])

# TODO: Think about the statistics to perform
# TODO: Using all reads for one sample try to do an average sequence ?
