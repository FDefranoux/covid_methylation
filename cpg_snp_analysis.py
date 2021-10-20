import pandas as pd
import glob
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.stats import zscore
import pingouin as pg


dir = 'nanopolish_grep_reads/*.csv'


def drop_datas(df, group_ls, thresh_zscore=5, rep=3):

    # Calculation of z-score for outliers
    df_zscore = df.select_dtypes(include=['float']).apply(zscore)
    outliers = df_zscore[df_zscore > thresh_zscore].dropna(
        how='all').index.tolist()
    # Table of replicate number per gouping var (can be list of var)
    replicate_df = df.drop(outliers, axis=0).groupby(group_ls).count()
    no_replicate_group = replicate_df[replicate_df < rep].dropna(
        how='all').index
    no_replicate_index = df[df[group_ls].isin(no_replicate_group)].index
    index_ls = set(list(no_replicate_index) + list(outliers))
    removed = df.loc[index_ls].copy()
    removed.loc[no_replicate_index, 'removed'] = 'Not enough replicates'
    removed.loc[outliers, 'removed'] = 'Outliers'
    removed.loc[set(no_replicate_index).intersection(
        outliers), 'removed'] = 'Both'
    print(removed['removed'].sort_index())
    # TODO: Modify the function to include automatic dropping of datas
    # per ex with Analysis of categorical datas (repartition of cat datas)
    # TODO: Do we need to remove the whole row when there is just one variable
    # not having appropriate number of replicates
    return df.drop(index_ls, axis=0), removed


def info_per_value(df, col, lim=2):
    expe_info = pd.DataFrame(columns=df.columns, index=df[col].unique())
    for val in df[col].unique():
        val_df = df[df[col] == val]
        try:
            expe_info.loc[val] = val_df.loc[val_df.index[0],
                                            val_df.nunique()[val_df.nunique() < lim].index]
        except IndexError:
            pass
    return expe_info.drop(col, axis=1).dropna(how='all').dropna(how='all', axis=1)


# Opening the file
all = pd.DataFrame()
for file in glob.glob(dir):
    inter = pd.read_csv(file)
    try:
        all = pd.concat([all, inter])
    except Exception as err:
        print(file, err)
all = pd.read_csv('nano_genotyped_5b.csv')
# Dropping empty rows and cols
all = all.dropna(how='all').dropna(how='all', axis=1).copy()  # Same shape
# Insert phenotype variable
all.loc[all['name'].str.contains('PROM1'), 'phenotype'] = 'Severe'
all['phenotype'].fillna('Mild', inplace=True)
all.loc[all['num_motifs'] == 1, 'distance_cpg_snp'] = abs(
    all['pos'] - all['start'])
# Data cleaning
# Deletions
all = all[all['ref'].str.len() < 2].copy()

# Filtering out the reads_names associated with several CHR
double_chr = all.groupby(
    ['read_name', 'CHR']).size().unstack().dropna(thresh=2).index
if len(double_chr) > 0:
    all.drop(double_chr, inplace=True)
try:
    nano_new, outliers = drop_datas(all, 'phenotype', thresh_zscore=5, rep=3)
except:
    print('Error in dropping outliers')
    nano_new = pd.DataFrame()

if len(nano_new) > 100:
    median_df = nano_new.groupby(
        ['Genotype', 'read_name', 'SNP', 'name', 'phenotype']).median().reset_index()
else:
    print('ERROR in droping outliers')
    median_df = nano_new.groupby(
        ['Genotype', 'read_name', 'SNP', 'name', 'phenotype']).median().reset_index()

# Removing the SNPs representing only one genotype
snp_index = median_df.groupby(['SNP', 'phenotype', 'Genotype']).size().unstack().filter(
    regex='^((?!2).)*$', axis=1).dropna(thresh=3, axis=0).reset_index()['SNP'].unique()

# Removing the 2 alleles ?
blou = info_per_value(median_df[median_df['SNP'].isin(snp_index)], 'Genotype')
blou

selected = median_df[(median_df['SNP'].isin(snp_index)) & (
    median_df['Genotype'].isin(['0/0', '0/1', '1/1']))].copy()
for chr in selected['CHR'].unique():
    pos_chr = selected.loc[selected['CHR'] == chr, 'pos'].copy()
    selected.loc[selected['CHR'] == chr, 'pos_bin'] = pd.qcut(
        pos_chr.astype(int), q=10, duplicates='drop')

selected.groupby(['CHR', 'pos_bin', 'Genotype']).size().unstack()

# When we categorize the SNPs per bins
g1 = sns.catplot(data=selected, kind='violin',
                 x='pos_bin', y='log_lik_ratio', hue='Genotype',
                 row='CHR', col='phenotype', aspect=4,
                 inner="quartile", linewidth=3,
                 sharex=False, sharey=False)

for chr in selected['CHR'].unique():
    n_snp = selected[selected['CHR'] == chr]['SNP'].nunique()
    print(int(chr), '--', n_snp)
    if n_snp < 5:
        g = sns.catplot(data=selected[selected['CHR'] == chr], kind='violin',
                        y='SNP', x='log_lik_ratio', hue='Genotype',
                        col='phenotype', height=4,
                        aspect=1,
                        inner="quartile", linewidth=2,
                        sharex=False, sharey=False)
    elif n_snp < 50:
        g = sns.catplot(data=selected[selected['CHR'] == chr], kind='violin',
                        y='SNP', x='log_lik_ratio', hue='Genotype',
                        col='phenotype', height=15,
                        aspect=.7,
                        inner="quartile", linewidth=2,
                        sharex=False, sharey=False)
    else:
        g = sns.catplot(data=selected[selected['CHR'] == chr], kind='violin',
                        y='SNP', x='log_lik_ratio', hue='Genotype',
                        col='phenotype', height=45,
                        aspect=.15,
                        inner="quartile", linewidth=1,
                        sharex=False, sharey=False)
    g.savefig(f'violinplot_median_all_SNPs_ratio_distribution_chr{chr}.png')
