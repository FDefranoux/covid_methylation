import pandas as pd
import seaborn as sns
from scipy.stats import zscore
import matplotlib.pyplot as plt
import pingouin as pg
from sklearn.linear_model import LinearRegression
from scipy.stats import spearmanr
import numpy as np
import sys

file = 'Filtered_nano_bam_files_all_samples.csv'
file_snp = 'significant_hits_COVID19_HGI_A2_ALL_leave_23andme_20210607.txt'


def select_SNP_per_pvalue(file, pval_col, dist_bp=500000):
    snp_df = pd.read_table(file, usecols=['#CHR', 'POS', 'SNP', pval_col])
    best_snp = []
    for chr in snp_df['#CHR'].unique():
        chr_df = snp_df[snp_df['#CHR'] == chr]
        while chr_df.shape[0] > 0:
            best_pval = chr_df[pval_col].min()
            best_snp += chr_df[chr_df[pval_col] == best_pval]['SNP'].unique().tolist()
            pos = chr_df[chr_df[pval_col] == best_pval]['POS'].unique().tolist()[-1]
            chr_df = chr_df[(chr_df['POS'] <  pos - dist_bp) | (chr_df['POS'] >  pos + dist_bp) ]
    return best_snp


def drop_datas(df, group_ls, thresh_zscore=3, rep=3):
    # Calculation of z-score for outliers
    df_zscore = df.select_dtypes(include=['float']).apply(zscore)
    outliers = df_zscore[abs(df_zscore) > thresh_zscore].dropna(
        how='all').index.tolist()

    # Table of replicate number per gouping var (can be list of var)
    replicate_df = df.drop(outliers, axis=0).groupby(group_ls).count()
    no_replicate_group = replicate_df[replicate_df < rep].dropna(
        how='all').index
    no_replicate_index = df[df[group_ls].isin(no_replicate_group)].index
    index_ls = set(list(no_replicate_index) + list(outliers))

    # Table for the summary
    removed = df.loc[index_ls].copy()
    removed['removed'] = ''
    if len(no_replicate_index) > 0:
        removed.loc[no_replicate_index, 'removed'] = 'Not enough replicates'
    if len(outliers) > 0:
        removed.loc[outliers, 'removed'] = 'Outliers'
    if len(set(no_replicate_index).intersection(outliers)) > 0:
        removed.loc[set(no_replicate_index).intersection(
            outliers), 'removed'] = 'Both'
    print(removed['removed'].value_counts().to_markdown())

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
    return expe_info.drop(col, axis=1)


def genotype_filter(df, n_genotypes=2):
    new_df = df.copy()

    # n genotypes represented  min
    snp_all_genotype = new_df.groupby(['SNP', 'Genotype']).size(
        ).unstack().dropna(thresh=n_genotypes).index.unique().tolist()
    cpg_all_genotype = new_df.groupby(['cpg', 'Genotype']).size(
        ).unstack().dropna(thresh=n_genotypes).index.unique().tolist()
    new_df = new_df[new_df['SNP'].isin(
        snp_all_genotype) & new_df['cpg'].isin(cpg_all_genotype)].copy()

    print(f'\nRows dropped: {df.shape[0]-new_df.shape[0]}\n\n',
          (new_df.nunique() - df.nunique()).to_markdown())
    return new_df


def count_filter(df, min_count=10, n_genotypes=2):
    new_df = df.copy()

    # minimal count
    snp_counts = new_df.groupby(['SNP', 'Genotype']).size(
        ).unstack()
    snp_ls = snp_counts[snp_counts > min_count].dropna(
        thresh=n_genotypes).index.unique().tolist()

    cpg_counts = new_df.groupby(['cpg', 'Genotype']).size(
        ).unstack()
    cpg_ls = cpg_counts[cpg_counts > min_count].dropna(
        thresh=n_genotypes).index.unique().tolist()

    new_df = new_df[new_df['SNP'].isin(
        snp_ls) & new_df['cpg'].isin(cpg_ls)].copy()

    print(f'\nRows dropped: {df.shape[0]-new_df.shape[0]}\n\n',
          (new_df.nunique() - df.nunique()).to_markdown())
    return new_df


def multiple_mann_whitney(df, col, measure, title=''):
    mwu = pd.DataFrame()
    for v in range(df[col].nunique() - 1):
        for v2 in range(v+1, df[col].nunique()):
            val = df[col].unique()[v]
            val2 = df[col].unique()[v2]
            try:
                mwu_int = pg.mwu(df[df[col] == val][measure],
                df[df[col] == val2][measure]
                    ).rename(index=lambda s: s +
                    f' {title} {"-".join(sorted([val, val2]))}')
                mwu = pd.concat([mwu, mwu_int], axis=0)
            except:
                pass
    return mwu


def run_stat(df, unit='', var='', measure='log_lik_ratio', suppl_title=''):
    df = df.copy()
    results = pd.DataFrame(index=df[unit].unique())
    for u in df[unit].unique():
        u_df = df[df[unit] == u]

        # Counts
        count_dict = u_df[var].value_counts().rename(
            index=lambda s: 'Counts ' + s).to_dict()

        # Diff between means
        diff_means = u_df[u_df['Gen'] == 'alt'][measure].mean(
            ) - u_df[u_df['Gen'] == 'ref'][measure].mean()

        # Mann Whitney
        mwu = multiple_mann_whitney(u_df, var, measure)

        # Linear regression
        # df['phenotype_dum'] = pd.get_dummies(df['phenotype']).iloc[:, 0]
        try:
            lm = pg.linear_regression(u_df[measure], u_df['Genotype_dum'])
        except:
            lm = pd.DataFrame()

        # Spearman correlation
        spear = spearmanr(u_df[measure], u_df[var])

        # Results
        results.loc[u,count_dict.keys()] = count_dict.values()
        results.loc[u, 'diff_means_altVSref'] = diff_means
        results.loc[u, f'Spearman correlation {var} p-value'] = spear[1]
        results.loc[u, f'Spearman correlation {var} rho'] = spear[0]
        if not lm.empty:
            try:
                results.loc[u, [f'LM {var} intercept',
                                f'LM {var} coef', 'r2']] = lm['coef'].tolist(
                                    ) + lm['r2'].tolist()[:1]
            except:
                pass
        if not mwu.empty:
            results.loc[u, mwu.index] = mwu['p-val']
    results.to_csv(f'Stat_Analysis_{measure}VS{var}_per_{unit}_{suppl_title}.csv')
    return results


### PLOTS ###
def violinplot(df, title_supp='', yvar='cpg', xvar='log_lik_ratio', huevar='Genotype', colvar='phenotype'):
    n_snp = df[yvar].nunique()
    # print(n_snp)
    if n_snp < 10:
        g_chr = sns.catplot(data=df, kind='violin',
                            y=yvar, x=xvar, hue=huevar,
                            col=colvar, height=6,
                            aspect=0.9, orient='h',
                            inner="quartile", linewidth=1.5,
                            sharex=False, sharey=False)
    elif n_snp < 50:
        g_chr = sns.catplot(data=df, kind='violin',
                            y=yvar, x=xvar, hue=huevar,
                            col=colvar, height=25,
                            aspect=.7, orient='h',
                            inner="quartile", linewidth=1.5,
                            sharex=False, sharey=False)
    elif n_snp < 500:
        g_chr = sns.catplot(data=df, kind='violin',
                            y=yvar, x=xvar, hue=huevar,
                            col=colvar, height=45,
                            aspect=.15, orient='h',
                            inner="quartile", linewidth=1,
                            sharex=False, sharey=False)
    else:
        g_chr = sns.catplot(data=df, kind='violin',
                            y=yvar, x=xvar, hue=huevar,
                            col=colvar, height=45,
                            aspect=.1, orient='h',
                            inner="quartile", linewidth=1,
                            sharex=False, sharey=False)
    g_chr.savefig(
        f'violinplot_median_all_{yvar}_{xvar}_distribution_{title_supp}.png')


def ridgeplot(df, chr='', rowvar='cpg', huevar='Genotype', colvar='phenotype', var='log_lik_ratio'):
    with sns.plotting_context('paper', font_scale=0.1):
        sns.set_theme(style="white", rc={
                      "axes.facecolor": (0, 0, 0, 0), "font.size": 16})
        g_ridge = sns.FacetGrid(df, col=colvar, row=rowvar,
                                hue=huevar,
                                aspect=10,
                                height=4,
                                palette="mako",
                                margin_titles=True)

        [plt.setp(ax.texts, text="") for ax in g_ridge.axes.flat]
        # Draw the densities in a few steps
        g_ridge.map(sns.kdeplot, var, bw_adjust=.2,
                    clip_on=False,
                    fill=True, alpha=0.4, linewidth=0,
                    legend=True)
        g_ridge.map(sns.kdeplot, var, clip_on=False, color="w",
                    lw=2, bw_adjust=.2, cut=0)

        # Set the subplots to overlap
        plt.subplots_adjust(hspace=-0.3, wspace=-0.5, top=5, bottom=4)

        # # Remove axes details that don't play well with overlap
        # g.set_titles(row_template = '{col_name}', col_template = '{row_name}', loc='left',
        #         fontdict = {'fontsize': 2, 'color': 'c'})
        g_ridge.set_titles("")
        g_ridge.set(yticks=[], ylabel="")
        g_ridge.despine(bottom=True, left=True)
        plt.legend(bbox_to_anchor=(0, 2),
                   fontsize='xx-large', facecolor='white')
        g_ridge.savefig(
            f'rigdeplot_median_all_{rowvar}_{var}_distribution_chr{chr}_color{huevar}.png')


def scatter_with_unique_cpg(median_df, huevar='Genotype', colvar='phenotype', xvar='cpg', yvar='log_lik_ratio'):
    n = 0
    median_df = median_df.copy()
    for chr in median_df['CHR'].sort_values().unique():
        var_df = median_df.loc[median_df['CHR'] == chr, xvar].sort_values()
        for val in var_df.unique():
            median_df.loc[(median_df[xvar] == val)
                          & (median_df['CHR'] == chr), f'{xvar}_unique'] = n
            n = n + 1
        n = n + 5
        median_df.loc[(median_df['CHR'] == chr), 'sep'] = n

    g = sns.relplot(kind='scatter', data=median_df, x=f'{xvar}_unique', y=yvar,
                    hue=huevar, col=colvar, aspect=4)
    g.set(xlabel="CHR", xticks=median_df['sep'].unique(
        ), xticklabels=median_df['CHR'].unique())
    g.savefig(
        f'scatterplot_median_all_{yvar}_unique{xvar}_color{huevar}.png')


def linear_reg_plot(df, var='', unit='', plot=False, title_suppl=''):
        # Linear regression
    g = sns.lmplot(data=u_df, x='Genotype',
                   y=var, hue='phenotype', legend=False)
    g.set(xticks=[0, 1, 2], xticklabels=[
          '0/0', '0/1', '1/1'])
    plt.legend(labels=['Mild', 'Severe'], loc='upper left')
    plt.title(unit + ' ' + u)


# MAIN
def main(file):
    # Selection of the SNPs
    snp_ls = select_SNP_per_pvalue(file_snp, pval_col='all_inv_var_meta_p',
        dist_bp=100000)
    gen_ls = ['0/0', '0/1', '1/1']

    # Reading in chunks
    # chunks = []
    # size=0
    # n=0
    # for chunk in pd.read_csv(file, chunksize=1000000, usecols=['cpg', 'SNP', 'phenotype', 'Genotype',
    #                                      'name', 'log_lik_ratio', 'CHR', 'Gen', 'read_name']):
    #     part = chunk[(chunk['Gen'] != 'other') & (chunk['SNP'].isin(snp_ls)) & (chunk['Genotype'].isin(gen_ls))]
    #     chunks.append(part)
    #     size = part.size + size
    #     n += 1
    #     print(n, size, flush=True)
    #
    # all_df = pd.concat(chunks)
    # all_df.to_csv('Filtered_nano_bam_files_all_samples.csv', index=False, mode='w')
    # del chunks
    all_df = pd.read_table(file, usecols=['cpg', 'SNP', 'phenotype', 'Genotype',
                                         'name', 'log_lik_ratio', 'CHR', 'Gen', 'read_name'])

    # QUESTION: Dropping outliers ?
    print('SHAPE all_df', all_df.shape, all_df.size)
    # Median over the read_name with same Allele calling
    median_df = all_df.groupby(['phenotype', 'name', 'CHR', 'cpg',
        'SNP', 'Genotype', 'Gen']).median().reset_index()
    del all_df
    # QUESTION: What should we do with the ALT calls in '0/0' and REF in '1/1'?
    print('SHAPE median_df', median_df.shape,  median_df.size, flush=True)

    # Filter per genotype
    print(median_df['Genotype'].value_counts())
    median_df = median_df[median_df['Gen'] != 'other']
    print('SHAPE median_df', median_df.shape,  median_df.size, flush=True)
    # median_df = genotype_filter(median_df, n_genotypes=2)
    median_df = count_filter(median_df, min_count=5, n_genotypes=2)

    # Filter per quantity
    median_df['Genotype_dum'] = median_df['Genotype'].replace({'0/0': 0, '0/1': 1, '1/1': 2})

    # STATS
    # for unit in ['cpg', 'SNP']:
    #     stat = run_stat(median_df, unit=unit, var='Genotype', measure='log_lik_ratio')
    #
    #     # Special plot (kinda MannHattan plot)
    #     stat = stat.reset_index()
    #     stat['-log10'] = np.log10(stat['Spearman correlation Genotype p-value'].astype(float)) * (-1)
    #     stat[['CHR', 'POS']] = stat['index'].str.split(':', expand=True)[[0, 1]].astype(int)
    #     stat = stat.sort_values('CHR')
    #     stat['CHR'] = stat['CHR'].astype('category')
    #     g = sns.relplot(kind='scatter', data=stat, x='index',
    #         y='-log10', aspect=4, hue='CHR', cmap='Spectral', legend=True)
    #     g.set(xlabel="CHR", xticks=stat.groupby(['CHR']).last()['index'].unique(),
    #         xticklabels=stat['CHR'].unique())
    #     g.savefig(f'-log10_Spearman_pvalue-{unit}.png')
    #
    #     g1 = sns.relplot(kind='scatter', data=stat, y='diff_means_altVSref', x='Spearman correlation Genotype rho')
    #     g1.savefig(f'Diff_meansVSrho_{unit}.png')
    #     del stat
    #
    #     # Stats heterozygotes
    #     stat_het = run_stat(median_df[median_df['Genotype'] == '0/1'], unit=unit,
    #         measure='log_lik_ratio', var='Gen', suppl_title='Het_only')
    #     g1 = sns.relplot(kind='scatter', data=stat_het, y='diff_means_altVSref', x='Spearman correlation Gen rho')
    #     g1.savefig(f'Diff_meansVSrho_{unit}_heterozygotes.png')
    #     del stat_het

    # PLOTS
    scatter_with_unique_cpg(median_df, huevar='Genotype',
                            colvar=None, xvar='SNP', yvar='log_lik_ratio')
    scatter_with_unique_cpg(median_df, huevar='Gen',
                            colvar=None, xvar='cpg', yvar='log_lik_ratio')

    # Violinplot only for the SNP
    violinplot(median_df[median_df['Genotype'] == '0/1'].sort_values(
        by=['CHR', 'SNP', 'Genotype']), title_supp='_heterozygotes', yvar='SNP',
        xvar='log_lik_ratio', huevar='Gen', colvar=None)
    violinplot(median_df.sort_values(by=['CHR', 'SNP', 'Genotype']), title_supp='_best_pval', yvar='SNP',
        xvar='log_lik_ratio', huevar='Genotype', colvar=None)
    violinplot(median_df.sort_values(by=['CHR','cpg', 'Genotype']), title_supp='_best_pval', yvar='cpg',
        xvar='log_lik_ratio', huevar='Genotype', colvar=None)


if __name__ == '__main__':
    main(file)

# TODO:
# Violin plot
# Fix ridge PLOTS not working with cpg, bigger font size
# Question:
# Arrange the script to have an easy analysis --> bin/all/clean_extreme ....
# Separation of the SNP heterozygote and analysis of both alleles separetely

def other_plots():
    # Plot all SNP, discarding the very different ones.
    for chr in median_df['CHR'].sort_values().unique():
        print('CHR = ', chr)
        cleaned_per_chr, _ = drop_datas(median_df[median_df['CHR'] == chr], 'name', thresh_zscore=2.5)
        try:
            violinplot(cleaned_per_chr, title_supp=str(chr)+'_removedextreme', yvar='cpg',
                       xvar='log_lik_ratio', huevar='Genotype', colvar=None)
        except ValueError:
            violinplot(median_df[median_df['CHR'] == chr], title_supp=str(chr)+'_all', yvar='cpg',
                       xvar='log_lik_ratio', huevar='Genotype', colvar=None)
        # violinplot(cleaned_per_chr, chr=chr, yvar='SNP',
        #            xvar='log_lik_ratio', huevar='Genotype', colvar='phenotype')

    # Plotting the most represented
    most_represented = median_df.groupby(['cpg', 'CHR']).count().reset_index(
                        ).groupby(['CHR']).first()['cpg']
    most_represented_df = median_df[median_df['cpg'].isin(most_represented)]
    violinplot(most_represented_df, title_supp=str(chr)+'_most_rep', yvar='cpg',
        xvar='log_lik_ratio', huevar='Genotype', colvar='phenotype')

    # Mean/Median over binned sites

    # binned = median_df.copy()
    binned, _ = drop_datas(median_df, 'name', thresh_zscore=3)
    binned['pos'] = binned['cpg'].str.split(':', expand=True)[1].astype(int)

    # binned['distance_to_median'] = binned['cpg'].str.split(':', expand=True)[1].astype(int)
    for chr in binned['CHR'].sort_values().unique():
        chr_index = binned[binned['CHR'] == chr].index
        binned.loc[chr_index, 'bin'] = (
            (binned.loc[chr_index, 'pos'] - binned.loc[chr_index, 'pos'].median()
            )/binned.loc[chr_index, 'pos'].median()).round(1)
        # print(f'CHR = {chr}\n\n', binned['bin'] .value_counts().to_markdown())
        # print(binned['pos'] .value_counts().to_markdown())
    binned['bin'] = binned['CHR'].astype(str) + ' : ' + binned['bin'].astype(str)
    violinplot(binned, title_supp=str(chr)+'_per_bin', yvar='bin',
        xvar='log_lik_ratio', huevar='Genotype', colvar='phenotype')

        # binned.loc[binned['CHR'] == chr, 'bin'] = pd.cut(
        #     x=binned[binned['CHR'] == chr]['cpg'].str.split(
        #         ':', expand=True)[1].astype(int), bins=4).astype(str)
        # binned.loc[binned['CHR'] == chr, 'bin'] = pd.qcut(
        #     x=binned[binned['CHR'] == chr]['cpg'].str.split(
        #         ':', expand=True)[1].astype(int), q=3, duplicates='drop').astype(str)
    binned['bin'] = binned['CHR'].astype(str) + ': ' +  binned['bin']
    binned.sort_values('cpg', inplace=True)
    violinplot(binned, title_supp=str(chr)+'_per_bin', yvar='bin',
        xvar='log_lik_ratio', huevar='Genotype', colvar=None)

    scatter_with_unique_cpg(median_df, huevar='phenotype',
                            colvar='Genotype', xvar='cpg', yvar='log_lik_ratio')
    scatter_with_unique_cpg(median_df, huevar='phenotype',
                            colvar='Genotype', xvar='SNP', yvar='log_lik_ratio')
    #
    for chr in median_df['CHR'].sort_values().unique():
        violinplot(median_df[median_df['CHR'] == chr], title_supp=chr, yvar='cpg',
                   xvar='log_lik_ratio', huevar='Genotype', colvar='phenotype')
        violinplot(median_df[median_df['CHR'] == chr], title_supp=chr, yvar='SNP',
                   xvar='log_lik_ratio', huevar='Genotype', colvar='phenotype')
        # ridgeplot(median_df[median_df['CHR'] == 1], chr=1, rowvar='SNP',
        #           huevar='Genotype', colvar='phenotype', var='log_lik_ratio')
