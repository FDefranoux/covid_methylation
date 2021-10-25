import pandas as pd
import glob
import seaborn as sns
from scipy.stats import zscore
import matplotlib.pyplot as plt
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
    return expe_info.drop(col, axis=1).dropna(how='all').dropna(how='all', axis=1)


def gather_dfs_fromdir(dir):
    df = pd.DataFrame()
    for file in glob.glob(dir):
        inter = pd.read_csv(file)
        try:
            df = pd.concat([df, inter])
        except Exception as err:
            print(file, err)
    return df


def filtering_datas(df, force=False):
    print(f'Initial dataframe shape: {df.shape}')
    # Dropping empty rows and cols
    new_df = df.dropna(how='all').dropna(how='all', axis=1).copy()

    # Deletions
    new_df = new_df[new_df['ref'].str.len() < 2].copy()
    print('Filtering deletions: ', new_df.shape)

    # Filtering out the reads_names associated with several CHR
    double_chr = new_df.groupby(
        ['read_name', 'CHR']).size().unstack().dropna(thresh=2).index
    if len(double_chr) > 0:
        new_df.drop(double_chr, inplace=True)
    print('Filtering read_name in several chr: ', new_df.shape)

    # Drop SNP with 'other' allele (non-ref non-alt)
    print(new_df[new_df['Genotype'].isin(['0/0', '0/1', '1/1'])
                 == False]['Genotype'].value_counts().to_markdown())
    new_df = new_df[new_df['Genotype'].isin(['0/0', '0/1', '1/1'])].copy()
    print('Filtering non-ref non-alt alleles: ', new_df.shape)

    # All genotypes represented
    # snp_all_genotype = new_df.groupby(['SNP', 'Genotype']).size(
    #     ).unstack().dropna(thresh=2).reset_index()['SNP'].unique().tolist()
    # print(len(snp_all_genotype))
    # new_df = new_df[new_df['SNP'].isin(snp_all_genotype)].copy()
    # print('Filtering SNP with one genotype represented: ', new_df.shape)

    # Drop outliers and samples missing replicates
    try:
        nano_new, outliers = drop_datas(new_df, 'name', thresh_zscore=3, rep=3)
    except:
        print('Error in dropping outliers')
    print(f'\nFiltered dataframe shape: {new_df.shape}\n')
    print('\nRemoval summary:\n',
          (df.nunique() - new_df.nunique()).T.to_markdown())
    if force:
        return nano_new
    else:
        if df.shape[0] - new_df.shape[0] > df.shape[0] / 3:

            print('More than one third of the datas would be removed !'
                  + 'Returning initial dataframe')
            return new_df
        else:
            return nano_new


def violinplot(df):
    n_snp = df['SNP'].nunique()
    if n_snp < 5:
        g_chr = sns.catplot(data=df, kind='violin',
                            y='SNP', x='log_lik_ratio', hue='Genotype',
                            col='phenotype', height=4,
                            aspect=1,
                            inner="quartile", linewidth=2,
                            sharex=False, sharey=False)
    elif n_snp < 50:
        g_chr = sns.catplot(data=df, kind='violin',
                            y='SNP', x='log_lik_ratio', hue='Genotype',
                            col='phenotype', height=15,
                            aspect=.7,
                            inner="quartile", linewidth=2,
                            sharex=False, sharey=False)
    else:
        g_chr = sns.catplot(data=df, kind='violin',
                            y='SNP', x='log_lik_ratio', hue='Genotype',
                            col='phenotype', height=45,
                            aspect=.15,
                            inner="quartile", linewidth=1,
                            sharex=False, sharey=False)
    g_chr.savefig(
        f'violinplot_median_all_SNPs_ratio_distribution_chr{chr}.png')


def ridgeplot(df):
    with sns.plotting_context('paper', font_scale=2):
        sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
        g_ridge = sns.FacetGrid(df, col="phenotype", row='pos',
                                hue="Genotype", aspect=10, height=4, palette="mako",
                                margin_titles=True,  legend_out=False)

        [plt.setp(ax.texts, text="") for ax in g_ridge.axes.flat]
        # Draw the densities in a few steps
        g_ridge.map(sns.kdeplot, "log_lik_ratio", bw_adjust=.2, clip_on=False,
                    fill=True, alpha=0.4, linewidth=0, legend=True, thresh=0.5)
        g_ridge.map(sns.kdeplot, "log_lik_ratio", clip_on=False, color="w",
                    lw=2, bw_adjust=.2, cut=0, thresh=0.5)

        # Set the subplots to overlap
        plt.subplots_adjust(hspace=-0.3, wspace=-0.5, top=5, bottom=4)

        # # Remove axes details that don't play well with overlap
        # g.set_titles(row_template = '{col_name}', col_template = '{row_name}', loc='left',
        #         fontdict = {'fontsize': 2, 'color': 'c'})
        g_ridge.set_titles("")
        g_ridge.set(yticks=[], ylabel="")
        g_ridge.despine(bottom=True, left=True)
        try:
            g_ridge.legend(bbox_to_anchor=(0.5, 1))
        except:
            pass
        g_ridge.savefig(
            f'rigdeplot_median_all_SNPs_ratio_distribution_chr{chr}.png')


def scatter_with_unique_pos(median_df):
    n = 0
    for chr in median_df['CHR'].sort_values().unique():
        pos = median_df.loc[median_df['CHR'] == chr, 'pos'].sort_values()
        for val in pos.unique():
            median_df.loc[(median_df['pos'] == val)
                          & (median_df['CHR'] == chr), 'pos_unique'] = n
            n = n + 1
        n = n + 10
        median_df.loc[(median_df['CHR'] == chr), 'sep'] = n

    g = sns.relplot(kind='scatter', data=median_df, x='pos_unique', y='log_lik_ratio',
                    hue='Genotype', col='phenotype', aspect=4)
    g.set(xlabel="CHR", xticks=median_df['sep'].unique(
    ), xticklabels=median_df['CHR'].unique())
    g.savefig('scatterplot_median_all_snp_ratio_uniquepos_colorg.png')


def scatter_with_unique_pos2(median_df):
    n = 0
    for chr in median_df['CHR'].sort_values().unique():
        pos = median_df.loc[median_df['CHR'] == chr, 'pos'].sort_values()
        for val in pos.unique():
            median_df.loc[(median_df['pos'] == val)
                          & (median_df['CHR'] == chr), 'pos_unique'] = n
            n = n + 1
        n = n + 10
        median_df.loc[(median_df['CHR'] == chr), 'sep'] = n

    g = sns.relplot(kind='scatter', data=median_df, x='pos_unique', y='log_lik_ratio',
                    hue='phenotype', col='Genotype', aspect=4)
    g.set(xlabel="CHR", xticks=median_df['sep'].unique(
    ), xticklabels=median_df['CHR'].unique())
    g.savefig('scatterplot_median_all_snp_ratio_uniquepos_colorph.png')


def stat_linear_reg(df):
    from sklearn.linear_model import LinearRegression
    # create linear regression object
    mlr = LinearRegression()
    df['phenotype'] = pd.get_dummies(df['phenotype']).iloc[:, 0]
    df['Genotype'].replace(
        {'0/0': 0, '0/1': 1, '1/1': 2}, inplace=True)

    # Linear model
    results = pd.DataFrame(
        index=['intercept', 'coeffs [Genotype, Phenotype]', 'r2'])
    for snp in df['SNP'].unique():
        snp_df = df[df['SNP'] == snp].copy()
        x, y = snp_df[['log_lik_ratio', 'phenotype']], snp_df['Genotype']
        mlr.fit(x, y)
        results[snp] = [mlr.intercept_, mlr.coef_, mlr.score(x, y)]
        try:
            aov = pg.anova(data=snp_df, dv='log_lik_ratio',
                           between=['Genotype', 'phenotype'])
            print(aov)
            results[snp] = pd.concat([results, aov.set_index(
                'Source')[['F']].rename(columns={'F': snp})])
        except:
            pass
        if (mlr.score(x, y) > 0.3) & (0 not in mlr.coef_):
            g = sns.lmplot(data=snp_df, x='Genotype',
                           y='log_lik_ratio', hue='phenotype')
            g.set(xticks=[0, 1, 2], xticklabels=[
                  '0/0', '0/1', '1/1'])

    return results.T


def main(dir):
    all = gather_dfs_fromdir(dir)
    # all = pd.read_csv('nano_genotyped_5b.csv')
    all = filtering_datas(all)

    # Insert phenotype variable
    all.loc[all['name'].str.contains('PROM1'), 'phenotype'] = 'Severe'
    all['phenotype'].fillna('Mild', inplace=True)
    all.loc[all['num_motifs'] == 1, 'distance_cpg_snp'] = abs(
        all['pos'] - all['start'])

    # Median over the read_name
    median_df = all.groupby(
        ['Genotype', 'SNP', 'name', 'phenotype', 'read_name']).median().reset_index()

    # STATS
    results_lm = stat_linear_reg(median_df)
    print('\n# Linear regression results\n')
    print(results_lm.sort_index().to_markdown())

    # PLOTS
    for chr in median_df['CHR'].sort_values().unique():
        violinplot(median_df[median_df['CHR'] == chr])
        ridgeplot(median_df[median_df['CHR'] == chr])
    # scatter_with_unique_pos(median_df)
    # scatter_with_unique_pos2(median_df)


if __name__ == '__main__':
    main(dir)

    # # Per Bins ?
    # for chr in selected['CHR'].unique():
    #     pos_chr = selected.loc[selected['CHR'] == chr, 'pos'].copy()
    #     selected.loc[selected['CHR'] == chr, 'pos_bin'] = pd.qcut(
    #         pos_chr.astype(int), q=10, duplicates='drop')
    #
    #
    # # When we categorize the SNPs per bins
    # g_bin = sns.catplot(data=selected, kind='violin',
    #                     x='pos_bin', y='log_lik_ratio', hue='Genotype',
    #                     row='CHR', col='phenotype', aspect=4,
    #                     inner="quartile", linewidth=3,
    #                     sharex=False, sharey=False)
    # g_bin.savefig('violinplot_median_bin_SNPs_ratio_distribution_chr.png')
