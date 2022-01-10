import pandas as pd


### PLOTS ###
def violinplot(df, title_supp='', yvar='cpg', xvar='log_lik_ratio', huevar='Genotype', colvar='phenotype', out_dir=''):
    n_snp = df[xvar].nunique()
    print(n_snp)
    if n_snp < 10:
        g_chr = sns.catplot(data=df, kind='violin',
                            y=yvar, x=xvar, hue=huevar,
                            col=colvar, height=6,
                            aspect=0.9, orient='v',
                            inner="quartile", linewidth=1.5,
                            sharex=False, sharey=False)
    elif n_snp < 50:
        g_chr = sns.catplot(data=df, kind='violin',
                            y=yvar, x=xvar, hue=huevar,
                            col=colvar, height=25,
                            aspect=.7, orient='v',
                            inner="quartile", linewidth=1.5,
                            sharex=False, sharey=False)
    elif n_snp < 500:
        g_chr = sns.catplot(data=df, kind='violin',
                            y=yvar, x=xvar, hue=huevar,
                            col=colvar, height=45,
                            aspect=.15, orient='v',
                            inner="quartile", linewidth=1,
                            sharex=False, sharey=False)
    else:
        g_chr = sns.catplot(data=df, kind='violin',
                            y=yvar, x=xvar, hue=huevar,
                            col=colvar, height=45,
                            aspect=.1, orient='v',
                            inner="quartile", linewidth=1,
                            sharex=False, sharey=False)
    plt.title(f'{title_supp}')
    g_chr.savefig(
        f'{out_dir}/violinplot_median_all_{yvar}_{xvar}_distribution_{title_supp}.png')


def ridgeplot(df, chr='', rowvar='cpg', huevar='Genotype', colvar='phenotype', var='log_lik_ratio', out_dir=''):
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
            f'{out_dir}/rigdeplot_median_all_{rowvar}_{var}_distribution_chr{chr}_color{huevar}.png')


def scatter_with_unique_cpg(median_df, huevar='Genotype', colvar='phenotype', xvar='cpg', yvar='log_lik_ratio', out_dir=''):
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
        f'{out_dir}/scatterplot_median_all_{yvar}_unique{xvar}_color{huevar}.png')


def linear_reg_plot(df, var='', unit='', plot=False, title_suppl=''):
        # Linear regression
    g = sns.lmplot(data=u_df, x='Genotype',
                   y=var, hue='phenotype', legend=False)
    g.set(xticks=[0, 1, 2], xticklabels=[
          '0/0', '0/1', '1/1'])
    plt.legend(labels=['Mild', 'Severe'], loc='upper left')
    plt.title(unit + ' ' + u)


def spearman_correlation_plot(stat_df, unit='cpg', n_site=2, out_dir=''):
        stat = stat_df.copy()
        stat[['CHR', 'POS']] = stat[unit].str.split(':', expand=True)[[
            0, 1]].astype(int)
        stat = stat.sort_values('CHR')
        stat['CHR'] = stat['CHR'].astype('category')
        best_cpg = stat_gen[stat_gen['cutoff'] < stat_gen['minus_log10']].sort_values('minus_log10', ascending=False).groupby('CHR').head(2)['cpg'].tolist()
        stat.loc[(stat[unit].isin(best_cpg)) &
            (stat['minus_log10'] > stat['cutoff']),#
            f'{unit}_best'] = stat.loc[(stat[unit].isin(best_cpg)) &
                        (stat['minus_log10'] > stat['cutoff']), 'cpg']
        stat.loc[stat[f'{unit}_best'].isna(), f'{unit}_best'] = ' '
        print(stat[(stat[unit].isin(best_cpg)) &
                    (stat['minus_log10'] > stat['cutoff'])])
        g = sns.FacetGrid(stat, aspect=4, height=4, palette='Spectral',
                          margin_titles=True)
        g.map(sns.lineplot,unit, 'cutoff', hue=None)
        g.map_dataframe(sns.scatterplot, unit, 'minus_log10', hue='CHR',
            legend=True)
        g.set(xlabel="CHR", xticks=stat.groupby(['CHR']).last()[unit].unique(),
            xticklabels=stat['CHR'].unique())
        try:
            for row in stat.itertuples():
                g.axes[0,0].text(row.cpg, row.minus_log10 + 0.5, row.cpg_best,
                    horizontalalignment='left')
        except:
            print('couldn print best points')
        g.savefig(f'{out_dir}/minuslog10_Spearman_pvalue_{unit}.png')


def spearmanRho_diffmean_plot(stat, unit='cpg', hue_var='CHR', col_var=None, out_dir=''):
    stat = stat.copy()
    stat = stat.reset_index().rename(columns={unit:unit})
    stat['CHR'] = stat[unit].str.split(':', expand=True)[0].astype('category')
    g1 = sns.relplot(kind='scatter', data=stat, y='diff_means_altVSref',
                     x='Spearman correlation Gen rho', hue=hue_var, col=col_var, col_wrap=3)
    g1.savefig(f'{out_dir}/Diff_meansVSrho_{unit}_heterozygotes_sepCHR.png')


def main():
    stat_gen = pd.read_csv('Stat_Analysis_log_lik_ratioVSGenotype_per_cpg_.csv', index_col=0)
    spearman_correlation_plot(stat[stat['minus_log10'].isna()==False], unit=unit, n_site=2, out_dir=dir_out)

    stat_het = pd.read_csv('Stat_Analysis_log_lik_ratioVSGen_per_cpg_Het_only.csv', index_col=0)
    stat_het['cut_log'] = pd.cut(x=stat_het['minus_log10'], bins=4)
    high_counts = stat_het.filter(
        regex='Counts*')[stat_het.filter(regex='Counts*')>5].dropna().index
    spearmanRho_diffmean_plot(stat_het.loc[high_counts], unit='cpg', col_var='SNP', out_dir=dir_out)


    # PLOTS
    # scatter_with_unique_cpg(median_df, huevar='Gen',
    #                         colvar=None, xvar='cpg', yvar='log_lik_ratio')

    # Violinplot only for the cpg
    # TODO: Fix Violin plot (all the same when coming back from the cluster)
    for snp in snp_ls:
        snp_df = median_new[median_new['SNP'] == snp].copy()
        print(snp_df.shape)
        try:
            g = sns.catplot(data=snp_df, y='log_lik_ratio',
                            x='Genotype', orient='v', kind= 'violin',
                            height=6, aspect=0.9, hue='phenotype',
                            sharex=False, sharey=False)
            g.savefig(f'{dir_out}/violinplot_median_all_ratio_GenPhen_distribution_{snp}.png')
            g1 = sns.catplot(data=snp_df,
                            y='log_lik_ratio', x='Gen', kind= 'swarm',
                            height=6, aspect=0.9, hue='phenotype',
                            sharex=False, sharey=False)
            g1.savefig(f'{dir_out}/swarmplot_median_all_ratio_GenPhen_distribution_{snp}.png')
        except Exception as err:
            print(f'ERROR WITH SNP {snp} ', err)
