import pandas as pd
import seaborn as sns
from scipy.stats import zscore
import matplotlib.pyplot as plt
import pingouin as pg
from sklearn.linear_model import LinearRegression
from scipy.stats import spearmanr
import numpy as np
import os
import sys
import matplotlib.patches as mpatches

file = 'Filtered_nano_bam_files_all_samples.csv'
file_snp = 'significant_hits_COVID19_HGI_A2_ALL_leave_23andme_20210607.txt'
dir_out='results_cpg_snp_analysis'
unit='cpg'


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


def outliers(df, thresh_zscore=3):
    # Calculation of z-score for outliers
    df_zscore = df.select_dtypes(include=['float']).apply(zscore)
    outliers = df_zscore[abs(df_zscore) > thresh_zscore].dropna(
        how='all').index.tolist()
    print(f'Ouliers rows: {len(outliers)}', flush=True)
    return df.drop(outliers, axis=0)


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


def count_filter(df, min_count=3, n_genotypes=2):
    new_df = df.copy()
    # minimal count
    snp_counts = new_df.groupby(['SNP', 'Genotype', 'Gen']).size(
                    ).unstack()
    snp_ls = snp_counts[snp_counts > min_count].dropna(
        thresh=n_genotypes).index.unique().tolist()
    cpg_counts = new_df.groupby(['cpg', 'Genotype', 'Gen']).size(
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


def run_stat(df, unit='', var='', measure='log_lik_ratio', suppl_title='', out_dir='', pval_spearman=0.01):
    # TODO: Add paired ttest/wilcoxon one
    # TODO: start analysis with phenotype
    df = df.copy()
    dict_dum = {x:i for i,x in enumerate(df[var].unique())}
    df[f'{var}_dum'] = df[var].replace(dict_dum)
    results = pd.DataFrame(index=df[unit].unique())
    for u in df[unit].unique():
        u_df = df[df[unit] == u].copy()

        # Counts
        count_dict = u_df[var].value_counts().rename(
            index=lambda s: 'Counts ' + s).to_dict()

        # Diff between means
        diff_means = u_df[u_df['Gen'] == 'ref'][measure].mean(
            ) - u_df[u_df['Gen'] == 'alt'][measure].mean()

        # Mann Whitney
        if isinstance(var, list):
            u_df['mwu_var'] = u_df[var].T.apply('&'.join)
            mwu = multiple_mann_whitney(u_df, 'mwu_var', measure)
        else:
            mwu = multiple_mann_whitney(u_df, var, measure)

        # Linear regression
        try:
            lm = pg.linear_regression(u_df[measure], u_df['Genotype_dum'])
        except:
            lm = pd.DataFrame()

        # Spearman correlation
        spear = spearmanr(u_df[measure], u_df[var])
        # spear = spearmanr(u_df[measure], u_df[var[0]])

        # Pairwise tests
        try:
            paired_t = pg.pairwise_ttests(dv=measure, between=var,
                data=u_df).set_index('PairTtest ' + paired_t['A'] + '-' + paired_t['B'])
            paired_w = pg.pairwise_ttests(dv=measure, between=var, data=u_df,
                parametric=False).set_index('PairWilcox ' + paired_t['A'] + '-' + paired_t['B'])
        except:
            paired_t = pd.DataFrame()
            paired_w = pd.DataFrame()

        # Results
        results.loc[u,'SNP'] =str(u_df['SNP'].unique())[2:-2]
        results.loc[u,count_dict.keys()] = count_dict.values()
        results.loc[u, 'diff_means_altVSref'] = diff_means
        results.loc[u, f'Spearman correlation p_value'] = spear[1]
        results.loc[u, f'Spearman correlation rho'] = spear[0]
        if not lm.empty:
            try:
                results.loc[u, [f'LM intercept',
                                f'LM coef', 'r2']] = lm['coef'].tolist(
                                    ) + lm['r2'].tolist()[:1]
            except:
                pass
        if not mwu.empty:
            results.loc[u, mwu.index] = mwu['p-val']
        if not paired_t.empty:
            results.loc[u, paired_t.index] = paired_t['p-unc']
        if not paired_w.empty:
            results.loc[u, paired_w.index] = paired_w['p-unc']
    results = results.dropna(how='all').dropna(how='all', axis=1)
    results['minus_log10'] = np.log10(results[
        f'Spearman correlation p_value'].astype(float)) * (-1)
    results = results.reset_index().rename(columns={'index':unit})
    results['cutoff'] = -1 * np.log10(pval_spearman/results.shape[0])
    print(f'\n SNP ABOVE CUTOFF {suppl_title}')
    print(results[results['minus_log10'] > results['cutoff']][unit].tolist(),
          flush=True)
    results.sort_values(f'Spearman correlation p_value').to_csv(
        f'{out_dir}/Stat_Analysis_{measure}VS{var}_per_{unit}_{suppl_title}.csv')
    return results


def description(df, title='Description.csv'):
    pd.DataFrame(df.dtypes).T.to_csv(title, mode='w', header=True)
    print('types', flush=True)
    pd.DataFrame(df.nunique()).T.to_csv(title, mode='a', header=False)
    print('nunique', flush=True)
    pd.DataFrame(df.count()).T.to_csv(title, mode='a', header=False)
    print('count', flush=True)
    pd.DataFrame(df.max()).T.to_csv(title, mode='a', header=False)
    print('max', flush=True)
    pd.DataFrame(df.min()).T.to_csv(title, mode='a', header=False)
    print('min', flush=True)
    pd.DataFrame(df.mean()).T.to_csv(title, mode='a', header=False)
    print('mean', flush=True)
    pd.DataFrame(df.median()).T.to_csv(title, mode='a', header=False)
    print('median', flush=True)
    return final_desc


def number_catvalues_per_var(df, var, lim_values=5):
    per_var = pd.DataFrame(columns=df[var].unique())
    for col in df.select_dtypes(exclude='number').drop(var, axis=1).columns:
        unique_val = df[col].sort_values().astype(str).unique().tolist()
        col_num = pd.DataFrame(df.groupby(var).nunique().unstack().loc[col], columns=[col]).T
        per_var = pd.concat([per_var, col_num])
        if df[col].nunique() < lim_values:
            per_var = pd.concat([per_var, df.groupby([col, var]).size().unstack(
                ).round(0).rename(index=lambda x: f'{col}: ' + x)])
    return per_var

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


def spearman_correlation_plot(stat_df, unit='cpg', n_site=2, out_dir='', title_supp=''):
        stat = stat_df.copy()
        stat[['CHR', 'POS']] = stat[unit].str.split(':', expand=True)[[
            0, 1]].astype(int)
        stat = stat.sort_values('CHR')
        stat['CHR'] = stat['CHR'].astype('category')
        best_cpg = stat[stat['cutoff'] < stat['minus_log10']].sort_values(
            'minus_log10', ascending=False).groupby('CHR').head(2)['cpg'].tolist()
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
        g.savefig(f'{out_dir}/minuslog10_Spearman_pvalue_{unit}_{title_supp}.png')


def boxplot_customized(df, x_var, y_var, hue_var=None, dict_colors='tab10', width_var=None, hatch_var=None, ax=None, replace_val={}):
    dict_colors = dict_colors.copy()
    if hue_var:
        if not isinstance(dict_colors, dict):
            dict_colors = {hue_val: sns.color_palette(dict_colors)[n] for n, hue_val in enumerate(df[hue_var].unique())}

    if hatch_var:
        hatch_list = ['/', '|', '-', '\ ', '+', 'x', 'o', 'O', '.', '*']
        assert len(hatch_list) > df[hatch_var].nunique(), 'Hatching list not big enough, recombination to do'
    else:
        hatch_list = ['', '', '']
    print(dict_colors)
    g2 = sns.boxplot(data=df, x=x_var, y=y_var, orient='v', hue=hue_var,
                     palette=dict_colors, ax=ax)

    dict_var = {'Hue': hue_var, 'Width': width_var, 'Hatch': hatch_var}
    blou = df.copy()
    for key, var in dict_var.items():
        n = 0
        if var:
            for val in blou[var].unique():
                blou.loc[blou[var] == val, key] = n
                n = n+1
        else:
            blou[key] = 0

    list_var = set([x_var] + list(dict_var.keys()) + [val for val in dict_var.values() if val is not None])
    art_df = blou[list_var].value_counts().reset_index()
    if hue_var:
        art_df.sort_values([x_var, hue_var], inplace=True)
    else:
        art_df.sort_values([x_var], inplace=True)

    # assert len(art_df.index) == len(g2.artists), 'One or several variables are affected to non-existent patches'
    art_df.index = g2.artists
    for art in g2.artists:
        art.set_hatch(hatch_list[int(art_df.loc[art, 'Hatch'])])
        art.set_linewidth((art_df.loc[art, 'Width']+1)*2)

    if replace_val:
        for key, val in replace_val.items():
            dict_colors[val] = dict_colors.pop(key)

    # print(art_df)
    contour_line = sns.color_palette('dark')[-3]
    print(dict_colors.items())
    if hue_var == hatch_var == width_var:
        legend_hue = g2.legend(handles=[mpatches.Patch(facecolor=val, label=key,
                               lw=(art_df.loc[art_df[hue_var] == key, 'Width'][0]+1)*2,
                               ec=contour_line,
                               hatch=hatch_list[int(art_df.loc[art_df[hue_var] == key, 'Hatch'][0])]
                               ) for key, val in dict_colors.items()],
                               title=hue_var, loc=2)

    elif hue_var == hatch_var:
        legend_hue = g2.legend(handles=[mpatches.Patch(facecolor=val, label=key,
                               lw=2, ec=contour_line,
                               hatch=hatch_list[int(art_df.loc[art_df[hue_var] == key, 'Hatch'][0])]
                               ) for key, val in dict_colors.items()],
                               title=hue_var, loc=2)
        g2.add_artist(legend_hue)
        try:
            legend_width = g2.legend(handles=[mpatches.Patch(facecolor='white', label=key,
                                     lw=(art_df.loc[art_df[width_var] == key, 'Width'][0]+1)*2,
                                     ec=contour_line
                                     ) for key, val in art_df[[width_var, 'Width']].value_counts().index.tolist()],
                                     title=width_var, loc=1)
        except:
            pass

    elif hue_var == width_var:

        legend_hue = g2.legend(handles=[mpatches.Patch(facecolor=val, label=key,
                               lw=(art_df.loc[art_df[hue_var] == key, 'Width'][0]+1)*2,
                               ec=contour_line
                               ) for key, val in dict_colors.items()],
                               title=hue_var, loc=2)
        g2.add_artist(legend_hue)
        try:
            legend_hatch = g2.legend(handles=[mpatches.Patch(facecolor='white',
                                     label=key, lw=2, ec=contour_line,
                                     hatch=hatch_list[int(art_df.loc[art_df[hatch_var] == key, 'Hatch'][0])]
                                     ) for key, val in art_df[[hatch_var, 'Hatch']].value_counts().index.tolist()],
                                     title=hatch_var, loc=1)
        except:
            pass
    elif hatch_var == width_var:
        legend_hue = g2.legend(handles=[mpatches.Patch(facecolor=val, label=key,
                               lw=2, ec=contour_line
                               ) for key, val in dict_colors.items()],
                               title=hue_var, loc=2)
        g2.add_artist(legend_hue)
        try:
            legend_hatch = g2.legend(handles=[mpatches.Patch(facecolor='white',
                                     label=key, ec=contour_line,
                                     lw=(art_df.loc[art_df[width_var] == key, 'Width'][0]+1)*2,
                                     hatch=hatch_list[int(art_df.loc[art_df[hatch_var] == key, 'Hatch'][0])]
                                     ) for key, val in art_df[[hatch_var, 'Hatch']].value_counts().index.tolist()],
                                     title=hatch_var, loc=1)
        except:
            pass
    else:
        legend_hue = g2.legend(handles=[mpatches.Patch(facecolor=val, label=key,
                               lw=2, ec=contour_line
                               ) for key, val in dict_colors.items()],
                               title=hue_var, loc=2)
        g2.add_artist(legend_hue)
        try:
            legend_hatch = g2.legend(handles=[mpatches.Patch(facecolor='white',
                                     label=key, lw=2, ec=contour_line,
                                     hatch=hatch_list[int(art_df.loc[art_df[hatch_var] == key, 'Hatch'][0])]
                                     ) for key, val in art_df[[hatch_var, 'Hatch']].value_counts().index.tolist()],
                                     title=hatch_var, loc=1)
            g2.add_artist(legend_hatch)
        except:
            pass
        try:
            legend_width = g2.legend(handles=[mpatches.Patch(facecolor='white',
                                     label=key, lw=(art_df.loc[art_df[width_var] == key, 'Width'][0]+1)*2,
                                     ec=contour_line
                                     ) for key, val in art_df[[width_var, 'Width']].value_counts().index.tolist()],
                                     title=width_var, loc=3)
        except:
            pass
    return g2


def setup_customizedboxplot_cpg_analysis(cpg_df, dir_out):
    cpg = str(cpg_df['cpg'].unique())[2:-2]
    snp = str(cpg_df['SNP'].unique())[2:-2]
    colors_hom = {'0/0':'#1678F5', '0/1':'#2aa69a', '1/1':'#3ED43E'} # Genotype colors
    colors_het = {0:'#1678F5', 1:'#3ED43E'} # Haplotype colors
    colors_phen = {'Mild':'#f0f0f5', 'Severe':'#c2c2d6'} # Phenotype colors
    cpg_df.sort_values(['phenotype', 'Genotype', 'Gen'], inplace=True)
    alls = snp.split(':')[-2:]
    replace_val = {'Genotype': {'0/0': alls[0]+'/'+alls[0],
                                '0/1': alls[0]+'/'+alls[1],
                                '1/1': alls[1]+'/'+alls[1]},
                   'Gen': {0: alls[0], 1: alls[1]}}
    fig, ax = plt.subplots(2, 3, figsize=(17,10))
    boxplot_customized(cpg_df, 'phenotype', 'log_lik_ratio', hue_var='phenotype',
                                 dict_colors=colors_phen,
                                 hatch_var='phenotype', ax=ax[0,0])
    ax[0,0].set(title='Phenotype')
    boxplot_customized(cpg_df, 'Genotype', 'log_lik_ratio', hue_var='Genotype',
                                 dict_colors=colors_hom,
                                 width_var=None, ax=ax[0,1], replace_val=replace_val['Genotype'])
    ax[0,1].set(title='Genotype', xticklabels=replace_val['Genotype'].values())
    boxplot_customized(cpg_df, 'phenotype', 'log_lik_ratio', hue_var='Genotype',
                                 dict_colors=colors_hom,
                                 hatch_var='phenotype', ax=ax[0,2], replace_val=replace_val['Genotype'])
    ax[0,2].set(title='Genotype X Phenotype')
    boxplot_customized(cpg_df[cpg_df['Genotype'] == '0/1'], 'phenotype', 'log_lik_ratio',
                                  hue_var='phenotype', dict_colors=colors_phen,
                                  hatch_var='phenotype', width_var=None, ax=ax[1,0])
    ax[1,0].set(title='Heterozygous Phenotype')
    boxplot_customized(cpg_df[cpg_df['Genotype'] == '0/1'], 'Gen', 'log_lik_ratio',
                                  hue_var='Gen', dict_colors=colors_het,
                                  hatch_var=None, width_var=None, ax=ax[1,1], replace_val=replace_val['Gen'])
    ax[1,1].set(title='Heterozygous Haplotype', xticklabels=replace_val['Gen'].values())
    boxplot_customized(cpg_df[cpg_df['Genotype'] == '0/1'], 'phenotype', 'log_lik_ratio',
                                  hue_var='Gen', dict_colors=colors_het,
                                  hatch_var='phenotype', width_var=None, ax=ax[1,2], replace_val=replace_val['Gen'])
    ax[1,2].set(title='Heterozygous Haplotype X Phenotype')
    plt.suptitle(f'CpG {cpg} associated with SNP {snp}')
    fig.savefig(f'{dir_out}/Multiplots_{cpg}.png')

# MAIN
def main(file, dir_out='FROZEN_results_cpg_snp_analysis/special_plots', unit='cpg'):
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)
    # TODO: Analysis with yaml software, moved in the result folder (with date and info for analysis)
    # Selection of the SNPs
    # TODO: Check selection of random/control SNPs with TOM
    snp_ls = select_SNP_per_pvalue(file_snp, pval_col='all_inv_var_meta_p',
        dist_bp=500000)
    gen_ls = ['0/0', '0/1', '1/1']

    # Opening file
    all_df = pd.read_csv(file)
    all_df = all_df[(all_df['SNP'].isin(snp_ls)) & (all_df['Gen'] != 'other')
                    & (all_df['Genotype'].isin(gen_ls))]

    # # READNAME
    # for var in ['read_name', 'cpg']:
    #     for unit in ['SNP', 'name']:
    #     # unit='cpg'
    #     # var='read_name'
    #         try:
    #             print(unit)
    #             count = all_df.groupby(['phenotype','CHR', unit]).count()[var].reset_index()
    #             count.to_csv(f'FROZEN_results_cpg_snp_analysis/special_plots/Count_{var}_per_{unit}.csv')
    #             nunique = all_df.groupby(['phenotype','CHR', unit]).nunique()[var].reset_index()
    #             nunique.to_csv(f'FROZEN_results_cpg_snp_analysis/special_plots/Nunique_{var}_per_{unit}.csv')
    #             unit_count = count[unit].nunique()
    #             if unit == 'SNP':
    #                 g_count = sns.catplot(data=count, y=unit, x=var, col='phenotype',
    #                                       kind='box', sharex=False, sharey=False)
    #                 g_nunique = sns.catplot(data=nunique, y=unit, x=var, col='phenotype',
    #                                       kind='box', sharex=False, sharey=False)
    #             else:
    #                 g_count = sns.catplot(data=count, y=unit, x=var, col='phenotype',
    #                                       kind='bar', aspect=1.5, height=15,
    #                                       sharex=False, sharey=False)
    #                 g_nunique = sns.catplot(data=nunique, y=unit, x=var, col='phenotype',
    #                                       kind='bar', aspect=1.5, height=15,
    #                                       sharex=False, sharey=False)
    #             g_count.savefig(f'FROZEN_results_cpg_snp_analysis/special_plots/Counts_{var}_per_{unit}.png')
    #             g_nunique.savefig(f'FROZEN_results_cpg_snp_analysis/special_plots/Nuniques_{var}_per_{unit}.png')
    #         except Exception as err:
    #             print(err)


    # Dataset description
    # number_catvalues_per_var(all_df.drop(['Allele', 'ref', 'alt'], axis=1),
    #     'name', lim_values=5).to_csv('Numbervalues_persamples.csv')

    # Median over the read_name with same Allele calling
    median_df = all_df.groupby(['phenotype', 'name', 'CHR', 'cpg',
                                'SNP', 'Genotype', 'Gen']).median().reset_index()
    median_df = median_df.sort_values(by=['Gen'], ascending=False)
    median_df = median_df.sort_values(by=['CHR', 'SNP', 'Genotype', 'phenotype'])
    #
    # # QUESTION: What should we do with the ALT calls in '0/0' and REF in '1/1'?
    #
    # # Filter
    # # median_new = outliers(median_df, thresh_zscore=3)
    # # median_new = count_filter(median_new, min_count=5, n_genotypes=2)
    #
    # # STATS ON THE OVERALL DF
    # # stat = run_stat(median_df, unit=unit, var='Genotype',
    # #                 measure='log_lik_ratio', out_dir=dir_out)
    # # spearman_correlation_plot(stat[stat['minus_log10'].isna()==False], unit=unit, n_site=2, out_dir=dir_out)
    # # print(stat.shape)
    #
    # # STATS FOR HETEROZYGOTES
    # # stat_het = run_stat(median_df[median_df['Genotype'] == '0/1'], unit=unit,
    # #     measure='log_lik_ratio', var='Gen', suppl_title='Het_only', out_dir=dir_out)
    # # # TODO: Check if you should not just NaN the rows with counts too low for one Genotype
    # # del stat_het
    #
    # # STATS PER PHENOTYPE
    # # stat_severe = run_stat(median_df[median_df['phenotype'] == 'Severe'],
    # #                         unit=unit, var='Genotype',
    # #                         measure='log_lik_ratio', out_dir=dir_out, suppl_title='Severe_phenotype')
    # # spearman_correlation_plot(stat_severe[stat_severe['minus_log10'].isna()==False], unit=unit, n_site=2, out_dir=dir_out, title_supp='Severe')
    #
    # # stat_mild = run_stat(median_df[median_df['phenotype'] == 'Mild'],
    # #                         unit=unit, var='Genotype',
    # #                         measure='log_lik_ratio', out_dir=dir_out, suppl_title='Mild_phenotype')
    # # spearman_correlation_plot(stat_mild[stat_mild['minus_log10'].isna()==False], unit=unit, n_site=2, out_dir=dir_out, title_supp='Mild')
    #
    # # ANALYSIS MILD VS SEVERE
    # # merge = pd.merge(stat_mild, stat_sev, on=['cpg', 'SNP'], suffixes=['_mild', '_sev'])
    # # index_counts = merge.filter(regex='Counts*')[merge.filter(regex='Counts*') < 3].dropna(thresh=6).index
    # # merge.drop(index_counts, inplace=True)
    # # merge.loc[merge['Spearman correlation p_value_mild'] < 1e-5, 'log_mild'] = 'High'
    # # merge.loc[merge['Spearman correlation p_value_sev'] < 1e-5, 'log_sev'] = 'High'
    # # merge.loc[merge['Spearman correlation p_value_sev'] > 1e-5, 'log_sev'] = 'Low'
    # # merge.loc[merge['Spearman correlation p_value_mild'] > 1e-5, 'log_mild'] = 'Low'
    # # cpg_highsev_lowmild = merge[(merge['log_sev'] == 'High') & (merge['log_mild'] == 'Low')]['cpg'].unique()
    # # cpg_lowsev_highmild = merge[(merge['log_sev'] == 'Low') & (merge['log_mild'] == 'High')]['cpg'].unique()
    #
    # # READING STAT FILES FOR PLOTS
    # # cols = ['Counts 0/0', 'Counts 1/1', 'Counts 0/1', 'Spearman correlation p_value', 'cpg']
    # # stat = pd.read_csv('FROZEN_results_cpg_snp_analysis/Stat_Analysis_log_lik_ratioVSGenotype_per_cpg_.csv', usecols=cols)
    # # stat_mild = pd.read_csv('FROZEN_results_cpg_snp_analysis/Stat_Analysis_log_lik_ratioVSGenotype_per_cpg_Mild_phenotype.csv')
    # # stat_sev = pd.read_csv('FROZEN_results_cpg_snp_analysis/Stat_Analysis_log_lik_ratioVSGenotype_per_cpg_Severe_phenotype.csv')
    #
    # # CPGs of interest
    # # cpgs_plot = stat[(stat['Counts 0/0'] > 3) & (stat['Counts 1/1'] > 3) & (stat['Counts 0/1'] > 3) & (stat['Spearman correlation p_value'] < 1e-5)]['cpg'].unique()
    # # highest p-val for highmildhighsev : ['17:46065976:1', '12:112942465:1', '21:33229986:3'] # highest rho

    dict_cpgs = {'INTEREST': ['17:46768336:3', '6:33069193:2'],
                 'HIGHmildHIGHsev': ['12:112925744:1', '17:46065976:1', '21:33229986:3'],
                 'EXTRA': ['3:45848456:1', '21:33226777:1', '21:33242527:1']}
    # dict_cpgs = median_df['cpg'][:10].to_dict('list')
    # TODO: get the number of bix tht we will have in each condition to set up the conditions automatically on i for pattern and line

    for supp_title, list in dict_cpgs.items():

    #     # df = median_df[median_df['cpg'].isin(list)]
    #     # GENERAL PLOT GENOTYPE
    #     # g = sns.catplot(data=df,
    #     #                 y='log_lik_ratio', col='cpg', col_wrap=5,
    #     #                 x='Genotype', orient='v', kind= 'box',
    #     #                 height=6, aspect=0.9, hue='phenotype',
    #     #                 sharex=False, sharey=False)
    #     # g.savefig(f'{dir_out}/Boxplot_median_ratio_GenPhen_all_{supp_title}.png')
    #     #
    #     # # GENERAL PLOT HETEROZYGOTES
    #     # g = sns.catplot(data=df[df['Genotype'] == '0/1')],
    #     #                 y='log_lik_ratio', col='cpg', col_wrap=5,
    #     #                 x='Gen', orient='v', kind= 'box',
    #     #                 height=6, aspect=0.9, hue='phenotype',
    #     #                 sharex=False, sharey=False, )
    #     # g.savefig(f'{dir_out}/Boxplot_median_ratio_GenPhen_all_{supp_title}_het.png')
    #
        # INDIVIDUAL PLOTS
        for cpg in list:
            cpg_df = median_df[median_df['cpg'] == cpg].copy()
            snp = str(cpg_df['SNP'].unique().tolist())[2:-2]
            cpg_df['Gen'].replace({'alt': 1, 'ref':0}, inplace=True)
            print(snp)
            setup_customizedboxplot_cpg_analysis(cpg_df, dir_out=dir_out)
            # try:
                # Individual per genotype
                # 3 way genotype
                # g1 = sns.catplot(data=cpg_df, y='log_lik_ratio',
                #                 x='Genotype', orient='v', kind= 'box',
                #                 height=6, aspect=0.9,
                #                 sharex=False, sharey=False)
                # g1.set(title=f'CpG {cpg} associated with SNP {snp}')
                # g1.savefig(f'{dir_out}/Boxplot_median_cpg_ratio_Gen_{supp_title}_{cpg}.png')
                # 2 way phenotype
                # g2 = sns.catplot(data=cpg_df, y='log_lik_ratio',
                #                 x='phenotype', orient='v', kind= 'bar',
                #                 height=6, aspect=0.9,
                #                 sharex=False, sharey=False, cmap='viridis')
                # g2.set(title=f'CpG {cpg} associated with SNP {snp}')
                # g2.savefig(f'{dir_out}/Boxplot_median_cpg_ratio_Phen_{supp_title}_{cpg}.png')
                # 6 way
                # g3 = sns.catplot(data=cpg_df, y='log_lik_ratio',
                #                 x='Genotype', orient='v', kind= 'box',
                #                 height=6, aspect=0.9, hue='phenotype',
                #                 sharex=False, sharey=False)
                # g3.set(title=f'CpG {cpg} associated with SNP {snp}')
                # g3.savefig(f'{dir_out}/Boxplot_median_cpg_ratio_GenPhen_{supp_title}_{cpg}.png')

                # Individual for heterozygotes
                # 4 way
                # g4 = sns.catplot(data=cpg_df[cpg_df['Genotype'] == '0/1'],
                #                 y='log_lik_ratio', x='Gen', kind= 'box',
                #                 height=6, aspect=0.9, hue='phenotype',
                #                 sharex=False, sharey=False)
                # g4.set(title=f'CpG {cpg} associated with SNP {snp}', xlabel='Allele')
                # g4.savefig(f'{dir_out}/Boxplot_median_cpg_ratio_GenPhen_{supp_title}_{cpg}_het.png')

                # # 2 way allele
                # g5 = sns.catplot(data=cpg_df[cpg_df['Genotype'] == '0/1'],
                #                  y='log_lik_ratio', x='Gen', kind= 'box',
                #                  height=6, aspect=0.9,
                #                  sharex=False, sharey=False)
                # g5.set(title=f'CpG {cpg} associated with SNP {snp}', xlabel='Allele')
                # g5.savefig(f'{dir_out}/Boxplot_median_cpg_ratio_Gen_{supp_title}_{cpg}_het.png')

                # 2 way phenotype
            #     g6 = sns.catplot(data=cpg_df[cpg_df['Genotype'] == '0/1'],
            #                     y='log_lik_ratio', x='phenotype', kind= 'box',
            #                     height=6, aspect=0.9,
            #                     sharex=False, sharey=False, cmap='viridis')
            #     g6.set(title=f'CpG {cpg} associated with SNP {snp}')
            #     g6.savefig(f'{dir_out}/Boxplot_median_cpg_ratio_Phen_{supp_title}_{cpg}_het.png')
            #
            # except Exception as err:
            #     print(f'ERROR WITH cpg {cpg} ', err)

# workflow slide and stat with number of loci cpg site in total average per loci

if __name__ == '__main__':
    main(file)

# descriptive statistics- number of cpg per samples (boxplot per chr?)
# slide with  which snp we selected first and all the parameters we selected for the analysis

# p-values et betas of the SNP on browser ???
# look at each locus that are above the cutoff
# is the cpg around promoter of a gene
# ensembl element features
# transform to normal distrib 9inverse norm)
# lm(medianlog ~ snp * clinical, data)
# lm(invnorm(medianlog) ~ snp + clinical + snp * clinical, data)
# lm(invnorm(median_log_hap) ~ allele * clinical , haplotype_data)

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
    # scatter_with_unique_cpg(median_df, huevar='Gen',
    #                         colvar=None, xvar='cpg', yvar='log_lik_ratio')
    for chr in median_df['CHR'].sort_values().unique():
        violinplot(median_df[median_df['CHR'] == chr], title_supp=chr, yvar='cpg',
                   xvar='log_lik_ratio', huevar='Genotype', colvar='phenotype')
        violinplot(median_df[median_df['CHR'] == chr], title_supp=chr, yvar='SNP',
                   xvar='log_lik_ratio', huevar='Genotype', colvar='phenotype')
        # ridgeplot(median_df[median_df['CHR'] == 1], chr=1, rowvar='SNP',
        #           huevar='Genotype', colvar='phenotype', var='log_lik_ratio')
# TODO: Fix ridge PLOTS not working with cpg, bigger font size
# TODO: Fix Violin plot (all the same when coming back from the cluster)
