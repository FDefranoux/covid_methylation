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
import socket

if 'Fanny' in socket.gethostname():
    abs_dir = '/home/fanny/Work/EBI/covid_nanopore'
else:
    abs_dir = '/nfs/research/birney/users/fanny/covid_nanopore'

file = f'{abs_dir}/FROZEN_Nov2021_cpg_snp_analysis/Filtered_nano_bam_files_all_samples.csv'
file_snp = f'{abs_dir}/FROZEN_Nov2021_cpg_snp_analysis/significant_hits_COVID19_HGI_A2_ALL_leave_23andme_20210607.txt'
dir_out= f'{abs_dir}/FROZEN_Nov2021_cpg_snp_analysis/results_cpg_snp_analysis'
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


def multiple_mann_whitney(df, col, measure):
    if isinstance(col, list):
        df['mwu_var'] = df[col].T.apply('_'.join)
        col.append('mwu_var')
    else:
        var = [col]
    mwu = pd.DataFrame()
    for c in col:
        for v in range(df[c].nunique() - 1):
            for v2 in range(v+1, df[c].nunique()):
                val = df[c].unique()[v]
                val2 = df[c].unique()[v2]
                eq_val = [(val.split('_')[n] == val2.split('_')[n])
                            for n in range(int(df[c].str.split(
                                '_').str.len().mean()))]
                if (c == 'mwu_var') & (eq_val.count(True) < len(eq_val) -1):
                    pass
                else:
                    try:
                        mwu_int = pg.mwu(df[df[c] == val][measure],
                                         df[df[c] == val2][measure]
                                         ).rename(index=lambda s: s +
                                         f' {"-".join(sorted([val, val2]))}')
                        mwu = mwu.append(mwu_int)
                    except:
                        pass
    return mwu


# def run_stat(df, unit='', var='', measure='log_lik_ratio', suppl_title='', out_dir=''):
#     # TODO: Joint analysis with phenotype
#     df = df.copy()
#
#     # Creation of dummy columns
#     if isinstance(var, list):
#         for v in var:
#             dict_dum = {x:i for i,x in enumerate(df[v].unique())}
#             df[f'{v}_dum'] = df[v].replace(dict_dum)
#             ls_dum = df.filter(regex='_dum').columns.tolist()
#     else:
#         dict_dum = {x:i for i,x in enumerate(df[var].unique())}
#         df[f'{var}_dum'] = df[var].replace(dict_dum)
#
#     results = pd.DataFrame(index=df[unit].unique())
#     for u in df[unit].unique():
#         u_df = df[df[unit] == u].copy()
#
#         # Counts
#         count_dict = u_df[var].value_counts().rename(
#             index=lambda s: 'Counts ' + s).to_dict()
#
#         # Diff between means
#         diff_means = u_df[u_df['Gen'] == 'ref'][measure].mean(
#             ) - u_df[u_df['Gen'] == 'alt'][measure].mean()
#
#         # Linear regression
#         try:
#             lm = pg.linear_regression(u_df[measure], u_df[f'{var}_dum'])
#         except:
#             lm = pd.DataFrame()
#
#     lm_dum = pd.DataFrame()
#     for dum in ls_dum:
#         pg.linear_regression(u_df[measure], u_df[dum]).set_index('names')
#         # Spearman correlation
#         if isinstance(var, list):
#             pg.pairwise_corr(u_df[[measure] + ls_dum], method='spearman')
#         data = pg.read_dataset('pairwise_corr').iloc[:, 1:]
#         spear = spearmanr(u_df[measure], u_df[var])
#
#         # spear = spearmanr(u_df[measure], u_df[var[0]])
#
#         # Pairwise tests
#         try:
#             paired_t = pg.pairwise_ttests(dv=measure, between=var,
#                 data=u_df).set_index('PairTtest ' + paired_t['A'] + '-' + paired_t['B'])
#             paired_w = pg.pairwise_ttests(dv=measure, between=var, data=u_df,
#                 parametric=False).set_index('PairWilcox ' + paired_t['A'] + '-' + paired_t['B'])
#         except:
#             paired_t = pd.DataFrame()
#             paired_w = pd.DataFrame()
#
#         # Results
#         results.loc[u,'SNP'] =str(u_df['SNP'].unique())[2:-2]
#         results.loc[u,count_dict.keys()] = count_dict.values()
#         results.loc[u, 'diff_means_altVSref'] = diff_means
#         results.loc[u, f'Spearman correlation p_value'] = spear[1]
#         results.loc[u, f'Spearman correlation rho'] = spear[0]
#         if not lm.empty:
#             try:
#                 results.loc[u, [f'LM intercept',
#                                 f'LM coef', 'r2']] = lm['coef'].tolist(
#                                     ) + lm['r2'].tolist()[:1]
#             except:
#                 pass
#         if not mwu.empty:
#             results.loc[u, mwu.index] = mwu['p-val']
#         if not paired_t.empty:
#             results.loc[u, paired_t.index] = paired_t['p-unc']
#         if not paired_w.empty:
#             results.loc[u, paired_w.index] = paired_w['p-unc']
#     results = results.dropna(how='all').dropna(how='all', axis=1)
#     results['minus_log10'] = np.log10(results[
#         f'Spearman correlation p_value'].astype(float)) * (-1)
#     results = results.reset_index().rename(columns={'index':unit})
#
#
#     results.sort_values(f'Spearman correlation p_value').to_csv(
#         f'{out_dir}/Stat_Analysis_{measure}VS{var}_per_{unit}_{suppl_title}.csv')
#     return results


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


def pval_plot(stat_df, xvar, pval_col, pval_cutoff=0.01, n_site=2, format_xaxis=True, out_dir='', title_supp=''):
        stat = stat_df.copy()
        stat[['CHR', 'POS']] = stat[xvar].str.split(':', expand=True)[[
            0, 1]].astype(int)
        stat = stat.sort_values('CHR')
        stat['CHR'] = stat['CHR'].astype('category')
        stat['minus_log10'] = -np.log10(stat[pval_col].astype(float))
        stat['cutoff'] = -np.log10(pval_cutoff/stat.shape[0])

        # Format X-axis
        n=0
        for chr in stat['CHR'].sort_values().unique():
            var_df = stat.loc[stat['CHR'] == chr, xvar].sort_values()
            for val in var_df.unique():
                stat.loc[(stat['cpg'] == val)
                    & (stat['CHR'] == chr), f'{xvar}_unique'] = n
                n = n + 1

        # Selection best snp
        best_cpg = stat[stat['cutoff'] < stat['minus_log10']].sort_values(
            'minus_log10', ascending=False).groupby('CHR').head(n_site)[
            'cpg'].tolist()
        stat.loc[(stat[xvar].isin(best_cpg)) &
            (stat['minus_log10'] > stat['cutoff']),#
            f'{xvar}_best'] = stat.loc[(stat[xvar].isin(best_cpg)) &
                        (stat['minus_log10'] > stat['cutoff']), 'cpg']
        stat.loc[stat[f'{xvar}_best'].isna(), f'{xvar}_best'] = ' '
        print(f'SNP ABOVE CUTOFF {title_supp}: {best_cpg}')

        # PLOT
        g = sns.FacetGrid(stat, aspect=4, height=4, palette='Spectral',
                          margin_titles=True)
        g.map(sns.lineplot, xvar, 'cutoff', hue=None)
        g.map_dataframe(sns.scatterplot, xvar, 'minus_log10', hue='CHR',
            legend=True)
        g.set(xlabel="CHR", xticks=stat.groupby(['CHR']).last()['POS'].unique(),
            xticklabels=stat['CHR'].unique())
        try:
            for row in stat.iterrows():
                row = row[1]
                g.axes[0,0].text(row[xvar], row.minus_log10 + 0.5, row.cpg_best,
                    horizontalalignment='left')
        except Exception as err:
            print('couldn\'t print best points', err)
        g.savefig(f'{out_dir}/minuslog10_{pval_col}_pvalue_{unit}_{title_supp}.png')


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
    snp = str(cpg_df['SNP'].unique().tolist())[2:-2]
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

    # Median over the read_name with same Allele calling
    median_df = all_df.groupby(['phenotype', 'name', 'CHR', 'cpg',
                                'SNP', 'Genotype', 'Gen']).median().reset_index()
    median_df = median_df.sort_values(by=['Gen'], ascending=False)
    median_df = median_df.sort_values(by=['CHR', 'SNP', 'Genotype', 'phenotype'])

    # QUESTION: What should we do with the ALT calls in '0/0' and REF in '1/1'?

    # STATS
    mann_whit = pd.DataFrame()
    # INDIVIDUAL PLOTS
    for cpg in median_df['cpg'].unique():
        cpg_df = median_df[median_df['cpg'] == cpg].copy()
        cpg_df['Gen'].replace({'alt': 1, 'ref':0}, inplace=True)
        mwu = multiple_mann_whitney(cpg_df, ['Genotype', 'phenotype'], 'log_lik_ratio')
        if not mwu.empty:
            mann_whit.loc[:, cpg] = mwu['p-val']
    # mann_whit[mann_whit < 0.05].dropna(how='all').dropna(how='all', axis=1)
    mann_whit = mann_whit.T.reset_index().rename(columns={'index': 'cpg'})
    mann_whit[['CHR', 'POS']] = mann_whit[unit].str.split(':', expand=True)[[
        0, 1]].astype(int)
    mann_whit.to_csv(f'{dir_out}/Mann_whitney_table.csv', index=False)
    pval_plot(mann_whit, 'cpg', 'MWU Mild-Severe', pval_cutoff=1, n_site=2, out_dir=dir_out, format_xaxis=True)

    stat = pd.read_csv('FROZEN_results_cpg_snp_analysis/Stat_Analysis_log_lik_ratioVSGenotype_per_cpg_.csv')
    pval_plot(stat, 'cpg', 'Spearman correlation p_value', pval_cutoff=1, n_site=2, title_supp='MannWhitney_GenPhen', out_dir=dir_out, format_xaxis=True)

    # Filter
    # median_new = outliers(median_df, thresh_zscore=3)
    # median_new = count_filter(median_new, min_count=5, n_genotypes=2)

    # STATS ON THE OVERALL DF
    # stat = run_stat(median_df, unit=unit, var='Genotype',
    #                 measure='log_lik_ratio', out_dir=dir_out)
    # spearman_correlation_plot(stat[stat['minus_log10'].isna()==False], unit=unit, n_site=2, out_dir=dir_out)
    # print(stat.shape)

    # STATS FOR HETEROZYGOTES
    # stat_het = run_stat(median_df[median_df['Genotype'] == '0/1'], unit=unit,
    #     measure='log_lik_ratio', var='Gen', suppl_title='Het_only', out_dir=dir_out)
    # # TODO: Check if you should not just NaN the rows with counts too low for one Genotype
    # del stat_het

    # CPGs of interest
    # cpgs_plot = stat[(stat['Counts 0/0'] > 3) & (stat['Counts 1/1'] > 3) & (stat['Counts 0/1'] > 3) & (stat['Spearman correlation p_value'] < 1e-5)]['cpg'].unique()
    # highest p-val for highmildhighsev : ['17:46065976:1', '12:112942465:1', '21:33229986:3'] # highest rho
    #
    dict_cpgs = {#'INTEREST': ['17:46768336:3', '6:33069193:2'],
                 # 'HIGHmildHIGHsev': ['12:112925744:1', '17:46065976:1', '21:33229986:3'],
                 # 'EXTRA': ['3:45848456:1', '21:33226777:1', '21:33242527:1'],
                 'Hsev_Lmild': ['9:133271878:1','1:155209089:1'],
                 'Lsev_Hmild': ['9:133271842:1']}

    for supp_title, list_cpg in dict_cpgs.items():
        # INDIVIDUAL PLOTS
        for cpg in list_cpg:
            cpg_df = median_df[median_df['cpg'] == cpg].copy()
            cpg_df['Gen'].replace({'alt': 1, 'ref':0}, inplace=True)
            setup_customizedboxplot_cpg_analysis(cpg_df)


if __name__ == '__main__':
    main(file)

# TODO: descriptive statistics- number of cpg per samples (boxplot per chr?)
# TODO: Joint model
# transform to normal distrib 9inverse norm)
# lm(medianlog ~ snp * clinical, data)
# lm(invnorm(medianlog) ~ snp + clinical + snp * clinical, data)
# lm(invnorm(median_log_hap) ~ allele * clinical , haplotype_data)
