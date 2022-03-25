import sys
import os
import glob
import pandas as pd
import socket
host = socket.gethostname()
if 'Fanny' in host:
    PATH_UTILS = '/home/fanny/Work/EBI/Utils'
    ABS_PATH = '/home/fanny/Work/EBI/covid_nanopore'
else:
    PATH_UTILS = '/nfs/research/birney/users/fanny/Utils'
    ABS_PATH = '/nfs/research/birney/users/fanny/covid_nanopore'
sys.path.insert(0, PATH_UTILS)
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.offsetbox import (AnchoredOffsetbox, DrawingArea, HPacker,
                                  TextArea)
import numpy as np
from utils import *
# from utils_plots import *
import re
import argparse
import subprocess


def boxplot_customized(df, x_var, y_var, hue_var=None, dict_colors='tab10', width_var=None, hatch_var=None, ax=None,  art_df=pd.DataFrame()):
    # df
    # x_var = 'phenotype'
    # y_var = 'log_lik_ratio'
    # hue_var='Genotype'
    # width_var=None
    # hatch_var='phenotype'

    dict_var = {'Hue': hue_var, 'Width': width_var, 'Hatch': hatch_var, 'X_order':x_var}
    # print('VARS', dict_var)
    art_df_bis = art_df[art_df != ''].dropna().copy()

    # if hue_var != xvar:
    #     art_df_bis['X_order'] =
    # else:
    #     art_df_bis['X_order'] =
    #
    None_art = [k for k,v in dict_var.items() if v == None]
    art_df_bis[None_art] = None
    art_df_bis = art_df_bis[art_df_bis['Variable'].isin([v for v in dict_var.values() if v])]
    for k, v in dict_var.items():
        art_df_bis.loc[art_df_bis['Variable'] != v, k] = None

    # print(art_df_bis)
    art_df_bis.sort_values(['Hatch'], inplace=True)
    # if hue_var:
    #     art_df_bis.sort_values([x_var, hue_var], inplace=True)
    # else:
    #     art_df_bis.sort_values([x_var], inplace=True)


    assert not art_df_bis.empty, print('ERROR, Function not implemented if art_df not given !')
    # if art_df_bis.empty:
    #     if hue_var:
    #         dict_colors = dict_colors.copy()
    #         if not isinstance(dict_colors, dict):
    #             dict_colors = {hue_val: sns.color_palette(dict_colors)[n] for n, hue_val in enumerate(df[hue_var].unique())}
    #     art_df.loc[hatch_var, 'Hatch'].nunique()
    #     if hatch_var:
    #         hatch_list = ['+', 'x', 'o', '/', '|', '-', '\ ', 'O', '.', '*']
    #         assert len(hatch_list) > df[hatch_var].nunique(), 'Hatching list not big enough, recombination to do'
    #     else:
    #         hatch_list = ['', '', '']

    # Creation of the plot
    contour_line = sns.color_palette('dark')[-3]
    g2 = sns.boxplot(data=df, x=x_var, y=y_var, orient='v', hue=hue_var, palette=dict_colors, ax=ax)

    # Customisation of the artists
    # print(g2.artists)
    art_df_bis.index = g2.artists
    for art in g2.artists:
        # print(art, art_df_bis.loc[art, 'Hatch'], art_df_bis.loc[art, 'Width'])
        art.set_hatch(art_df_bis.loc[art, 'Hatch'])
        art.set_linewidth(art_df_bis.loc[art, 'Width'])

    print('FIX HATCHING FOR DOUBLE HUE XVAR!!!!!!!!!')
    # TODO: FIX HATCHING FOR DOUBLE HUE XVAR!!!!!!!!!

    # blou = df.copy()
    # for key, var in dict_var.items():
    #     n = 0
    #     if var:
    #         for val in blou[var].unique():
    #             blou.loc[blou[var] == val, key] = n
    #             n = n+1
    #     else:
    #         blou[key] = 0

    # list_var = set([x_var] + list(dict_var.keys()) + [val for val in dict_var.values() if val is not None])
    # art_df = blou[list_var].value_counts().reset_index()

    # TODO: reduce art_df to the current axis
    # assert len(art_df.index) == len(g2.artists), 'One or several variables are affected to non-existent patches'

    #
    # if hatch_var:
    #     g2.legend(handles=[mpatches.Patch(facecolor='white',
    #                              label=key, lw=2, ec=contour_line,
    #                              hatch=hatch_list[int(art_df.loc[art_df[hatch_var] == key, 'Hatch'][0])]
    #                              ) for key, val in art_df[[hatch_var, 'Hatch']].value_counts().index.tolist()],
    #                              title=hatch_var, loc=3)

    return g2


def setup_customizedboxplot_cpg_analysis(cpg_df, unit='control_snp', dir_out='', dict_colors=None):
    # General modification of the datas
    cpg_df = df.copy()
    cpg = str(cpg_df['cpg'].unique())[2:-2]
    snp = str(cpg_df[unit].unique().tolist())[2:-2]
    ref, alt = snp.split(':')[-2:]
    replace_dict = {'haplotype': {0: ref, 1: alt}, 'Genotype': {'0/0': f'{ref}/{ref}', '0/1': f'{ref}/{alt}', '1/1': f'{alt}/{alt}'}}
    cpg_df.replace(replace_dict, inplace=True)

    # Aesthetic parameters
    hatch_list = ['+', 'O', '.', '*', '-', '|', '/', '\ ', 'x', 'o',]
    if not dict_colors:
        dict_colors = {'Genotype': {'0/0':'#1678F5', '0/1':'#2aa69a', '1/1':'#3ED43E'},
                       'haplotype': {0:'#1678F5', 1:'#3ED43E'},
                       'phenotype': {'Mild':'#f0f0f5', 'Severe':'#c2c2d6'}}

    # Modification of the color dataframe to fit the replacement
    repl_colors, flat_dict = {}, {}
    for key, dict in dict_colors.items():
        if replace_dict.get(key):
            flat_dict.update({replace_dict.get(key).get(k): dict_colors.get(key).get(k) for k in dict.keys() if replace_dict.get(key).get(k) != None})
            repl_colors.update({key: {replace_dict.get(key).get(k): dict_colors.get(key).get(k) for k in dict.keys() if replace_dict.get(key).get(k) != None}})
        else:
            repl_colors.update({key:{k: dict_colors.get(key).get(k) for k in dict.keys() }})
            flat_dict.update({k: dict_colors.get(key).get(k) for k in dict.keys() })

    # Creation of the dataframe containing the aesthetic parameters
    art_df = pd.DataFrame(flat_dict.values(), index=flat_dict.keys(), columns=['Hue'])
    art_df['Hatch'] = hatch_list[:art_df.shape[0]]
    art_df['Width'] = range(1, art_df.shape[0]*2, 2)
    val='G/G'
    art_df = art_df.reset_index().rename(columns={'index' : 'Value'})

    art_df['Variable'] = [cpg_df[cpg_df == val].dropna(how='all',axis=1).columns[0]  if not cpg_df[cpg_df == val].dropna(how='all',axis=1).empty else None for val in art_df['Value']]
    art_df.dropna(inplace=True)

    # Creation of the plot
    fig, ax = plt.subplots(2, 3, figsize=(17,10))
    boxplot_customized(cpg_df, 'phenotype', 'log_lik_ratio', hue_var='phenotype',
                                 dict_colors=repl_colors['phenotype'], width_var=None,
                                 hatch_var='phenotype', ax=ax[0,0], art_df=art_df)
    ax[0,0].set(title='Phenotype')
    boxplot_customized(cpg_df, 'Genotype', 'log_lik_ratio', hue_var='Genotype',
                                 dict_colors=repl_colors['Genotype'],
                                 width_var=None, ax=ax[0,1],  art_df=art_df)
    ax[0,1].set(title='Genotype')#, xticklabels=replace_val['Genotype'].values())
    boxplot_customized(cpg_df, 'phenotype', 'log_lik_ratio', hue_var='Genotype',
                                 dict_colors=repl_colors['Genotype'], width_var=None,
                                 hatch_var='phenotype', ax=ax[0,2], art_df=art_df)
    ax[0,2].set(title='Genotype X Phenotype')
    print(cpg_df[cpg_df['Genotype'] == '0/1'], cpg_df[cpg_df['Genotype'] == '0/1'].shape)
    boxplot_customized(cpg_df[cpg_df['Genotype'] == replace_dict['Genotype']['0/1']], 'phenotype',
                                   'log_lik_ratio', hue_var='phenotype',
                                   dict_colors=repl_colors['phenotype'],
                                   hatch_var='phenotype', width_var=None,
                                   ax=ax[1,0], art_df=art_df)
    ax[1,0].set(title='Heterozygous Phenotype')
    boxplot_customized(cpg_df[cpg_df['Genotype'] == replace_dict['Genotype']['0/1']], 'haplotype', 'log_lik_ratio',
                                  hue_var='haplotype', dict_colors=repl_colors['haplotype'],
                                  hatch_var=None, width_var=None, ax=ax[1,1],
                                  art_df=art_df)
    ax[1,1].set(title='Heterozygous Haplotype') #, xticklabels=replace_val['haplotype'].values())
    boxplot_customized(cpg_df[cpg_df['Genotype'] == replace_dict['Genotype']['0/1']], 'phenotype', 'log_lik_ratio',
                                  hue_var='haplotype', dict_colors=repl_colors['haplotype'],
                                  hatch_var='phenotype', width_var=None, ax=ax[1,2],
                                  art_df=art_df)
    ax[1,2].set(title='Heterozygous Haplotype X Phenotype')
    plt.suptitle(f'CpG {cpg} associated with SNP {snp}')

    # Creation of one common legend for all
    lines, labels = [], []
    x = 0.6
    for a in fig.axes:
        Line, Label = a.get_legend_handles_labels()
        a.get_legend().remove()
        if (set(Label) & set(labels)) == set():
            lines.extend(Line)
            labels.extend(Label)
            Label_title =  art_df.loc[art_df['Value'].isin(Label), 'Variable'].tolist()[0]
            print(Label, Label_title)
            if Label_title == 'phenotype':
                print('TOFIX')
                # fig.legend(handles=[mpatches.Patch(facecolor='white',
                #                          label=key, lw=2, ec=contour_line,
                #                          hatch=hatch_list[int(art_df.loc[art_df[hatch_var] == key, 'Hatch'][0])]
                #                          ) for key, val in art_df[[hatch_var, 'Hatch']].value_counts().index.tolist()],
                #                          title=hatch_var, loc=3)
            else:
                fig.legend(Line, Label, loc='center left', bbox_to_anchor=(0.9, x), title=Label_title)
            x = x - 0.08


    # save_plot_report(f'{dir_out}/Multiplots_{cpg}.png', fig, file=None)
    # fig.savefig(f'{dir_out}/Multiplots_{cpg}.png')

# setup_customizedboxplot_cpg_analysis(df)

# TODO change pvalplot to be as inclusive as possible -> tranfert in utils_plots
def pval_plot_new(stat, xvar, pval_col, pval_cutoff=0.01, n_site=2, format_xaxis=True, out_dir='', title_supp=''):
    # stat = mw[(mw['index'] == val) & (mw['data'] == data)].copy()
    # xvar='cpg'
    # pval_col='p-val'
    # pval_cutoff=0.01
    # n_site=2
    # title_supp=f"MannWhitney_{val}_{data.replace('/', '-')}"
    # out_dir=''
    # format_xaxis=True

    stat[['CHR', 'POS']] = stat[xvar].str.split(':', expand=True)[[0, 1]].astype(int)
    stat['CHR'] = stat['CHR'].astype('category')
    stat['minus_log10'] = -np.log10(stat[pval_col].astype(float))
    stat['cutoff'] = pval_cutoff
    stat.sort_values('CHR', inplace=True)

    # Selection best snp
    best_cpg = stat[stat['cutoff'] < stat['minus_log10']].sort_values(
        'minus_log10', ascending=False).groupby('CHR').head(n_site)[
        'cpg'].tolist()
    stat.loc[(stat[xvar].isin(best_cpg)) &
        (stat['minus_log10'] > stat['cutoff']),#
        f'{xvar}_best'] = stat.loc[(stat[xvar].isin(best_cpg)) &
                    (stat['minus_log10'] > stat['cutoff']), 'cpg']
    stat.loc[stat[f'{xvar}_best'].isna(), f'{xvar}_best'] = ' '
    print(f'ABOVE CUTOFF {title_supp}: {best_cpg}')

    # Format X-axis
    if format_xaxis:
        last_pos = 0
        for chr in stat.sort_values(['CHR', 'POS'])['CHR'].astype(int).unique():
            stat.loc[stat['CHR'] == chr, 'new_xvar'] = stat.loc[stat['CHR'] == chr, 'POS'] - stat.loc[stat['CHR'] == chr, 'POS'].min() + last_pos + 1
            last_pos = stat.loc[stat['CHR'] == chr, 'new_xvar'].max()
            # print(last_pos)
        xvar = 'new_xvar'

    # PLOT
    g = sns.FacetGrid(stat, aspect=4, height=4, palette='Spectral',
                      margin_titles=True)
    g.map(sns.lineplot, xvar, 'cutoff', hue=None)
    g.map_dataframe(sns.scatterplot, xvar, 'minus_log10', hue='CHR',
        legend=True)
    g.set(xlabel="CHR", xticks=stat.groupby(['CHR']).last()[xvar].unique(),
        xticklabels=stat['CHR'].unique())
    try:
        for row in stat.iterrows():
            row = row[1]
            g.axes[0,0].text(row[xvar], row.minus_log10 + 0.2, row.cpg_best,
                horizontalalignment='left')
    except Exception as err:
        print('couldn\'t print best points', err)
    plt.suptitle(title_supp)
    save_plot_report(f'minuslog10_{pval_col}_pvalue_{title_supp}', g, output=out_dir, file=sys.stdout)

    return stat[stat['cutoff'] < stat['minus_log10']].sort_values(
        'minus_log10', ascending=False)['cpg'].tolist()


def pval_plot_old(stat, xvar, pval_col, pval_cutoff=0.01, n_site=2, format_xaxis=True, out_dir='', title_supp=''):
    stat = stat.copy()
    stat[['CHR', 'POS']] = stat[xvar].str.split(':', expand=True)[[
        0, 1]].astype(int)
    stat = stat.sort_values('CHR')
    stat['CHR'] = stat['CHR'].astype('category')
    stat['minus_log10'] = -np.log10(stat[pval_col].astype(float))
    stat['cutoff'] = -np.log10(pval_cutoff/stat.shape[0])
    stat.sort_values('CHR', inplace=True)

    # Format X-axis
    if format_xaxis:
        for chr in stat.sort_values(['CHR', 'POS'])['CHR'].astype(int).unique():
            stat.loc[stat['CHR'] == chr, 'new_xvar'] = stat.loc[stat['CHR'] == chr, 'POS'] - stat.loc[stat['CHR'] == chr, 'POS'].min() + 1

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
    xvar = 'new_xvar'
    g = sns.FacetGrid(stat, aspect=4, height=4, palette='Spectral',
                      margin_titles=True)
    g.map(sns.lineplot, xvar, 'cutoff', hue=None)
    g.map_dataframe(sns.scatterplot, xvar, 'minus_log10', hue='CHR',
        legend=True)
    g.set(xlabel="CHR", xticks=stat.groupby(['CHR']).last()[xvar].unique(),
        xticklabels=stat['CHR'].unique())
    try:
        for row in stat.iterrows():
            row = row[1]
            g.axes[0,0].text(row[xvar], row.minus_log10 + 0.2, row.cpg_best,
                horizontalalignment='left')
    except Exception as err:
        print('couldn\'t print best points', err)
    save_plot_report(f'minuslog10_{pval_col}_pvalue_{title_supp}', g, output=out_dir, file=sys.stdout)
    # g.savefig(f'{out_dir}minuslog10_{pval_col}_pvalue_{title_supp}.png')
    return stat[stat['cutoff'] < stat['minus_log10']].sort_values(
        'minus_log10', ascending=False).groupby('CHR')['cpg'].tolist()


def save_results_count_significant_cpg(res_cpg):
    res_cpg['ratio_count'] = res_cpg['cpg_significant'].astype(float) * 100 / res_cpg['total_cpg'].astype(float)
    res_cpg['ratio_plot'] = res_cpg['plot_cpg_number'].astype(float) * 100 / res_cpg['total_cpg'].astype(float)
    res_cpg.reset_index(inplace=True)
    res_cpg.dropna(inplace=True)
    res_cpg['test'] = res_cpg['test'].str[4:-4]
    res_cpg['value'] = res_cpg['value'].str.replace('0/1', 'HET')
    res_cpg['data'] = res_cpg['dataset'].str.replace('0/1', 'HET')
    blou = res_cpg.pivot(values=['ratio_count', 'ratio_plot'], columns=['type'], index=['test', 'value', 'dataset']).reset_index().round(2)
    # blou.columns = ['index', 'covid', 'random']
    blou.to_csv('Ratio_hits.csv', index=False)


def main(files, snp_type, output_dir='plots', pval_cutoff=0.01):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    file_ls = files.split('--')
    print(file_ls)

    # QUESTION: do we need to run it in separate job ?
    if snp_type == 'covid_snp':
        cutoff = -np.log10(pval_cutoff/df['cpg'].nunique())
        cutoff = round(cutoff, 1)
    else:
        cutoff = pval_cutoff
    print(cutoff)

    for file in file_ls:
        df = pd.read_csv(file, usecols=['index', 'data']).drop_duplicates()
        list_data = df['data'].unique()
        if 'Mann_Whitney' in file:
            list_val = ['Mild-Severe', 'alt-ref']
        else:
            list_val = df['index'].unique()
        if (list_val != []) & (list_data != []):
            df = df[(df['index'].isin(list_val)) & (df['data'].isin(list_data))]
        res_df = pd.DataFrame(columns=['test', 'value', 'dataset', 'type', 'total_cpg', 'cpg_significant', 'plot_cpg_number'])
        for i in df.index:
            val, data = df.loc[i].tolist()
            grep_out = subprocess.Popen([f"grep {data} {file} | grep {val}"], stdout=subprocess.PIPE,shell=True)
            df_data = pd.DataFrame([n.split(',') for n in str(grep_out.communicate()[0]).split('\\n')], columns=['index', 'p-val', 'cpg', 'data'])
            df_data = df_data[(df_data['index'] == val) & (df_data['data'] == data)]
            # Pval plot
            list_cpg = pval_plot_new(df_data, xvar='cpg', pval_col='p-val', pval_cutoff=cutoff, n_site=2,
                      title_supp=f"{file}_{val}_{data.replace('/', '-')}_{snp_type}", out_dir=output_dir,
                      format_xaxis=False)
            res = pd.DataFrame([file, val, data, snp_type, df_data['cpg'].nunique(), df_data[df_data['p-val'].astype(float) < pval_cutoff]['cpg'].nunique(), len(list_cpg)]).T
            res.columns = res_df.columns
            res_df = res_df.append(res, sort=False)
            print(f'\n{file} {data} {val}:\n', df_data[df_data['p-val'].astype(float) < 0.01]['cpg'].unique())

    save_results_count_significant_cpg(res_df)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='STEP2 - Pipeline for mQTLs'
        + ' - Statisques (Spearman correlation and MannWhitney test) on specific datasets')
    parser.add_argument('files', type=str, help='files to plot in concatenated string')
    parser.add_argument('-o', '--output_dir', type=str, default='plots',
                        help='directory to save output files')
    parser.add_argument('-s', '--snp_type', type=str, default='snp',
                        help='unit to perform analysis')
    parser.add_argument('-p', '--pval_cutoff', type=float, default=0.01,
                        help='pvalue to define significant cpg')
    args = parser.parse_args()
    main(**vars(args))



################################################################################
#############     Functions to detect CPG of interest !!     ###################
################################################################################
def define_effect(df):
    df = df.copy()
    df[['index1', 'index2']] = df['index'].str.split('-', expand=True)
    #todo find common substring between index1 and 2, if its contains '/' --> genotype analysis
    df.loc[(df['index1'].str.split('_', expand=True).astype(str)[0] == df['index2'].str.split('_', expand=True).astype(str)[0])
                | (df['index1'].str.contains('/') == False), 'effect'] = 'Phenotype'
    df.loc[(df['index1'].str.split('_', expand=True).astype(str)[1] == df['index2'].str.split('_', expand=True).astype(str)[1])
                & (df['index1'].str.contains('/')), 'effect'] = 'Genotype'

    return df.drop(['index1', 'index2'], axis=1)


def select_SNP_per_pvalue(df, pval_col, dist_bp=500000):
    df = df.copy()
    df[['CHR', 'POS']] = df['cpg'].str.split(':', expand=True)[[0, 1]].astype(int)
    best_snp = []
    for chrom in df['CHR'].unique():
        chr_df = df[df['CHR'] == chrom]
        while chr_df.shape[0] > 0:
            best_pval = chr_df[pval_col].min()
            best_snp += chr_df[chr_df[pval_col]
                               == best_pval]['cpg'].unique().tolist()
            pos = chr_df[chr_df[pval_col]
                         == best_pval]['POS'].unique().tolist()[-1]
            chr_df = chr_df[(chr_df['POS'] < pos - dist_bp)
                            | (chr_df['POS'] > pos + dist_bp)]
    return best_snp


# mw = define_effect(mw[mw['p-val'] < 0.01])
# dict_cpg = {}
# for effect in mw['effect'].unique():
#     blou = mw[mw['effect'] == effect]
#     dict_cpg[effect] = {}
#     for index in blou[blou['index'].str.contains('_') == False]['index'].unique():
#         dict_cpg[effect][f'{index}'] = select_SNP_per_pvalue(blou[blou['index'] == index], 'p-val', dist_bp=500000000)
#         print(effect, index, len(dict_cpg[effect][f'{index}'] ))
#
# effect= 'Genotype'
# blou = dict_cpg[effect]
# len(blou.values())
#
# # INTERSECTION OF THE FOUR
# set_list = [set(list_cpg) for list_cpg in blou.values()]
# set.intersection(*set_list)
#
# # INTERSECTION OF 3
# for index in blou.keys():
#     set_list = [set(list_cpg) for key, list_cpg in blou.items() if key != index]
#     print(len(set_list))
#     print(index.upper(), set.intersection(*set_list))
#
# # INTERSECTION OF 2
# for index in blou.keys():
#     bli = blou.copy()
#     del bli[index]
#     for i in bli.keys():
#         set_list = [set(list_cpg) for key, list_cpg in bli.items() if key != i]
#         print([k for k in bli.keys() if k != i])
#         print(set.intersection(*set_list))
#
# cpg_ls =  ['4:125732000:1', '4:125726363:1', '4:125723523:1', '4:125733776:1', '4:125725982:1','4:125729588:1', '11:118924979:1', '19:835833:3', '21:46133825:1']
# #### Funtion to extract the cpg of interest from MWU and info in ensembl




################################################################################
#                           Test code for information                          #
################################################################################
def get_info_cpg(cpg_ls):
    import ensembl_rest
    import pyensembl
    from bioservices.kegg import KEGG
    ensembl = pyensembl.EnsemblRelease()
    k = KEGG()
    cpg_df = pd.DataFrame(sorted(cpg_ls))
    cpg_df[['chr', 'pos']] = cpg_df[0].str.split(':', expand=True)[[0,1]]
    cpg_df['pos'] = cpg_df['pos'].astype(int)
    results = {}
    for row in cpg_df.index:
        cpg = list(cpg_df.loc[row])
        gene_ls = ensembl.gene_names_at_locus(contig=cpg[1], position=int(cpg[2]))
        for gene in gene_ls:
            a, b, c = '', '', ''
            if gene:
                results[cpg[0]] = {gene: []}
                try:
                    a = ensembl_rest.symbol_lookup(species='homo sapiens', symbol=gene)
                    if (a != ''):
                        results[cpg[0]][gene].append(a.get('display_name'))
                        results[cpg[0]][gene].append(re.sub("\[.*\]"," ", a.get('description')))
                except:
                    pass
                try:
                    b = k.find("hsa", gene).split('\t')[0]
                    if (b != '') & (b != '\n'): c = k.get(b)
                    if (c != '') & isinstance(c, str):
                            c = c[c.find('ORTHOLOGY')+10: c.find('ORGANISM')-1]
                            results[cpg[0]][gene].append(c)
                except:
                    pass
    return results

# import os
# import glob
# list_cpg = glob.glob('/home/fanny/Work/EBI/covid_nanopore/covid_snp old results/Multiplots*')
# list_cpg = [os.path.basename(l).split('_')[1][:-4] for l in list_cpg]
# list_cpg = set(sorted(list_cpg))
# cpg_man = ['22:22497096:2', '14:94038739:2.0', '13:41390703:2.0', '13:41390313:1.0', '14:101290108:1', '16:84186751:6', '10:20510924:1', '15:40718252:1.0', '22:22498857:2', '19:861278:1', '17:38574125:1', '1:54438691:1.0', '20:8831166:1', '16:28454050:1', '11:72839989:1', '4:35507549:1', '21:33275223:1', '20:8832013:1', '5:172265703:1', '3:86481684:1', '4:35509327:1', '12:70206114:1', '5:140194463:1', '11:2283191:1', '3:21420995:1', '2:141164217:1', '2:36568940:1', '7:17479501:1', '1:83739112:3', '6:29874488:1', '6:29871769:1', '19:27942785:1', '17:64844277:1', '7:64304879:1', '12:25490551:1', '9:19189676:1', '8:95417340:1', '10:20570117:1', '15:40733257:1', '18:55845130:2', '18:62375583:1', '21:33276218:1', '9:100366071:1', '8:95416218:1']
# dict_man = get_info_cpg(list_cpg)
# dict_man
#
#
# cpg_spear = ['14:77321123:3', '14:103946111:1', '16:919128:1', '17:64788918:1']
#
# dict_spear = get_info_cpg(cpg_spear)
# dict_spear


################################################################################
#                      Test code for individual plots                          #
################################################################################

# file = 'Filtered_finemapped.csv10.csv'
# df = pd.read_csv(file, header=None, names=['chromosome', 'strand', 'start', 'end', 'read_name', 'log_lik_ratio', 'log_lik_methylated', 'log_lik_unmethylated', 'num_calling_strands', 'num_motifs', 'sequence', 'sample_id', 'control_snp', 'covid_snp', 'control_snp_rs', 'covid_snp_rs', 'base_called', 'pos', 'ref', 'alt', 'haplotype', 'Allele1', 'Allele2', 'Genotype', 'phenotype', 'cpg', 'distance_cpg_snp'])
# median_df = df.groupby(['phenotype', 'sample_id', 'chromosome', 'cpg',
#                                 'control_snp', 'Genotype', 'haplotype']).median().reset_index()
# median_df.groupby('cpg').count()['log_lik_ratio'].sort_values()
# cpg = '10:65168344:1'
# df = median_df[(median_df['cpg'] == cpg) & (median_df['haplotype'] != 'other')].copy()
# df['phenotype'].unique()
# df['haplotype'].replace({'alt': 1, 'ref':0}, inplace=True)
# setup_customizedboxplot_cpg_analysis(df)
#


################################################################################
#                            COUNT-HITS  (PVAL < 0.01)                         #
################################################################################

# import pandas as pd
# res_cpg = pd.DataFrame()
#
# # COVID SNPS
# mw_file = '/home/fanny/Work/EBI/covid_nanopore/covid_snp new results/All_Mann_Whitney.csv'
# sp_file = '/home/fanny/Work/EBI/covid_nanopore/covid_snp new results/All_Spearmann_corr.csv'
# mw = pd.read_csv(mw_file)
# sp = pd.read_csv(sp_file)
# snp='covid'
#
# for val in ['Mild-Severe', 'alt-ref']:
#     for data in mw[mw['index'] == val]['data'].unique():
#         df = mw[(mw['index'] == val) & (mw['data'] == data)]
#         res_cpg.loc[f'MW-{val}-{data}-{snp}', ['tot_cpg', 'cpg_sign', 'snp']] = df['cpg'].nunique(), df[df['p-val'] < 0.01]['cpg'].nunique(), snp
#
# for val in sp['index'].unique():
#     for data in sp[sp['index'] == val]['data'].unique():
#         df = sp[(sp['index'] == val) & (sp['data'] == data)]
#         res_cpg.loc[f'SP-{val}-{data}-{snp}', ['tot_cpg', 'cpg_sign', 'snp']] = df['cpg'].nunique(), df[df['p-val'] < 0.01]['cpg'].nunique(), snp
#
# #Random SNPS
# mw_file = '/home/fanny/Work/EBI/covid_nanopore/random-SNPS_first_results/final_Mann_Whitney.csv'
# mw = pd.read_csv(mw_file, header=None, names= ['index', 'p-val', 'cpg', 'data'])
# sp_file = '/home/fanny/Work/EBI/covid_nanopore/random-SNPS_first_results/final_sper.csv'
# sp = pd.read_csv(sp_file)
# snp='random'
#
# for val in ['Mild-Severe', 'alt-ref']:
#     for data in mw[mw['index'] == val]['data'].unique():
#         df = mw[(mw['index'] == val) & (mw['data'] == data)]
#         res_cpg.loc[f'MW-{val}-{data}-{snp}', ['tot_cpg', 'cpg_sign', 'snp']] = df['cpg'].nunique(), df[df['p-val'] < 0.01]['cpg'].nunique(), snp
#
# for val in sp['index'].unique():
#     for data in sp[sp['index'] == val]['data'].unique():
#         df = sp[(sp['index'] == val) & (sp['data'] == data)]
#         res_cpg.loc[f'SP-{val}-{data}-{snp}', ['tot_cpg', 'cpg_sign', 'snp']] = df['cpg'].nunique(), df[df['p-val'] < 0.01]['cpg'].nunique(), snp
#
#
# res_cpg['ratio'] = res_cpg['cpg_sign'] * 100 / res_cpg['tot_cpg']
# res_cpg.index = res_cpg.reset_index()['index'].str.rsplit('-', n=1, expand=True)[0]
# res_cpg.reset_index(inplace=True)
# res_cpg.dropna(inplace=True)
# res_cpg['data'] = res_cpg[0].str.rsplit('-', n=1, expand=True)[1]
# res_cpg['test'] = res_cpg[0].str.split('-', n=1, expand=True)[0]
# res_cpg=res_cpg[res_cpg['data'].isin(['Mild', 'Severe'])==False]
# res_cpg[0] = res_cpg[0].str.replace('MW', 'MannWhitney')
# res_cpg[0] = res_cpg[0].str.replace('SP', 'SpearmanCorr')
# res_cpg[0] = res_cpg[0].str.replace('0/1', 'HET')
# res_cpg['data'] = res_cpg['data'].str.replace('0/1', 'HET')
# res_cpg.rename(columns={0: 'stat_test'}, inplace=True)
#
# blou = res_cpg.pivot(values=['ratio'], columns=['snp'], index=['stat_test']).reset_index().round(2)
# blou.columns = ['index', 'covid', 'random']
# blou.to_csv('Ratio_hits.csv')
#
# import seaborn as sns
# g = sns.catplot(data=res_cpg, y='stat_test', x='ratio', hue='snp', col='data', kind='bar', sharey=False, aspect=1.5)
# save_plot_report('random-SNPS_first_results/hits_per_test.jpg', g)

################################################################################
############################    P-VAL  PLOTS    ################################
################################################################################

# # COVID SNPS
# mw_file = '/home/fanny/Work/EBI/covid_nanopore/covid_snp new results/All_Mann_Whitney.csv'
# sp_file = '/home/fanny/Work/EBI/covid_nanopore/covid_snp new results/All_Spearmann_corr.csv'
# mw = pd.read_csv(mw_file)
#
# #Random SNPS
# mw_file = '/home/fanny/Work/EBI/covid_nanopore/random-SNPS_first_results/final_Mann_Whitney.csv'
# mw = pd.read_csv(mw_file, header=None, names= ['index', 'p-val', 'cpg', 'data'])
# sp_file = '/home/fanny/Work/EBI/covid_nanopore/random-SNPS_first_results/final_sper.csv'
#
#
# ##RUN
# mw = mw[(mw == mw.columns) == False].dropna(how='all')
# dict_cpg = {}
# for val in ['Mild-Severe', 'alt-ref']:
#     for data in mw[mw['index'] == val]['data'].unique():
#         print(val.upper(), data.upper())
#         # print(mw[(mw['index'] == val) & (mw['data'] == data)].head(3))
#         tot_mw = mw[(mw['index'] == val) & (mw['data'] == data)]['cpg'].nunique()
#         # list_cpg = info_pval_plot(mw[(mw['index'] == val) & (mw['data'] == data)],
#         #           xvar='cpg', pval_col='p-val', pval_cutoff=0.01, n_site=2,
#         #           title_supp=f"MannWhitney_{val}_{data.replace('/', '-')}", out_dir='',
#         #           format_xaxis=False)
#         list_cpg = pval_plot_new(mw[(mw['index'] == val) & (mw['data'] == data)],
#                   xvar='cpg', pval_col='p-val', pval_cutoff=0.01, n_site=2,
#                   title_supp=f"MannWhitney_{val}_{data.replace('/', '-')}", out_dir='',
#                   format_xaxis=False)
#         dict_cpg[f'MW-{val}-{data}'] = sorted(set(list_cpg)), tot_mw
#
#
# sp = pd.read_csv(sp_file)
# sp = sp[(sp == sp.columns) == False].dropna(how='all')
# sp.nunique()
# for val in sp['index'].unique():
#     for data in sp[sp['index'] == val]['data'].unique():
#         print(val.upper(), data.upper())
#         # print(sp[(sp['index'] == val) & (sp['data'] == data)].head(3))
#         tot_sp = sp[(sp['index'] == val) & (sp['data'] == data)]['cpg'].nunique()
#         # list_cpg = info_pval_plot(sp[(sp['index'] == val) & (sp['data'] == data)],
#         #           xvar='cpg', pval_col='p-val', pval_cutoff=0.01, n_site=2,
#         #           title_supp=f"Spearman_{val}_{data.replace('/', '-')}", out_dir='',
#         #           format_xaxis=False)
#         list_cpg = pval_plot_new(sp[(sp['index'] == val) & (sp['data'] == data)],
#                   xvar='cpg', pval_col='p-val', pval_cutoff=0.01, n_site=2,
#                   title_supp=f"Spearman_{val}_{data.replace('/', '-')}", out_dir='',
#                   format_xaxis=False)
#         dict_cpg[f'SP-{val}-{data}'] = sorted(set(list_cpg)), tot_sp
#
#
# with open('Cpg_list_covid.txt', 'w') as f:
#     print('key len(list) tot ratio ', tot_sp, file=f)
#     for key, val in dict_cpg.items():
#         print('\n', key, ' ', len(val[0]), val[1], len(val)/val[1], '\n', val, '\n', file=f)
#
#


################################################################################
##           !!!                   TODO                     !!!               ##
################################################################################
# - RHO vs Means plots
# - RHO mild vs Sev with pval --> We do not have the dataset to do only heterozygotes
# - Distance from SNP plot
# - plot p-val MWU phenotype vs p-val Spearman genotype
# - individuals (TODO: search how to arrang the legend)
#   NB : detection of the good CPGS ; good enough ?  !!!!!    SELECT BY COUNTS
#   NB : Search the ensembl database for info : MORE ?

# QUESTION: Could it be interesting to have a look at the ref alt differences in homozygotes ?
# QUESTION: Do we need to remove the called 'alt' in '0/0' and 'ref' in '1/1' ??

# CPG selection: the most represented ???
# most_represented = median_df.groupby(['cpg', 'CHR']).count().reset_index(
#                     ).groupby(['CHR']).first()['cpg']

# # TODO Simplify the calling of the function and the diferent created dataframes for running the script on the cluster

#     sp = pd.read_csv('Spearmann_corr.csv')
#     # QUESTION: When should we reduce by count number
#     mw[mw['index'] == 'Mild-Severe']
#
#     # COMBAK: pval plot modifies so df needs to include thos modifications:


#     # Rho-means plots
#     count = pd.read_csv('Counts_Diff_means.csv')
#     count.rename(columns={'index' : 'cpg', 'means_ref-alt': 'means_refalt'}, inplace=True)
#     sp_rho = pd.merge(sp[(sp['variable'] == 'Rho')], count[['cpg', 'means_refalt']], on='cpg')
#     for data in ['ALL', '0/1']:
#         for val in sp['index'].unique():
#             g = sns.relplot(data=sp_rho[(sp_rho['data'] == data) & (sp_rho['index'] == val)],
#                             kind='scatter', y='means_refalt', x='value')
#             save_plot_report(f"Rho_{data.replace('/', '-')}_{val}", g, file=None)
#
#     # RHo vs RHO
#     for val in sp['index'].unique():
#         try:
#             g = sns.relplot(data=sp_rho[(sp_rho['index']==val)].pivot_table(index=['index', 'cpg'],
#                                        columns='data')['value'].reset_index(),
#                 x='Mild', y='Severe')
#             save_plot_report(f"Rho_MildSevere_{val}", g, file=None)
#         except ValueError:
#             pass

#     # Spearman VS Mann Whitney
#     mw.loc[mw['index'].str.contains('/'), 'type'] = 'Genotype'
#     mw.loc[mw['index'] == 'alt-ref', 'type'] = 'Gen'
#     mw.loc[mw['index'] == 'Mild-Severe', 'type'] = 'phenotype'
#     mw.loc[mw['index'].str.contains('_'), 'type'] = 'mix'
#
#     for data in mw['data'].unique():
#         print(data)
#         spmw = pd.merge(mw[mw['data']== data], sp[sp['data']== data], right_on=['index', 'cpg', 'data'],
#                  left_on=['type', 'cpg', 'data'], copy=False, how='inner')
#         for val in ['phenotype', 'Gen']:
#             sns.relplot(data=spmw[(spmw['index_y'] == val) & (spmw['variable']=='p-val')], x='value', y='p-val', hue='variable')
#
#     # Distance from SNP
#     ### Log lik ratio plotting
#     for col in ['cpg', 'SNP']:
#         median_df[f'POS_{col}'] = median_df[col].str.split(':', expand=True)[1]
#     median_df['distance'] = median_df['POS_SNP'].astype(int) - median_df['POS_cpg'].astype(int)
#     sp = pd.merge(sp, median_df[['cpg', 'SNP', 'distance', 'Genotype', 'Gen']], on = 'cpg')
#
#     # QUESTION: Is this the way for the distance plots ?
#     sp[['data', 'index']].value_counts()
#     for data in sp['data'].unique():
#         for index in sp['index'].unique():
#             data = 'ALL'
#             index = 'Genotype'
#             print(data, index)
#             sns.relplot(data=sp[(sp['index'] == index) & (sp['data']==data) & (sp['variable'] == 'p-val')],
#                         x='distance', y='value', hue='Gen', col='Genotype')



################################################################################
################################################################################

    # Mean/Median over binned sites

    # binned = median_df.copy()
    # binned, _ = drop_datas(median_df, 'name', thresh_zscore=3)
    # binned['pos'] = binned['cpg'].str.split(':', expand=True)[1].astype(int)
    #
    # # binned['distance_to_median'] = binned['cpg'].str.split(':', expand=True)[1].astype(int)
    # for chr in binned['CHR'].sort_values().unique():
    #     chr_index = binned[binned['CHR'] == chr].index
    #
    #     # Manual binning
    #     binned.loc[chr_index, 'bin'] = (
    #         (binned.loc[chr_index, 'pos'] - binned.loc[chr_index, 'pos'].median()
    #         )/binned.loc[chr_index, 'pos'].median()).round(1)
    #     # print(f'CHR = {chr}\n\n', binned['bin'] .value_counts().to_markdown())
    #     # print(binned['pos'].value_counts().to_markdown())
    #
    #     # Binning per cut
    #     binned.loc[binned['CHR'] == chr, 'bin'] = pd.cut(
    #         x=binned[binned['CHR'] == chr]['cpg'].str.split(
    #         ':', expand=True)[1].astype(int), bins=4).astype(str)
    #     binned.loc[binned['CHR'] == chr, 'bin'] = pd.qcut(
    #         x=binned[binned['CHR'] == chr]['cpg'].str.split(
    #         ':', expand=True)[1].astype(int), q=3, duplicates='drop').astype(str)
    #
    # binned['bin'] = binned['CHR'].astype(str) + ':' +  binned['bin']
    # binned.sort_values('cpg', inplace=True)
