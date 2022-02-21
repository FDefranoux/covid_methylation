import sys
import os
import glob
import pandas as pd
sys.path.insert(0, "../Utils")
from utils import *
import seaborn as sns
import matplotlib.pyplot as plt
from utils_plot import *


### PLOTS already in utils_plots but to be tested !###
def violinplot(df, xvar, yvar, huevar=None, colvar=None, title_supp='', out_dir='',, file=None):
    n_snp = df[xvar].nunique()
    print(n_snp)
    if n_snp < 10:
        g = sns.catplot(data=df, kind='violin',
                            y=yvar, x=xvar, hue=huevar,
                            col=colvar, height=6,
                            aspect=0.9, orient='v',
                            inner="quartile", linewidth=1.5,
                            sharex=False, sharey=False)
    elif n_snp < 50:
        g = sns.catplot(data=df, kind='violin',
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
        g = sns.catplot(data=df, kind='violin',
                            y=yvar, x=xvar, hue=huevar,
                            col=colvar, height=45,
                            aspect=.1, orient='v',
                            inner="quartile", linewidth=1,
                            sharex=False, sharey=False)
    plt.title(f'{title_supp}')
    # g.savefig(
    #     f'{out_dir}violinplot_{yvar}_{xvar}_distribution_{title_supp}.png')
    save_plot_report(f"violinplot_{yvar}_{xvar}{huevar}{colvar}_distribution_{title_supp}", g, output=out_dir, file=file)


def ridgeplot(df, var, rowvar=None, huevar=None, colvar=None,  out_dir='',, file=None):
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
        save_plot_report(f'rigdeplot_{rowvar}_{var}_color{huevar}',
                                 g_ridge, output=out_dir, file=file)


def linear_reg_plot(df, xvar, yvar, huevar=None,  xticks_labels={}, leg_lab=[], title_suppl='', out_dir='', file=None):
        # L,, file=Noneinear regression
    g = sns.lmplot(data=df, x=xvar,
                   y=var, hue=huevar, legend=False)
    if xticks_labels:
        g.set(xticks=list(xticks_labels.keys()), xticklabels=list(xticks_labels.values()))
    if leg_lab: plt.legend(labels=leg_lab, loc='upper left')
    plt.title(title_suppl)
    save_plot_report(f'LinearRegPlot_{xvar}_{yvar}{huevar}', g, output=out_dir, file=file)



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


    # TODO change legend !!!!
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



def main():


    # CPG selection: the most represented
    # most_represented = median_df.groupby(['cpg', 'CHR']).count().reset_index(
    #                     ).groupby(['CHR']).first()['cpg']


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
