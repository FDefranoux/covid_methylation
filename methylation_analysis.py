import pandas as pd


def spearmanRho_diffmean_plot(stat, unit='cpg', hue_var='CHR', col_var=None, out_dir=''):
    stat = stat.copy()
    stat = stat.reset_index().rename(columns={unit:unit})
    stat['CHR'] = stat[unit].str.split(':', expand=True)[0].astype('category')
    g1 = sns.relplot(kind='scatter', data=stat, y='diff_means_altVSref',
                     x='Spearman correlation Gen rho', hue=hue_var, col=col_var, col_wrap=3)
    g1.savefig(f'{out_dir}/Diff_meansVSrho_{unit}_heterozygotes_sepCHR.png')


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
        stat.sort_values('CHR', inplace=True)

        # Format X-axis
        if format_xaxis:
            n=0
            for chr in stat['CHR'].sort_values().unique():
                var_df = stat.loc[stat['CHR'] == chr, xvar].sort_values()
                for val in var_df.unique():
                    stat.loc[(stat['cpg'] == val)
                        & (stat['CHR'] == chr), f'{xvar}'] = n
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
        g.set(xlabel="CHR", xticks=stat.groupby(['CHR']).last()[xvar].unique(),
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


def setup_customizedboxplot_cpg_analysis(cpg_df, dir_out=''):
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
    # mann_whit = pd.read_csv('FROZEN_Nov2021_cpg_snp_analysis/Mann_whitney_table.csv')
    # pval_plot(mann_whit, 'cpg', 'MWU Mild-Severe', pval_cutoff=1, n_site=2, title_supp='MannWhitney_GenPhen_sorted', out_dir=dir_out, format_xaxis=True)
    # pval_plot(mann_whit, 'cpg', 'MWU Mild-Severe', pval_cutoff=1, n_site=2, title_supp='MannWhitney_GenPhen', out_dir=dir_out, format_xaxis=False)
    # pval_plot(stat, 'cpg', 'Spearman correlation p_value', pval_cutoff=1, n_site=2, title_supp='Spearman_sorted', out_dir=dir_out, format_xaxis=True)
    # pval_plot(stat, 'cpg', 'Spearman correlation p_value', pval_cutoff=1, n_site=2, title_supp='Spearman', out_dir=dir_out, format_xaxis=False)
    # spearman_correlation_plot(stat[stat['minus_log10'].isna()==False], unit=unit, n_site=2, out_dir=dir_out)

    # INDIVIDUAL PLOTS FOR CPGs of interest
    # highest p-val for highmildhighsev : ['17:46065976:1', '12:112942465:1', '21:33229986:3'] # highest rho
    # dict_cpgs = {#'INTEREST': ['17:46768336:3', '6:33069193:2'],
    #              # 'HIGHmildHIGHsev': ['12:112925744:1', '17:46065976:1', '21:33229986:3'],
    #              # 'EXTRA': ['3:45848456:1', '21:33226777:1', '21:33242527:1'],
    #              'Hsev_Lmild': ['9:133271878:1','1:155209089:1'],
    #              'Lsev_Hmild': ['9:133271842:1']}
    #
    # for supp_title, list_cpg in dict_cpgs.items():
    #     for cpg in list_cpg:
    #         cpg_df = median_df[median_df['cpg'] == cpg].copy()
    #         cpg_df['Gen'].replace({'alt': 1, 'ref':0}, inplace=True)
    #         setup_customizedboxplot_cpg_analysis(cpg_df, dir_out=dir_out)




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

# PLOT CPG INTEREST
# dict_cpgs = {'INTEREST': ['17:46768336:3', '6:33069193:2'],
#              'HIGHmildHIGHsev': ['12:112925744:1', '17:46065976:1', '21:33229986:3'],
#              'EXTRA': ['3:45848456:1', '21:33226777:1', '21:33242527:1']}
# for supp_title, list in dict_cpgs.items():
#
#     df = median_df[median_df['cpg'].isin(list)]
#     GENERAL PLOT GENOTYPE
#     g = sns.catplot(data=df,
#                     y='log_lik_ratio', col='cpg', col_wrap=5,
#                     x='Genotype', orient='v', kind= 'box',
#                     height=6, aspect=0.9, hue='phenotype',
#                     sharex=False, sharey=False)
#     g.savefig(f'{dir_out}/Boxplot_median_ratio_GenPhen_all_{supp_title}.png')
#
#     # GENERAL PLOT HETEROZYGOTES
#     g = sns.catplot(data=df[df['Genotype'] == '0/1')],
#                     y='log_lik_ratio', col='cpg', col_wrap=5,
#                     x='Gen', orient='v', kind= 'box',
#                     height=6, aspect=0.9, hue='phenotype',
#                     sharex=False, sharey=False, )
#     g.savefig(f'{dir_out}/Boxplot_median_ratio_GenPhen_all_{supp_title}_het.png')

# def other_plots():
#
#             # Individual per genotype
#             # 3 way genotype
#             # g1 = sns.catplot(data=cpg_df, y='log_lik_ratio',
#             #                 x='Genotype', orient='v', kind= 'box',
#             #                 height=6, aspect=0.9,
#             #                 sharex=False, sharey=False)
#             # g1.set(title=f'CpG {cpg} associated with SNP {snp}')
#             # g1.savefig(f'{dir_out}/Boxplot_median_cpg_ratio_Gen_{supp_title}_{cpg}.png')
#             # 2 way phenotype
#             # g2 = sns.catplot(data=cpg_df, y='log_lik_ratio',
#             #                 x='phenotype', orient='v', kind= 'bar',
#             #                 height=6, aspect=0.9,
#             #                 sharex=False, sharey=False, cmap='viridis')
#             # g2.set(title=f'CpG {cpg} associated with SNP {snp}')
#             # g2.savefig(f'{dir_out}/Boxplot_median_cpg_ratio_Phen_{supp_title}_{cpg}.png')
#             # 6 way
#             # g3 = sns.catplot(data=cpg_df, y='log_lik_ratio',
#             #                 x='Genotype', orient='v', kind= 'box',
#             #                 height=6, aspect=0.9, hue='phenotype',
#             #                 sharex=False, sharey=False)
#             # g3.set(title=f'CpG {cpg} associated with SNP {snp}')
#             # g3.savefig(f'{dir_out}/Boxplot_median_cpg_ratio_GenPhen_{supp_title}_{cpg}.png')
#
#             # Individual for heterozygotes
#             # 4 way
#             # g4 = sns.catplot(data=cpg_df[cpg_df['Genotype'] == '0/1'],
#             #                 y='log_lik_ratio', x='Gen', kind= 'box',
#             #                 height=6, aspect=0.9, hue='phenotype',
#             #                 sharex=False, sharey=False)
#             # g4.set(title=f'CpG {cpg} associated with SNP {snp}', xlabel='Allele')
#             # g4.savefig(f'{dir_out}/Boxplot_median_cpg_ratio_GenPhen_{supp_title}_{cpg}_het.png')
#
#             # # 2 way allele
#             # g5 = sns.catplot(data=cpg_df[cpg_df['Genotype'] == '0/1'],
#             #                  y='log_lik_ratio', x='Gen', kind= 'box',
#             #                  height=6, aspect=0.9,
#             #                  sharex=False, sharey=False)
#             # g5.set(title=f'CpG {cpg} associated with SNP {snp}', xlabel='Allele')
#             # g5.savefig(f'{dir_out}/Boxplot_median_cpg_ratio_Gen_{supp_title}_{cpg}_het.png')
#
#             # 2 way phenotype
#         #     g6 = sns.catplot(data=cpg_df[cpg_df['Genotype'] == '0/1'],
#         #                     y='log_lik_ratio', x='phenotype', kind= 'box',
#         #                     height=6, aspect=0.9,
#         #                     sharex=False, sharey=False, cmap='viridis')
#         #     g6.set(title=f'CpG {cpg} associated with SNP {snp}')
#         #     g6.savefig(f'{dir_out}/Boxplot_median_cpg_ratio_Phen_{supp_title}_{cpg}_het.png')
#         #
#         # except Exception as err:
#         #     print(f'ERROR WITH cpg {cpg} ', err)
#
#
#     # Plot all SNP, discarding the very different ones.
#     for chr in median_df['CHR'].sort_values().unique():
#         print('CHR = ', chr)
#         cleaned_per_chr, _ = drop_datas(median_df[median_df['CHR'] == chr], 'name', thresh_zscore=2.5)
#         try:
#             violinplot(cleaned_per_chr, title_supp=str(chr)+'_removedextreme', yvar='cpg',
#                        xvar='log_lik_ratio', huevar='Genotype', colvar=None)
#         except ValueError:
#             violinplot(median_df[median_df['CHR'] == chr], title_supp=str(chr)+'_all', yvar='cpg',
#                        xvar='log_lik_ratio', huevar='Genotype', colvar=None)
#         # violinplot(cleaned_per_chr, chr=chr, yvar='SNP',
#         #            xvar='log_lik_ratio', huevar='Genotype', colvar='phenotype')
#
#     # Plotting the most represented
#     most_represented = median_df.groupby(['cpg', 'CHR']).count().reset_index(
#                         ).groupby(['CHR']).first()['cpg']
#     most_represented_df = median_df[median_df['cpg'].isin(most_represented)]
#     violinplot(most_represented_df, title_supp=str(chr)+'_most_rep', yvar='cpg',
#         xvar='log_lik_ratio', huevar='Genotype', colvar='phenotype')
#
#     # Mean/Median over binned sites
#
#     # binned = median_df.copy()
#     binned, _ = drop_datas(median_df, 'name', thresh_zscore=3)
#     binned['pos'] = binned['cpg'].str.split(':', expand=True)[1].astype(int)
#
#     # binned['distance_to_median'] = binned['cpg'].str.split(':', expand=True)[1].astype(int)
#     for chr in binned['CHR'].sort_values().unique():
#         chr_index = binned[binned['CHR'] == chr].index
#         binned.loc[chr_index, 'bin'] = (
#             (binned.loc[chr_index, 'pos'] - binned.loc[chr_index, 'pos'].median()
#             )/binned.loc[chr_index, 'pos'].median()).round(1)
#         # print(f'CHR = {chr}\n\n', binned['bin'] .value_counts().to_markdown())
#         # print(binned['pos'] .value_counts().to_markdown())
#     binned['bin'] = binned['CHR'].astype(str) + ' : ' + binned['bin'].astype(str)
#     violinplot(binned, title_supp=str(chr)+'_per_bin', yvar='bin',
#         xvar='log_lik_ratio', huevar='Genotype', colvar='phenotype')
#
#         # binned.loc[binned['CHR'] == chr, 'bin'] = pd.cut(
#         #     x=binned[binned['CHR'] == chr]['cpg'].str.split(
#         #         ':', expand=True)[1].astype(int), bins=4).astype(str)
#         # binned.loc[binned['CHR'] == chr, 'bin'] = pd.qcut(
#         #     x=binned[binned['CHR'] == chr]['cpg'].str.split(
#         #         ':', expand=True)[1].astype(int), q=3, duplicates='drop').astype(str)
#     binned['bin'] = binned['CHR'].astype(str) + ': ' +  binned['bin']
#     binned.sort_values('cpg', inplace=True)
#     violinplot(binned, title_supp=str(chr)+'_per_bin', yvar='bin',
#         xvar='log_lik_ratio', huevar='Genotype', colvar=None)
#
#     scatter_with_unique_cpg(median_df, huevar='phenotype',
#                             colvar='Genotype', xvar='cpg', yvar='log_lik_ratio')
#     scatter_with_unique_cpg(median_df, huevar='phenotype',
#                             colvar='Genotype', xvar='SNP', yvar='log_lik_ratio')
#     # scatter_with_unique_cpg(median_df, huevar='Gen',
#     #                         colvar=None, xvar='cpg', yvar='log_lik_ratio')
#     for chr in median_df['CHR'].sort_values().unique():
#         violinplot(median_df[median_df['CHR'] == chr], title_supp=chr, yvar='cpg',
#                    xvar='log_lik_ratio', huevar='Genotype', colvar='phenotype')
#         violinplot(median_df[median_df['CHR'] == chr], title_supp=chr, yvar='SNP',
#                    xvar='log_lik_ratio', huevar='Genotype', colvar='phenotype')
#         # ridgeplot(median_df[median_df['CHR'] == 1], chr=1, rowvar='SNP',
#         #           huevar='Genotype', colvar='phenotype', var='log_lik_ratio')
