import sys
import os
import pandas as pd
import socket
import numpy as np
host = socket.gethostname()
if 'Fanny' in host:
    PATH_UTILS = '/home/fanny/Work/EBI/Utils'
    ABS_PATH = '/home/fanny/Work/EBI/covid_nanopore'
else:
    PATH_UTILS = '/nfs/research/birney/users/fanny/Utils'
    ABS_PATH = '/nfs/research/birney/users/fanny/covid_nanopore'
sys.path.insert(0, PATH_UTILS)
from utils import *
import argparse

def boxplot_customized(df, x_var, y_var, hue_var=None, dict_colors='tab10', width_var=None, width_dict={}, hatch_var=None, hatch_dict={}, ax=None):
    df = df[df[hue_var].isin(dict_colors.keys())]

    # Creation of the plot
    # plt.rcParams['patch.edgecolor'] = 'black'
    g2 = sns.boxplot(data=df, x=x_var, y=y_var, orient='v', hue=hue_var, palette=dict_colors, ax=ax, color='black')

    var_ls = {x for x in [x_var, hue_var, hatch_var, width_var] if x}
    value_df = df[var_ls].drop_duplicates()

    if (hatch_var != None) & (hatch_dict == {}):
        hatch_list = ['+', 'O', '.', '*', '-', '|', '/', '\ ', 'x', 'o',]
        hatch_dict = {val: hatch_list[n] for n, val in enumerate(value_df[hatch_var].unique())}
    if (width_var != None) & (width_dict == {}):
        width_dict = {val: n*2 + 1 for n, val in enumerate(value_df[width_var].unique())}

    value_df['Hue'] = value_df[hue_var].replace(dict_colors)

    value_df.index = g2.artists
    for art in g2.artists:
        if hatch_var:
            value_df['Hatch'] = value_df[hatch_var].replace(hatch_dict)
            art.set_hatch(value_df.loc[art, 'Hatch'])
        if width_var:
            value_df['Width'] = value_df[width_var].replace(width_dict)
            art.set_linewidth(value_df.loc[art, 'Width'])
    return g2


def setup_customizedboxplot_cpg_analysis(cpg_df, unit='control_snp', dir_out='.', dict_colors=None):
    # General modification of the datas
    cpg_df = cpg_df.copy()
    unit='covid_snp'
    cpg = str(cpg_df['cpg'].unique())[2:-2]
    snp = str(cpg_df[unit].unique().tolist())[2:-2]
    ref, alt = snp.split(':')[-2:]
    replace_dict = {'haplotype': {'ref': ref, 'alt': alt}, 'Genotype': {'0/0': f'{ref}/{ref}', '0/1': f'{ref}/{alt}', '1/1': f'{alt}/{alt}'}}
    cpg_df.replace(replace_dict, inplace=True)

    # Aesthetic parameters
    dict_hatch = {'Severe': '||', 'Mild': '/'}
    if not dict_colors:
        dict_colors = {'Genotype': {'0/0':'#1678F5', '0/1':'#2aa69a', '1/1':'#3ED43E'},
                       'haplotype': {'ref':'#1678F5', 'alt':'#3ED43E'},
                       'phenotype': {'Mild':'#f0f0f5', 'Severe':'#c2c2d6'}}

    # Modification of the color dataframe to fit the replacement
    repl_colors = {}
    for key, dict in dict_colors.items():
        if replace_dict.get(key):
            repl_colors.update({key: {replace_dict.get(key).get(k): dict_colors.get(key).get(k) for k in dict.keys() if replace_dict.get(key).get(k) != None}})
        else:
            repl_colors.update({key:{k: dict_colors.get(key).get(k) for k in dict.keys() }})

    # Creation of the plot
    fig, ax = plt.subplots(2, 3, figsize=(17,10))
    boxplot_customized(cpg_df, 'phenotype', 'log_lik_ratio', hue_var='phenotype',
                                 dict_colors=repl_colors['phenotype'], width_var=None,
                                 hatch_var='phenotype', ax=ax[0,0], hatch_dict=dict_hatch)
    ax[0,0].set(title='Phenotype')
    boxplot_customized(cpg_df, 'Genotype', 'log_lik_ratio', hue_var='Genotype',
                                 dict_colors=repl_colors['Genotype'],
                                 width_var=None, ax=ax[0,1], hatch_dict=dict_hatch)
    ax[0,1].set(title='Genotype')#, xticklabels=replace_val['Genotype'].values())
    boxplot_customized(cpg_df, 'phenotype', 'log_lik_ratio', hue_var='Genotype',
                                 dict_colors=repl_colors['Genotype'], width_var=None,
                                 hatch_var='phenotype', ax=ax[0,2], hatch_dict=dict_hatch)
    ax[0,2].set(title='Genotype X Phenotype')
    if not cpg_df[cpg_df['Genotype'] == replace_dict['Genotype']['0/1']].empty:
        boxplot_customized(cpg_df[cpg_df['Genotype'] == replace_dict['Genotype']['0/1']], 'phenotype',
                                       'log_lik_ratio', hue_var='phenotype',
                                       dict_colors=repl_colors['phenotype'],
                                       hatch_var='phenotype', width_var=None,
                                       ax=ax[1,0], hatch_dict=dict_hatch)
        ax[1,0].set(title='Heterozygous Phenotype')
        boxplot_customized(cpg_df[cpg_df['Genotype'] == replace_dict['Genotype']['0/1']], 'haplotype', 'log_lik_ratio',
                                      hue_var='haplotype', dict_colors=repl_colors['haplotype'],
                                      hatch_var=None, width_var=None, ax=ax[1,1],
                                      hatch_dict=dict_hatch)
        ax[1,1].set(title='Heterozygous Haplotype') #, xticklabels=replace_val['haplotype'].values())
        boxplot_customized(cpg_df[cpg_df['Genotype'] == replace_dict['Genotype']['0/1']], 'phenotype', 'log_lik_ratio',
                                      hue_var='haplotype', dict_colors=repl_colors['haplotype'],
                                      hatch_var='phenotype', width_var=None, ax=ax[1,2],
                                      hatch_dict=dict_hatch)
        ax[1,2].set(title='Heterozygous Haplotype X Phenotype')
    else:
        fig.delaxes(ax[1, 0])
        fig.delaxes(ax[1, 1])
        fig.delaxes(ax[1, 2])
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
            print(Label)
            label_title = df[df == Label[0]].dropna(how='all', axis=1).columns.tolist()[0]
            print(label_title)
            if label_title == 'phenotype':
                fig.legend(handles=[mpatches.Patch(facecolor=val, label=key, lw=1, ec=sns.color_palette('dark')[-3],
                                         hatch=dict_hatch[key]) for key, val in dict_colors['phenotype'].items()],
                                         title='Phenotype', loc='center left', bbox_to_anchor=(0.9, x))
            elif label_title == 'base_called':
                label_title = 'Haplotype'
                fig.legend(Line, Label, loc='center left', bbox_to_anchor=(0.9, x), title=label_title)
            else:
                fig.legend(Line, Label, loc='center left', bbox_to_anchor=(0.9, x), title=label_title)
            x = x - 0.08


    save_plot_report(f'{dir_out}/Multiplots_{cpg}.png', fig, file=None)
    # fig.savefig(f'{dir_out}/Multiplots_{cpg}.png')


def MannWhitney_Spearman_stats(df, measure, vars,  output='', add_col={}, pval=0.05):
    if not df.empty:
        stat = Stats(df)
        try:
            res_mw = stat.multiple_mann_whitney(measure=measure, var=vars)
            if not res_mw.empty:
                if add_col: res_mw[list(add_col.keys())] = list(add_col.values())
                res_mw[['index', 'p-val', 'cpg', 'data']].dropna(
                    subset=['p-val']).to_csv(output + 'Mann_Whitney.csv',
                                             mode='a', header=False, index=False)
        except Exception as err:
            print('MWU', err)
        try:
            res_sp = stat.Spearman_correlation(measure=measure, var=vars)
            if add_col: res_sp[list(add_col.keys())] = list(add_col.values())
            res_sp.dropna(subset=['p-val']).reset_index().to_csv(output+'Spearmann_corr.csv',
                                                mode='a', header=False, index=False)
        except Exception as err:
            print('SP', err)


def count_and_mean_values(df, vars, mean_measure):
    count_mean = pd.Series(dtype='int')
    for col in vars:
        count_mean = count_mean.append(df[col].value_counts())
    alt_vs_ref = df[df['haplotype'] == 'ref'][mean_measure].mean() - df[df['haplotype'] == 'alt'][mean_measure].mean()
    count_mean = count_mean.append(pd.Series(alt_vs_ref, index=['means_ref-alt']))
    return count_mean


def Loop_stats(df, unit, output='', phen_ls=['Severe', 'Mild'], gen_ls=['0/1'], cpg_ls=[]):

    if not cpg_ls:
        cpg_ls = df['cpg'].unique()
    count_cols = ['cpg', '0/0', '0/1', '1/1', 'alt', 'ref', 'Mild', 'Severe', 'means_ref-alt']
    phen_ls=['Severe', 'Mild']
    gen_ls=['0/1']
    print(phen_ls, gen_ls)
    for cpg in cpg_ls:

        if cpg == cpg_ls[0]:
            pd.DataFrame(columns=['index', 'p-val', 'cpg', 'data']).to_csv(output + 'Mann_Whitney.csv', mode='w', index=False)
            pd.DataFrame(columns=['index', 'rho', 'p-val', 'cpg', 'data']).to_csv(output + 'Spearmann_corr.csv', mode='w', index=False)
            pd.DataFrame(columns=count_cols).to_csv(f'{output}Counts_Diff_means.csv', mode='w', header=True, index=False)

        samp = df[df['cpg'] == cpg].copy()
        samp = samp[samp.drop(unit, axis=1).duplicated() == False]
        if (samp.shape[0] > 1):
            if not np.isnan(samp['log_lik_ratio'].std()) or (samp['log_lik_ratio'].std() == 0):

                # Counts
                count_mean = count_and_mean_values(samp, mean_measure='log_lik_ratio', vars=['Genotype', 'phenotype', 'haplotype'])
                count_mean = pd.DataFrame(count_mean, columns=[cpg], index=count_cols).T
                count_mean.reset_index().to_csv(f'{output}Counts_Diff_means.csv', mode='a', index=False, header=False)

                # Mann Whitney & Spearman tests
                MannWhitney_Spearman_stats(samp, measure='log_lik_ratio' ,
                    vars=['Genotype', 'phenotype'], output=output, add_col={'cpg': cpg,
                                                                            'data': 'ALL'})
                for phen in phen_ls:
                    MannWhitney_Spearman_stats(samp[samp['phenotype'] == phen], measure='log_lik_ratio',
                        vars=['haplotype', 'Genotype'], output=output, add_col={'cpg': cpg,
                                                                          'data': phen})
                for gen in gen_ls:
                    MannWhitney_Spearman_stats(samp[samp['Genotype'] == gen], measure='log_lik_ratio',
                        vars=['haplotype', 'phenotype'], output=output, add_col={'cpg': cpg,
                                                                           'data': gen})


# MAIN
def main(file, unit, output_dir='.', gen_ls=['0/1'], phen_ls=['Mild', 'Severe'], list_snp_file=''):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Opening file
    all_df = pd.read_csv(file, dtype='object')

    # TODO ERASE ???:
    print('Colnames in df: ', all_df[(all_df == all_df.columns)].dropna(how='all').shape[0], flush=True)
    all_df = all_df[(all_df == all_df.columns) == False].dropna(how='all')
    print('\nNumber of rows with wrong genotype: ', all_df[all_df['Genotype'].isin(['0/0', '0/1', '1/1']) == False].shape[0], flush=True)
    all_df = all_df[all_df['Genotype'].isin(['0/0', '0/1', '1/1'])]
    if list_snp_file:
        snp_ls = pd.read_table(list_snp_file, header=None)[0].tolist()
        print('\nNumber of rows with wrong SNPs: ', all_df[(all_df[unit].isin(snp_ls) == False)].shape[0], flush=True)
        all_df = all_df[(all_df[unit].isin(snp_ls))]
    print('\nNumber of duplicated lines: ', all_df[all_df.duplicated(keep=False)].shape[0], flush=True)
    del snp_ls
    all_df = all_df[all_df.duplicated() == False]
    all_df['log_lik_ratio'] = all_df['log_lik_ratio'].astype(float)

    # QUESTION: Where should we do the outliers detection ?
    stat_raw = Stats(all_df[['log_lik_ratio']])
    print('\nOutliers:', stat_raw.outliers.shape[0]/all_df[['log_lik_ratio']].shape[0], flush=True)
    all_df = all_df.loc[stat_raw.no_outliers.index]

    # Median over the read_name with same Allele calling
    # QUESTION: Before doing the median, should we delete the cpgs/snps with not enough counts for one file ?
    median_df = all_df.groupby(['phenotype', 'sample_id', 'chromosome', 'cpg',
                                    unit, 'Genotype', 'haplotype']).median().reset_index()
    del all_df
    median_df = median_df.sort_values(by=['haplotype'], ascending=False)
    median_df = median_df.sort_values(by=['chromosome', unit, 'Genotype', 'phenotype'])

    # Filter
    # snp_counts = median_df.groupby([unit, 'Genotype']).size().unstack()
    # cpg_counts = median_df.groupby(['cpg', unit, 'Genotype']).size().unstack()
    # cpg_counts_10 = cpg_counts[cpg_counts > 10].dropna(how= 'all').index.levels[0]

    # Loop_stats(median_df, unit, gen_ls=gen_ls, phen_ls=phen_ls, output=f'{os.path.join(output_dir, os.path.basename(file)[:-4])}_')
    cpg_ls = ['21:33243685:1', '21:33247116:1', '17:45501201:2', '11:16996188:1',
              '1:85220632:1', '15:43591259:1', '12:112925195:1', '19:10354080:3',
              '10:79535851:1', '6:31133406:2', '9:133262493:1', '13:27151562:1',
              '19:50314700:1', '15:43691102:3', '3:45847444:1', '18:62362967:1',
              '1:155220731:1', '13:41393886:2', '12:112935779:1', '18:70635452:1',
              '17:81597736:2', '10:79491025:1', '11:72864523:1', '6:31142922:2',
              '9:133286263:1', '14:62083099:1', '3:45847362:1', '14:103504667:4',]
    cpg_ls += ['21:33242527:1', '17:46065976:1', '21:33229986:3', '3:45848456:1',
               '12:112925744:1', '17:46768336:3', '9:133271878:1', '1:155209089:1',
               '9:133271842:1', '6:33069193:2', '21:33226777:1']

    for cpg in set(cpg_ls):
        cpg_df = median_df[median_df['cpg'] == cpg]
        if not cpg_df.empty:
            setup_customizedboxplot_cpg_analysis(cpg_df, unit=unit, dir_out='plots', dict_colors=None)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='STEP2 - Pipeline for mQTLs'
        + ' - Statisques (Spearman correlation and MannWhitney test) on specific datasets')
    parser.add_argument('file', type=str, help='filtered working file')
    parser.add_argument('-o', '--output_dir', type=str, default='',
                        help='directory to save output files')
    parser.add_argument('-u', '--unit', type=str, default='snp',
                        help='unit to perform analysis')
    parser.add_argument('-g', '--gen_ls', type=str, default='0/1',
                        help='list of genotypes to perform specific analysis on')
    parser.add_argument('-p', '--phen_ls', type=str, default='Mild-Severe',
                        help='list of phenotypes to perform specific analysis on')
    parser.add_argument('-f', '--list_snp_file', type=str, default=f'{ABS_PATH}/finemapped',
                        help='list of snp to verify validity of snp list we got')
    args = parser.parse_args()
    main(**vars(args))


def joint_model():
    pass
    # TODO: Joint model
    # transform to normal distrib inverse norm)
    # INVERSE NORMALISATION
    # JOIN MODEL
    # > anova(vsx2_linear,vsx2_interaction, test = "F")
    # Analysis of Variance Table
    #
    # Model 1: invnorm(ONL) ~ standing_height_0_0 + weight_0_0 + age_when_attended_assessment_centre_0_0 +
    # genetic_sexuses_datacoding_9_0_0 + as.factor(opticalcoherence_tomography_device_id_0_0) +
    # PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
    # PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 +
    # PC20 + (vsx2_locus) + rs375435
    # Model 2: invnorm(ONL) ~ standing_height_0_0 + weight_0_0 + age_when_attended_assessment_centre_0_0 +
    # genetic_sexuses_datacoding_9_0_0 + as.factor(opticalcoherence_tomography_device_id_0_0) +
    # PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
    # PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 +
    # PC20 + (vsx2_locus) * rs375435
    # Res.Df RSS Df Sum of Sq F Pr(>F)
    # 1 31061 29341
    # 2 31054 29314 7 27.003 4.0866 0.0001714 ***
    # ---
    # Signif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

    # lm(medianlog ~ snp * clinical, data)
    # lm(invnorm(medianlog) ~ snp + clinical + snp * clinical, data)
    # lm(invnorm(median_log_hap) ~ allele * clinical , haplotype_data)





# CPGs of interest --> save the median table just for those
# TODO: Automate the finding of cpgs of interest ! eg :
# cpgs_interest = stat[(stat['Counts 0/0'] > 3) & (stat['Counts 1/1'] > 3) & (stat['Counts 0/1'] > 3) & (stat['Spearman correlation p_value'] < 1e-5)]['cpg'].unique()
# highest p-val for highmildhighsev : ['17:46065976:1', '12:112942465:1', '21:33229986:3'] # highest rho
# dict_cpgs = {#'INTEREST': ['17:46768336:3', '6:33069193:2'],
#              # 'HIGHmildHIGHsev': ['12:112925744:1', '17:46065976:1', '21:33229986:3'],
#              # 'EXTRA': ['3:45848456:1', '21:33226777:1', '21:33242527:1'],
#              'Hsev_Lmild': ['9:133271878:1','1:155209089:1'],
#              'Lsev_Hmild': ['9:133271842:1']}

# CPGs interesting from MILD vs SEVERE analysis
# merge = pd.merge(stat_mild, stat_sev, on=['cpg', 'SNP'], suffixes=['_mild', '_sev'])
# index_counts = merge.filter(regex='Counts*')[merge.filter(regex='Counts*') < 3].dropna(thresh=6).index
# merge.drop(index_counts, inplace=True)
# merge.loc[merge['Spearman correlation p_value_mild'] < 1e-5, 'log_mild'] = 'High'
# merge.loc[merge['Spearman correlation p_value_sev'] < 1e-5, 'log_sev'] = 'High'
# merge.loc[merge['Spearman correlation p_value_sev'] > 1e-5, 'log_sev'] = 'Low'
# merge.loc[merge['Spearman correlation p_value_mild'] > 1e-5, 'log_mild'] = 'Low'
# cpg_highsev_lowmild = merge[(merge['log_sev'] == 'High') & (merge['log_mild'] == 'Low')]['cpg'].unique()
# cpg_lowsev_highmild = merge[(merge['log_sev'] == 'Low') & (merge['log_mild'] == 'High')]['cpg'].unique()
# merge[(merge['log_sev'] == 'High') & (merge['log_mild'] == 'Low') & (merge['Spearman correlation rho_mild'] > -0.25) &(merge['Spearman correlation rho_mild'] < 0) & (merge['Spearman correlation rho_sev'] < -0.5) & (merge['Spearman correlation rho_sev'] < 0)]['cpg'].unique()
# merge[(merge['log_sev'] == 'Low') & (merge['log_mild'] == 'High') & (merge['Spearman correlation rho_sev'] > -0.25) &(merge['Spearman correlation rho_sev'] < 0) & (merge['Spearman correlation rho_mild'] < -0.5) & (merge['Spearman correlation rho_mild'] < 0)]['cpg'].unique()
# merge['CHR'] = merge['SNP'].str.split(':', expand=True)[0].astype('category')
# merge[(merge['log_sev'] == 'High') & (merge['log_mild'] == 'High')].sort_values(by=['Spearman correlation rho_sev', 'Spearman correlation rho_mild']).groupby('CHR').head(1)['cpg'].tolist()
# merge[(merge['log_sev'] == 'High') | (merge['log_mild'] == 'High')].sort_values(by=['Spearman correlation rho_sev', 'Spearman correlation rho_mild']).groupby('CHR').head(1)['cpg'].tolist()
# index_rhoneg = cpg_highsev_lowmild[(cpg_highsev_lowmild['Spearman correlation rho_sev'] < -0.25)
#                         & (cpg_highsev_lowmild['Spearman correlation rho_mild'] > -0.35)
#                         ].dropna(how='all').dropna(how='all', axis=1).index
# cpg_highsev_lowmild.loc[index_rhoneg]['cpg'].unique().tolist()
