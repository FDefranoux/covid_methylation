import sys
import os
import pandas as pd
import socket
import numpy as np
import sqlite3
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
# import seaborn as sns
# import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches
# from matplotlib.offsetbox import (AnchoredOffsetbox, DrawingArea, HPacker,
#                                   TextArea)
import pingouin as pg
from tqdm import tqdm
OUTDATED_IGNORE=1

# TODO: Modify the stat function to have the double ttests...DOUBLE ???


def stat_test(df, x, y, test, n_min=5, data='', **kwargs):
    dict_rename = {'corr': [pg.corr, {'r':'coeff'}],
                   'mann_whitney': [pg.mwu, {'U-val': 'coeff'}],
                   'ttest': [pg.ttest, {'T':'coeff'}]}

    if (df.empty) or (df[x].std() == 0) or (df[x].nunique() < 2) or (df[y].std() == 0) or (df[y].nunique() < 2) or (df.shape[0] < n_min):
       res = pd.DataFrame(columns=['cpg', 'data', 'coeff', 'p-val', 'x', 'y'])
    else:
        try:
                res = dict_rename[test][0](df[x], df[y],**kwargs
                    ).rename(columns=dict_rename[test][1])
                res['x'] = x
                res['y'] = y
                res['cpg'] = list(df['cpg'].unique())[0]
                res['data'] = data
        except Exception as err:
            print(err)
            print(df)
            res = pd.DataFrame(columns=['cpg', 'data', 'coeff', 'p-val', 'x', 'y'])

    return res[['cpg', 'data', 'coeff', 'p-val', 'x', 'y']].reset_index()


def Stats_MultiTests(df, measure, target_snp, cpg_ls=[], output='Multi_Stats', n_min=5, pval=0.05):
    df = df.copy()
    df = df.sort_values(by=['Genotype', 'phenotype'])
    df = df.sort_values(by=['haplotype'], ascending=False)
    for col in ['Genotype', 'phenotype', 'haplotype']:
        df[f'{col}_test'] = pd.factorize(df[col])[0]
    if not cpg_ls:
        cpg_ls = df['cpg'].unique()
    results = pd.DataFrame()
    print('Running through the CpGs\n')
    for cpg in tqdm(cpg_ls):
        samp = df[df['cpg'] == cpg].drop(target_snp, axis=1).drop_duplicates()
        if (samp.shape[0] > 1):
            results = results.append(stat_test(samp, 'Genotype_test', measure, test='corr', data='all', method='spearman'))
            results = results.append(stat_test(samp, 'phenotype_test', measure, test='mann_whitney', data='all'))

        # Mann_Whitney tests on heterozygotes
        samp_het = samp[samp['Genotype'] == '0/1']
        if (samp_het.shape[0] > 1):
            results = results.append(stat_test(samp_het, 'phenotype_test', measure, test='mann_whitney', data='het'))
            results = results.append(stat_test(samp_het, 'haplotype_test', measure, test='mann_whitney', data='het'))
        del samp_het

        # T-test between alt and ref for Mild and Severe specific
        for val in samp['phenotype'].unique():
            samp_ph = samp[samp['phenotype'] == val].drop_duplicates()
            results = results.append(stat_test(samp_ph, 'Genotype_test', measure, test='corr', data=val, method='spearman'))
            samp_hetph = samp_ph[samp_ph['Genotype'] == '0/1'].drop_duplicates()
            samp_hetph = samp_hetph.set_index(['cpg', 'sample_id', 'phenotype', 'Genotype', 'haplotype'])[measure].unstack().dropna(thresh=2).reset_index()
            res = stat_test(samp_hetph, 'alt', 'ref', test='ttest', data=f'het-{val}', n_min=n_min)
            results = results.append(res)

    # Formatting forsaving
    results.rename(columns={'index':'stat', 'x':'variable', 'y':'measure'}, inplace=True)
    results.loc[results['variable'] == 'alt', ['variable', 'measure']] = ['alt-ref', measure]
    results.to_csv(output + '.csv', index=False)

    # Prinout of significative results
    print(f'Number of significative cpgs at pvalue threshold {pval} (non-corrected)\n')
    print(results[results['p-val'] < pval].astype(str).groupby(['variable', 'stat', 'data']).size().unstack().sort_index(axis=1, key=lambda x: x.str.lower()).reset_index(level='stat').reindex(['phenotype_test', 'Genotype_test', 'haplotype_test', 'alt-ref']).reset_index().fillna('').to_markdown())


### OLD FUNCTIONS
def MannWhitney_Spearman_stats(df, measure, vars,  output='', add_col={}, pval=0.05):
    df = df.copy()
    if not df.empty:
        df.sort_values('Genotype', inplace=True)
        stat = Stats(df)
        stat.define_dummies(vars)
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
            for var in vars:
                print(stat.dum[var].unique(), df[var].unique())
                res_sp = pg.corr(stat.dum[var], df[measure], method="spearman").drop('CI95%', axis=1)
                add_col['var'] = var
                # res_sp = stat.Spearman_correlation(measure=measure, var=vars)
                if add_col: res_sp[list(add_col.keys())] = list(add_col.values())
                print(res_sp.head())
                res_sp.dropna(subset=['p-val']).reset_index().to_csv(output+'Spearmann_corr.csv',
                                                    mode='a', header=False, index=False)
        except Exception as err:
            print('ERR SP', err)


def count_and_mean_values(df, vars, mean_measure):
    count_mean = pd.Series(dtype='int')
    for col in vars:
        count_mean = count_mean.append(df[col].value_counts())
    alt_vs_ref = df[df['haplotype'] == 'ref'][mean_measure].mean() - df[df['haplotype'] == 'alt'][mean_measure].mean()
    count_mean = count_mean.append(pd.Series(alt_vs_ref, index=['means_ref-alt']))
    return count_mean


def Loop_stats(df, target_snp, output='', phen_ls=['Severe', 'Mild'], gen_ls=['0/1'], cpg_ls=[]):

    if not cpg_ls:
        cpg_ls = df['cpg'].unique()
    count_cols = ['cpg', '0/0', '0/1', '1/1', 'alt', 'ref', 'Mild', 'Severe', 'means_ref-alt']
    phen_ls=['Severe', 'Mild']
    gen_ls=['0/1']
    print(phen_ls, gen_ls)
    for cpg in cpg_ls:

        if cpg == cpg_ls[0]:
            pd.DataFrame(columns=['test', 'pval', 'cpg', 'data']).to_csv(output + 'Mann_Whitney.csv', mode='w', index=False)
            pd.DataFrame(columns=['n', 'rho', 'pval', 'power', 'cpg', 'data', 'test']).to_csv(output + 'Spearmann_corr.csv', mode='w', index=False)
            pd.DataFrame(columns=count_cols).to_csv(f'{output}Counts_Diff_means.csv', mode='w', header=True, index=False)

        samp = df[df['cpg'] == cpg].copy()
        samp = samp[samp.drop(target_snp, axis=1).duplicated() == False]
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
def main(file_db, target_snp='covid_snp', chr=None, output_dir='.', n_min=5, pval=0.05, zscore=3):
    # Create output path
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Opening file
    # db = sqlite3.connect(file_db)
    # if chr:
    #     all_df = pd.read_sql(f"SELECT * FROM median_datas WHERE chromosome IS {chr}", con=db)
    # else:
    #     all_df = pd.read_sql("SELECT * FROM median_datas", con=db)


#### Temporary solution
    file_db = '/home/fanny/Work/EBI/covid_nanopore/covid_snp_March2022/Filtered_nano_bam_files.csv'
    target_snp='covid_snp'
    chr = 12
    all_df = pd.read_csv(file_db)
    all_df.groupby('chromosome').size()
    all_df=all_df[all_df['chromosome']== chr]
    print(all_df.head(2).to_markdown())
    # all_df.drop('Unnamed: 0', axis=1, inplace=True)

    all_df = all_df.groupby(['phenotype', 'sample_id', 'chromosome', 'cpg', target_snp, 'Genotype', 'haplotype']).median().reset_index()

    # OUTLIERS DETECTION
    # QUESTION: Where should we do the outliers detection ?
    stat_raw = Stats(all_df[['log_lik_ratio']])
    stat_raw.no_outliers(zscore_thresh=zscore)
    print('\nOutliers:', round(stat_raw.outliers.shape[0]/all_df[['log_lik_ratio']].shape[0], 3), '%', flush=True)
    all_df = all_df.loc[stat_raw.no_outliers.index]
    Stats_MultiTests(all_df, 'log_lik_ratio', 'covid_snp', cpg_ls=[], output=os.path.join(output_dir, f'Multi_Stats_{chr}'), n_min=n_min, pval=pval)

    # Old function
    # Loop_stats(all_df, target_snp, gen_ls=['0/1'], phen_ls=['Mild', 'Severe'], output=f'{os.path.join(output_dir, os.path.basename(file)[:-4])}_')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='STEP2 - Pipeline for mQTLs'
        + ' - Statisques (Spearman correlation and MannWhitney test) on specific datasets')
    parser.add_argument('file_db', type=str, help='filtered db')
    parser.add_argument('-t', '--target_snp', type=str, default='covid_snp',
        help='target snp on which to perform analysis')
    parser.add_argument('-c', '--chr', type=int, default=None,
                        help='chromosome selection')
    parser.add_argument('-o', '--output_dir', type=str, default='.',
                        help='directory to save output files')
    parser.add_argument('-n', '--n_min', type=int, default=5,
                        help='number of minimum replicates')
    parser.add_argument('-p', '--pval', type=float, default=0.05,
                        help='pval threshold')
    parser.add_argument('-z', '--zscore', type=int, default=3,
        help='zscore value for outlier detection')
    args = parser.parse_args()
    print('## Statistical analysis. \nThe arguments are: ', vars(args), '\n')
    main(**vars(args))


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
