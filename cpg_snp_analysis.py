import sys
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

file = 'Filtered_finemapped.csv'
unit='covid_snp'


def MannWhitney_Spearman_stats(df, measure, vars,  output='', add_col={}, pval=0.05):
    if not df.empty:
        stat = Stats(df)
        # try:
        #     res_mw = stat.multiple_mann_whitney(measure=measure, var=vars)
        #     if not res_mw.empty:
        #         if add_col: res_mw[list(add_col.keys())] = list(add_col.values())
        #         res_mw[['index', 'p-val', 'cpg', 'data']].dropna(
        #             subset=['p-val']).to_csv(output + 'Mann_Whitney.csv',
        #                                      mode='a', header=False, index=False)
        # except Exception as err:
        #     print('MWU', err)
        try:
            res_sp = stat.Spearman_correlation(measure=measure, var=vars)
            if add_col: res_sp[list(add_col.keys())] = list(add_col.values())
            res_sp.dropna(subset=['p-val']).to_csv(output+'Spearmann_corr.csv',
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


def Loop_stats(df, output='', phen_ls=['Severe', 'Mild'], gen_ls=['0/1'], cpg_ls=[]):
    # TODO: Include automatic selection of the cpg of interest?
    # df = median_df.copy()
    # output=''

    if not cpg_ls:
        cpg_ls = df['cpg'].unique()
    count_cols = ['0/0', '0/1', '1/1', 'alt', 'ref', 'Mild', 'Severe', 'means_ref-alt']
    for cpg in cpg_ls:

        if cpg == cpg_ls[0]:
            pd.DataFrame(columns=['index', 'p-val', 'cpg', 'data']).to_csv(output + 'Mann_Whitney.csv', mode='w', index=False)
            pd.DataFrame(columns=['index', 'variable', 'value', 'cpg', 'data']).to_csv(output + 'Spearmann_corr.csv', mode='w', index=False)
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
def main(file, unit):

    # TODO: Automation of the path management (in, out, frozen etc)
    # TODO: Include title in the MAIN
    # TODO: Analysis with yaml software, moved in the result folder (with date and info for analysis)
    # if not os.path.exists(dir_out):
    #     os.makedirs(dir_out)

    # Opening file
    all_df = pd.read_csv(file, dtype='object')
    # TODO ERASE:
    all_df = all_df[(all_df == all_df.columns) == False].dropna(how='all')

    # TODO: check if we still need filtering at this point
    snp_ls = pd.read_table('base_called_from_bam_files/finemapped', header=None)[0].tolist()
    print('\nNumber of rows with wrong genotype: ', all_df[all_df['Genotype'].isin(['0/0', '0/1', '1/1']) == False].shape[0], flush=True)
    print('\nNumber of rows with wrong SNPs: ', all_df[(all_df['covid_snp'].isin(snp_ls) == False)].shape[0], flush=True)
    all_df = all_df[(all_df['covid_snp'].isin(snp_ls))]
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
    Loop_stats(median_df, gen_ls=['0/1'])


if __name__ == '__main__':
    main(file, unit)
    # parser = argparse.ArgumentParser(description='STEP2 - Pipeline for mQTLs'
    #     + ' - Statisques (Spearman correlation and MannWhitney test) on specific datasets')
    # parser.add_argument('filtered_file', type=str, help='basecalling file')
    # parser.add_argument('-u', '--unit', type=str, default='cpg',
    #                     help='unit to perform analysis')
    # parser.add_argument('-g', '--gen_ls', type=str, default='0/1',
    #                     help='list of genotypes to perform specific analysis on')
    # parser.add_argument('-p', '--phen_ls', type=str, default='Mild-Severe',
    #                     help='list of phenotypes to perform specific analysis on')
    # args = parser.parse_args()
    # main(**vars(args))


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
