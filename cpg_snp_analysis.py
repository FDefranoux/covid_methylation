import sys
import os
import pandas as pd
import socket
import numpy as np
import sqlite3
host = socket.gethostname()
import subprocess
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
import multiprocessing as mp
import warnings
OUTDATED_IGNORE=1
import csv
import time

# TODO: Modify the stat function to have the double ttests...DOUBLE ???


def stat_test(df, cpg, snp, measure, test_args, target_snp='covid_snp', n_min=5, **kwargs):
    start = time.time()
    # print(cpg, snp, df.shape, flush=True)
    dict_rename = {'corr': [pg.corr, {'r':'coeff'}, dict(method='spearman')],
                   'mwu': [pg.mwu, {'U-val': 'coeff'}],
                   'ttest': [pg.ttest, {'T':'coeff'}]}

    # Getting the stat test parameters
    df = df[(df['cpg'] == cpg) & (df[target_snp] == snp)].copy()
    stat_test, var = test_args.pop('test'), test_args.pop('var')
    mask = (df[var] == df[var])
    for key, val in test_args.items():
        mask &= (df[key] == val)
    df = df[mask].copy()
    assert stat_test in dict_rename.keys(), 'Statistical test not recognized'


    # Setting the variables
    if (df[var].nunique() < 2) or (df.shape[0] < n_min) :
        pass
    else:
        if stat_test in ['mwu', 'ttest']:
            vars = list(df[var].unique())
            X = df.loc[df[var] == vars[0], measure].astype(float)
            Y = df.loc[df[var] == vars[1], measure].astype(float)
            if len(vars) > 2:
                warnings.warn(f'WARNING: More than two values for variable {var} ({vars})!')
        elif stat_test == 'corr':
            X = df[var].astype(int)
            Y = df[measure].astype(float)

        # Computing results
        if (X.std() != 0) & (Y.std() != 0):
            try:
                if len(dict_rename[stat_test]) > 2:
                    res = dict_rename[stat_test][0](X, Y, **dict_rename[stat_test][2]).rename(columns=dict_rename[stat_test][1])
                else:
                    res = dict_rename[stat_test][0](X, Y).rename(columns=dict_rename[stat_test][1])
                n_samples = df[var].value_counts().to_dict()
                # print('Statistical function time', time.time() - start, '\n', flush=True)
                return [cpg, snp, stat_test, var, n_samples, '-'.join(test_args.values()),
                    res['coeff'][0], res['p-val'][0]]

            except Exception as err:
                print(cpg, snp, err)
                # print(df)

# MAIN
def main(file, target_snp='covid_snp', select=None, measure='log_lik_ratio', output_dir='.', n_min=5, pval=0.05, zscore=3, multiprocess=False):
    start_time = time.time()

    # List of tests to do:
    #   - test-mwu_var-phenotype
    #   - test-mwu_var-phenotype_Genotype-0/1
    #   - test-mwu_var-haplotype
    #   - test-mwu_var-haplotype_Genotype-0/1
    #   - test-ttest_var-haplotype_Genotype-0/1_phenotype-Severe
    #   - test-ttest_var-haplotype_Genotype-0/1_phenotype-Mild
    #   - test-corr_var-genotypefact
    #   - test-corr_var-genotypefact_phenotype-Mild
    #   - test-corr_var-genotypefact_phenotype-Severe

    test_ls = [
        dict(test='mwu', var='phenotype'),
        dict(test='mwu', var='phenotype', Genotype='0/1'),
        dict(test='mwu', var='haplotype'),
        dict(test='mwu', var='haplotype', Genotype='0/1'),
        dict(test='ttest', var='haplotype', Genotype='0/1', phenotype='Severe'),
        dict(test='ttest', var='haplotype', Genotype='0/1', phenotype='Mild'),
        dict(test='corr', var='genotypefact'),
        dict(test='corr', var='genotypefact', phenotype='Severe'),
        dict(test='corr', var='genotypefact', phenotype='Mild'),
    ]

    # Create output path
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Opening file
    if select:
        # Getting dictionnary or arguments
        select = select.replace('=', '-')
        select_dict = {val.split('-')[0]:val.split('-')[1] for val in select.split('_')}
        test_type = {key:select_dict.pop(key) for key in ['test', 'var'] if key in select_dict.keys()}
        sql_query = "SELECT * FROM median_datas"
        if select_dict:
            sql_query += " WHERE "
            value_ls = [str(key) + " IS '" + str(val) + "'" if not '%' in val else str(key) + " LIKE \"" + str(val) + "\"" for key, val in select_dict.items()]
            sql_query += " AND ".join(value_ls)
            sql_query = sql_query.replace('snp', target_snp)
            sql_query += ";"
            print(f'selection: {sql_query}')

    # Loading the DF
    if file[-3:] == '.db':
        print('DB MODE')
        db = sqlite3.connect(file)
        all_df = pd.read_sql(sql_query, con=db)
        db.close()
    elif file[-4:] == '.csv':
        print('CSV MODE')
        command = f"sqlite3 :memory: -cmd '.mode csv' -cmd '.import {file} median_datas' \ '{sql_query}'"
        print(command)
        res = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        res_out = res.stdout.read().decode('utf8')
        res_out = res_out.split('\n')
        res_out = [row.split(',') for row in res_out if row != '']
        all_df = pd.DataFrame(res_out)
        print(all_df.shape)

        # Adding columns
        if not all_df.empty:
            header = subprocess.Popen(f"head -n1 {file}", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            header_ls = header.stdout.read().decode('utf8')[:-1].split(',')
            print(all_df.head(2).to_markdown())
            all_df.columns = header_ls
        else:
            print('DF EMPTY')
    else:
        all_df = pd.DataFrame()

    print(all_df.head(2))

    if not all_df.empty:
        # OUTLIERS DETECTION
        # QUESTION: Where should we do the outliers detection ?
        stat_raw = Stats(all_df[[measure]])
        stat_raw.no_outliers(zscore_thresh=zscore)
        print('\nOutliers:', round(stat_raw.outliers.shape[0]/all_df[[measure]].shape[0], 3), '% -- n_lines = ', stat_raw.outliers.shape[0], flush=True)
        all_df = all_df.loc[stat_raw.no_outliers.index]
        all_df['genotypefact'] = all_df['Genotype'].replace({'0/0':0, '0/1':1, '1/1': 2})
        cpg_snp = all_df[['cpg', target_snp]].drop_duplicates()
        print('DF shapes', all_df.shape, cpg_snp.shape)
        read_time = time.time()
        print ('Read the DF and remove the outliers time', read_time - start_time, flush=True)

        if test_type:
            # Initiation of multiprocess
            if multiprocess:
                pool = mp.Pool(mp.cpu_count())
                print ('CPU count', mp.cpu_count())
                result_ls = [pool.apply(stat_test, args=(all_df[(all_df['cpg'] == cpg) & (all_df[target_snp] == snp)].drop_duplicates(), cpg, snp, measure, test_type.copy()), kwds=dict(target_snp=target_snp, n_min=n_min, method='spearman')) for cpg, snp in zip(cpg_snp['cpg'], cpg_snp[target_snp])]
                pool.close()
            else:
                result_ls = []
                for _, (cpg, snp) in cpg_snp.iterrows():
                    result_ls.append(stat_test(all_df[(all_df['cpg'] == cpg) & (all_df[target_snp] == snp)].drop_duplicates(), cpg, snp, measure, test_type.copy(), target_snp=target_snp, n_min=n_min, method='spearman'))
        else:
            print('No test defined !')
            result_ls = []
            for test_type in test_ls:
                for _, (cpg, snp) in cpg_snp.iterrows():
                    result_ls.append(stat_test(all_df[(all_df['cpg'] == cpg) & (all_df[target_snp] == snp)].drop_duplicates(), cpg, snp, measure, test_type.copy(), target_snp=target_snp, n_min=n_min, method='spearman'))
            pass

        # Save to csv
        end_process = time.time()
        print(f'Process time (MP = {multiprocess})', end_process - read_time )
        fields = ['cpg', 'snp', 'stat_test', 'var', 'n_samples', 'data',
        'coeff', 'pval']
        result_ls = [ res for res in result_ls if res != None]
        # select = '_'.join([str(select) + 'select' for select in [test_select, chr_select, snp_select] if select != None])
        select = select.replace('/', '')
        with open(os.path.join(output_dir, f'Results_stats_select_{select}.csv'), 'w') as f:
            write = csv.writer(f)
            write.writerow(fields)
            write.writerows(result_ls)
        print('Saving time', time.time() - end_process)
        print('Complete time', time.time() - start_time)

    else:
        print(f'Dataframe is empty ! {file} -- {select}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='STEP2 - Pipeline for mQTLs'
        + ' - Statisques (Spearman correlation and MannWhitney test) on specific datasets')
    parser.add_argument('file', type=str, help='filtered db')
    parser.add_argument('-t', '--target_snp', type=str, default='covid_snp',
        help='target snp on which to perform analysis')
    parser.add_argument('-s', '--select', type=str, default=None,
                        help='selection to do on data/test')
    parser.add_argument('-m', '--measure', type=str, default='log_lik_ratio',
                        help='name of the measurement')
    parser.add_argument('-o', '--output_dir', type=str, default='.',
                        help='directory to save output files')
    parser.add_argument('-n', '--n_min', type=int, default=5,
                        help='number of minimum replicates')
    parser.add_argument('-p', '--pval', type=float, default=0.05,
                        help='pval threshold')
    parser.add_argument('-z', '--zscore', type=int, default=3,
        help='zscore value for outlier detection')
    parser.add_argument('-f', '--multiprocess', action='store_true', default=False,
        help='Run with multiprocess')
    args = parser.parse_args()
    print('Arguments (cpg_snp_analysis.py): ', ' - '.join([str(x)+': '+str(y) for x,y in vars(args).items()]), '\n', flush=True)
    main(**vars(args))

# Old function
# Stats_MultiTests(all_df, 'log_lik_ratio', 'covid_snp', cpg_ls=[], output=os.path.join(output_dir, f'Multi_Stats_{chr}'), n_min=n_min, pval=pval)
# Loop_stats(all_df, target_snp, gen_ls=['0/1'], phen_ls=['Mild', 'Severe'], output=f'{os.path.join(output_dir, os.path.basename(file)[:-4])}_')

## To test functions
# file_db = '/home/fanny/Work/EBI/covid_nanopore/covidVScontrol_March2022/covid_snp_March2022/Filtered_nano_bam_files.csv'
# target_snp='covid_snp'
# chr = 12
# all_df = pd.read_csv(file_db)
# all_df.groupby('chromosome').size()
# all_df=all_df[all_df['chromosome']== chr]
# # print(all_df.head(2).to_markdown())
# # all_df.drop('Unnamed: 0', axis=1, inplace=True)
#
# all_df = all_df.groupby(['phenotype', 'sample_id', 'chromosome', 'cpg', target_snp, 'Genotype', 'haplotype']).median().reset_index()


### OLD FUNCTIONS
def Stats_MultiTests(df, measure, target_snp, cpg_ls=[], output='Multi_Stats', n_min=5, pval=0.05):
    df = df.copy()
    print(df.head(2).to_markdown())
    df = df.sort_values(by=['Genotype', 'phenotype'])
    df = df.sort_values(by=['haplotype'], ascending=False)
    for col in ['Genotype', 'phenotype', 'haplotype']:
        df[f'{col}_test'] = pd.factorize(df[col])[0]
    if not cpg_ls:
        cpg_ls = df['cpg'].unique()
    print('Running through the CpGs\n')
    for cpg in tqdm(cpg_ls):
        results = pd.DataFrame()
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
            del samp_ph
            if not samp_hetph.empty:
                try:
                    samp_hetph = samp_hetph.set_index(['cpg', 'sample_id', 'phenotype', 'Genotype', 'haplotype'])[measure].unstack().dropna(thresh=2).reset_index()
                    res = stat_test(samp_hetph, 'alt', 'ref', test='ttest', data=f'het-{val}', n_min=n_min)
                    results = results.append(res)
                except:
                    print(f'error with cpg {cpg} for het samples stats')
                    print(samp_hetph.head(3).to_markdown())
            del samp_hetph

        # Formatting forsaving
        results.rename(columns={'index':'stat', 'x':'variable', 'y':'measure'}, inplace=True)
        results.loc[results['variable'] == 'alt', ['variable', 'measure']] = ['alt-ref', measure]
        if cpg == cpg_ls[0]:
            results.to_csv(output + '.csv', index=False, header=True, mode='w')
        else:
            results.to_csv(output + '.csv', index=False, header=False, mode='a')

    # # Printout of significative results
    # print(f'Number of significative cpgs at pvalue threshold {pval} (non-corrected)\n')
    # print(results[results['p-val'] < pval].astype(str).groupby(['variable', 'stat', 'data']).size().unstack().sort_index(axis=1, key=lambda x: x.str.lower()).reset_index(level='stat').reindex(['phenotype_test', 'Genotype_test', 'haplotype_test', 'alt-ref']).reset_index().fillna('').to_markdown())


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
