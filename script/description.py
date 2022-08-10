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

file = 'Filtered_nano_bam_files_all_samples.csv'
file_snp = 'significant_hits_COVID19_HGI_A2_ALL_leave_23andme_20210607.txt'
dir_out='FROZEN_results_cpg_snp_analysis'
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


def description(df, title='Description.csv'):
    try:
        pd.DataFrame(df.dtypes, columns=['types']).T.to_csv(title, mode='w', header=True)
        print('types', flush=True)
        pd.DataFrame(df.nunique(), columns=['nunique']).T.to_csv(title, mode='a', header=False)
        print('nunique', flush=True)
        pd.DataFrame(df.count(), columns=['count']).T.to_csv(title, mode='a', header=False)
        print('count', flush=True)
        pd.DataFrame(df.mean(), columns=['mean']).T.to_csv(title, mode='a', header=False)
        print('mean', flush=True)
        pd.DataFrame(df.median(), columns=['median']).T.to_csv(title, mode='a', header=False)
        print('median', flush=True)
        pd.DataFrame(df.max(), columns=['max']).T.to_csv(title, mode='a', header=False)
        print('max', flush=True)
        pd.DataFrame(df.min(), columns=['min']).T.to_csv(title, mode='a', header=False)
    except:
        pass


# MAIN
snp_ls = select_SNP_per_pvalue(file_snp, pval_col='all_inv_var_meta_p',
    dist_bp=500000)
gen_ls = ['0/0', '0/1', '1/1']

# Opening file
all_df = pd.read_csv(file)
all_df = all_df[(all_df['SNP'].isin(snp_ls)) & (all_df['Gen'] != 'other')
                & (all_df['Genotype'].isin(gen_ls))]

# Dataset description
description(all_df)
