import sys
import numpy as np
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
from utils_plots import save_plot_report
import seaborn as sns
import matplotlib.pyplot as plt

def main():
    df = pd.read_csv('Filtered_nano_bam_files.csv')
    df['distance_cpg_snp'] = df['start'].astype(float) - df['pos'].astype(float)
    df['log_distance'] = np.log10(abs(df['distance_cpg_snp']))
    # g = sns.relplot(data=df, x='log_distance', y='log_lik_ratio', col='Genotype', hue='haplotype')
    # save_plot_report('Log_Distance_plot_haplotype_genotype.jpg', g,  bbox_inches='tight', dpi=100)
    for n, gen in enumerate(df['Genotype'].unique()):
        g = sns.jointplot(data=df[df['Genotype'] == gen], x='log_distance', y='log_lik_ratio', hue='haplotype')
        plt.suptitle(gen, pad=5)
        save_plot_report(f'Log_Distance_jointplot_haplotype_genotype_{n}.jpg', g)

if __name__ == '__main__':
    main()
