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
from bioinfokit import visuz
import sqlite3


def genomic_info_plot(chr, min, max, db, mark_ls=[], legend_type='color', df_measure=pd.DataFrame(), xmeasure='', ymeasure=''):

    # Recuperation of the genomic info
    info_df = pd.read_sql(f"SELECT seqid, source, type, start, end, attributes FROM hsa_ensembl_annot WHERE seqid IS {chr} AND start > {min} AND end < {max}", con=db)
    info_df = info_df[info_df['source'] != 'GRCh38'].copy()
    print(info_df.shape)

    # Initiation of the plot params
    fig, ax = plt.subplots()
    dict_color = {type:color for color, type in zip(sns.color_palette("tab10",n_colors=info_df['type'].nunique() ), info_df['type'].unique())}
    dict_height = {type:n for n, type in zip(range(info_df['type'].nunique()), info_df['type'].unique())}
    ax.set_xlabel('Genomic coordinates')
    ax2 = ax.twinx()
    ax2.set_ylabel('Genomic regions')

    # adding a measure as background (x axis corresponding to genomic coordinate)
    if not df_measure.empty:
        df_measure[xmeasure] = df_measure[xmeasure].astype(int)
        sns.histplot(df_measure, x=xmeasure, y=ymeasure, alpha=0.5, ax=ax)
    else:
        ax.get_yaxis().set_visible(False)

    # Plots the genomic regions
    for i in info_df.index:
        ax.plot([int(info_df.loc[i,'start']),int(info_df.loc[i,'end'])], [dict_height[info_df.loc[i,'type']],dict_height[info_df.loc[i,'type']]], lw=4, alpha=0.8, label=info_df.loc[i,'type'],  color=dict_color[info_df.loc[i,'type']])

    # Get a legend by color
    if legend_type == 'color':
        ax2.get_yaxis().set_visible(False)
        ax.legend()
        Lines, Labels = ax.get_legend_handles_labels()
        ax.get_legend().remove()
        unique_labels, unique_lines = [], []
        for line, label in zip(Lines, Labels):
            if not label in unique_labels:
                unique_lines.append(line)
                unique_labels.append(label)
        unique_labels.reverse()
        unique_lines.reverse()
        plt.legend(unique_lines, unique_labels, bbox_to_anchor=(1.9, 1.05))

    # # get legend on y-axis2
    # elif legend_type=='yaxis':
    #
    #     ax2.set_yticks(list(dict_height.values()))
    #     ax2.set_yticklabels(dict_height.keys())

    # Add positional marker as red vertical lines
    for mark in mark_ls:
        plt.axvline(int(mark), lw=2, color='red', alpha=0.5)
    return fig


def finemapping(df, pval_col, x='cpg', dist_bp=500000, group='CHR'):
    df = df.copy()
    df[['CHR', 'POS']] = df[x].str.split(':', expand=True)[[0, 1]].astype(int)
    best_snp_df = pd.DataFrame()
    for val in df[group].unique():
        chr_df = df[df[group] == val]
        while chr_df.shape[0] > 0:
            best_pval = chr_df[pval_col].min()
            chr_df[pval_col].min()
            best_snp_df = best_snp_df.append(chr_df[chr_df[pval_col]
                               == best_pval].drop_duplicates())
            pos = chr_df[chr_df[pval_col]
                         == best_pval]['POS'].unique().tolist()[-1]
            chr_df = chr_df[(chr_df['POS'] < pos - dist_bp)
                            | (chr_df['POS'] > pos + dist_bp)]
    return best_snp_df


def info_formatting(df):
    # DONE?: Find a way to select minimal feature to describe position of cpg
    df = df[df['type'] != 'chromosome'].copy()
    attributes_ls = set(df['attributes'].str.split(pat = '=.*?;', expand=True).dropna().values.flatten())
    att_ls = list(sorted(attributes_ls))[:4]
    # att_ls = st.multiselect('Which attributes ?', attributes_ls)
    for var in att_ls:
        df[var] = df['attributes'].str.extract(f'({var}=.*;)')[0].str.split(';', expand=True)[0].str[len(var)+1:]
    df.drop('attributes', axis=1, inplace=True)
    dict_id_to_name = df[['ID', 'Name']].dropna().set_index('ID').sort_index().to_dict()['Name']
    low_annot = lower_annotation_cpg(df)
    df = df.loc[low_annot]
    df['Parent_name'] = df['Parent'].str.split(':', expand=True)[1].replace(dict_id_to_name)
    df['Name_simplified'] = df['Name'].str.rsplit('-', n=1, expand=True)[0]
    df.loc[df['Name'].astype(str).str.contains('ENS'), 'Name_simplified'] = df.loc[df['Name'].astype(str).str.contains('ENS'), 'Parent_name'].str.rsplit('-', n=1, expand=True)[0]
    return df[['type', 'biotype', 'Name_simplified', 'ID', 'Parent']].drop_duplicates()


def lower_annotation_cpg(df):
    df = df.copy()
    # df['len'] = df['end'] - df['start']
    index_list = []
    for cpg in df['cpg'].unique():
        cpg_df = df[df['cpg'] == cpg].copy()
        if cpg_df.shape[0] != 1:
            for n, row in cpg_df.iterrows():
                if not str(row['Parent']) == 'nan':
                    index_list += cpg_df[cpg_df['ID'] == row['Parent'].split(':')[1]].index.tolist()
                elif row['type'] == 'chromosome':
                    print('chr')
                    index_list.append(n)
    return df.drop(index_list).index



def pval_plot_new(stat, xvar, pval_col, pval_cutoff=0.01, annot=True, finemap_dist=500000, snp_list=[], n_site=2, format_xaxis=True, save=False, out_dir='', title_supp=''):

    stat = stat.copy()
    stat[['CHR', 'POS']] = stat[xvar].str.split(':', expand=True)[[0, 1]].astype(int)
    stat['CHR'] = stat['CHR'].astype('category')
    stat['minus_log10'] = -np.log10(stat[pval_col].astype(float))
    pval_cutoff = -np.log10(pval_cutoff)
    stat.sort_values('CHR', inplace=True)

    # Selection best snp
    if annot == 'top':  # By top values
        best_cpg = stat[pval_cutoff < stat['minus_log10']].sort_values('minus_log10', ascending=False).groupby('CHR').head(n_site)['cpg'].tolist()
    else:   # by finemapping
        best_cpg = finemapping(stat[pval_cutoff < stat['minus_log10']], pval_col=pval_col, dist_bp=finemap_dist, group='CHR')['cpg'].unique()

    print(f'{xvar} ABOVE CUTOFF ({-np.log10(pval_cutoff/stat["cpg"].nunique())}) {title_supp}: {best_cpg}\n', stat.loc[stat[xvar].isin(best_cpg)].to_markdown())

    # If position annotation (by vertical lines)
    if snp_list:
        snp_df = pd.DataFrame(snp_list, columns=[xvar])
        snp_df[['CHR', 'POS']] = snp_df[xvar].str.split(':', expand=True)[[0, 1]].astype(int)
    else:
        snp_df = pd.DataFrame()

    # Format X-axis
    stat.sort_values(['CHR', 'POS'], inplace=True)
    if format_xaxis:
        last_pos = 0
        for chr in stat.sort_values(['CHR', 'POS'])['CHR'].astype(int).unique():
            stat.loc[stat['CHR'] == chr, 'new_xvar'] = stat.loc[stat['CHR'] == chr, 'POS'] - stat.loc[stat['CHR'] == chr, 'POS'].min() + last_pos + 1
            if not snp_df.empty:
                snp_df.loc[snp_df['CHR'] == chr, 'new_xvar'] = snp_df.loc[snp_df['CHR'] == chr, 'POS'] - stat.loc[stat['CHR'] == chr, 'POS'].min() + last_pos + 1
            last_pos = stat.loc[stat['CHR'] == chr, 'new_xvar'].max()
        print(stat[stat['new_xvar'] <= 0])
        stat['new_xvar'] = np.log10(stat['new_xvar'])

        if not snp_df.empty:
            print(snp_df[snp_df['new_xvar'] <= 0])
            snp_df['new_xvar'] = pd.Series(np.where(snp_df['new_xvar'] > 0, np.log10(snp_df['new_xvar']), 1), index=snp_df.index)
            print(snp_df[snp_df['new_xvar'] <= 0])
            # print(last_pos)
    else:
        if not snp_df.empty:
            blou = stat[['CHR', 'POS', xvar]].append(snp_df[['CHR', 'POS', xvar]]).sort_values(['CHR', 'POS']).reset_index(drop=True)
            snp_df['new_xvar'] = blou[blou[xvar].isin(snp_df[xvar])].index
            stat['new_xvar'] = blou[blou[xvar].isin(stat[xvar])].index
        else:
            stat['new_xvar'] = stat.sort_values(['CHR', 'POS']).reset_index(drop=True).index

    # PLOT
    xvar = 'new_xvar'
    stat.sort_values(['CHR', xvar], inplace=True)
    g = sns.FacetGrid(stat, aspect=2, height=4, palette='Spectral',
                      margin_titles=True)
    g.map(plt.axhline, y=pval_cutoff, color='r')
    g.map_dataframe(sns.scatterplot, xvar, 'minus_log10', hue='CHR',
        legend=True)
    if snp_list:
        for snp in snp_df.sort_values(xvar)[xvar].unique():
            plt.axvline(x=snp, alpha=0.5)
    if format_xaxis:
        g.set(xlabel="CHR", xticks=stat.groupby(['CHR']).median()[xvar].unique(),
            xticklabels=stat['CHR'].unique())
    else:
        g.set(xlabel="CHR", xticks=stat.groupby(['CHR']).last()[xvar].unique(),
        xticklabels=stat['CHR'].unique())

    if annot:
        try:
            for row in stat[stat['cpg'].isin(best_cpg)].iterrows():
                row = row[1]
                g.axes[0,0].text(row[xvar], row.minus_log10 + 0.2, row.cpg,
                    horizontalalignment='left')
        except Exception as err:
            print('couldn\'t print best points', err)
    if save:
        save_plot_report(f'minuslog10_{pval_col}_{title_supp}', g, output=out_dir, file=sys.stdout)

    return g


def stats_plots(db):
    # TODO: Add SQL calling
    df = pd.read_csv('All_files_temp.csv') # TEMP:  Example file
    # Basic plots to print out
    df_cpg_lists = pd.DataFrame()
    for _, row in df[['stat', 'variable', 'data', ]].drop_duplicates().T.iteritems():
        df_stat = df[(df['stat'] == row['stat']) & (df['variable'] == row['variable']) & (df['data'] == row['data'])]
        cutoff = 0.05 / df_stat['cpg'].nunique()

        # Save list of cpfs above the cutoff
        df_cpg_lists['_'.join(row)] = df_stat[df_stat['p-val'] < cutoff].reset_index()['cpg']
        finemapped_cpgs = finemapping(df_stat, 'p-val', group='CHR', dist_bp=10000000)[['cpg', 'p-val']]
        print(finemapped_cpgs)
        df_cpg_lists.to_csv('CpGs_significantly_associated.csv', index=False)


        # Finemapping per SNP
        # QUESTION: Is double finemapping to select the right SNPs to analyse is useful?
        # QUESTION: How else?
        list_snp = list(finemapping(finemapping(df_stat, 'p-val', group='snp_bis', dist_bp=10000000),'p-val',  x='snp_bis', group='CHR', dist_bp=10000000)['snp_bis'].unique())

        # Plots
        pval_plot_new(df_stat, 'cpg', 'p-val', snp_list=list_snp, pval_cutoff=cutoff,format_xaxis=False, out_dir='', title_supp='_'.join(row), save=True)
        pval_plot_new(df_stat, 'cpg', 'p-val', pval_cutoff=cutoff,format_xaxis=True, out_dir='', title_supp='_'.join(row), snp_list=list_snp, save=True)
        visuz.marker.mhat(df=df_stat.copy(), chr='chr',pv='p-val', show=False, gwas_sign_line=True, gwasp=0.05 / df_stat['cpg'].nunique(), markernames=False, markeridcol='cpg')

    # Other way of having the finemapped cpgs in a table
    # finemapped_cpgs_df = df.groupby(['stat', 'data', 'variable']).apply(lambda x: list(finemapping(x, 'p-val', group='CHR')['cpg'].unique())).reset_index()

    return list_snp


def distance_plots(list_snp, db):
    # TEMP:  Example file2:
    df = pd.read_csv('All_files_temp.csv')
    # TODO: Add SQL calling
    # TODO: Join the SNPs? df = join_snps(df, snp_file)
    try:
        snp_col = df.filter(regex='[Ss][Nn][Pp]').columns[0]   # TEMP
        df[f'pos_{snp_col}'] = df[snp_col].str.split(':', expand=True)[1]
        df['pos_cpg'] = df['cpg'].str.split(':', expand=True)[1]
        df['distance_cpg_snp'] = df['pos_cpg'].astype(float) - df[f'pos_{snp_col}'].astype(float)
        df['log_distance'] = np.log10(abs(df['distance_cpg_snp']))
    except:
        pass

    # TODO: Distance plot with SNPs of interest
    # TODO: Same with p-val?
    # for snp in list_snp:
    #     g = sns.relplot(data=df[df[snp_col] == snp], x='log_distance', y='log_lik_ratio', col='Genotype', hue='haplotype', row='phenotype')
    #     save_plot_report(f'Log_Distance_plot_{snp}.jpg', g,  bbox_inches='tight', dpi=100)

    # Plots the region info plot
    # TODO: Change selection of the regions?
    try:
        markers_ls = [int(str.split(':')[1]) for str in list_snp]
        chr_snp = [int(str.split(':')[0]) for str in list_snp]
        for chr in set(chr_snp):
            max_pos, min_pos = max(markers_ls), min(markers_ls)
            fig_region = genomic_info_plot(15, min_pos, max_pos, db, mark_ls=markers_ls, legend_type='color', df_measure=df[(df['start'] >= min_pos - 100000) & (df['start'] <= (max_pos + 100000))], xmeasure='start', ymeasure='minuslog10')
            save_plot_report(f'GenomicInfo_plot_{chr}', fig_region)
    except:
        pass


def join_snps(df, snp_file):
    ## TODO: find the SNPS whose genotype is most associated
    # Method to join the snp to the table
    # Or to create a table and then pick up the merge/join by sql ?
    snp_file = '/hps/nobackup/birney/projects/gel_methylation/covid_snps_current/snps/merged_significant_snp_list_genmoicc_covid_hg_hg38_20072022.txt'
    snp_file = 'merged_significant_snp_list_genmoicc_covid_hg_hg38_20072022.txt'
    snp_df = pd.read_table(snp_file, header=None, names=['snp'])
    snp_df[['chr', 'pos', 'ref', 'alt']] = snp_df['snp'].str.split(':', expand=True)

    df[['chr','start', 'n_cpg']] = df['cpg'].str.split(':', expand=True)
    join = pd.merge(df[['chr', 'start', 'cpg']], snp_df[['chr', 'pos', 'snp']], on='chr')
    join = join[(join['pos'].astype(int) - 500000 < join['start'].astype(int)) & (join['pos'].astype(int) + 500000 > join['start'].astype(int))]
    join['snp_bis'] = (join['snp'].str.rsplit(':', n=2, expand=True)[0].str[:-5] + '00000')
    join = pd.merge(df, join[['snp', 'snp_bis', 'cpg']], on='cpg')
    return join


def main(db_file):
    db_file = '/home/fanny/Work/EBI/covid_nanopore/new_covid_snp/covid_snp_March2022.db' # TEMP
    db = sqlite3.connect(db_file)
    list_snp = stats_plots(db)
    distance_plots(list_snp, db)


# TODO: Add table with SNP - CPG link --> join_snps(df, snp_file)
# TODO: ADD genomic ref table to SQL database !
# db_file = '/home/fanny/Work/EBI/covid_nanopore/new_covid_snp/covid_snp_March2022.db' # TEMP
# db = sqlite3.connect(db_file)
# # # Add Annotation file to the Database
# col_names = ['seqid', 'source', 'type', 'start', 'end', 'score',
#     'strand', 'phase', 'attributes']
# df = pd.read_csv('/home/fanny/Work/EBI/Homo_sapiens.GRCh38.85.gff3.gz', compression='gzip', sep='\t',
#                  comment='#', low_memory=False, header=None, names=col_names)
# df.to_sql('hsa_ensembl_annot', con=db)


if __name__ == '__main__':
    main()
