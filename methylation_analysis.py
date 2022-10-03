import sys
import os
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
from utils import *
from utils_plots import *
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import streamlit as st
import sqlite3
import plotly.express as px

################################################################################
############################ GENERAL TODOs #####################################
################################################################################

# TODO: Change the previous scriptS #with_updated_snplist
#   - Decide which table are weel suited
#   - Report of the output for each steps
#   - # QUESTION: What should we do with controls ?
#   - Implement read_count function anad save table in sql
#   - Add read_number in the median table?
#   - add counts for the individual boxplot (table next to the plot or n= over each box)
# QUESTION: Do we discard completely the cpg containing less than 5reads ?
# TODO: Remove the big plots from this script

# TODO: Change the access of datas to the median table instead ! #ALL
# TODO: Verify that all sql qery pass by the function (how to join table ?) #ALL
# TODO: Change the way the counts are selected (instead of Mild+severe > count_min) #ALL
# TODO: Add function for saving each plot ?  #ALL
# MAYBE: Possibility to run the main plots without streamlit ? #ALL
# TODO: Other info: Info on the nearest gene? pathway GO annotation? #distanceplot #ing_analysis #info
# TODO: Distance plots = zoomplots (around one SNP) #distanceplot
# TODO: Anotations plots around one SNP #distanceplot #info
# TODO: Info on the nearest gene? #info

# QUESTION: Could it be interesting to have a look at the ref alt meth differences in homozygotes ? #ALL
# QUESTION: Do we need to remove the called 'alt' in '0/0' and 'ref' in '1/1' ?? #ALL
# QUESTION: When do we discard the cpgs with not enough count (reads or patients) #ALL


###################### CPG SELECTION ###########################################
def finemapping(df, pval_col, dist_bp=500000):
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


def interactive_param(db, main_tests_dict={}):
    dict_sql = {}
    dict_sql['p'] = st.sidebar.number_input('P-value', value=0.01000, min_value=0.0000001, max_value=1.000, step=0.00001, format="%.8f")
    if main_tests_dict:
        main_test = st.sidebar.radio('Test selection', main_tests_dict.keys())
        dict_sql = main_tests_dict[main_test]
    else:
        dict_sql['table'] = st.sidebar.radio('Statistical test', ['mann whitney', 'spearman correlation']).replace(' ', '_')
        datacol, subtest = st.sidebar.columns(2)
        datas_ls = pd.read_sql(f"SELECT DISTINCT data FROM {dict_sql['table']}", con=db)['data'].tolist()
        if 'ALL' in datas_ls: default='ALL'
        else: default='0/1'
        dict_sql['data'] = datacol.radio('Datas', sorted(datas_ls), index=sorted(datas_ls).index(default))
        tests_ls = pd.read_sql(f"SELECT DISTINCT test FROM {dict_sql['table']} WHERE data IS \'{dict_sql['data']}\'", con=db)['test'].tolist()
        dict_sql['test'] = subtest.radio('Tests', tests_ls)
    bonf_corr = st.sidebar.checkbox('Bonferoni correction')
    dist_bp = st.sidebar.number_input('Finemapping: (minimal distance in kb)', 0, 1000000, step=5000)
    return dict_sql, bonf_corr, dist_bp


def dict_to_sql_query(db, dict_sql, count=0, count_cols=[], list_cpgs=[]):
    dict_sql= dict_sql.copy()
    # Creation of sql query from dict_sql
    table = dict_sql.pop('table')
    if dict_sql.get('cols'):
        cols = ','.join([f'{table}.{col}' for col in dict_sql.pop('cols')])
    else:
        cols= '*'
    sql_query = f"SELECT {cols} FROM {table}"

    if count != 0:
        sql_query += f" JOIN counts_diff_means ON counts_diff_means.cpg={table}.cpg"
        if not count_cols:
            count_cols = ['Mild', 'Severe', 'read_count']
        if 'read_count' in count_cols:
            sql_query += f" JOIN reads_samples ON reads_samples.sample_id || '_' || reads_samples.cpg={table}.sample_id || '_' || {table}.cpg"
        sql_query += " WHERE " + ' AND '.join([f"{col} > {count}" for col in count_cols])
    if (dict_sql != {}) & (count == 0):
        sql_query += " WHERE "
    elif dict_sql != {}:
        sql_query += " AND "
    try:
        p = dict_sql.pop('p')
    except KeyError:
        p = None
    # Recup other parameters
    param = " AND ".join([key + f" IN ({str(val)[1:-1]})" if isinstance(val, list) else key + f" IS '{val}'" for key, val in dict_sql.items()])
    sql_query += f" {param}"
    if p != None: sql_query += f" AND pval < {p}"
    if list_cpgs:
        sql_query += f" AND cpg IN ({str(list_cpgs)[1:-1]})"
    # st.write(sql_query)
    return sql_query


def cpg_selection_pval(db, interactive=True, dict_sql={}, bonf_corr=True, dist_bp=0, count=0, count_cols=[], list_cpgs=[]):
    if interactive:
        dict_sql, bonf_corr, dist_bp = interactive_param(db, main_tests_dict=dict_sql)

    if (bonf_corr) & (dict_sql.get('p') != None):
        tot_query = dict_to_sql_query(db, dict_sql, count=count, count_cols=count_cols, list_cpgs=list_cpgs)
        tot = pd.read_sql(tot_query, con=db).loc[0, 'tot']
        # tot = pd.read_sql(f"SELECT COUNT (DISTINCT cpg) AS tot FROM {dict_sql['table']}", con=db).loc[0, 'tot']
        dict_sql['p'] = dict_sql['p'] / tot

    sql_query = dict_to_sql_query(db, dict_sql, count=count, count_cols=count_cols, list_cpgs=list_cpgs)
    df = pd.read_sql(sql_query, con=db)


    if (dist_bp != 0) & (not df.empty):
        pval_cpg_ls = finemapping(df, 'pval', dist_bp=dist_bp * 1000)
    else:
        # st.dataframe(pd.DataFrame({"Min": df.set_index('cpg')['pval'].idxmin(),
        #                            "Max": df.set_index(['cpg'])['pval'].idxmax()}))
        pval_cpg_ls = list(df['cpg'].unique())

    # st.dataframe(df[df['cpg'].isin(pval_cpg_ls)].drop_duplicates().style.format({"pval": "{:.2E}",}).hide_index(), 700, 300)
    # st.write(df['cpg'].nunique())
    return pval_cpg_ls


# CHECK: Mild Vs Severe for diff in %meth alt VS severe (in alt ???) #mildsev
# TODO: Change to to plot the 3 possibily interesting plots (stop giving possibilities for too much plots) #mildsev
def mild_severe_selection(db, count=0, count_cols=[]):
    dict_test = {'Spearman': {'cols': ['cpg', 'data'],
                              'table': 'spearman_correlation',
                              'test': ['Genotype', 'haplotype'],
                              'data':['Mild', 'Severe'],
                              'vars': ['pval', 'rho'],
                              'count_cols' : ['Mild', 'Severe']},
                 'Mann Whitney': {'cols': ['cpg', 'data'],
                                  'table': 'mann_whitney',
                                  'test': ['alt-ref', '0/0-1/1'],
                                  'data':['Mild', 'Severe'],
                                  'vars' : ['pval'],
                                  'count_cols' : ['Mild', 'Severe']},
                 'TTest': {'cols': ['cpg', 'data'],
                            'table': 'ttest',
                            'test': ['ref-alt'],
                            'data': ['Mild', 'Severe'],
                            'vars': ['pval', 'T'],
                            'count_cols' : ['Mild', 'Severe', 'alt', 'ref']},
                 'Percent log-likelyhood': {'cols': ['cpg', 'phenotype', 'Genotype',
                                                     'haplotype', 'sample_id'],
                                            'table': 'median',
                                            'vars': ['log_lik_ratio', '\"%meth\"']}}
    cols = st.columns(4)
    test = cols[0].radio('Which test?', dict_test.keys())
    dict_sql = dict_test[test]
    if dict_sql.get('test'):
        dict_sql['test'] = cols[1].radio('Test type: ', dict_sql.get('test'))
    vars_ls = dict_sql.pop('vars')
    measure = cols[2].radio('Measure', vars_ls)
    vars_ls = set([None] + dict_sql['cols'] + vars_ls) - set(['data', 'cpg'])
    facet = cols[3].radio('Facet', vars_ls)
    vars_ls.remove(None)
    dict_sql['cols'] = set(dict_sql['cols'] + list(vars_ls))
    count_cols = dict_sql.pop('count_cols')
    sql_query = dict_to_sql_query(db, dict_sql, count=count, count_cols=count_cols, list_cpgs=[])
    st.write(sql_query)
    df_mild_severe = pd.read_sql(sql_query, con=db, index_col=None).rename(columns={'phenotype':'data'})
    if not df_mild_severe.empty:
        val_ls = [val for val in [measure, facet] if val != None]
        # list_cols = df_mild_severe.drop(['data'] + val_ls, axis=1).columns
        pivot = df_mild_severe.pivot_table(columns=['data'], index=['cpg'], values=val_ls)
        pivot.columns = [ a +'_'+ b for a,b in pivot.columns]
        pivot.reset_index(inplace=True)
        pivot['chromosome'] = pivot['cpg'].str.split(':', expand=True)[0].astype(int)
        pivot.sort_values(['chromosome', 'cpg'], inplace=True)
        pivot['chromosome'] = pivot['chromosome'].astype('category')
        st.dataframe(pivot.dropna())
        # plotly
        x, y = pivot.filter(regex=measure).columns.tolist()
        mildsevere_kws = {'color':'chromosome',
                          'labels':{'x': x, 'y': y,
                                    # 'facet_row': facet_x, 'facet_col': facet_y
                                    },
                          'hover_name':'cpg', #'trendline': 'ols',
                          'opacity':0.6,
                          #'trendline_color_override':'darkblue',
                          # 'hover_data': {'pval_Mild':True, 'pval_Severe':True,
                                         # 'rho_Mild':True, 'rho_Severe':True,
                           # 'pval_Mild_cut':False, 'pval_Severe_cut':False,
                           # 'chromosome': False}
                          }

        if facet:
            facet_x, facet_y = pivot.filter(regex=facet).columns.tolist()
            pivot[facet_x + '_cut'] = pd.qcut(pivot[facet_x].astype(float), q=2)
            #     ).astype(str).replace({f'({bin}, 1.0]': 'high',
            #         f'(0.0, {bin}]': 'low', 'nan': None})
            pivot[facet_y + '_cut'] = pd.qcut(pivot[facet_y].astype(float), q=2)
            #     ).astype(str).replace({f'({bin}, 1.0]': 'high',
            #         f'(0.0, {bin}]': 'low', 'nan': None})
            mildsevere_kws.update({'facet_row':facet_x +'_cut', 'facet_col':facet_y+'_cut'})
        pivot.dropna(inplace=True)

        df = dynamic_pointselect(pivot, x=x, y=y, **mildsevere_kws)
        st.dataframe(df.set_index('cpg'))
        # pivot.sort_values(pivot.filter(regex='cut').columns.tolist())
        if not df.empty:
            st.dataframe(df.set_index('cpg'))
            return list(df['cpg'].unique())
        else:
            st.error('No CpG selected !')
            return []


######################## INDIVIDUAL CPG INFO AND PLOTS #########################
def individual_cpg_analysis(db, count=0, count_cols=[]):
    st.subheader("Let's discover CpGs of interest")
    selection_dict = {'Significatives cpgs selection':  cpg_selection_pval,
                      'Mild VS Severe response': mild_severe_selection}

    choice = st.sidebar.selectbox('How would you like to select CpGs?', selection_dict)
    list_cpgs = selection_dict[choice](db, count=count, count_cols=['Mild', 'Severe'])
    new_list = st.multiselect('CpGs to plot', sorted(list_cpgs))
    if new_list:
        list_cpgs = new_list
    if len(list_cpgs) == 0:
        st.error('No data selected')
    elif len(list_cpgs) > 20:
        st.warning(f'You have selected {len(list_cpgs)} CpGs. This might take a while... Continue?')
        cont_button = st.checkbox('Continue')
    else:
        st.info(f'You have pre-selected {len(list_cpgs)} CpGs')
    if (len(list_cpgs) <= 20) or cont_button:
        col1, col2, col3 = st.columns(3)
        space_info = st.empty()
        boxplot=col1.checkbox('Boxplots')
        genom_feat=col2.checkbox('Genomic Features')

        # Plotting the cpgs
        if boxplot:
            path = st.text_input('Enter images path: ', value=os.path.join(ABS_PATH, 'plots'))
            if not os.path.exists(path):
                os.makedirs(path)
            # TODO: Change sql query!
            df = pd.read_sql(f'SELECT * FROM log_methylation WHERE cpg IN ({str(list_cpgs)[1:-1]})', con=db)
            if not df.empty:
                df = df.groupby(['phenotype', 'sample_id', 'chromosome', 'cpg',
                    'covid_snp', 'Genotype', 'haplotype']).median().reset_index()


        # Info about cpgs OLD WAY # TODO: DELETE #oldinfofunction
        # if col3.checkbox('Retrieve informations'):
            # df_info = get_info_cpg(list_cpgs)
            # space_info.dataframe(df_info)
            # space_info.dataframe(df_info.groupby(['Gene', 'Description']).apply(lambda x: str(x.index.tolist())).reset_index())

        df_info2 = pd.DataFrame()
        if list_cpgs:
            for cpg in list_cpgs:
                if boxplot:
                    st.subheader(cpg)
                    boxplot_st(df[df['cpg'] == cpg], path=path)

                if genom_feat:
                    # TODO: Put in function callable also for distance analysis #info
                    chr, pos, _ = cpg.split(':')
                    blou = pd.read_sql(f"SELECT seqid, source, type, start, end, attributes FROM hsa_ensembl_annot WHERE start < {pos} AND end > {pos} AND seqid IS {chr}",
                        con=db)
                    blou['cpg'] = cpg
                    df_info2 = df_info2.append(blou)
        if not df_info2.empty:
            info = info_formatting(df_info2)
            space_info.dataframe(info)
            # df_info2['len'] = df_info2['end'] - df_info2['start']
            # df_info2.to_csv('test_cpg_info.csv', index=False)
            # st.success('saving successful')


def setup_customizedboxplot_cpg_analysis(cpg_df, unit='control_snp', dir_out='.', dict_colors=None):
    # General modification of the datas
    cpg_df = cpg_df.copy()
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
    cpg_df.sort_values(['phenotype', 'Genotype', 'haplotype'], inplace=True)
    # Creation of the plot
    fig, ax = plt.subplots(2, 3, figsize=(17,10))
    boxplot_customized(cpg_df, 'phenotype', 'log_lik_ratio', hue_var='phenotype',
                                 dict_colors=repl_colors['phenotype'], width_var=None,
                                 hatch_var='phenotype', ax=ax[0,0], hatch_dict=dict_hatch)
    ax[0,0].set(title='Symptom severity')
    boxplot_customized(cpg_df, 'Genotype', 'log_lik_ratio', hue_var='Genotype',
                                 dict_colors=repl_colors['Genotype'],
                                 width_var=None, ax=ax[0,1], hatch_dict=dict_hatch, dodge=False)
    ax[0,1].set(title='Genotype correlation')#, xticklabels=replace_val['Genotype'].values())
    boxplot_customized(cpg_df, 'phenotype', 'log_lik_ratio', hue_var='Genotype',
                                 dict_colors=repl_colors['Genotype'], width_var=None,
                                 hatch_var='phenotype', ax=ax[0,2], hatch_dict=dict_hatch)
    ax[0,2].set(title='Genotype correlation X Symptom severity')
    if not cpg_df[cpg_df['Genotype'] == replace_dict['Genotype']['0/1']].empty:
        boxplot_customized(cpg_df[cpg_df['Genotype'] == replace_dict['Genotype']['0/1']], 'phenotype',
                                       'log_lik_ratio', hue_var='phenotype',
                                       dict_colors=repl_colors['phenotype'],
                                       hatch_var='phenotype', width_var=None,
                                       ax=ax[1,0], hatch_dict=dict_hatch)
        ax[1,0].set(title='Heterozygous Symptom Severity')
        boxplot_customized(cpg_df[cpg_df['Genotype'] == replace_dict['Genotype']['0/1']], 'haplotype', 'log_lik_ratio',
                                      hue_var='haplotype', dict_colors=repl_colors['haplotype'],
                                      hatch_var=None, width_var=None, ax=ax[1,1],
                                      hatch_dict=dict_hatch, dodge=False)
        ax[1,1].set(title='Heterozygous Allele difference') #, xticklabels=replace_val['haplotype'].values())
        boxplot_customized(cpg_df[cpg_df['Genotype'] == replace_dict['Genotype']['0/1']], 'phenotype', 'log_lik_ratio',
                                      hue_var='haplotype', dict_colors=repl_colors['haplotype'],
                                      hatch_var='phenotype', width_var=None, ax=ax[1,2],
                                      hatch_dict=dict_hatch)
        ax[1,2].set(title='Het. Allele difference X Symptom severity')
    else:
        fig.delaxes(ax[1, 0])
        fig.delaxes(ax[1, 1])
        fig.delaxes(ax[1, 2])
    plt.suptitle(f'CpG {cpg} associated with SNP {snp}', weight='bold')

    # Creation of one common legend for all
    lines, labels = [], []
    x = 0.6
    for a in fig.axes:
        Line, Label = a.get_legend_handles_labels()
        a.get_legend().remove()
        if (set(Label) & set(labels)) == set():
            lines.extend(Line)
            labels.extend(Label)
            label_title = cpg_df[cpg_df == Label[0]].dropna(how='all', axis=1).columns.tolist()[0]
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
    return fig


def boxplot_st(cpg_df, path=os.path.join(ABS_PATH, 'plots')):
    cpg = cpg_df['cpg'].unique()[0]
    title_boxplot = os.path.join(path, f'Multiplots_{cpg}')
    # st.write(title_boxplot)
    if os.path.exists(title_boxplot + '.jpg'):
        st.write(title_boxplot + '.jpg exist already. Retrieved saved plot.')
        st.image(title_boxplot + '.jpg')
    if (os.path.exists(title_boxplot + '.jpg') == False) or (st.button('Rerun')):
        # st.dataframe(df[df['cpg']== cpg])
        cpg_boxplot = setup_customizedboxplot_cpg_analysis(cpg_df, unit='covid_snp', dir_out='.', dict_colors=None)
        st.pyplot(cpg_boxplot)
        # cpg_boxplot = sns.catplot(kind='box', data=df[(df['cpg']== cpg) & (df['Genotype'] == '0/1') & (df['haplotype'] != 'other')], x='phenotype', y='log_lik_ratio',
        #                               hue='haplotype')
        # st.pyplot(cpg_boxplot)

        blou = cpg_df[(cpg_df['Genotype'] == '0/1')].groupby(['cpg', 'sample_id', 'phenotype', 'Genotype', 'haplotype']).median()['log_lik_ratio'].unstack().reset_index()
        blou['diff_alt_ref'] = blou['ref'] - blou['alt']
        blou.sort_values(['phenotype', 'Genotype'], inplace=True)
        dict_colors = {'Genotype': {'0/0':'#1678F5', '0/1':'#2aa69a', '1/1':'#3ED43E'},
                       'haplotype': {'ref':'#1678F5', 'alt':'#3ED43E'},
                       'phenotype': {'Mild':'#f0f0f5', 'Severe':'#c2c2d6'}}
        new_box = sns.catplot(kind='box', data=blou[blou['Genotype'] == '0/1'],
            x='phenotype', hue='phenotype', y='diff_alt_ref', palette=dict_colors['phenotype'])

        col1, col2, col3  = st.columns(3)
        col2.pyplot(new_box)
        if st.button(f'Save {cpg}'):
            # path = st.text_input('Path for saving', value=os.path.abspath(title_boxplot))
            if not os.path.exists(os.path.dirname(title_boxplot)):
                st.error(f'Path {os.path.abspath(title_boxplot)} not found')
            else:
                save_plot_report(os.path.basename(title_boxplot), cpg_boxplot, output=os.path.dirname(title_boxplot))
                st.success(f'File {os.path.abspath(title_boxplot)} saved successfully')


# TODO: Delete function ? Or replace by new way ? #info #oldinfofunction
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
        if gene_ls:
            results[cpg[0]] = {'genes': gene_ls, 'symbol':[], 'pathway':[], 'description':[], 'other':[]}
            for gene in gene_ls:
                a, b, c = '', '', ''
                if gene:
                    try:
                        a = ensembl_rest.symbol_lookup(species='homo sapiens', symbol=gene)
                        if a:
                            results[cpg[0]]['symbol'].append(a.get('display_name'))
                            results[cpg[0]]['description'].append(re.sub("\[.*\]"," ", a.get('description')))
                    except:
                        pass
                    try:
                        b = k.find("hsa", gene).split('\t')[0]
                        if (b != '') & (b != '\n'):
                            c = k.get(b)
                            results[cpg[0]]['pathway'].append(b)
                        if (c != '') & isinstance(c, str):
                                c = c[c.find('ORTHOLOGY')+10: c.find('ORGANISM')-1]
                                results[cpg[0]]['other'].append(c)
                    except:
                        pass

    # df_info = pd.DataFrame()
    # st.write(pd.DataFrame(results))
    # for key, dict_val in results.items():
    #     for n, gene in enumerate(dict_val.keys()):
    #         df_info.loc[key, f'Gene_{n+1}'] = str(dict_val[gene])
    # if df_info.shape[1] == 1:
    #     return df_info['Gene_1'].str[1:-1].str.split(', ', n=3, expand=True).rename(columns={0:'Gene', 1:'Description', 2:'Other'}
    #         ).replace('\'', '')
    # else:
    return pd.DataFrame(results).replace('[]', '').dropna(how='all')


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
    return df[['cpg', 'type', 'biotype', 'Name_simplified', 'ID', 'Parent']].drop_duplicates()


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



############ STATS PLOTS  & ASSOCIATION TESTING ANALYSIS #######################
def stat_plots(db, count=0, count_cols=[]):
    # TODO: spearman VS mann whitney (plot p-val MWU phenotype vs p-val Spearman genotype) #statplots
    main_tests_dict = {'Genotype correlation': {'table': 'spearman_correlation',
                            'test': 'Genotype', 'data':'ALL'},
                       'Haplotype difference': {'table': 'mann_whitney',
                            'test': 'alt-ref', 'data':'0/1'},
                       'Symptom severity':  {'table': 'mann_whitney',
                            'test': 'Mild-Severe', 'data':'ALL'}}
    res_df = pd.DataFrame(columns=main_tests_dict.keys(), index=['covid_snp', 'control_snp'])
    dict_cpgs = {}

    test_ls = st.multiselect('Mahnattan plots', main_tests_dict.keys())
    cols = st.columns(2)
    bonf_corr = cols[0].checkbox('Bonferroni Correction')
    norm = cols[1].checkbox('Normalisation')
    plots = st.button('Plots')

    # For table
    st.subheader('Final counts')
    pval = st.radio('Choose p-val', [0.1, 0.05, 0.01, 0.001])
    for test in test_ls:
        total_cpgs = cpg_selection_pval(db, interactive=False, dict_sql=main_tests_dict.get(test), bonf_corr=bonf_corr, dist_bp=0)
        # TODO: Implement the pval plots with possibility to zoom in on some SNPs #statplots
        # df = pd.read_sql("", con=db)
        # pval_plot = pval_plot_new(df, xvar, pval_col, pval_cutoff=0.01, n_site=2, format_xaxis=True, out_dir='', title_supp='')
        # st.pyplot(pval_plot)
        new_sql = main_tests_dict.get(test).copy()
        new_sql['p'] = pval
        dict_cpgs['covid_snp' + test] = cpg_selection_pval(db, interactive=False, dict_sql=new_sql, bonf_corr=bonf_corr, dist_bp=0)
        tot = len(total_cpgs)
        res_df.loc['covid_tot', test] = tot
        if (plots == False) & (norm != False):
            tot = len(cpg_selection_pval(db, interactive=False, dict_sql=main_tests_dict.get(test), bonf_corr=bonf_corr, dist_bp=0))
        else:
            tot=1
        res_df.loc['covid_snp', test] = len(dict_cpgs['covid_snp' + test])/tot
        st.write(len(set(dict_cpgs['covid_snp' + 'Genotype correlation']) & set(dict_cpgs['covid_snp' + 'Haplotype difference'])))
        # res_df = res_df.reset_index().melt(id_vars=['index'])
        st.write(res_df)
        plotcol = st.columns(3)
        st.pyplot(sns.catplot(kind='bar', data=res_df.reset_index().melt(id_vars=['index']), x='index', y='value', hue='variable'))
        for n, key in enumerate(main_tests_dict.keys()):
            plotcol[n].pyplot(sns.catplot(kind='bar', data=res_df.iloc[:, n].reset_index(), x='index', y=key))


    # TODO: ADD possibility to put the SNPs in vertical line #statplots
    # TODO: Bonferroni correction in ---- and one other line (bit less that 7/8?) #statplots


def pval_plot_new(stat, xvar, pval_col, pval_cutoff=0.01, n_site=2, format_xaxis=True, out_dir='', title_supp=''):
    # stat = mw[(mw['index'] == val) & (mw['data'] == data)].copy()
    # xvar='cpg'
    # pval_col='p-val'
    # pval_cutoff=0.01
    # n_site=2
    # title_supp=f"MannWhitney_{val}_{data.replace('/', '-')}"
    # out_dir=''
    # format_xaxis=True
    stat = stat.copy()
    stat[['CHR', 'POS']] = stat[xvar].str.split(':', expand=True)[[0, 1]].astype(int)
    stat['CHR'] = stat['CHR'].astype('category')
    stat['minus_log10'] = -np.log10(stat[pval_col].astype(float))
    stat['cutoff'] = -np.log10(pval_cutoff/stat['cpg'].nunique())
    print(stat['cpg'].nunique())
    stat.sort_values('CHR', inplace=True)

    # Selection best snp
    best_cpg = stat[stat['cutoff'] < stat['minus_log10']].sort_values(
        'minus_log10').groupby('CHR').head(n_site)[
        'cpg'].tolist()
    stat.loc[(stat[xvar].isin(best_cpg)) &
        (stat['minus_log10'] > stat['cutoff']),#
        f'{xvar}_best'] = stat.loc[(stat[xvar].isin(best_cpg)) &
                    (stat['minus_log10'] > stat['cutoff']), 'cpg']
    stat.loc[stat[f'{xvar}_best'].isna(), f'{xvar}_best'] = ' '
    print(f'{xvar} ABOVE CUTOFF {title_supp}: {best_cpg}')

    # Format X-axis
    if format_xaxis:
        last_pos = 0
        for chr in stat.sort_values(['CHR', 'POS'])['CHR'].astype(int).unique():
            stat.loc[stat['CHR'] == chr, 'new_xvar'] = stat.loc[stat['CHR'] == chr, 'POS'] - stat.loc[stat['CHR'] == chr, 'POS'].min() + last_pos + 1
            last_pos = stat.loc[stat['CHR'] == chr, 'new_xvar'].max()
            # print(last_pos)

    # PLOT
    if format_xaxis:
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

    return g


# QUESTION: PUT in the stats analysis part ? #powerplot
# DONE: x = diff in %meth alt VS severe & y= rho haplotype test #powerplot
def power_analysis_plot(db, list_cpgs=[], count=0, count_cols=[]):
    # Spearman correlation haplotype
    # diff alt - means
    sql_sp = dict_to_sql_query(db, {'table': 'spearman_correlation',
                                    'cols': ['pval', 'rho', 'cpg'],
                                    'test': 'haplotype', 'data': '0/1'}, count=count, count_cols=count_cols)
    df_sp = pd.read_sql(sql_sp, con=db)
    sql_data = dict_to_sql_query(db, {'table': 'median',
                                      'cols': ['log_lik_ratio', '\"%meth\"',
                                             'sample_id', 'haplotype', 'cpg',
                                             'covid_snp', 'phenotype'],
                                      'Genotype': '0/1'}, count=count, count_cols=count_cols)
    df_data = pd.read_sql(sql_data, con=db)
    df_data = df_data[df_data['haplotype'] != 'other'].copy()
    df_data.dropna()
    st.dataframe(df_data.head())
    pivot = df_data.pivot_table(columns=['haplotype'], index=['sample_id', 'phenotype', 'covid_snp',
        'cpg'], values=['%meth'], aggfunc='median')
    st.dataframe(pivot.head(2))
    pivot.columns = [ '_'.join([a,b]) for a,b in pivot.columns]
    # pivot = pivot.reset_index().dropna()
    pivot['meth_delta'] = pivot['%meth_alt'] - pivot['%meth_ref']
    # pivot['scaled_minmax_haplo_delta'] = pivot['scaled_minmax_alt'] - pivot['scaled_minmax_ref']
    # pivot['scaled_median_haplo_delta'] = pivot['scaled_median_alt'] - pivot['scaled_median_ref']
    # log = px.scatter(pd.merge(df_sp, pivot, on='cpg'), x='rho', y='log_lik_ratio_haplo_delta', opacity=0.65,
    #     trendline='ols', trendline_color_override='darkblue')
    # scaled = px.scatter(pd.merge(df_sp, pivot, on='cpg'), x='rho', y='scaled_minmax_haplo_delta', opacity=0.65,
    #     trendline='ols', trendline_color_override='darkblue')
    # median = px.scatter(pd.merge(df_sp, pivot, on='cpg'), x='rho', y='scaled_median_haplo_delta', opacity=0.65,
    #     trendline='ols', trendline_color_override='darkblue')
    # st.plotly_chart(log)
    scaled = px.scatter(pd.merge(df_sp, pivot, on='cpg').drop_duplicates(), x='rho', y='meth_delta', opacity=0.65,
    trendline='ols', trendline_color_override='darkblue')
    st.plotly_chart(scaled)
    # st.plotly_chart(median)

############### DISTANCE AND INFO PLOTS PER SNP ################################
def distance_plot(db, list_cpgs=[], count=0, count_cols=[]):
    # TODO: CHANGE SQL QUERY !
    # Distance plot
    snp_ls = pd.read_sql("SELECT DISTINCT covid_snp FROM datas", con=db)['covid_snp'].tolist()
    snp_select = st.multiselect('SNP choice', sorted(snp_ls), snp_ls[0])

    if snp_select:
        df = pd.read_sql(f"SELECT * FROM datas WHERE covid_snp IN ({str(snp_select)[1:-1]}) AND haplotype <> 'other'", con=db)  # AND Genotype IS "0/1"
        df = df.groupby(['phenotype', 'sample_id', 'chromosome', 'cpg', 'covid_snp', 'Genotype', 'haplotype']).median().reset_index()
        df['distance_cpg_snp'] = abs(df['start'].astype(float) - df['pos'].astype(float))
        df['log_distance'] = np.log10(abs(df['distance_cpg_snp']))
        df.sort_values(by="distance_cpg_snp", inplace=True)

    if st.button('Distance plots'):
        fig1 = sns.relplot(kind='scatter', alpha=0.3, data=df, x='distance_cpg_snp', y='log_lik_ratio', hue='haplotype', col='phenotype', row='covid_snp', facet_kws={'sharey': True, 'sharex': False})
        fig1.set(xscale="log")
        fig1.map(plt.axvline, x=0, color='red')
        fig1.map(sns.rugplot)
        st.pyplot(fig1)

    if st.button('Info plots'):
        df_info = pd.DataFrame()
        progress_bar = st.empty()
        for n,snp in enumerate(snp_select):
            progress_bar.progress(n/df['snp'].nunique())
            chr, pos, _ = snp.split(':')
            blou = pd.read_sql(f"SELECT seqid, source, type, start, end, attributes FROM hsa_ensembl_annot WHERE start < {pos - 500000}  AND end > {pos + 500000} AND seqid IS {chr}",
                con=db)
            blou = info_formatting(blou)
            # TODO: Continue here to associate SNP with cpgs #info_distance_plots
            blou[['cpg', 'snp']] = snp
            df_info = df_info.append(blou)
        if not df_info.empty:
            st.dataframe(info)
            info.sort_values('snp', inplace=True)
            type = sns.countplot(data=info, y='type', row='snp')
            gene = sns.countplot(data=info, y='Name_simplified', row='snp')
            plots = st.columns(2)
            plots[0].pyplot(type)
            plots[1].pyplot(gene)
    ## All SNPs
    # fig = sns.relplot(kind='scatter', alpha=0.3, data=df, x='distance_cpg_snp', y='log_lik_ratio', hue='haplotype')
    # fig = sns.rugplot(data=df, x='distance_cpg_snp', y='log_lik_ratio', hue='haplotype', legend=False)
    # plt.suptitle(f'CpGs methylation likelyhood against distance all SNPs', y=1.05, va='top')
    # fig.set(xscale="log")
    # fig = plt.axvline(x=0, color='red')
    # st.pyplot(fig.get_figure())

        for snp in snp_ls:
            blou = df[df['covid_snp'] == snp]
            if blou['cpg'].nunique() > 3:
                st.write(snp, ' Number of cpg: ', blou['cpg'].nunique())
                fig1 = sns.relplot(kind='scatter', alpha=0.3, data=blou, x='log_distance', y='log_lik_ratio', hue='Genotype',
                                 col='phenotype')
                fig1 = sns.rugplot(data=blou, x='log_distance', y='log_lik_ratio', hue='Genotype')
                plt.suptitle(f'CpGs methylation likelyhood against distance to {snp}', y=1.05, va='top')
                fig1.set(xscale="log")
                st.pyplot(fig1.get_figure())
                fig2 = sns.relplot(kind='scatter', alpha=0.3, data=blou[blou['Genotype']=='0/1'], x='log_distance', y='log_lik_ratio', hue='haplotype',
                                 row='Genotype', col='phenotype')
                fig2.set(xscale="log")
                plt.suptitle(f'CpGs methylation likelyhood against distance to {snp} HET', y=1.05, va='top')
                st.pyplot(fig2)

        # snp = snp_ls[0]
        # chr, pos, _, _ = snp.split(':')
        # min = df["distance_cpg_snp"].min()
        # max = df["distance_cpg_snp"].max()
        # blou = pd.read_sql(f"SELECT seqid, source, type, start, end, attributes FROM hsa_ensembl_annot WHERE start > {int(pos) - min} AND end < {int(pos) + max} AND seqid IS {chr}", con=db)
        # yu = blou[blou['source'] != 'GRCh38'].copy()
        # st.write(min, max)
        # dict_color = {type:color for color, type in zip(sns.color_palette("tab10",n_colors=yu['type'].nunique() ), yu['type'].unique())}
        # dict_height = {type:n for n, type in zip(range(yu['type'].nunique()), yu['type'].unique())}
        # st.write(dict_color)
        # st.dataframe(yu)
        # fig, ax = plt.subplots()
        # plt.axhline(0.5, lw=1, alpha=0.5)
        # plt.axvline(int(pos), lw=2, color='red')
        # ax.set_xlim([yu['start'].min(), yu['end'].max()])
        # # ax.set_ylim([0, 1])
        # for i in yu.index:
        #     # st.write(str(yu.loc[i]))
        #     plt.plot( [int(yu.loc[i,'start']),int(yu.loc[i,'end'])], [dict_height[yu.loc[i,'type']],dict_height[yu.loc[i,'type']]], lw=4, alpha=0.2, label=yu.loc[i, 'type'], color=dict_color[yu.loc[i,'type']])
        # # ax.get_yaxis().set_visible(False)
        # plt.legend(bbox_to_anchor =(1.75, 0.5))
        # st.pyplot(fig)


############## First page to appear ############################################
def readme(db, list_cpgs=[], count=0, count_cols=[]):
    st.write('# README')
    df = pd.read_sql("SELECT DISTINCT * FROM cpg_snp", con=db)
    snp_ls = list(df['snp'].unique())
    info_snp = pd.read_table('/home/fanny/Work/EBI/covid_nanopore/significant_hits_COVID19_HGI_A2_ALL_leave_23andme_20210607.txt')

    st.write('## Covid SNP informations')
    st.dataframe(info_snp[info_snp['SNP'].isin(snp_ls)])
    if st.button('CpG & SNP couples'):
        st.dataframe(df)
    st.write('## Infos dataset')
    st.dataframe(pd.read_sql("SELECT chromosome, COUNT(DISTINCT cpg), COUNT(DISTINCT read_name), MAX(log_lik_ratio), MIN(log_lik_ratio) FROM log_methylation GROUP BY chromosome", con=db))

    st.dataframe(pd.read_sql("SELECT covid_snp, COUNT(DISTINCT cpg), COUNT(DISTINCT read_name), MAX(log_lik_ratio), MIN(log_lik_ratio) FROM log_methylation GROUP BY covid_snp", con=db))
    st.expander('cpgs').write(' - ' + '\n - '.join(sorted(list_cpgs)))
    for table in ['spearman_correlation', 'mann_whitney']:
        count_df = pd.read_sql(f"SELECT DISTINCT {table}.cpg,test,data FROM {table} JOIN counts_diff_means ON counts_diff_means.cpg={table}.cpg WHERE Mild+Severe > {count}", con=db)
        st.write(table)
        st.write(count_df.groupby(['test', 'data']).size().unstack())


def main(database_name):
    database_name=os.path.join(ABS_PATH, 'new_covid_snp', 'covid_snp_March2022.db')
    db = sqlite3.connect(database_name)

    dict_analysis = { 'Intro' : readme,
                     'Individual CpGs': individual_cpg_analysis,
                     'Association testing': stat_plots,
                     'power analysis': power_analysis_plot,
                     'distance': distance_plot}

    st.title('mQTLs analysis')
    analysis = st.sidebar.selectbox('Analysis type: ', dict_analysis.keys())

    # Run the appropriate function
    count = st.sidebar.slider('Minimal count for each cpg stats ', value=5, max_value=50)
    exp = st.sidebar.expander('More precise selection for counts')
    count_cols = exp.multiselect('Columns on which to apply count selection', ['read_count', 'Mild', 'Severe', 'alt', 'ref', 'other', '0/0', '0/1', '1/1'])
    dict_analysis[analysis](db, count=count, count_cols=count_cols)


if __name__ == '__main__':
    main(database_name=os.path.join(ABS_PATH, 'new_covid_snp', 'covid_snp_March2022.db'))
################################################################################
    ############################    NOTES   ###############################
################################################################################

# HOSTING WEBSITE ?
# https://www.ebi.ac.uk/birney-srv/medaka-ref-panel/

# FOR DISTANCE PLOTS EXAMPLE ADRIAN
# fig 5 MIKK panel
# https://github.com/birneylab/MIKK_genome_companion_paper/blob/master/docs/DNA_methylation/code/Interactive_comp_report.ipynb

# IMPORTANT NB: log_lik_ratio = np.log10(10**(log_lik_methylated) / 10**(log_lik_unmethylated))

################################################################################
######################### DATABASE MODIFICATION ################################
################################################################################
# import sqlite3
# import pandas as pd
# db = sqlite3.connect('new_covid_snp/covid_snp_March2022.db')
# c=db.cursor()
# pd.read_csv('Spearmann_corr_power.csv').to_sql('spearman_power', con=db)
# pd.read_csv('new_covid_snp/All_Spearmann_corr.csv').rename( columns = {'index': 'test', 'p-val': 'pval'}).to_sql('spearman_correlation', con=db, if_exists='replace')
# pd.read_csv('new_covid_snp/All_Counts_Diff_means.csv').to_sql('counts', con=db, if_exists='replace')
# pd.read_csv('Filtered_nano_bam_files.csv').to_sql('log_methylation', con=db, if_exists='replace')
# c.execute('SELECT * FROM datas LIMIT 10')
# c.fetchall()
# c.execute('ALTER datas RENAME TO log_methylation')
# c.close()


################# TO ADD IN THE DB #############################################
# TODO: Calculate percent methylation from number of reads
# import numpy as np
# import seaborn as sns
# # count = pd.read_sql("SELECT * FROM reads_samples", con=db)
# df = pd.read_sql("SELECT * FROM log_methylation", con=db)
# df.columns
# df['methylated'] = None
# thr = 0
# df.loc[df['log_lik_ratio'] > thr, 'methylated'] = 'methylated'
# df.loc[df['log_lik_ratio'] < -thr, 'methylated'] = 'unmethylated'
# df.loc[abs(df['log_lik_ratio']) < thr, 'methylated'] = 'unsure'
# g = sns.displot(data=df[abs(df['log_lik_ratio']) < 20], x='log_lik_ratio', hue='methylated', kde=True)
# save_plot_report(f'distribution_threshold_{thr}', g)
# px.histogram(df[abs(df['log_lik_ratio']) < 20 ], x="log_lik_ratio")


# perc = pd.read_sql('SELECT cpg, sample_id, haplotype, \"%meth\" FROM perc_meth', con=db)
#
# df.shape
# df.head()
# df = pd.merge(df, perc, on=['haplotype', 'sample_id', 'cpg'], how= 'outer')
# df.to_sql('median', if_exists='replace', con=db)

# TODO: Mix median_df + percent-meth to have the plots ewan wanted all along.

# DONE: Run t-tests between alt and ref in Mild-Severe
# cpg_ls = df['cpg'].unique()
# # count_cols = ['cpg', '0/0', '0/1', '1/1', 'alt', 'ref', 'Mild', 'Severe', 'means_ref-alt']
# phen_ls=['Severe', 'Mild']
# phen = phen_ls[1]
# gen_ls=['0/1']
# print(phen_ls, gen_ls)
# unit='covid_snp'
# cpg = '12:112919316:1'
# import pingouin as pg
# for cpg in cpg_ls:
#     add_col = {'cpg': cpg, 'var': 'ref-alt'}
#     if cpg == cpg_ls[0]:
#         pd.DataFrame(columns=['T', 'dof', 'tail', 'pval',  'CI95%', 'cohend', 'BF10', 'power', 'cpg', 'var', 'data']).to_csv('TTest.csv', mode='w', index=False)
#     samp = df[(df['cpg'] == cpg) & (df['Genotype'] == '0/1')].set_index(['cpg', 'sample_id', 'phenotype', 'Genotype', 'haplotype'])['log_lik_ratio'].unstack().reset_index().dropna()
#     for phen in phen_ls:
#         try:
#             res_sp = pg.ttest(samp[(samp['phenotype'] == phen) ]['ref'].tolist(), samp[(samp['phenotype'] == phen)]['alt'])
#             add_col['data'] = phen
#             # res_sp = stat.Spearman_correlation(measure=measure, var=vars)
#             if add_col: res_sp[list(add_col.keys())] = list(add_col.values())
#             res_sp.dropna(subset=['p-val']).reset_index().to_csv('TTest.csv',
#                                                 mode='a', header=False, index=False)
#         except Exception as err:
#             print(err)
#
# samp1 = df[(df['cpg'] == '12:112919316:1' ) & (df['Genotype'] == '0/1')].set_index(['cpg', 'sample_id', 'phenotype', 'Genotype', 'haplotype'])['log_lik_ratio'].unstack().reset_index().dropna()








# #### MANN WHITNEY CHANGE INDEXING :
# mw = pd.read_csv('new_covid_snp/All_Mann_Whitney.csv')
# mw = pd.concat([mw, mw['index'].str.split('_',expand=True)], axis=1)
# try:
#     mask = (mw[1].str.split('-', expand=True)[0] == mw[0])
#     mw.loc[mask, 'new_index'] = mw.loc[mask, 0] + ':' + (mw.loc[mask, 1].str.split('-', expand=True)[1] + '-' + mw.loc[mask, 2])
# except:
#     pass
# try:
#     mask = (mw[1].str.split('-', expand=True)[0] == mw[2])
#     mw.loc[mask, 'new_index'] = mw.loc[mask, 2] + ':' + mw.loc[mask, 0] + '-' + mw.loc[mask, 1].str.split('-', expand=True)[1]
# except:
#     pass
# try:
#     mask = (mw[1].str.split('-', expand=True)[1] == mw[0])
#     mw.loc[mask, 'new_index'] = mw.loc[mask, 0] + ':' + (mw.loc[mask, 1].str.split('-', expand=True)[0] + '-'+ mw.loc[mask, 2])
# except:
#     pass
# try:
#     mask = (mw[1].str.split('-', expand=True)[1] == mw[2])
#     mw.loc[mask, 'new_index'] = mw.loc[mask, 2] + ':' + mw.loc[mask, 0] + '-' + mw.loc[mask, 1].str.split('-', expand=True)[0]
# except:
#     pass
#
# mw.loc[(mw[1].astype(str).str.contains('-') == False) , 'new_index'] = mw.loc[(mw[1].astype(str).str.contains('-') == False), 'index']
#
# mw['detail_test'] = mw['data'] + ':' + mw['index']
#
# mw.loc[mw['new_index'].str.contains(':'), 'data'] = mw.loc[mw['new_index'].str.contains(':'), 'new_index'].str.split(':', expand=True)[0]
# mw['new_index'] = mw['new_index'].str.split(':', expand=True)[1]
# mw.loc[mw['new_index'].isnull(), 'new_index'] = mw.loc[mw['new_index'].isnull(), 'detail_test'].str.split(':', expand=True)[1]
#
#
# mw[['new_index', 'detail_test', 'p-val', 'data', 'cpg']].rename( columns = {'new_index': 'test', 'p-val': 'pval'}).to_sql('mann_whitney', con=db, if_exists='replace')


# # Add Annotation file to the Database
# col_names = ['seqid', 'source', 'type', 'start', 'end', 'score',
#     'strand', 'phase', 'attributes']
# df = pd.read_csv('/home/fanny/Work/EBI/Homo_sapiens.GRCh38.85.gff3.gz', compression='gzip', sep='\t',
#                  comment='#', low_memory=False, header=None, names=col_names)
# df.to_sql('hsa_ensembl_annot', con=db)
