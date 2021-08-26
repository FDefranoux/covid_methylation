import pandas as pd
import pysam
import os

file_list = 'all_files.txt'
file_list = '../covid_nanopore/all_files.txt'
os.chdir('../covid_nanopore/')
#in /nfs/research/birney/users/fanny/Covid
hits_table = 'significant_hits_COVID19_HGI_A2_ALL_leave_23andme_20210607.txt'


def region_select_fromSNP(pos_df):
    pos_df = pos_df.astype(int).copy()
    pos_df['new_POS'] = (pos_df['POS'] - 11).astype(str).str[:-6]
    pos_df['POS_end'] = (pos_df['new_POS'].astype(int) + 1) * 1000000
    pos_df['new_POS'] = pos_df['new_POS'] + '000000'
    pos_df['#CHR'] = pos_df['#CHR'].astype(str)
    pos_df[['new_POS', 'POS', 'POS_end']] = pos_df[['new_POS', 'POS', 'POS_end']].astype(int)
    pos_df.set_index('POS', inplace=True)
    pos_df['tuple'] = pos_df.set_index(['#CHR', 'new_POS', 'POS_end']).index.tolist()
    return pos_df['tuple']


def tabix_listregion_listfiles(list_region, file, cols=[]):
    df_final = pd.DataFrame(columns = cols)
    error_region = []
    tbx = pysam.TabixFile(file)
    for n, region_n in enumerate(list_region):
        try:
            rows_n = [x for x in tbx.fetch(*region_n)]
        except ValueError:
            rows_n = []
            error_region.append((file, region_n))
        if rows_n:
            try:
                df_reg = pd.DataFrame([r.split('\t') for r in [x for x in rows_n]], columns=cols)
            except ValueError as err:
                print(f'\n{err} with file {file} at region {region_n}')
                print(pd.DataFrame([r.split('\t') for r in [x for x in rows_n]]).head(5).to_markdown())
            df_reg['region'] = str(region_n)
            df_final = pd.concat([df_final, df_reg], axis=0)
    df_final['#CHR'] = df_final['#CHR'].astype(int)
    print(f'\n{file} with region error : {error_region}')
    return df_final


def associate_snp_to_sequence(df, snps):
    df = df[df.duplicated() == False].copy()
    # If working with regions, how to associate SNP value to region asked
    df['start_sequence'] = df['start'].astype(int) - 5
    df['end_sequence'] = df['start_sequence'].astype(int) + df['sequence'].str.len()
    df_out = pd.DataFrame()
    for snp in snps:
        df_int = df[(df['start_sequence'] < snp) & (df['end_sequence'] > snp)].copy()
        df_int['SNP'] = snp
        df_out = pd.concat([df_out, df_int], axis=0)
    return df_out


def associate_snp_alleles(df, snps_refs):
    df['rel_snp'] = df['SNP'] - df['start_sequence']
    df['allele'] = '0'
    df['#CHR'] = df['#CHR'].astype(int)
    for n in df['rel_snp'].unique():
        df.loc[df['rel_snp'] == n, 'allele'] = df['sequence'].str[n]
    allele_df = pd.merge(snps_refs, df, right_on=['SNP', '#CHR'], left_on=['POS', '#CHR'])
    allele_df.loc[allele_df['REF'] == allele_df['allele'], 'Genotype'] = '0'
    allele_df.loc[allele_df['ALT'] == allele_df['allele'], 'Genotype'] = '1'
    # Is this normal to have some allele different from our count ?
    allele_df.loc[(allele_df['ALT'] != allele_df['allele']) & (allele_df['REF'] != allele_df['allele']), 'Genotype'] = '2'
    print(allele_df['Genotype'].value_counts())
    return allele_df


def main(file_list, hits_table):
    if 'significant_SNPs_across_nanopore.csv' in os.listdir():
        print('Removing previous file')
        os.remove('significant_SNPs_across_nanopore.csv')
    file_ls = pd.read_table(file_list, header=None)[0].tolist()
    hits_df = pd.read_table(hits_table)
    # os.chdir('../covid_nanopore')
    # hits_df['#CHR'] = 1  # TO REMOVE
    # hits_df = hits_df[hits_df['#CHR'] == 1]  # TO REMOVE
    cols = ['#CHR', 'strand', 'start', 'end', 'read_name', 'log_lik_ratio',
            'log_lik_methylated', 'log_lik_unmethylated', 'num_calling_strands',
            'num_motifs', 'sequence'] # find a way to select columns directly from the file
    list_region = set(region_select_fromSNP(hits_df[['#CHR', 'POS']]).values)
    error_dict = {}
    for file in file_ls:
        try:
            df_final = tabix_listregion_listfiles(list_region, file, cols=cols)
            df_snps = associate_snp_to_sequence(df_final[['sequence', 'start', '#CHR']], hits_df['POS'])
            df_final = pd.merge(df_final, df_snps, on=['sequence', 'start', '#CHR'])
            snp_genotype = associate_snp_alleles(df_final[['start_sequence', 'SNP', 'sequence', '#CHR']].copy(),
                                                 hits_df[['#CHR', 'POS', 'REF', 'ALT']].copy())
            df_final = pd.merge(df_final, snp_genotype, on=['#CHR', 'start_sequence', 'sequence', 'SNP'])
            # To find the other motifs
            df_final['distance_snp_motif'] = df_final['start'].astype(int) - df_final['SNP']
            # df_final[df_final['num_motifs'].astype(int) > 1]['sequence'].str.find('GC')
            df_final['File'] = file
            if 'PROM1' in file:
                df_final['Phenotype'] = 'Severe'
            else:
                df_final['Phenotype'] = 'Mild'
            df_final.to_csv('significant_SNPs_across_nanopore.csv', mode='a')
        except Exception as err:
            error_dict[file] = err
        print(error_dict)

def GC_motifs(df):

    cols = pd.read_csv('significant_SNPs_across_nanopore.csv').columns
    df = pd.read_csv('test.txt', header=None, index=None)
    df.head(2)
    df.columns = cols.tolist() + ['SNP', '?1', '?2', 'REF', 'ALT', 'len', 'Allele', '?3', 'distance_snp_motif', 'Phenotype']
    df = df.copy()


if __name__ == '__main__':
    main(file_list, hits_table)
# df_final.columns
# import seaborn as sns
# sns.boxplot(x = df_final['Phenotype'], y=df_final['log_lik_methylated'].astype(float), hue=df_final['Genotype'])

# TODO: Add gestion of the number of GC motifs to handle
# TODO: Measure the distance between GC motif and SNP
# TODO: Plot the right things
# TODO: Think about the statistics to perform
# TODO: genotype homozygote ref or homozygote alt --> table of counts
# TODO: Using all reads for one sample try to do an average sequence ?
