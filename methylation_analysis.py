import pandas as pd
from cigar import Cigar


hits_table = 'significant_hits_COVID19_HGI_A2_ALL_leave_23andme_20210607.txt'
bam_file = 'Bam_test.csv'


def region_select_fromSNP(pos_df):
    pos_df = pos_df.astype(int).copy()
    pos_df['new_POS'] = (pos_df['SNP_POS'] - 1)
    pos_df['POS_end'] = (pos_df['SNP_POS'] + 1)
    pos_df['#CHR'] = pos_df['#CHR'].astype(str)
    pos_df['tuple'] = pos_df.set_index(
        ['#CHR', 'new_POS', 'POS_end']).astype(str).index.tolist()
    pos_df['SNP_POS'] = pos_df['#CHR'] + '_' + pos_df['SNP_POS'].astype(str)
    return pos_df[['#CHR', 'tuple', 'SNP_POS']]


def associate_snp_alleles(df_final, snps_refs):
    df = df_final[df_final.duplicated() == False].copy()
    df['rel_snp'] = df['SNP'].astype(int) - df['start'].astype(int)
    df['allele'] = '0'
    df[['SNP', '#CHR']] = df[['SNP', '#CHR']].astype(int)
    for n in df['rel_snp'].unique():
        df.loc[df['rel_snp'] == n, 'allele'] = df['sequence'].str[n]
        allele_df = pd.merge(snps_refs, df, right_on=[
            'SNP', '#CHR'], left_on=['POS', '#CHR'])
        return allele_df


def cigar2SNP(row, pos=1):
    seq = row['SEQ']
    cigar_seq = row['CIGAR']
    snp = int(row['rel_SNP']) - 1
    # Transforming CIGAR sea in list of tuples
    c = Cigar(cigar_seq)
    cigar_list = list(c.items())
    assert len(c) == len(seq), 'Sequence and cigar are not same lenght'
    new_seq = ''
    i, d, sm = 0, 0, 0  # Serves as verification of right lenghts
    for item in cigar_list:
        if item[1] == 'S' or item[1] == 'M':
            # recuperation of the read sequence for the matches and mismatches
            new_seq += seq[int(pos): int(pos) + int(item[0])]
            pos += int(item[0])
            sm += int(item[0])
            # Replacement of the deletions by empty characters (indexing similar to template)
        elif item[1] == 'D':
            new_seq += int(item[0]) * '-'
            pos += int(item[0])
            d += int(item[0])
        # We do not take in account the insertion as we want to have same indexing as template
        elif item[1] == 'I':
            i += int(item[0])
        else:
            print('ERROOOOR letter not known ', item)
    # assert len(new_seq) == (
    #     sm+d), 'New sequence do not have same lenght has template one  '
    # assert len(Cigar(cigar_seq)) == (sm+i)
    if len(new_seq) <= snp:
        print('LEN ERROR')
        return '??'
    else:
        return new_seq[snp]


def bam_filtering(bam, mapq=20, full=False):
    index_flag = bam[bam['FLAG'].isin(['0', '16']) == False].index.tolist()
    index_pos = bam[bam['POS'].astype(int) <= 0].index.tolist()
    index_rname = bam[bam['RNAME'] == '0'].index.tolist()
    index_relsnp = bam[bam['rel_SNP'].astype(int) < 0].index.tolist()
    # Check here why some SNPs are not in the position
    index_tlen_snp = bam[bam['TLEN'].astype(
        int) < bam['rel_SNP'].astype(int)].index.tolist()
    index_mapq = bam[bam['MAPQ'].astype(int) < mapq].index.tolist()
    all_index_drop = set(index_flag + index_pos + index_rname + index_relsnp
                         + index_tlen_snp + index_mapq)
    bam.drop(all_index_drop, axis=0, inplace=True)
    if full:
        print(
            f'FLAG not 0 or 16: {len(index_flag)}',
            f'\nPOS negative: {len(index_pos)}',
            f'\nRNAME is 0: {len(index_rname)}',
            f'\nrel_SNP negative: {len(index_relsnp)}',
            f'\nTLEN < SNP position: {len(index_tlen_snp)}',
            f'\nMAPQ < {mapq}: {len(index_mapq)}',
            f'\nTotal index filtered: {len(all_index_drop)}'
            )
    print(f'Filtered out {len(all_index_drop)} row(s)')


def reverse_seq(seq):
    trans = 'ATGC'.maketrans("ATGC", "TACG")
    new_seq = seq.str.translate(trans)
    return new_seq


def main(hits_table, bam_file):
    # SNPs table
    hits_df = pd.read_table(hits_table)
    hits_df[['#CHR', 'POS']] = hits_df[['#CHR', 'POS']].astype(int)
    hits_df.rename(columns={'POS': 'SNP_POS'}, inplace=True)

    # Opening table to associate region with SNP
    region_file = region_select_fromSNP(
        hits_df[['#CHR', 'SNP_POS']]).astype(str)

    # Opening bam file summary
    bam = pd.read_table(bam_file, header=0, dtype=object, sep=',')
    bam.replace(region_file.set_index(
        'tuple').to_dict()['SNP_POS'], inplace=True)
    bam[['#CHR', 'SNP_POS']] = bam['region'].str.split('_', expand=True)
    bam['rel_SNP'] = bam['SNP_POS'].astype(int) - bam['POS'].astype(int)

    # Filtering the reads
    bam_filtering(bam, mapq=20, full=True)

    # Reversing the sequences corresponding to flag=16
    flag16 = bam.loc[bam['FLAG'] == '16', ['SEQ']].copy(deep=False)
    flag16 = flag16.apply(lambda x: reverse_seq(x))
    bam.loc[bam['FLAG'] == '16', ['SEQ']] = flag16

    # Genotyping
    bam['Allele'] = bam.apply(lambda x: cigar2SNP(x), axis=1)
    bam = pd.merge(hits_df[['#CHR', 'SNP_POS', 'REF', 'ALT']].astype(str), bam.astype(
        str), right_on=['#CHR', 'SNP_POS'], left_on=['#CHR', 'SNP_POS'])
    bam.loc[bam['Allele'] == bam['REF'], 'Genotype'] = 'REF'
    bam.loc[bam['Allele'] == bam['ALT'], 'Genotype'] = 'ALT'
    bam.loc[bam['Genotype'].isna(), 'Genotype'] = 'OTHER'
    bam['Genotype'].value_counts()


# pd.options.display.html.table_schema = True
# pd.options.display.max_rows = None


# TABLE COUNTS
# TODO: genotype homozygote ref or homozygote alt --> table of counts
# df.loc[df['strand'] == '+', ['File', '#CHR', 'SNP',
#                              'sequence', 'start', 'end', 'REF', 'ALT', 'allele']]
# df[['File', '#CHR', 'SNP', 'allele']].sort_values(
#     ['#CHR', 'SNP']).value_counts().sort_index()
# df.groupby(['File', 'SNP', 'Genotype']).count()
# # Until you solve the problem of the columsn names
# df = df[df['SNP'] == df['POS']]


#  GC MOTIFS
# TODO: Add gestion of the number of GC motifs to handle
# TODO: Measure the distance between GC motif and SNP
# # To find the other motifs
# df_final['distance_snp_motif'] = df_final['start'].astype(int) - df_final['SNP']
# df_final[df_final['num_motifs'].astype(int) > 1]['sequence'].str.find('GC')
#     df[df['num_motifs'] != '1']['sequence']
#     df['num_motifs'].astype(int).max()
#     df['sequence'].str.rfind('GC').max()
#     GC_motifs = pd.DataFrame()
#     for motif_n in range(1, df['num_motifs'].astype(int).max()):
#     df['sequence'].str.rfind('GC')
#     blou.value_counts()
#     df = df.copy()

# PLOTS
# TODO: Plot the right things
# import seaborn as sns
# sns.boxplot(x = df_final['Phenotype'], y=df_final['log_lik_methylated'].astype(float), hue=df_final['Genotype'])

# TODO: Think about the statistics to perform
# TODO: Using all reads for one sample try to do an average sequence ?
