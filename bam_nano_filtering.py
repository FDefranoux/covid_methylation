import pandas as pd
import glob

dir = 'nanopolish_grep_reads/*.csv'
# dir = 'nano*5b.csv'


def gather_dfs_fromdir(dir):
    info_per_file = pd.DataFrame()
    filtered_summary = pd.DataFrame(columns=glob.glob(dir),
                                    index=['initial shape', 'deletions',
                                           'read_name in several chr',
                                           'Non ALT or REF alleles',
                                           'final shape'])
    for file in glob.glob(dir):
        inter = pd.read_csv(file)
        # if 'PROM1' in file:
        #     inter['phenotype'] = 'Severe'
        # else:
        #     inter['phenotype'] = 'Mild'
        inter.loc[inter['name'].str.contains('PROM1'), 'phenotype'] = 'Severe'
        inter.loc[inter['name'].str.contains('PROM1') == False, 'phenotype'] = 'Mild'
        inter['cpg'] = inter['CHR'].astype(
            str) + ':' + inter['start'].astype(
            str) + ':' + inter['num_motifs'].astype(str)
        inter.loc[inter['num_motifs'] == 1, 'distance_cpg_snp'] = abs(
            inter['pos'].astype(int) - inter['start'].astype(int))
        info_per_file[file + '0'] = inter.nunique()
        inter = filtering_datas(inter, filtered_summary[file])
        info_per_file[file + '1'] = inter.nunique()
        if file == glob.glob(dir)[0]:
            inter.to_csv('Filtered_nano_bam_files_all_samples.csv', mode='w',
                         header=True, index=False)
        else:
            inter.to_csv('Filtered_nano_bam_files_all_samples.csv', mode='a',
                         header=False, index=False)
    print(filtered_summary.to_markdown())
    print('\n\n')
    print(info_per_file.T.to_markdown())


def filtering_datas(df, results):
    results.loc['initial shape'] = str(df.shape)

    # Deletions
    new_df = df[df['ref'].str.len() < 2].copy()
    results.loc['deletions'] = int(
        df.shape[0] - new_df.shape[0])

    # Filtering out the reads_names associated with several CHR
    double_chr = new_df.groupby(
        ['read_name', 'CHR']).size().unstack().dropna(thresh=2).index
    if len(double_chr) > 0:
        new_df.drop(double_chr, inplace=True)
    results.loc['read_name in several chr'] = int(
        df.shape[0] - new_df.shape[0])

    # Drop SNP with 'other' allele (non-ref non-alt)
    new_df = new_df[new_df['Gen'] == 'other'].copy()
    results.loc['Non ALT or REF alleles'] = int(
        df.shape[0] - new_df.shape[0])
    results.loc['final shape'] = str(new_df.shape)

    return new_df


def main(dir):
    gather_dfs_fromdir(dir)


if __name__ == '__main__':
    main(dir)
