import pandas as pd
import pysam
import os
import sys

# Variable definition
file_list = 'file_list.txt'
hits_table = 'significant_hits_COVID19_HGI_A2_ALL_leave_23andme_20210607.txt'
output = {'datas': 'Bam_per_significant_SNPs.csv',
          'reference': 'Ref_names_bam_significant_SNPs.csv'}


class SamFiles:
    # def __init__(self):
    #     self.error_files = {}
    def region(file):
        return file.fetch

    def reads(file):
        return file.find

    def sam_iterators(func, args, cols=None):
        try:
            args = tuple(args)
            rows = [x for x in func(*args)]
            array = [str(r).split('\t') for r in rows]
            df = pd.DataFrame(array, columns=cols)
        except Exception as err:
            print(f'Error with {args}\n{err}')
            # error_files[file].append(f'Error with {args}\n{err}')
            df = pd.DataFrame()
        return df

    def bam_ref_names(bam):
        # Bam file header saving
        try:
            bam_head = pd.DataFrame([el.split('\t')
                                     for el in str(bam.header).split('\n')])
            bam_head = bam_head[bam_head[0] == '@SQ'][[1, 2]]
            ref_namelenght = bam_head[1].str.cat(bam_head[2], '_')
            return ref_namelenght.T
        except:
            pass

    def open(file):
        try:
            if '.bam' in file:
                sam_file = pysam.AlignmentFile(file, 'rb')
            else:
                sam_file = pysam.TabixFile(file)
            return sam_file
        except OSError as err:
            # error_files[file] = f'{type} file {file} or its index were not found'
            return err


def region_select_fromSNP(pos_df):
    pos_df = pos_df.copy()
    pos_df['#CHR'] = pos_df['#CHR'].astype(str)
    # Creating temporary columns
    pos_df['new_POS'] = (pos_df['POS'] - 1)
    pos_df['POS_end'] = (pos_df['POS'] + 1)
    # Regions to look up
    pos_df['tuple'] = pos_df.set_index(
        ['#CHR', 'new_POS', 'POS_end']).astype(str).index.tolist()
    # Concat of CHR and POS as potentiel index
    pos_df['SNP_POS'] = pos_df['#CHR'] + '_' + pos_df['POS'].astype(str)
    # Removing temporary columns
    pos_df.drop(['new_POS', 'POS_end', 'POS'], axis=1)
    return pos_df


def recup_bam_perregion(file_ls, bam_list, method='',
                        output_names={'datas': 'Bam_summary.csv',
                                      'reference': 'Ref_table.csv'}):
    bam_cols = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ',
                'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', '?']
    df_final = pd.DataFrame()
    for file in file_ls:
        bam_file = SamFiles.open(file)
        SamFiles.bam_ref_names(bam_file).to_csv(output_names['reference'],
                                                mode='a', header=False)

        for region in bam_list:
            if method == 'bam_regions':
                bam_df = SamFiles.sam_iterators(
                    SamFiles.region(bam_file), region, cols=bam_cols)
                bam_df['SNP_hit'] = region[1] + 1
                bam_df['file'] = file
                df_final = pd.concat([df_final, bam_df], axis=0)

        if method == 'bam_reads':
            name_indexed = pysam.IndexedReads(bam_file)
            name_indexed.build()
            for read in bam_list:
                bam_df = SamFiles.sam_iterators(
                    name_indexed, SamFiles.reads, read, cols=bam_cols)
                bam_df['file'] = file
                df_final = pd.concat([df_final, bam_df], axis=0)

        if file == file_list[0]:
            df_final.to_csv(output_names['datas'], mode='a',
                            header=True, index=False)
        else:
            df_final.to_csv(output_names['datas'], mode='a',
                            header=False, index=False)


def main(file_list, hits_table):
    # Remove existing file
    for name in output:
        if name in os.listdir():
            print('Removing previous file')
            os.remove(name)

    # Reading Hits table
    file_ls = pd.read_table(file_list, header=0).iloc[:, 0].tolist()
    hits_df = pd.read_table(hits_table)
    assert hits_df[hits_df[['#CHR', 'POS']].duplicated(
        )].empty, 'Chromosome and SNP Position not enough to create unique indexes'
    table_region = region_select_fromSNP(
        hits_df[['#CHR', 'POS', 'REF', 'ALT']])
    list_region = set(table_region['tuple'])

    # Look up bam files
    recup_bam_perregion(file_ls, list_region,
                        method='bam_regions', output_names=output)


if __name__ == '__main__':
    main(file_list, hits_table)
