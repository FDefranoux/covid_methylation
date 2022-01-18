import pandas as pd
import pysam
import os
import argparse

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
            ref_namelenght = pd.DataFrame(
                bam_head[1].str.cat(bam_head[2], '_'))
            return ref_namelenght
        except Exception as err:
            print(err)
            return pd.DataFrame()

    def open(file):
        try:
            if '.bam' in file:
                sam_file = pysam.AlignmentFile(file, 'rb')
            else:
                sam_file = pysam.TabixFile(file)
            return sam_file
        except OSError as err:
            print(file)
            print(err)
            # error_files[file] = f'{type} file {file} or its index were not found'


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


def recup_bam_perregion(file_ls, bam_list, output_names={'datas': 'Bam_summary.csv',
                                                         'reference': 'Ref_table.csv'}):
    bam_cols = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ',
                'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', '?']
    for file in file_ls:
        bam_file = SamFiles.open(file)
        ref_names = SamFiles.bam_ref_names(bam_file)
        try:
            ref_names.columns = [file]
            ref_names.T.to_csv(output_names['reference'],
                               mode='a', header=False)
        except:
            print(f'No reference name for file {file}')

        for region in bam_list:
            bam_df = pd.DataFrame()
            try:
                bam_df = SamFiles.sam_iterators(
                    SamFiles.region(bam_file), region, cols=bam_cols)
                bam_df['SNP_hit'] = str(
                    region[0]) + '_' + str(region[1] + 1)  # WORKING?
                bam_df['file'] = file
            except:
                print(f'Error with iterating over file {file}')

            try:
                if file == file_ls[0]:
                    print(bam_df.head(3))
                    bam_df.to_csv(output_names['datas'], mode='a',
                                  header=True, index=False)
                else:
                    bam_df.to_csv(output_names['datas'], mode='a',
                                  header=False, index=False)
            except:
                print(f'unexpected error during the saving of file {file}')


def recup_bam_perregion_pileup(file_ls, bam_list, output_names={'datas': 'Bam_summary.csv',
                                                                'reference': 'Ref_table.csv'}):
    bam_cols = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ',
                'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', '?']
    for file in file_ls:
        bam_file = SamFiles.open(file)
        ref_names = SamFiles.bam_ref_names(bam_file)
        try:
            ref_names.columns = [file]
            ref_names.T.to_csv(output_names['reference'],
                               mode='a', header=False)
        except:
            print(f'No reference name for file {file}')

        bases_and_read_ids = {}
        samfile = pysam.Samfile(file, "rb")
        for region in bam_list:
            for pup in samfile.pileup(*region):
                for read in pup.pileups:
                    if not read.is_del:
                        if pup.pos == region[2]:
                            bases_and_read_ids[read.alignment.query_name] = read.alignment.query_sequence[read.query_position]
                            print('\tbase in read %s = %s' %
                                  (read.alignment.query_name,
                                   read.alignment.query_sequence[read.query_position]))
            bam_df = pd.DataFrame(bases_and_read_ids)
            try:
                if file == file_ls[0]:
                    bam_df.to_csv(output_names['datas'], mode='a',
                                  header=True, index=False)
                else:
                    bam_df.to_csv(output_names['datas'], mode='a',
                                  header=False, index=False)
            except:
                print(f'unexpected error during the saving of file {file}')


def main(file_list, hits_table, output, method='pileup'):
    # Remove existing file
    for name in output.values():
        if name in os.listdir():
            print(f'Removing previous file {name}')
            os.remove(name)

    # Reading Hits table
    file_ls = pd.read_table(file_list, header=None).iloc[:, 0].tolist()
    hits_df = pd.read_table(hits_table)
    assert hits_df[hits_df[['#CHR', 'POS']].duplicated(
        )].empty, 'Chromosome and SNP Position not enough to create unique indexes'
    table_region = region_select_fromSNP(
        hits_df[['#CHR', 'POS', 'REF', 'ALT']])
    list_region = set(table_region['tuple'])
    print('list region len:', len(list_region))

    # Look up bam files
    print(method)
    if method == 'pandas':
        recup_bam_perregion(file_ls, list_region, output_names=output)
    else:
        recup_bam_perregion_pileup(file_ls, list_region, output_names=output)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Recuperation of BAM datas')
    parser.add_argument('-f', '--file_list', type=str, default=file_list)
    parser.add_argument('-t', '--hits_table', type=str, default=hits_table)
    parser.add_argument('-o', '--output', default=output)
    parser.add_argument('-m', '--bam_method', type=str, default='pileup')
    args = parser.parse_args()
    main(**vars(args))
