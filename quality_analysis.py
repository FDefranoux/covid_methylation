import pandas as pd
import seaborn as sns
import os
import pysam
import glob

dir = 'nanopolish_indexed/*.tsv.gz'
# file = 'Filtered_nano_bam_files_test.csv'

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


def main(dir):
    region_list = [(str(chr), pos, pos+1000000) for chr in range(1, 24) for pos in range(1,250000000, 1000000)]
    len(region_list)
    # nanocols = ['file', 'region', 'chromosome', 'strand', 'start', 'end', 'read_name',
                # 'log_lik_ratio', 'log_lik_methylated', 'log_lik_unmethylated',
                # 'num_calling_strands', 'num_motifs', 'sequence']
    for file in glob.glob(dir):
        print(file, flush=True)
        # pd.DataFrame(nanocols).T.to_csv(f'Nunique_nanopolish_indexed_{file}.csv', mode='w', header=False)
        # pd.DataFrame(nanocols).T.to_csv('Size_nanopolish_indexed.csv', mode='w', header=False)
        # Opening the allele_table
        nano_file = SamFiles.open(file)
        n = 0
        for region in region_list:
            print(region, flush=True)
            nano_df = pd.DataFrame()
            try:
                nano_df = SamFiles.sam_iterators(SamFiles.region(nano_file), region)
                print(nano_df.memory_usage(deep=True), nano_df.shape, flush=True)
                # nano_df.columns = nanocols
                nano_df[['file', 'region', 'N']] = file, str(region), n
                nano_df[['region']].head(1).to_csv(f'Nunique_nanopolish_indexed_{os.path.basename(file)}.csv', mode='a', header=False)
                nano_df[['N']].head(1).to_csv(f'Nunique_nanopolish_indexed_{os.path.basename(file)}.csv', mode='a', header=False)
                nano_df.nunique().to_csv(f'Nunique_nanopolish_indexed_{os.path.basename(file)}.csv', mode='a', header=False)
                n = n + 1
                # pd.DataFrame(nano_df.groupby(['file', 'region']).size()).T.to_csv('Size_nanopolish_indexed.csv', mode='a', header=False)
            except Exception as err:
                print(err, flush=True)
                print(f'Error with iterating over file {file}-{region}', flush=True)

    # df = pd.read_csv(file)
    # os.system('mkdir quality')
    # for name in df['sample_id'].unique()[0]:
    #     print(file, flush=True)
    #     df[df['sample_id'] == name].groupby('chromosome').nunique().to_csv(f'quality_test/Nunique_{name}.csv')
    #     df[df['sample_id'] == name].groupby('chromosome').size().to_csv(f'quality_test/Size_{name}.csv')



if __name__ == '__main__':
    main(dir)
