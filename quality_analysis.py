import pandas as pd
import os
import pysam
import glob
import matplotlib.pyplot as plt
import seaborn as sns

# TODO: descriptive statistics after grepping ??
# TODO: descriptive statistics number of cpg per samples (boxplot per chr?) in Filtered dataset

dir = 'lsf_out.txt'
lsb = True

def lsf_arrray(file_list):
    lsb_index = int(os.environ['LSB_JOBINDEX'])-1
    file_array = file_list[lsb_index]
    return file_array


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


def quality_analysis(file, region_list, nanocols):
    nano_file = SamFiles.open(file)
    for region in region_list:
        nano_df = pd.DataFrame()
        try:
            nano_df = SamFiles.sam_iterators(SamFiles.region(nano_file), region, cols=nanocols)
            pd.DataFrame(nano_df[['strand', 'start', 'end', 'read_name']].nunique(), columns=[os.path.basename(file)[:-3] + '_' + str(region)]).T.to_csv(f'Nunique_nanopolish_indexed_{os.path.basename(file)}.csv', mode='a', header=False)
        except Exception as err:
            print(f'Error with iterating over file {os.path.basename(file)}-{region}', flush=True)
            print(err, flush=True)
            print(region, nano_df.shape, flush=True)
            print(nano_df.head(5).to_markdown(), flush=True)
        del nano_df


def main(dir, lsb=False):
    region_list = [(str(chr), pos, pos+1000000) for chr in range(1, 24) for pos in range(1,250000000, 1000000)]
    nanocols = ['chromosome', 'strand', 'start', 'end', 'read_name',
                'log_lik_ratio', 'log_lik_methylated', 'log_lik_unmethylated',
                'num_calling_strands', 'num_motifs', 'sequence']

    if os.path.isdir(dir):
        ls_files = os.listdir(dir)
    elif '*' in dir:
        ls_files = glob.glob(dir)
    else:
        with open(dir, 'r') as f:
            ls_files = f.readlines()
        ls_files = [file_str[:-1] for file_str in ls_files]

    assert isinstance(ls_files, list), 'List of files is not of the right type !'

    if not lsb:
        for file in ls_files:
            print(file, flush=True)
            quality_analysis(file, region_list, nanocols)
    else:
        file = lsf_arrray(ls_files)
        print(file, flush=True)
        quality_analysis(file, region_list, nanocols)


if __name__ == '__main__':
    main(dir, lsb=lsb)


# ## Analysis part
# error_df = pd.read_table('Errors_quality_file.txt', sep=':', header=None )
# error_df['index'] = error_df[0].str.split('-', expand=True)[0]
# dict_df = error_df.dropna().copy()
# dict_df[[0, 1]] = dict_df[1].str.split('-', expand=True)
# dict_df[0] = dict_df[0].str[31:]
# dict_index = dict_df.set_index('index')[0].to_dict()
# dict_region = dict_df[1].to_dict()
# len(dict_index)
# blou = error_df.drop(error_df.dropna().index).dropna(axis=1, how='all').copy()
# blou.drop(blou[blou[0] == '--'].index, inplace=True)
# blou[0] = blou[0].str.split('-', n=1, expand=True)[1]
# blou['index'].replace(dict_index, inplace=True)
# blou.groupby('index').nunique().sort_values(by=0)
# blou[0].value_counts()
#
#
#
# quality_df = pd.read_csv('Nunique_nanopolish_indexed_all.csv').reset_index()
# quality_df.columns = ['region', 'strand', 'start', 'end', 'read_name']
# quality_df[['file', 'region']] = quality_df['region'].str.rsplit('_', n=1, expand=True)
# quality_df['file'] = quality_df['file'].str[:-4]
# quality_df['chr'] = quality_df['region'].str.split(',', expand=True)[0].str[2:-1].astype(int)
# sns.catplot(data=quality_df[quality_df['chr'] == 1], x='region', y='read_name', row='chr', hue='file')
#
#
# row_to_drop = quality_df[(quality_df['strand'] == 0)].index.tolist()
#
# # Rows without 2strands
# row_to_drop += blou[blou['strand'] != 2]  # TO DROP ???
#
# blou.groupby(['file']).nunique()['chr'].value_counts().sort_index()
