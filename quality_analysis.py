import pandas as pd
import seaborn as sns
import os
import pysam
import glob

# TODO: descriptive statistics after grepping ??
# TODO: descriptive statistics number of cpg per samples (boxplot per chr?) in Filtered dataset

dir = 'nanopolish_indexed/*.tsv.gz'
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


def quality_analysis(nanopolish_file, region_list, nanocols):
    nano_file = SamFiles.open(nanopolish_file)
    for region in region_list:
        nano_df = pd.DataFrame()
        try:
            nano_df = SamFiles.sam_iterators(SamFiles.region(nano_file), region, cols=nanocols)
            pd.DataFrame(nano_df[['strand', 'start', 'end', 'read_name']].nunique(), columns=[os.path.basename(file)[:-3] + '_' + str(region)]).T.to_csv(f'Nunique_nanopolish_indexed.csv', mode='a', header=False)
        except Exception as err:
            print(f'Error with iterating over file {file}-{region}', flush=True)
            print(err, flush=True)
            print(region, nano_df.shape, flush=True)
            print(nano_df.head(5).to_markdown(), flush=True)
        del nano_df


def main(dir, lsb=False):
    region_list = [(str(chr), pos, pos+1000000) for chr in range(1, 24) for pos in range(1,250000000, 1000000)]
    nanocols = ['chromosome', 'strand', 'start', 'end', 'read_name',
                'log_lik_ratio', 'log_lik_methylated', 'log_lik_unmethylated',
                'num_calling_strands', 'num_motifs', 'sequence']

    if not lsb:
        for file in glob.glob(dir):
            print(file, flush=True)
            quality_analysis(file, region_list, nanocols)
    else:
        file = lsf_arrray(glob.glob(dir))
        print(file, flush=True)
        quality_analysis(file, region_list, nanocols)


if __name__ == '__main__':
    main(dir, lsb=lsb)
