N=${LSB_JOBINDEX}'q;d'
reads_file=$(sed ${N} read_name_files.txt)
echo $reads_file
file_basename=$(basename $reads_file)
nanopolish_file='nanopolish_indexed/'${file_basename%.*}'.tsv.gz'
echo $nanopolish_file

zcat $nanopolish_file | grep -f $reads_file >> nanopolish_grep_reads/${file_basename}
