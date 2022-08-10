#!/usr/bin/env nextflow

params.outdir = "$baseDir/covid_snp_all_updated_Summer_2022"
// params.basecall = "$baseDir/base_called_from_bam_files/covid_snps_finemapped/gcc000*.txt"
params.basecall = "$baseDir/base_called_from_bam_files/covid_snps_updated_summer2022/gcc008*.txt"
params.nanopolish= "/hps/nobackup/birney/projects/gel_methylation/nanopolish/gcc008*.tsv.gz"
params.snp_type = "covid_snp"
// params.snp_list = "$baseDir/base_called_from_bam_files/control_snp_finemapped"

// Channel creation
// nano_files = Channel.fromPath(params.nanopolish, checkIfExists:true)
// basecall_files = Channel.fromPath(params.basecall, checkIfExists:true)

Channel.fromFilePairs([params.basecall, params.nanopolish], flat:true).set {file_pairs}
// .ifEmpty( error 'ERROR: No corresponding files !')
println file_pairs

process bamnanoFiltering {
  tag "Filtering"
  publishDir "$params.outdir/filtering"
  executor 'lsf'
  memory { 10.GB * task.attempt }
  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 3

  input:
    set name, val(nano), val(bam) from file_pairs

  output:
    path "Filtered_nano_bam_files_*" into all_filt
    // path "gcc*_readnames.temp"
    // path "gcc*_greped.txt"

  """
  echo $name
  cut -f3 $bam > "${name}_readnames.temp"
  bash $baseDir/clean_nanopolish.sh $nano $name
  python3 $baseDir/bam_filtering.py $bam ${name}_greped.txt -n $name -s $params.snp_type
  echo $task.memory
  """
}

process gatherFiltered {
  tag "Concat_files"
  publishDir "$params.outdir"
  executor 'lsf'
  memory 1.GB

  input:
    file('*') from all_filt.collect()

  output:
    path "Filtered_nano_bam_files_all.csv"

  """
  grep . * | head -n1 > Filtered_nano_bam_files_all.csv
  tail -n+2 -q Filtered_nano_bam_files_* >> Filtered_nano_bam_files_all.csv
  """
}

process SplitChr{
  tag "SplitChr"
  publishDir "$params.outdir/per_chr"
  executor 'lsf'
  memory { 1.GB * task.attempt }
  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 3

  input:
    "Filtered_nano_bam_files_all.csv"

  output:
    "Filtered_nano_bam_files_chr_*.csv"

  """
  bash split_filter.sh "Filtered_nano_bam_files_all.csv" .
  """
}


// workflow {
//   Channel.fromFilePairs([params.basecall, params.nanopolish], flat:true)
//     .set {file_pairs}
//     .ifEmpty( error 'ERROR: No corresponding files !')
//   bamnanoFiltering(file_pairs) | gatherFiltered(file_all)
//
// }
// cut -f3 $bam > "${name}_readnames.temp"
// zcat $nano | grep -f "${name}_readnames.temp" > "${name}_greped.txt"
// echo "Not recognized read_name number:"
// grep -v -f "${name}_greped.txt" "${name}_readnames.temp" | uniq | wc -l
// cut -f12,13,14,15,16,17,18,19,20 "${name}_greped.txt" > "${name}_extracols1.txt"
// echo "cut good"
// sed "/^[[:space:]]*\$/d" "${name}_extracols.txt" > "${name}_extracols2.txt"
// echo "blou"
// grep -v -f "${name}_extracols2.txt" "${name}_greped.txt" > "${name}_new.txt"
// #mv "${name}_new.txt" "${name}_greped.txt"
// echo "number of row with extra col: "
// wc -l ${name}_extracols2.txt
