#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// General parameters
params.nanopolish = "/hps/nobackup/birney/projects/gel_methylation/nanopolish/gcc*.tsv.gz"

// Parameters list of snp specific
params.title = "covid_snp_all_updated_Summer_2022"
params.basecall = "$baseDir/base_called_from_bam_files/covid_snps_updated_summer2022/gcc*.txt"
params.snp_type = "covid_snp"
params.outdir = "$baseDir/$params.title"
params.list_regions = "List_region.txt"

// Parameters step specific
params.jump_filtering = "True"
// TODO: Make work the split per chr (currently going directly to analysis)
params.splitCHR = "False"
params.medianfile_basename = "$baseDir/covid_snp_all_updated_Summer_2022/All_filtered_nanobamfiles.csv"

workflow {
          // Recuperation and filtering of the bam/nanopolish files
          file_pairs = Channel.fromFilePairs([params.basecall, params.nanopolish], flat:true)
          (file_filt, out1) = bamnanoFiltering(file_pairs)
          file_filt = file_filt.collect()
          gatherFile(file_filt)

          // Adjustment of the regions on which to perform stat analysis (more jobs less time)
          // chr_files = SplitPerCHR(params.medianfile_basename)
          // TODO: Find a way to create the channel with precedent files (with splitchr files output) #important
          create_channel = {params.medianfile_basename.replaceAll(/\..+/, "") + '_chr_' + it.replaceFirst(/snp-/, "").replaceAll(/:.+/, "") + '.csv'}
          selected_region = file(params.list_regions)
          region_list = selected_region.readLines()
          // println chr_files


          // Statistical analysis
          (region_files, out2) = Analysis_cpg_snp( Channel.from(region_list.collect(create_channel)), Channel.from(region_list))
          region_files = region_files.collect()
          GatherStats(region_files)

          // Recuperation of output on the various pipeline works
          // SaveLog(out1.collect(), out2.collect())
          }

process bamnanoFiltering {
  tag "Filtering $name"
  publishDir "$params.outdir/filtering", mode: 'copy', overwrite: false
  executor 'lsf'
  memory { 100.GB * task.attempt }
  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 3
  when:
    params.jump_filtering == "False"

  input:
    tuple val(name), val(nano), val(bam)

  output:
    file("Filtered_nano_bam_files_*") optional true
    stdout()

  """
  echo $name
  cut -f3 $bam > "${name}_readnames.temp"
  bash $baseDir/clean_nanopolish.sh $nano $name
  python3 $baseDir/bam_filtering.py $bam ${name}_greped.txt -n $name -s $params.snp_type
  echo $task.memory
  """
}

process gatherFile{
  tag "Concatenating files"
  publishDir "$params.outdir/per_chr", mode: 'copy', overwrite: false
  executor 'lsf'
  memory 10.GB
  // when:
  //   params.jump_filtering == "False"

  input:
    file("Filtered_nano_bam_files_*")

  output:
    file("$params.medianfile_basename")

  """
  grep . -h Filtered_nano_bam_files_* | head -n1 > $params.medianfile_basename
  tail -n+2 -q * >> $params.medianfile_basename
  bash $baseDir/split_filter.sh $params.medianfile_basename ./
  """
  // (echo .separator ,; echo .import All_files_temp.csv $table_name) | sqlite3 $params.db_name
}

process SplitPerCHR{
  tag "Split per chr"
  publishDir "$params.outdir", mode: 'copy', overwrite: false
  executor 'lsf'
  memory { 5.GB }
  when:
    params.splitCHR == "True"

  input:
    file("$params.medianfile_basename")

  output:
    file("*.csv")

  """
  bash $baseDir/split_filter.sh $params.medianfile_basename ./
  """
}

process Analysis_cpg_snp{
  tag "Analysis of cpg/snp association $file_select $select"
  publishDir "$params.outdir/stats_nodb", mode: 'copy', overwrite: false
  executor 'lsf'
  memory { 200.GB * task.attempt }
  errorStrategy { task.exitStatus in [130] ? 'retry' : 'terminate' } // ,255
  maxRetries 3

  input:
    val(file_select)
    val(select)

  output:
    file("Results_stats_select_*.csv") optional true
    stdout()

  // readarray -d: -t splitNoIFS<<<"$select"
  // chr="\${splitNoIFS[0]:4:3}"
  """
  python3 $baseDir/cpg_snp_analysis.py $file_select -s $select
  """
}

process GatherStats {
  tag "Concatenating files"
  publishDir "$params.outdir", mode: 'copy', overwrite: false
  executor 'lsf'
  memory 1.GB

  input:
    file("*")

  output:
    file("Results_stats_all.csv")

  """
  grep . -h * | head -n1 > Results_stats_all.csv
  tail -n+2 -q * >> Results_stats_all.csv
  """
  // (echo .separator ,; echo .import Results_stats_all.csv $table_name) | sqlite3 $params.db_name
}

// TODO Add a split per statistique type ??
// TODO Add ploting script ?
// TODO Add copy of the config - region ls etc into the dir ?


process SaveLog {
  tag "Saving Logs"
  publishDir "$params.outdir"
  memory 1.GB

  input:
    val(res1)
    val(res2)

  output:
    file("Log_file.txt") optional true

  exec:
    log_file = file("$params.outdir/Log_file.txt")
    log_file.text = res1
    log_file.append(res2)
}
