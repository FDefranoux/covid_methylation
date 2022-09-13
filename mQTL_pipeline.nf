#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
params.title = "covid_snp_all_updated_Summer_2022"
// params.basecall = "$baseDir/base_called_from_bam_files/covid_snps_finemapped/gcc008*.txt"
params.basecall = "$baseDir/base_called_from_bam_files/covid_snps_updated_summer2022/gcc*.txt"
params.nanopolish= "/hps/nobackup/birney/projects/gel_methylation/nanopolish/gcc*.tsv.gz"
params.snp_type = "covid_snp"
params.db_name = "$params.title" + ".db"
params.outdir = "$baseDir/$params.title"


workflow {
        file_pairs = Channel.fromFilePairs([params.basecall, params.nanopolish], flat:true)
        (file_filt, res1) = bamnanoFiltering(file_pairs)
        file_filt = file_filt.collect()
        db = gatherFiletoDB(file_filt, "median_datas")
        // region_ch = Channel.of(1..24)
        region_list = file('List_region.txt').readLines()
        region_ch = Channel.from(region_list)
        (region_files, res2) = Analysis_cpg_snp(db, region_ch)
        region_files = region_files.collect()
        UpdateDB(region_files, "stat_datas")
        SaveLog(res1.collect(), res2.collect())
}

process bamnanoFiltering {
  tag "Filtering $name"
  publishDir "$params.outdir/filtering", mode: 'copy', overwrite: false
  executor 'lsf'
  memory { 100.GB * task.attempt }
  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 3

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

process gatherFiletoDB {
  tag "Concatenating files"
  publishDir "$params.outdir", mode: 'copy', overwrite: false
  executor 'lsf'
  memory 10.GB

  input:
    file("Filtered_nano_bam_files_*")
    val(table_name)

  output:
    file("*.db")

  """
  grep . -h Filtered_nano_bam_files_* | head -n1 > All_files_temp.csv
  tail -n+2 -q * >> All_files_temp.csv
  (echo .separator ,; echo .import All_files_temp.csv $table_name) | sqlite3 $params.db_name
  """
}

process Analysis_cpg_snp{
  tag "Analysis of cpg/snp association $select"
  publishDir "$params.outdir/stats", mode: 'copy', overwrite: false
  executor 'lsf'
  memory { 200.GB * task.attempt }
  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 3

  input:
    file("$params.db_name")
    val(select)

  output:
    file("Results_stats_*.csv") optional true
    stdout()

  """
  python3 $baseDir/cpg_snp_analysis.py $params.outdir/$params.db_name -s $select
  """

}

process UpdateDB {
  tag "Concatenating files"
  publishDir "$params.outdir", mode: 'copy', overwrite: false
  executor 'lsf'
  memory 1.GB

  input:
    file("Results_stats_*.csv")
    val(table_name)

  output:
    file("$params.db_name") optional true

  """
  grep . -h Results_stats_* | head -n1 > All_files_temp.csv
  tail -n+2 -q * >> All_files_temp.csv
  (echo .separator ,; echo .import All_files_temp.csv $table_name) | sqlite3 $params.db_name
  """
}

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
