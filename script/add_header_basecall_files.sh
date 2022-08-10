for i in $( ls /hps/nobackup/birney/projects/gel_methylation/covid_snps_current/reads/*.txt ); do sed '1 
s/^/sample_id\tcovid_snp\tread_name\tbase_called\n/' $i > base_called_from_bam_files/covid_snps_updated_summer2022/$(basename $i); done
