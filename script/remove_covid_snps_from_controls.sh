N=${LSB_JOBINDEX}'q;d'
file=$(sed ${N} ${1})
echo $file
new_file=$(basename $file)
grep -f nocovid_in_control_snps $file > 'control_snps_finemapped_clean/'${new_file}
