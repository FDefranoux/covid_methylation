
filenames_file="file_list3.txt"
outputdir=output2

file_str="${LSB_JOBINDEX}p"
fileselect="sed -n $file_str ${filenames_file}"
filename="$(${fileselect})"
fbname=$(basename -a $filename)

outfile=${outputdir}/${fbname}

echo ${filename} '-->' ${outfile}

gunzip -c ${filename} | grep -v MT | grep -v GL | grep -v KI | sed s'/X/23/' | sed s'/Y/24/' | sort -n -k 1 -k 3 -k 4 | /hps/software/users/birney/tomas/tabix/tabix-0.2.6/bgzip -c > ${outfile}.temp

echo 'sorting done'

mv ${outfile}.temp ${outfile}

/hps/software/users/birney/tomas/tabix/tabix-0.2.6/tabix -s 1 -b 3 -e 4 ${outfile}

echo 'tabix done'
