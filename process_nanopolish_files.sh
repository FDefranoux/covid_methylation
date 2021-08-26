outputdir=output
filename=$1
sed=$2

fbname=$(basename $filename)
outfile='new_'${fbname}

echo ${filename} > /dev/stderr
echo ${filename} '-->' ${outfile}

gunzip -c ${filename} | grep -v MT | grep -v GL | grep -v KI | sed s'/X/23/' | sed s'/Y/24/' | sed $sed | sort -n -k 1 -k 3 -k 4 -T /nfs/research/birney/users/fanny/Covid/temp | /hps/software/users/birney/tomas/tabix/tabix-0.2.6/bgzip -c > ${outfile}.temp

echo 'sorting done'

mv ${outfile}.temp ${outfile}

/hps/software/users/birney/tomas/tabix/tabix-0.2.6/tabix -s 1 -b 3 -e 4 ${outfile}

echo 'tabix done'
