for chr in `seq 1 24`; do
	blou="^${chr},"
	echo $chr
	head -n1 $1 > 'temp/'$1'_chr_'$chr'.csv'
	grep "${blou}" $1 >> 'temp/'$1'_chr_'$chr'.csv'
done
