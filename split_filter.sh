for chr in `seq 1 24`; do
	blou="^${chr},"
	base=$(basename $1)
	outdir=$2
	echo $outdir$base
	head -n1 $1 > $outdir$base'_chr_'$chr'.csv'
	grep "${blou}" $1 >> $outdir$base'_chr_'$chr'.csv'
	echo $chr 'Success'

done
