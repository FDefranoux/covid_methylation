for chr in `seq 1 24`; do
	blou="^${chr},"
	echo $chr
	grep "${blou}" $1 > $1$chr'.csv'
done

