import os
import gzip

file = "finemapped_readnames_chr.txt"
outputfile = "finemapped_readchr_count.txt"

chr_count = {}
with open(file,'rb') as f:
		for line in f:
			values = line.decode("utf-8").rstrip().split("\t")
			if str(values[0:2]) not in chr_count:
				chr_count[str(values[0:2])] = 1
			else:
				chr_count[str(values[0:2])] = chr_count[str(values[0:2])]+1


with open(outputfile, "w") as out:
	for c in chr_count:
		out.write(str(c) + "\t" + str(chr_count[c]) + "\n")

out.close()


