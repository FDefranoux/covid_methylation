import os
import gzip

file_list = "file_list_all.txt"
with open(file_list) as file:
	nanos = file.readlines()

file.close()
index = int(os.environ['LSB_JOBINDEX'])-1
nano = nanos[index].rstrip()

outputfile = str(index) + ".txt"

chr_count = {}
with open(nano,'rb') as f:
		for line in f:
			values = line.decode("utf-8").rstrip().split("\t")
			if values[0] not in chr_count:
				chr_count[values[0]] = 1
			else:
				chr_count[values[0]] = chr_count[values[0]]+1


with open(outputfile, "w") as out:
	for c in chr_count:
		out.write(nano +"\t" + str(c) + "\t" + str(chr_count[c]) + "\n")

out.close()


