import os

# This filters the gziped fastq files for only those which have primer bindig flanking sequences.
for i in range(1,19):
	os.system("zcat ~/Desktop/amp130-148/amp{}_S{}_L001_R2_001.fastq.gz | seqkit grep -s -i -p GTGGCTAA | seqkit grep -s -i -p CGGAAGAG > ~/Desktop/filtered/amp{}_S{}_L001_R2_001_filtered.fastq.gz".format(130 + i, 1 + i, 130 + i, 1 + i))
  
  
  
  
  
  
