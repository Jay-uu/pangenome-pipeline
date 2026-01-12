#!/bin/bash
PROJNR="" #SRA project accession
mkdir prefetch_reads
mkdir raw_reads
#Create file with all accesion nr associated with a specific project
esearch -db sra -query $PROJNR | efetch -format runinfo | cut -d "," -f 1 > SRR.numbers

#Loop through file to download each accession to a gzipped fastq file sequentially.
cat SRR.numbers | while read line || [[ -n $line ]];
do
  echo "======Pre-fetching $line======"
  prefetch $line -O prefetch_reads/
  echo "======Validating $line======"
  vdb-validate prefetch_reads/$line
  echo "======Retrieving $line======"
  fastq-dump prefetch_reads/$line -v --split-3 --gzip --outdir raw_reads
done
 echo "Finished!"
