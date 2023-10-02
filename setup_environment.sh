#!/bin/bash -l
echo "=====Creating squeezemeta environment====="
#mamba create -n squeezemeta -c conda-forge -c bioconda -c anaconda -c fpusan squeezemeta-dev --no-channel-priority
#mamba activate squeezemeta

echo "=====Configuring database====="
#configure_nodb.pl /home/fer/Data/ssd/SQM/db
configure_nodb.pl /proj/fume/databases/SqueezeMeta/220723/db

echo "=====Installing nextflow====="
mamba install -c bioconda nextflow 
mamba update nextflow
echo "nextflow version is:"
nextflow -v

echo "=====Instaling nf-core tools====="
mamba install nf-core
mamba update nf-core
echo "nf-core version:"
nf-core --version

echo "=====Installing seqtk====="
#git clone https://github.com/lh3/seqtk.git;
#cd seqtk; make
mamba install -c bioconda seqtk

echo "=====Done====="
echo "nf-validation should be automatically installed after running nextflow"

#echo "=====Testing nextflow====="
#nextflow run test_params.nf --project sq_dir --samples /data/fer/SQMtestData/freshmock/mock.samples --fastq /data/fer/SQMtestData/freshmock/raw
