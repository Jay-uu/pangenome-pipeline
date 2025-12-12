#!/bin/bash -l
echo "=====Creating squeezemeta environment====="
#mamba create -n squeezemeta -c conda-forge -c bioconda -c anaconda -c fpusan squeezemeta-dev --no-channel-priority
#mamba activate squeezemeta
#conda create -n sqm_alpha4 -c conda-forge -c bioconda -c fpusan --no-channel-priority --override-channels squeezemeta-dev=1.7.3.alpha4

echo "=====Configuring database====="
#configure_nodb.pl /home/fer/Data/ssd/SQM/db #Cass
#/home/fer/Data/ssd/gtdb_release220 #gtdb is symlinked in the above
configure_nodb.pl /proj/fume/databases/SqueezeMeta/230903/db/ #Uppmax
#configure_nodb.pl /cfs/klemming/projects/supr/fume/databases/SqueezeMeta/230903/db #Dardel'
#configure_nodb.pl /scratch/slu_dasa_1/slu_dasa_1_2/dbs/SQM #SCAYLE
#configure_nodb.pl /data/jay/databases/SQM/db



#these steps are now unnecessary

echo "=====Configuring Checkm2 database====="
export CHECKM2DB="/data/jay/databases/CheckM2_database/CheckM2_database" #Cass
export CHECKM2DB="/cfs/klemming/projects/supr/fume/databases/checkm2/CheckM2_database" #Dardel
export CHECKM2DB="/crex/proj/fume/nobackup/private/jay/dbs/CheckM2_database" #Uppmax


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
