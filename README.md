# pangenome-pipeline
A Nextflow pipeline for Pangenome analysis using [SuperPang](https://github.com/fpusan/SuperPang) and [SqueezeMeta](https://github.com/jtamames/SqueezeMeta) and many other programs.

# Installation
The pipeline is still in development, so changes may happen.
For now to use the pipeline you need to clone this repository to a location on the machine you wish to run it on:
`git clone https://github.com/Jay-uu/pangenome-pipeline.git`

Then you need to install and set up the SqueezeMeta dev Conda environment, I recomment that you use Mamba for this:
`mamba create -n SqueezeMeta -c conda-forge -c bioconda -c anaconda -c fpusan squeezemeta-dev --no-channel-priority`
When this is done, you can continue setting up the environment following the instructions on [SqueezeMeta's GitHub](https://github.com/jtamames/SqueezeMeta?tab=readme-ov-file#3-downloading-or-building-databases)

# Running the pipeline
The minimum things you need to run the pipeline is:
1. A directory with raw reads in fastq.gz format.
2. A tab-delimited file with sample names in the first column, the name of a read file in the second column, and "pair1" or "pair2" in the third column, to signify whether the read file has forward/single reads (pair1) or reverse reads (pair2).
3. A project name. This is where your output will be stored. Make this descriptive for your own sake!

# Simplest way of running the pipeline
`nextflow run <path/to/modular/main.nf> --project <project_name> --samples <tsv.samples> --fastq <path/to/dir> --threads <nr>`

# Other Options
To get a full list of available commands run:
`nextflow run <path/to/modular/main.nf> --help` 
