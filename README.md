# pangenome-pipeline
A Nextflow pipeline for Pangenome analysis using [SuperPang](https://github.com/fpusan/SuperPang) and [SqueezeMeta](https://github.com/jtamames/SqueezeMeta) and many other programs.
Need a description here about what the pipeline does specifically.

# Installation
The pipeline is still in development, so changes may happen.
For now to use the pipeline you need to clone this repository to a location on the machine you wish to run it on:
```git clone https://github.com/Jay-uu/pangenome-pipeline.git```.

Then you need to install and set up the SqueezeMeta dev Conda environment, I recommend that you use Mamba for this:
```mamba create -n SqueezeMeta -c conda-forge -c bioconda -c anaconda -c fpusan squeezemeta-dev --no-channel-priority``` and then activate it using ```mamba activate SqueezeMeta```
When this is done, you can continue setting up the environment following the instructions on [SqueezeMeta's GitHub](https://github.com/jtamames/SqueezeMeta?tab=readme-ov-file#3-downloading-or-building-databases) under title 3. Downloading or building databases.

# Running the pipeline
The minimum things you need to run the pipeline is:
1. A directory with raw reads in fastq.gz format.
2. A tab-delimited file with sample names in the first column, the name of a read file in the second column, and "pair1" or "pair2" in the third column, to signify whether the read file has forward/single reads (pair1) or reverse reads (pair2). This is referred to as tsv.samples, or samples file. Each individual Sample (column 1) may either have only paired reads or only single reads.
   Example:

 | ----------- | ----------- | ----------- |
 | Sample1 | Sample1.fwd.fq.gz | pair1 |
 | Sample1 | Sample1.rev.fq.gz | pair2 |
 | Sample2 | Sample2.R1.fastq.gz | pair1 |
 | Sample2 | abcdefg.R1.fastq.gz | pair1 |

3. A project name. This is where your output will be stored. Make this descriptive for your own sake!

## Run command

```nextflow run <path/to/modular/main.nf> --project <project_name> --samples <tsv.samples> --fastq <path/to/dir> --threads <nr>```

## Other Options
To get a full list of available commands run:

```nextflow run <path/to/modular/main.nf> --help```

or if you're in the modular dir:

```nextflow run main.nf --help```

## Alternative entrypoints and exitpoints:

You might not always want to run the whole pipeline depending on your purposes. There are several options to control which parts of the pipeline you run.
 - Stopping the pipeline after assembly and binning:
  If you only want to run assembly from raw reads to bins.

```nextflow run <path/to/modular/main.nf> --project <path/project_name> --samples <tsv.samples> --fastq <path/to/dir> --threads <nr> -entry raw_to_bins```

 - Existing bins entrypoint:
   If you already have bins that you want to assemble into pangenomes you need to provide a directory with the fasta files for those bins.

   ```nextflow run <path/to/modular/main.nf> --project <path/project_name> --samples <tsv.samples> --fastq <path/to/dir> --bins <path/to/dir/with/bins/fastas>```

   For more info see the -extbins flag for [SqueezeMeta](https://github.com/jtamames/SqueezeMeta?tab=readme-ov-file#5-execution-restart-and-running-scripts)

 - Reference genome entrypoint:
   If you have a pangenome/reference genome that you want to do some diversity analysis with provide a fasta file.

   ```nextflow run <path/to/modular/main.nf> --project <path/project_name> --samples <tsv.samples> --fastq <path/to/dir> --ref_genome <path/to/file>``` it is also recommended that you provide a file with the contigs of interest using ```--contigs <file>```. The contigs file has the name of a contig on each line, make sure it matches the names in the provided fasta file.

 - Skipping variant calling:
   This can be used for the regular run command, or with the existing bins entrypoint.

   ```nextflow run <path/to/modular/main.nf> --project <path/project_name> --samples <tsv.samples> --fastq <path/to/dir> --run_VC false```

## Subsampling


# Configurations

# Results structure
