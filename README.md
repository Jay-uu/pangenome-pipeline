# pangenome-pipeline
A Nextflow pipeline for Pangenome analysis using [SuperPang](https://github.com/fpusan/SuperPang) and [SqueezeMeta](https://github.com/jtamames/SqueezeMeta) and many other programs.
It allows an automatic way to get data for intra-species diveristy analysis on a large amount of samples.
The pipeline has 3 parts which can be run separately or together.<br>
Step 1: Assembles raw reads into bins using [SqueezeMeta](https://github.com/jtamames/SqueezeMeta) which means that the bins also will be evaluated for quality using Checkm2 and taxonomially classified using GTDB-Tk (and more! Look at the Squeezemeta link for more info about these results and what you can do with them!)<br>
Step 2: Clusters and assembles bins into mOTUs which are then assembled into pangenomes. This allows you to investigate the core and accessory genomes of species.<br>
Step 3: Maps reads towards the previously constructed pangenomes or a reference genome and runs [POGENOM](https://pogenom.readthedocs.io/en/latest/) for population genomic analysis.

# Installation
The pipeline is still in development, so changes may happen.
For now to use the pipeline you need to clone this repository to a location on the machine you wish to run it on:
```git clone https://github.com/Jay-uu/pangenome-pipeline.git```.

Then you need to install and set up the SqueezeMeta dev Conda environment, I recommend that you use Mamba for this. The current version the pipeline works with is SqueezeMeta-dev 1.7.2.
```mamba create -n SqueezeMeta -c conda-forge -c bioconda -c anaconda -c fpusan squeezemeta-dev --no-channel-priority``` and then activate it using ```mamba activate SqueezeMeta```
When this is done, you can continue setting up the environment following the instructions on [SqueezeMeta's GitHub](https://github.com/jtamames/SqueezeMeta?tab=readme-ov-file#3-downloading-or-building-databases) under title 3. Downloading or building databases.

# Running the pipeline
The minimum things you need to run the pipeline is:
1. A directory with raw reads in fastq.gz format.
2. A tab-delimited file with sample names in the first column, the name of a read file in the second column, and "pair1" or "pair2" in the third column, to signify whether the read file has forward/single reads (pair1) or reverse reads (pair2). This is referred to as tsv.samples, or samples file. Each individual Sample (column 1) may either have only paired reads or only single reads.<br>
Example:

    Sample1&nbsp;&nbsp;&nbsp;&nbsp;Sample1.fwd.fq.gz&nbsp;&nbsp;&nbsp;&nbsp;pair1&nbsp;  
    Sample1&nbsp;&nbsp;&nbsp;&nbsp;Sample1.rev.fq.gz&nbsp;&nbsp;&nbsp;&nbsp;pair2&nbsp;  
    Sample2&nbsp;&nbsp;&nbsp;&nbsp;Sample2.R1.fastq.gz&nbsp;&nbsp;&nbsp;&nbsp;pair1&nbsp;  
    Sample2&nbsp;&nbsp;&nbsp;&nbsp;abcdefg.R1.fastq.gz&nbsp;&nbsp;&nbsp;&nbsp;pair1&nbsp; 

3. A project name. This is where your output will be stored. Make this descriptive for your own sake!

## Run command

```nextflow run <path/to/pipeline/main.nf> --project <project_name> --samples <tsv.samples> --fastq <path/to/dir> --threads <nr>```

## Other Options
To get a full list of available commands run:

```nextflow run <path/to/pipeline/main.nf> --help```

or if you're in the pipeline dir:

```nextflow run main.nf --help```

## Alternative entrypoints and exitpoints:

You might not always want to run the whole pipeline depending on your purposes. There are several options to control which parts of the pipeline you run.
 - Stopping the pipeline after assembly and binning:
  If you only want to run assembly from raw reads to bins.

```nextflow run <path/to/pipeline/main.nf> --project <path/project_name> --samples <tsv.samples> --fastq <path/to/dir> --threads <nr> --only_bins true```

 - Existing bins entrypoint:
   If you already have bins that you want to assemble into pangenomes you need to provide a directory with the fasta files for those bins.
   <For this entry you are assumed to have already subsampled your reads, and therefore need to provide a readcount file.> (DOUBLE CHECK)

   ```nextflow run <path/to/pipeline/main.nf> --project <path/project_name> --samples <tsv.samples> --fastq <path/to/dir> --bins <path/to/dir/with/bins/fastas> --readcount <path/to/readcount.txt>```

   For more info see the -extbins flag for [SqueezeMeta](https://github.com/jtamames/SqueezeMeta?tab=readme-ov-file#5-execution-restart-and-running-scripts)

 - Reference genome entrypoint:
   If you have a pangenome/reference genome that you want to do some diversity analysis with you provide a fasta file.

   ```nextflow run <path/to/pipeline/main.nf> --project <path/project_name> --samples <tsv.samples> --fastq <path/to/dir> --ref_genome <path/to/file>``` it is also recommended that you provide a file with the contigs of interest using ```--contigs <file>```. The contigs file has the name of a contig on each line, make sure it matches the names in the provided fasta file.

 - Skipping variant calling:
   This can be used for the regular run command, or with the existing bins entrypoint.

   ```nextflow run <path/to/pipeline/main.nf> --project <path/project_name> --samples <tsv.samples> --fastq <path/to/dir> --run_VC false```

## Subsampling
By default the pipeline uses subsampling of the raw reads to map them to the pangenomes/reference genomes to get an estimate of expected coverage for a sample to a genome. This saves on computation time by not needing to map all reads to all genomes to determine which samples can contribute to the species diversity. Unless you have very few samples this is the recommended way to run the pipeline.

- Skipping subsampling:
If you despite this want to skip subsampling you can add ```--subsample false``` to the run command. This option is not allowed when starting from Entrypoint 1, using just raw reads, but works for Entrypoint 2 and 3.

- Skipping subsampling with pre-existing bins:
If you've previously run the pipeline and stopped after the binning step you can start from the existing bins entrypoint and skip subsampling. This presumes that you already have subsampled reads and therefore you need to provide a file with the total number of original reads in each sample using ```--readcount```. The pipeline provides you with the files to make this easy.

```nextflow run <path/to/pipeline/main.nf> --project <path/new_project_name> --samples <path/old_project_name/subsamples/old_project_name.subsampled.samples> --fastq <path/old_project_name/subsamples/fastqs> --subsample false --readcount <path/old_project_name/subsamples/original_readcounts.tsv --run_VC false```

# Configurations
Instead of writing all your parameters on the command line you can write a parameter file and provide it using ```-params-file <parameter_file.yaml>```.
The file needs to be either in YAML or JSON format.
YAML example:<br>

    project: 'My_project_date'
    samples: 'path/to/my_project.samples'
    fastq: 'path/to/fastqs'
    threads: 12
    taxSort: 'class'
    nr_samps_threshold: 2
    min_cov: 10
    min_breadth: 50

JSON example:<br>
```
{
"project": "My_project_date",
"samples": "path/to/my_project.samples",
"fastq": "path/to/fastqs",
"threads": 12
}
```
    
NB! There is a [bug](https://github.com/nextflow-io/nextflow/issues/2662) related to defining parameters within a config file, meaning that it is much safer to provide parameters on the command line or in a separate parameter file.
 
The same can be done with specific configurations. If you're running the pipeline on a cluster you might want to optimize the use of resources. There are general labels for this that you can use, or use the name of the individual process to determine how many resources you request for it and it is allowed to use.
The command for using a config file is `-c <config-file> `. Nf-core has some for different clusters, but here's an example on how I ran it on KTH's cluster Dardel:  
```
//General settings for the pipeline execution:
executor.queueSize = 100

process {
    //These settings will apply to all processes (except the ones with other configs using withLabel).
    executor = 'slurm'
        scratch = true //each process execution uses the scratch dir so temporary files won't take up disk space. (CHECK THIS)
        clusterOptions = { '-A <project_name>' }
        cpus = params.threads
        queue = 'shared' //This is the partition directive
        time = '1d'
        memory = '10 GB' //total mem for each task
        
    //This process is the most resource intensive, so I booked a full node for it.
    withLabel: 'fastq_to_bins' {
        time = { 24.h }
        queue = 'shared'
        memory = '111 GB'
        cpus = 24
    }
    
    withLabel: low_cpu {
        executor = 'slurm'
        clusterOptions = { '-A <project>' }
        cpus = { 1 }
        queue = 'shared'
        time = '1d 10h 30m'
        scratch = true
        memory = { 4.GB * task.attempt }
        time = { 10.h * task.attempt }
        errorStrategy = { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
        maxRetries = 2
    }
   }
   ```

The more general labels you can use to configure the pipeline are low_cpu (these processes use one cpu effectively), high_mem and the individual process names.
If you want to know more about how Nextflow uses configurations you can [read the docs](https://www.nextflow.io/docs/latest/config.html).

# Results structure
After running the pipeline you might wonder where your results are and what they mean. My best recommendation is to explore, but here's a tree showing the structure of the results output from a full run of the pipeline. Words within brackets [] are variable, depending on your input, sample names, taxonomic identification, tool etc. The .zip file in <Project/results/[taxonomic_group_Y]_mOTU_[N]/pangenome/[taxonomic_group_Y]_mOTU_[N].zip> can be explored using [SQMtools in R](https://github.com/jtamames/SqueezeMeta/wiki/Using-R-to-analyze-your-SQM-results).

[Project]/<br>
├─ bins/<br>
│  ├─ bintables/<br>
│  │  ├─ [SampX].bintable<br>
│  ├─ fastas/<br>
│  │  ├─ [SampX].[binner].[binNr].[tool_suffix].contigs.fa<br>
│  ├─ summarized_bintable.tsv<br>
├─ mOTUs/<br>
│  ├─ mOTUlizer/<br>
│  │  ├─ [taxonomic_group_Y]/<br>
│  │  │  ├─ [tax_group_Y]_mOTUs.tsv<br>
│  │  c  ├─ [tax_group_Y]_similarities.txt<br>
│  ├─ results/<br>
│  │  ├─ [taxonomic_group_Y]_mOTU_[N]/<br>
│  │  │  ├─ pang_bins.txt<br>
│  │  │  ├─ pangenome/<br>
│  │  │  │  ├─ pogenom/<br>
|  |  |  |  |  ├─ results/<br>
|  |  |  |  |  |  ├─ lots of files, will update later. FST, allele frequencies, pNpS etc. <br>
|  |  |  |  |  |─ [taxonomic_group_Y]_mOTU_[N]_unfiltered.vcf<br>
│  │  │  │  ├─ superpang/<br>
|  |  |  |  |  ├─ params.tsv<br>
|  |  |  |  |  ├─ Superpang output. Lots of files. Update later. Find core and accessory fasta files here.<br>
│  │  │  │  ├─ cov_breadth.txt<br>
│  │  │  │  ├─ [taxonomic_group_Y]_mOTU_[N]_cM2_summary.txt<br>
│  │  │  │  ├─ [taxonomic_group_Y]_mOTU_[N].samples<br>
│  │  │  │  ├─ [taxonomic_group_Y]_mOTU_[N].cov<br>
│  │  │  │  ├─ [taxonomic_group_Y]_mOTU_[N].zip/<br>
│  │  │  │  ├─ [taxonomic_group_Y]_mOTU_[N].cpm<br>
│  ├─ [Project].cov.tsv<br>
│  ├─ [Project].cpm.tsv<br>
├─ subsamples/<br>
│  ├─ fastqs/<br>
│  │  ├─ sub_[SampX]_R1.fq.gz<br>
│  │  ├─ sub_[SampX]_R2.fq.gz<br>
│  ├─ [Project].subsampled.samples<br>
│  ├─ original_readcounts.tsv<br>
├─ README.md<br>

