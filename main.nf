#!/usr/bin/env nextflow
/*
========================================================================================
Pipeline for running SqueezeMeta from start
========================================================================================
With only raw reads and samples file as input
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from "plugin/nf-validation"

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run main.nf --project <project_name> --samples <tsv.samples> --fastq <path/to/dir> --threads <nr> \n\n  To get more info about a specific parameter write nextflow run main.nf --help <parameter_name>")
   log.info """\
   You can supply a config file with the parameters using -c. For more info see the usr_template.config or
   ${params.manifest.homePage}
   """

   exit 0
}

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

//Check project parameter
def badChars = ["^","(",")","+", " ", "|"]
if ( params.project.findAll { a -> badChars.any { a.contains(it) } } ) {
	throw new Exception("Invalid project name. Special characters and whitespaces not allowed.")
}

if (workflow.resume == false) {
	//Workflow was not resumed, checking project dir
	Path projDir = new File(params.project).toPath()
	if (projDir.exists() == true) {
		throw new Exception("Project directory $params.project already exists, choose a new name or use the -resume flag. WARNING: Note that if you resume the wrong job, this might overwrite previous results.")
	}
}

workflow {

sam_chan = Channel.fromPath(params.samples, type: "file", checkIfExists: true)
/*The fastq_dir is needed for:
	- Formating the individual sample files
	- Running SqueezeMeta for each sample.
	- Subsampling
This creates three channels for the three different processes.
*/
Channel.fromPath(params.fastq, type: "dir", checkIfExists: true)
	    .multiMap { dir -> format: to_bins: subsample: dir }.set { fastq_chan }

/*Runs the process that creates individual samples files and creates two output channels:
	- For Squeezemeta (fastq_to_bins process)
	- For the subsampling
*/
format_samples(sam_chan, fastq_chan.format)
format_samples.out.flatten().multiMap { samp -> to_bins: subsample: samp }.set { single_samps }

/* Runs the fastq_to_bins process.
	flatten() to handle the created singles samples form format_samples individually
	first() to supply the fastq dir for each sample
*/
fastq_to_bins(single_samps.to_bins, fastq_chan.to_bins.first())

//Concatenating fastqs and subsampling for later mapping for each singles sample
subsample_fastqs(single_samps.subsample, fastq_chan.subsample.first())

/*
To create mOTUs, all the bins need to be gathered into one process run, and all the checkm files need to be gathered
into one and provided to the same process.
*/
//all_checkm = fastq_to_bins.out.checkm.collectFile(name: "all.checkM", newLine: true) //concatenating checkm-files
//all_checkm.multiMap { file -> to_mOTUlizer: to_mOTU_dirs: file }.set { all_checkm }
fastq_to_bins.out.checkm.collectFile(name: "all.checkM", newLine: true)
	.multiMap { file -> to_mOTUlizer: to_SuperPang: file }.set { all_checkm }

//all_bins = fastq_to_bins.out.bins.collect()
fastq_to_bins.out.bins.collect().multiMap { bins -> to_mOTUlizer: to_mOTU_dirs: bins }.set { all_bins }
bins_to_mOTUs(all_bins.to_mOTUlizer, all_checkm.to_mOTUlizer)

/*
Creating dirs for the mOTUs
*/
create_mOTU_dirs(bins_to_mOTUs.out.mOTUs_file, all_bins.to_mOTU_dirs)

/*
Running SuperPang
*/
mOTUs_to_pangenome(create_mOTU_dirs.out.flatten(), all_checkm.to_SuperPang.first())

//=======END OF WORKFLOW=============
}

process format_samples {
    label "short_time"
    cache "deep"
    input:
    path(samples_file)
    path(fastq_dir)
    output:
    path("singles/*")
    shell:
    $/
    #!/usr/bin/env python
    import os
    from collections import defaultdict
    from collections import Counter
    
    os.makedirs("singles", exist_ok=True)
    fastq_files = os.listdir("!{fastq_dir}")
    
    samples2lines = defaultdict(list)
    with open("!{samples_file}") as infile:
        for line in open("!{samples_file}"):
            fields = line.strip().split("\t")
            if len (fields) < 3:
                raise Exception(f"Missing columns or wrong delimiter on line: {line}") # this halts execution with exit code -1 and emits a message
            sample, filename, pair, *_ = fields
            if pair not in ("pair1", "pair2"):
                raise Exception("Error with provided samples file. Make sure there's no header and that the third column says 'pair1' or 'pair2'") #a header, or spelled wrong
            if filename not in fastq_files:
                 raise Exception(f"{filename} not found in !{fastq_dir}")
            if not filename.endswith(".fastq.gz") and not filename.endswith(".fq.gz"):
                 raise Exception(f"Wrong file format for {filename}. Only .fastq.gz and .fq.gz are accepted.")
            samples2lines[sample].append(line)
            
    paired_flag = set()
    for sample, lines in samples2lines.items():
        pairs = []
        for row in lines:
            pairs.append(row.strip().split("\t")[2])
        count = Counter(pairs)
        if count["pair1"] == count["pair2"]:
            paired_flag.add("paired")
        elif "pair2" not in count:
            paired_flag.add("unpaired")
        else:
            raise Exception(f"Sample {sample} has a paired file mismatch in the samples file.\n Make sure that there are only single reads or only paired end reads.")
        with open(f"singles/{sample}.samples", "w") as outfile:
            for line in lines:
                outfile.write(line)
    /$       
}

process fastq_to_bins {
    cache "deep"
    publishDir "${params.project}/sqm_res", mode: "copy", pattern: "${sample.baseName}"
    input:
    path(sample)
    path(fastq_dir)
    output:
    path("${sample.baseName}", emit: sample_dir)
    path("${sample.baseName}/results/bins/*.contigs.fa", emit: bins)
    path("${sample.baseName}/intermediate/17.*.checkM", emit: checkm)
    shell:
    '''
    echo "The sample file is !{sample.baseName} and the fastq dir is !{fastq_dir}"
    SAMPLE_ID="!{sample.baseName}"
    SqueezeMeta.pl -m coassembly -f !{fastq_dir} -s !{sample} -p $SAMPLE_ID -binners maxbin,metabat2,concoct -t !{params.threads} -test 1 -contigid $SAMPLE_ID
    10.mapsamples.pl $SAMPLE_ID
    14.runbinning.pl $SAMPLE_ID
    15.dastool.pl $SAMPLE_ID
    17.checkM_batch.pl $SAMPLE_ID
    '''

}

process subsample_fastqs {
    cache "deep"
    label "medium_time"
    input:
    path(sample)
    path(fastq_dir)
    output:
    tuple val("${sample.baseName}"), path("sub_*.fq.gz")
    shell:
    $/
    #!/usr/bin/env python
    import os
    from pathlib import Path
    from subprocess import call
    import gzip
    import shutil
    
    """
    A function which takes a list of fastq files and subsamples a million read from the combined files.
    Result is a compressed fq file.
    """
    def concat_subtk_compress(file_list, direction):
        #concatenating
        with open(f"concat_{direction}.fq.gz", "wb") as outfile:
            for f in file_list:
                outfile.write(open("!{fastq_dir}"+"/"+f,"rb").read())
        #Subsampling
        with open(f"sub_{sample_ID}_{direction}.fq", "w") as subout:
            call(["seqtk", "sample", "-s100", f"concat_{direction}.fq.gz", "1000000"], stdout = subout)
        #Compressing
        with open(f"sub_{sample_ID}_{direction}.fq", "rb") as in_f, gzip.open(f"sub_{sample_ID}_{direction}.fq.gz", "wb") as out_f:
            shutil.copyfileobj(in_f, out_f)
        #Removing intermediate files
        os.remove(f"concat_{direction}.fq.gz")
        os.remove(f"sub_{sample_ID}_{direction}.fq")

    fastq_files = os.listdir("!{fastq_dir}")
    sample_ID = Path("!{sample}").stem
    
    fwds = []
    revs = []

    with open("!{sample}") as infile:
        for line in open("!{sample}"):
            fields = line.strip().split("\t")
            if len(fields) < 3:
                raise Exception(f"Missing columns or wrong delimiter on line: {line}")
            sample, filename, pair, *_ = fields
            if filename not in fastq_files:
                raise Exception(f"{filename} not found in !{fastq_dir}")
            if pair == "pair1":
                fwds.append(filename)
            elif pair == "pair2":
                revs.append(filename)
            else:
                raise Exception("Error with provided samples file. Make sure there's no header and that the third column says 'pair1' or 'pair2'")

    fwds = sorted(fwds)
    revs = sorted(revs)

    concat_subtk_compress(fwds, "R1")
    if len(revs) > 0:
        concat_subtk_compress(revs, "R2")
    
    /$
}

process bins_to_mOTUs {
    //the conda part might be changed later. If for example SuperPang gets updated to run with the newest version
    //of mOTUlizer this process doesn't need a separate environment
    conda 'bioconda::mOTUlizer=0.3.2'
    publishDir "${params.project}/", mode: "copy"
    input:
    path(bins)
    path(all_checkm)
    output:
    path("mOTUs.tsv", emit: mOTUs_file)
    path("MAG_similarities.txt", emit: simi_file)
    shell:
    """
    mOTUlize.py --fnas *.fa --checkm !{all_checkm} --MAG-completeness !{params.MAGcomplete} --MAG-contamination !{params.MAGcontam} --threads !{params.threads} --keep-simi-file MAG_similarities.txt -o mOTUs.tsv
    """
}

process create_mOTU_dirs {
    label "short_time"
    publishDir "${params.project}/mOTUs", mode: "copy"
    input:
    path(motus_file)
    path(bins)
    output:
    path("mOTU_*", type: "dir")
    shell:
    """
    #!/usr/bin/env python
    import os
    import shutil
    min_genomes = !{params.min_mOTU_MAGs} ##add nextflow param
    with open("!{motus_file}") as infile:
        for line in open("!{motus_file}"):
            if line.startswith("mOTU_"):
                fields = line.strip().split("\t")
                mOTU, rep, mean_ANI, min_ANI, missing_edges, nb_MAGs, nb_SUBS, MAGs, *SUBs = fields
                MAGs = MAGs.split(";")
                if int(nb_MAGs) >= min_genomes:
                    os.mkdir(mOTU)
                    #will add here to make subs optional
                    if int(nb_SUBS) > 0:
                        SUBs = SUBs[0].split(";")
                        MAGs.extend(SUBs)
                    for genome in MAGs:
                        #move file to mOTU directory
                        shutil.copy2(genome + ".fa", mOTU + "/")
                        print(f"{genome} being copied into {mOTU} directory")
    """
}

process mOTUs_to_pangenome {
    publishDir "${params.project}/pangenomes", mode: "copy"
    input:
    path(mOTU_dir)
    path(checkm_file)
    output:
    path("pan_*", type: "dir")
    shell:
    """
    #!/bin/bash -ue
    #nextflow didn't interact well with the dollar-sign for command substitution, so I had to use the deprecated backquote method
    if (( `ls !{mOTU_dir}/* | wc -l` > 1 )); then #checking number of fasta files in mOTU_dir, no need for SuperPang if only one
        echo "Enough genomes to run pangenome computation"
        SuperPang.py --fasta !{mOTU_dir}/* --checkm !{checkm_file} --output-dir pan_!{mOTU_dir}
    else
        echo "Only one genome in mOTU. Copying to pangenome dir."
        mkdir pan_!{mOTU_dir}
        cp !{mOTU_dir}/* pan_!{mOTU_dir}/
    fi
    """
}
