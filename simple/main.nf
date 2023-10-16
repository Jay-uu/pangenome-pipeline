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
Before running mOTUlizer, the checkM and GTDB-Tk outputs (bintables) need to be parsed.
All bintables and all bins from different samples need to be collected so the taxonomy_parser process can run once with all data.
*/

fastq_to_bins.out.bintable.collect().multiMap { bintables -> to_tax_parser: to_SuperPang: bintables }.set { all_bintables }
fastq_to_bins.out.bins.collect().multiMap { bins -> to_tax_parser: to_mOTU_dirs: bins }.set { all_bins }

parse_taxonomies(all_bins.to_tax_parser, all_bintables.to_tax_parser) 

/*
Clustering of bins, if they've been presorted to lower taxonomic ranks this can spawn parallell processes
*/
bins_to_mOTUs(parse_taxonomies.out.tax_bin_dirs.flatten())


/*
Creating dirs for the mOTUs by sorting based on the mOTUlizer output, so each mOTU directory has the correct bins.
*/
create_mOTU_dirs(bins_to_mOTUs.out.mOTUs_file, all_bins.to_mOTU_dirs)

/*
Running SuperPang, creating pangenomes. Transpose makes it so that each mOTU from the same grouping within
the taxonomy selection will also be sent individually to the process together with the matching bintable.
*/

mOTUs_to_pangenome(create_mOTU_dirs.out.transpose())

/*
create several channels for pangenomes since they are going to multiple processes
*/
mOTUs_to_pangenome.out.pangenome_dir.multiMap { dir -> to_map: to_checkm: dir }.set { pang_dirs }

/*
Index pangenomes for read mapping
*/
index_pangenomes(pang_dirs.to_map)

/*
map subset reads to pangenome
*/
map_subset(index_pangenomes.out.pang_index.combine(subsample_fastqs.out.sub_reads))

//pang_dirs.to_map.combine(subsample_fastqs.out.sub_reads).view()
//map_subset(pang_dirs.to_map.combine(subsample_fastqs.out.sub_reads))


//=======END OF WORKFLOW=============
}
//=======START OF PROCESS DEFINITIONS=============

/*
This process takes a tab-delimited samples file and converts it to individual files per sample.
It also checks that the provided fastq-dir has the files specified in the samples file.
Output is a list of all the individual samples files.
*/
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

/*
Takes raw reads and runs them through SqueezeMeta, resulting in bins.
Output is the dir with all SqueezeMeta results, the bins, and the combined checkM and GTDB-Tk results.
*/
process fastq_to_bins {
    publishDir "${params.project}/sqm_res", mode: "copy", pattern: "${sample.baseName}"
    input:
    path(sample)
    path(fastq_dir)
    output:
    path("${sample.baseName}", emit: sample_dir)
    path("${sample.baseName}/results/bins/*.contigs.fa", emit: bins)
    path("${sample.baseName}/results/18.*.bintable", emit: bintable)
    shell:
    '''
    echo "The sample file is !{sample.baseName} and the fastq dir is !{fastq_dir}"
    SAMPLE_ID="!{sample.baseName}"
    SqueezeMeta.pl -m coassembly -f !{fastq_dir} -s !{sample} -p $SAMPLE_ID -binners maxbin,metabat2,concoct -t !{params.threads} -contigid $SAMPLE_ID --only-bins --gtdbtk
    '''

}
/*
A process that subsamples a million reads from each raw reads file for a sample,
and then concatenates paired subsampled reads.
Input is a tab delimited samples file and path to the directory with the raw reads.
Output is a tuple of the sample name and the two resulting concatenated subsample files.
*/
process subsample_fastqs {
    label "medium_time"
    input:
    path(sample)
    path(fastq_dir)
    output:
    tuple(val("${sample.baseName}"), path("sub_*.fq.gz"), emit: sub_reads)
    shell:
    $/
    #!/usr/bin/env python
    import os
    from pathlib import Path
    from subprocess import call
    import gzip
    import shutil
    
    """
    A function which takes a list of fastq files and subsamples a million reads from the combined files.
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
    #sorting the lists to make sure that matching fwds and revs are in the same index
    fwds = sorted(fwds)
    revs = sorted(revs)

    concat_subtk_compress(fwds, "R1")
    if len(revs) > 0:
        concat_subtk_compress(revs, "R2")
    
    /$
}
/*
Reads the GTDB-Tk results from the bintable files then sorts bins based on the group selected with params.taxSort.
This can for example lead to species-directories, with all bins that were identified a specific species within the same directory.
Input is all bins, and all bintables.
Output is the new directories which contains bins and a tsv with completeness and contamination.
*/
process parse_taxonomies {
    label "medium_time" //maybe test if short_time could be reasonable for this
    input:
    path(bins)
    path(all_bintables)
    output:
    path("*_bins", emit: tax_bin_dirs)
    shell:
    $/
    #!/usr/bin/env python
    import os
    import glob
    import pandas as pd
    import shutil
    ranks = ["root","auto","domain", "phylum","class","order", "family","genus","species"]
    rank_param = "!{params.taxSort}" #NF PARAM
    rank_param = rank_param.lower() #allow the user to spell with big letters
    if rank_param not in ranks:
        raise Exception(f"The provided sorting rank is not valid. Please use a taxonomic rank from the following list {ranks}")
    #then read files to df
    print("Reading bintables to dataframe")
    all_bintable_files = glob.glob("*.bintable")
    all_bintables = pd.DataFrame()
    for file in all_bintable_files:
        #read next file to be added, remove unneeded columns to save memory
        next_table = pd.read_csv(file, sep='\t', header=1)
        c_list = next_table.columns.to_list()
        c_list = [e for e in c_list if e not in 
                  ("Bin ID", "Tax GTDB-Tk", "Completeness", "Contamination")]
        next_table = next_table.drop(labels=c_list, axis=1)
        all_bintables = pd.concat((next_table, all_bintables), ignore_index=True)
    
    #rename bin id an taxonomy cols
    print("Renaming columns")
    all_bintables.rename(columns={"Bin ID": "Bin Id"}, inplace=True)
    all_bintables.rename(columns={"Tax GTDB-Tk": "root"}, inplace=True)
    #Convert to categorical data to save memory
    print("Convert taxonomy column to categorical")
    all_bintables["root"] = all_bintables["root"].astype("category")
    #too small genomes won't give any completeness and contamination data, which is necessary for mOTUlizer
    
    all_bintables.fillna(value={"Completeness" : 1, "Contamination" : 1}, inplace = True)
    #split Tax GTDB-Tk/root column so I can easily search ranks.
    print("Splitting taxonomy, and fixing Unclassified data.")
    all_bintables[ranks[2:]] = all_bintables["root"].str.split(";", expand=True)
    #making sure each rank column has either a classification or is unclassified
    rank_cols = (dict([[e,"Unclassified"] for e in ranks[2:]]))
    all_bintables.fillna(value=rank_cols, inplace=True)
    #it's aslo possible for bins that were only classified at lower levels
    #to have example s__ instead of a species. This should be unclassified instead.
    tax_short = ["d__", "p__","c__","o__","f__","g__","s__"]
    all_bintables.replace(to_replace=tax_short, value="Unclassified", inplace=True)
    if rank_param == "auto":
        #find lowest rank with more than 90% classified bins and change rank_param to that
        print("Selected rank: auto. Finding lowest well classified rank.")
        rank_not_selected = True
        i = 8
        tot_bins = len(all_bintables)
        while rank_not_selected:
            print(f"Checking unclassified proportion in {ranks[i]}.")
            if i == 2:
                rank_param = "root" #no lower classification was good enough quality
                rank_not_selected = False
            elif all_bintables[ranks[i]].value_counts()["Unclassified"]/tot_bins < 0.1:
                rank_param = ranks[i]
                rank_not_selected = False
            else:
                i = i-1 
    #start sorting based on rank_param
    #create dirs and taxonomic-sorted bintables 
    print(f"Sorting based on {rank_param}.")
    if rank_param == "root":
        all_bintables["root"] = "root"
    for group in all_bintables[rank_param].unique():
        print(f"Writing {group} to dirs")
        g = all_bintables[all_bintables[rank_param] == group]
        #check that the group has at least one bin of good enough quality for mOTUlizer. Otherwise skip.
        #MIGHT BE REMOVED IF ANOTHER SOLUTION IS DECIDED
        if len(g.loc[(g["Completeness"] > !{params.MAGcomplete}) & (g["Contamination"] < !{params.MAGcontam})]) == 0:
            print(f"No {group} bin is good enough quality for clustering. Exlcuding from further analysis.")
        else:
            os.makedirs(group+"_bins")
            g.to_csv(f"{group}_bins/{group}.bintable", index=None,
                     sep='\t', columns = ["Bin Id","Completeness","Contamination"])
            for bin_id in g['Bin Id']:
                shutil.copy2(bin_id+".fa", group+"_bins")
    /$
    
}


/*
Clusters bins based on similarity.
Input is a directory containing bins and a file with bin-name, completeness and contamination.
Output is the name of the taxonomic classification of the bins (unless taxSort = root),
a tsv with which bins belong to which mOTU, and the bintable file with quality data for the bins.
*/
process bins_to_mOTUs {
    /*the conda part might be removed later. If for example SuperPang gets updated to run with the newest version
    of mOTUlizer this process doesn't need a separate environment */
    conda 'bioconda::mOTUlizer=0.3.2'
    publishDir "${params.project}/mOTUs", mode: "copy"
    input:
    path(tax_dir)
    output:
    tuple(env(group), path("*_mOTUs.tsv"), path("${tax_dir}/*.bintable"), emit: mOTUs_file) //maybe change name to better represent content
    path("*_similarities.txt", emit: simi_file)
    shell:
    '''
    #!/bin/bash
    group="!{tax_dir}"
    group=${group%"_bins"}
    echo $group
    mOTUlize.py --fnas !{tax_dir}/*.fa --checkm !{tax_dir}/*.bintable --MAG-completeness !{params.MAGcomplete} --MAG-contamination !{params.MAGcontam} --threads !{params.threads} --keep-simi-file ${group}_similarities.txt -o ${group}_mOTUs.tsv
    '''
}

/*
Takes the output from mOTUlizer and the bins and sorts them into new directories based on which mOTU the belong to.
Input: Tuple with: group, name of the taxonomic classification/prefix for filenames. motus_file, has information about which bins belong to which
mOTU, and the bintable with bins quality data. Also takes the bins as input.
Output is a tuple with the new mOTU directory and the bintable.
*/
process create_mOTU_dirs {
    label "short_time"
    publishDir "${params.project}/mOTUs", mode: "copy", pattern: "${group}_mOTU_*"
    input:
    tuple(val(group), path(motus_file), path(bintable))
    path(bins)
    output:
    tuple(path("${group}_mOTU_*", type: "dir"), path("${bintable}"))
    shell:
    '''
    #!/usr/bin/env python
    import os
    import shutil
    min_genomes = !{params.min_mOTU_MAGs} #nextflow param
    with open("!{motus_file}") as infile:
        for line in open("!{motus_file}"):
            if line.startswith("mOTU_"):
                fields = line.strip().split("\t")
                mOTU, rep, mean_ANI, min_ANI, missing_edges, nb_MAGs, nb_SUBS, MAGs, *SUBs = fields
                MAGs = MAGs.split(";")
                if int(nb_MAGs) >= min_genomes:
                    os.mkdir("!{group}_" + mOTU)
                    #will add here to make subs optional
                    if int(nb_SUBS) > 0:
                        SUBs = SUBs[0].split(";")
                        MAGs.extend(SUBs)
                    for genome in MAGs:
                        #move file to mOTU directory
                        shutil.copy2(genome + ".fa", "!{group}_" + mOTU + "/")
                        print(f"{genome} being copied into !{group}_{mOTU} directory")
    '''
}
/*
Creates pangenomes from a directory with bins and a tsv with bin completeness and contamination by running SuperPang.
If there's not enough genomes to a mOTU, it will just copy the genome to results.
Input is a tuple containing a directory with bins (preferably belonging to the same mOTU) and a bintable with completeness and contamination.
Output is a directory with a subdirectory for the mOTU, containing the pangenome fasta file.
*/
process mOTUs_to_pangenome {
    publishDir "${params.project}/", mode: "copy"
    input:
    tuple(path(mOTU_dir), path(bintable))
    output:
    path("pangenomes/${mOTU_dir}", type: "dir", emit: pangenome_dir)
    shell:
    '''
    #!/bin/bash -ue
    #nextflow didn't interact well with the dollar-sign for command substitution, so I had to use the deprecated backquote method
    mkdir pangenomes
    if (( `ls !{mOTU_dir}/* | wc -l` > 1 )); then #checking number of fasta files in mOTU_dir, no need for SuperPang if only one
        echo "Enough genomes to run pangenome computation"
        SuperPang.py --fasta !{mOTU_dir}/* --checkm !{bintable} --output-dir pangenomes/!{mOTU_dir} --header-prefix !{mOTU_dir} --output-as-file-prefix --nice-headers --debug #====REMOVE DEBUG LATER======
    else
        echo "Only one genome in mOTU. Copying to pangenome dir."
        mkdir pangenomes/!{mOTU_dir}/ 
        #change name in pangenome dir to end with singlemOTU.core.fasta
        #also add to change header names?
        cp !{mOTU_dir}/*.fa pangenomes/!{mOTU_dir}/!{mOTU_dir}_singlemOTU.core.fasta
    fi
    '''
}

/*
Skeleton for running checkm on pangenomes
*/

process checkm_pangenomes {
    input:
    path(pangenome_dir)
    shell:
    '''
    #run checkm
    #run checkm2
    echo "nothing"
    '''
}

/*
indexing
*/
process index_pangenomes {
    input:
    path(pangenome_dir)
    output:
    tuple(path(pangenome_dir),path("index*"), env(pang_id), emit: pang_index) 
    shell:
    '''
    pang_file=!{pangenome_dir}/*.core.fasta
    pang_id=!{pangenome_dir.baseName}
    #*/ remove comment
    #Checking core genome existence
    if [ -s ${pang_file} ]; then
        echo "Mapping samples to core genome"
    else
        echo "The core pangenome file is empty. Using the consensus assembly instead."
        echo "This should not happen unless you're using mock communities."
        pang_file=!{pangenome_dir}/*.NBPs.fasta
        pang_id=$(basename $pang_file .NBPS.fasta)
        pang_id=${pang_id}_consensus
    fi
    
    echo "Building index"
    bowtie2-build $pang_file index
    '''
}

/*
Skeleton for mapping the subsets of the raw reads on the pangenomes
*/
process map_subset {
    input:
    tuple(path(pangenome_dir),path(index), val(pang_id), val(sample_ID), path(sub_reads)) //use the combine operator on the channels in the workflow.
    output:
    path("*_coverage.tsv")
    shell:
    '''
    #run bowtie2
    reads_id=$(basename sub_*_R1.fq.gz _R1.fq.gz)
    #check if there are sub*R2 reads, if yes:
    if stat --printf='' sub_*_R2.fq.gz 2>/dev/null; then
        echo "Running paired-end mode"
        bowtie2 -x index -1 sub_*_R1.fq.gz -2 sub_*_R2.fq.gz | samtools view -bS > tmp_alignment.bam
    else
        echo "Running unpaired reads mode"
        bowtie2 -x index -U sub_*_R1.fq.gz | samtools view -bS > tmp_alignment.bam
    fi
    
    echo "Sorting bam files"
    samtools sort tmp_alignment.bam -O BAM -o !{pang_id}_${reads_id}_alignment.bam --threads !{params.threads}
    
    echo "Computing coverage"
    samtools coverage !{pang_id}_${reads_id}_alignment.bam -o !{pang_id}_${reads_id}_coverage.tsv
    
    echo "Removing .bam files" #to save space
    rm *.bam #the sorted bams are named despite being deleted in case we decide a downstream process needs them.
    
    '''
}

