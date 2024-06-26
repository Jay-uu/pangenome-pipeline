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
//File with which fastq files belong to which samples. Tab delimited with sample-name, fastq file name and pair.
sam_chan = Channel.fromPath(params.samples, type: "file", checkIfExists: true)

/*The fastq_dir is needed for:
	- Formating the individual sample files
	- Running SqueezeMeta for each sample.
	- Subsampling
	- Mapping reads to the pangenomes
This creates four channels for the four different processes.
*/
Channel.fromPath(params.fastq, type: "dir", checkIfExists: true)
	    .multiMap { dir -> format: to_bins: subsample: to_pang_to_bams: dir }.set { fastq_chan }

/*Runs the process that creates individual samples files and creates three output channels:
	- For Squeezemeta (fastq_to_bins process)
	- For the subsampling
	- For creating samples files to the pangenomes
*/
format_samples(sam_chan, fastq_chan.format)
format_samples.out.flatten().multiMap { samp -> to_bins: subsample: to_covpang: samp }.set { single_samps }

/* Runs the fastq_to_bins process.
	flatten() to handle the created singles samples from format_samples individually
	first() to supply the fastq dir for each sample
*/
fastq_to_bins(single_samps.to_bins, fastq_chan.to_bins.first())

//Concatenating fastqs and subsampling for later mapping for each singles sample
subsample_fastqs(single_samps.subsample, fastq_chan.subsample.first())

/*
Before running mOTUlizer, the checkM and GTDB-Tk outputs (bintables) need to be parsed.
All bintables and all bins from different samples need to be collected so the taxonomy_parser process can run once with all data.
*/

fastq_to_bins.out.bintable.collect().multiMap { bintables -> to_tax_parser: bintables }.set { all_bintables } //SHOULD EDIT HERE TO REMOVE MULTIMAP
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
mOTUs_to_pangenome.out.pangenome_dir.multiMap { dir -> to_map: to_checkm: to_checkm2: dir }.set { pang_dirs }

/*
Run checkM on pangenomes for validation
*/
checkm_pangenomes(pang_dirs.to_checkm)
checkm2_pangenomes(pang_dirs.to_checkm2)
/*
Index pangenomes for read mapping
*/
index_pangenomes(pang_dirs.to_map)

/*
map subset reads to pangenome and get coverage information
*/
map_subset(index_pangenomes.out.pang_index.combine(subsample_fastqs.out.sub_reads))

/*
Using the coverage from the mapping, decides which reads "belong" to which pangenome and creates new .samples files
*/
cov_to_pang_samples(map_subset.out.coverage.collect(),single_samps.to_covpang.collect(), subsample_fastqs.out.readcount.collect())

/*
Create keys to match right samples file to right NBPs fasta (from same pangenome) as input to pang_to_bams.
*/
cov_to_pang_samples.out.pang_samples
		.flatten()
		.map { [it.getSimpleName(), it] }
		.set { pang_samples }

mOTUs_to_pangenome.out.NBPs_fasta
		.map { [it.getSimpleName(), it] }
		.set { NBPs_fasta }

/*
Using the generated samples files for the pangenome, the raw reads and the pangenome assembly to map reads using SqueezeMeta.
*/
pang_to_bams(pang_samples.combine(NBPs_fasta, by: 0),
				fastq_chan.to_pang_to_bams.first())

/*
Checking the breadth and the coverage of bams on the pangenome/ref-genome. Downsampling to even coverage and merging into one bam-file.
*/
downsample_bams_merge(pang_to_bams.out.pang_sqm)

/*
Running freebayes on the merged bam to get a filtered vcf file.
*/
detect_variants(downsample_bams_merge.out.ref_merged)

/*

*/
//pogenom(detect_variants.out.vcf, pang_sqm.to_pogenom)

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
    SqueezeMeta.pl -m coassembly -f !{fastq_dir} -s !{sample} -p $SAMPLE_ID -binners maxbin,metabat2,concoct -t !{params.threads} --onlybins --gtdbtk
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
    path("*_readcount.txt", emit: readcount)
    shell:
    $/
    #!/usr/bin/env python
    import os
    from pathlib import Path
    from subprocess import call
    from subprocess import check_output
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
    
    fq_count = 0 #number of read files for this sample
    with open("!{sample}") as infile:
        for line in open("!{sample}"):
            fq_count = fq_count+1
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
    
    tot_reads = 0
    for fq in (fwds + revs):
        print(f"Counting reads in {fq}")
        reads_bases = check_output(["seqtk", "size",f"!{fastq_dir}/{fq}"], text=True)
        reads = int(reads_bases.split()[0]) #0 is nr reads, 1 is nr nucleotides
        tot_reads = tot_reads + reads
        
    #write file with samp_name, nr fqs and tot_reads
    with open(f"{sample_ID}_readcount.txt", "w") as out:
        out.write('Sample\tNr_fastqs\tTotal_reads\n')
        out.write('\t'.join([f"{sample_ID}", str(fq_count), str(tot_reads)]))    
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
Output is a directory with a subdirectory for the mOTU, containing the pangenome fasta file and the NBPs.fasta in an individual channel.
This could possibly be changed to have different script parts, which would mean that I can have the python code for changing the name of the file and the fasta headers (the else part) directly in the process instead of in a separate script.
*/
process mOTUs_to_pangenome {
    publishDir "${params.project}/", mode: "copyNoFollow"
    input:
    tuple(path(mOTU_dir), path(bintable))
    output:
    path("pangenomes/${mOTU_dir}", type: "dir", emit: pangenome_dir)
    path("pangenomes/${mOTU_dir}/*.NBPs.fasta", emit: NBPs_fasta)
    shell:
    $/
    
    #!/usr/bin/env python
    from Bio import SeqIO
    from subprocess import call
    from pathlib import Path
    import os
    import glob
    
    genomes = glob.glob("!{mOTU_dir}/*.fa")
    nr_genomes = len(genomes)
    pg_dir_name = "pangenomes"
    os.makedirs(pg_dir_name)
    
    if nr_genomes > 1:
        print("Enough genomes to run pangenome computation")
        with open("input.fa", "w") as fastas:
            fastas.write("\n".join(genomes))
        call(["SuperPang.py", "--fasta", "input.fa", "--output-dir", f"{pg_dir_name}/!{mOTU_dir}", "--header-prefix", f"!{mOTU_dir}",
        "--output-as-file-prefix", "--nice-headers", "--debug"]) #====REMOVE DEBUG LATER======
    
    elif nr_genomes == 1:
        print("Only one genome in mOTU. Renaming headers and copying to pangenome dir.")

        outfile= f"{pg_dir_name}/!{mOTU_dir}/" + "!{mOTU_dir}" + ".singlemOTU.core.fasta"
        symfile = f"{pg_dir_name}/!{mOTU_dir}/" + "!{mOTU_dir}" + ".singlemOTU.NBPs.fasta"
        pangenome_file=f"{genomes[0]}"
        os.makedirs(pg_dir_name+"/!{mOTU_dir}")
        with open(pangenome_file, "rt") as pg:
            all_contigs = list(SeqIO.parse(pg, "fasta")) 
            pg_name = pangenome_file.split('.',1)[0]
            for i,s in enumerate(all_contigs):
                s.id =f"{pg_name}_NODE_Sc{i}-0-noinfo_length_{len(s.seq)}_cov_1.tag_noinfo" #renames contigs
                s.description="" #remove the old name
            with open(outfile, "w") as output_handle:
                SeqIO.write(all_contigs, output_handle, "fasta")
        #symlinking so that it can be sent as NBPs_fasta output
        Path(symfile).symlink_to("!{mOTU_dir}" + ".singlemOTU.core.fasta")
    else:
        raise Exception(f"No fastas found in !{mOTU_dir}")
    /$
}

/*
Running checkm on pangenomes
*/
process checkm_pangenomes {
    publishDir "${params.project}/checkm_pangenomes", mode: "copy" //change after deciding whether to use checkm or checkm2
    input:
    path(pangenome_dir)
    output:
    path("*_cM1_summary.txt")
    shell:
    '''
    pang_id=!{pangenome_dir.baseName}
    #using the SqueezeMeta installation of checkm
    installpath=$CONDA_PREFIX/SqueezeMeta/bin

    echo "Running checkM on all fastas in $pang_id"
    PATH=$installpath:$installpath/pplacer:$installpath/hmmer:$PATH $installpath/checkm lineage_wf -t !{params.threads} -x fasta !{pangenome_dir} ${pang_id}_cM1
    PATH=$installpath:$installpath/pplacer:$installpath/hmmer:$PATH $installpath/checkm qa ${pang_id}_cM1/lineage.ms ${pang_id}_cM1 > ${pang_id}_cM1_summary.txt
    
    '''
}

/*
This process might be removed if we decide to use checkM instead
*/
process checkm2_pangenomes {
    publishDir "${params.project}/checkm_pangenomes", mode: "copy" //change after deciding whether to use checkm or checkm2
    conda '/home/jay/mambaforge/envs/checkm2' //temporary until we decide which checkM version to use
    input:
    path(pangenome_dir)
    output:
    path("*_cM2", type: "dir")
    shell:
    '''
    pang_id=!{pangenome_dir.baseName}
    
    echo "Running checkM2 on core genome."
    checkm2 predict --threads !{params.threads} --input !{pangenome_dir} -x fasta --output-directory ${pang_id}_core_cM2 --output-directory ${pang_id}_cM2

    '''
}

/*
Indexing the pangenomes to use for read mapping.
Input is the directory with SuperPang output for a pangenome.
Output is the same directory, all of the index files and the environment variable set to either the base pangenome name,
or the base pangenome name + "_consensus" if there was no core genome.
This is separate from map_subset, since each pangenome only needs to be indexed once but will likely have multiple samples mapped to it.
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
This process aligns raw reads to previously indexed genomes.
It uses the names of the read files to determine if it should run in paired-end mode or not.
Input is somewhat complicated. It takes a tuple of a directory with SuperPang output, the index files for the pangenome,
a string with the pangenome name, a string with the name of the sample and the raw reads for that sample.
Output is the coverage information of how well the reads mapped to the genome.
*/
process map_subset {
    input:
    tuple(path(pangenome_dir),path(index), val(pang_id), val(sample_ID), path(sub_reads)) //use the combine operator on the channels in the workflow.
    output:
    path("*_coverage.tsv", emit: coverage)
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

/*
This process calculates which pangenomes have enough samples that pass the expected average coverage threshold to create new .samples files
for the pangenomes that will then be used to align the reads of those samples to the pangenomes for further analysis. Input should be the 
collected output from map_subset coverage, single_samps and subsample_fastqs readcount, so that all samples are processed with one process run.
Input: 
       coverage: A file with samtools coverage output, which has the depth of how reads from a sample mapped to a pangenome/reference genome.
       singles: The previously generated .samples files with which reads belong to which sample.
       readcounts: A file named {sample}_readcount.txt which conatains how many fastq files with reads belong to the sample, and the total number of reads 
       in those fastqs.
Output:
       pang_cpm_cov: Two files, one with the calculated results for CPM and the other with coverage.
       pang_samples: The new .samples files for the pangenomes.
*/
process cov_to_pang_samples {
    label "medium_time"
    publishDir "${params.project}/pangenomes", mode: "copy"
    input:
    path(coverage)
    path(singles)
    path(readcounts)
    output:
    path("samples/*.tsv", emit: pang_cpm_cov)
    path("samples/*.samples", emit: pang_samples)
    shell:
    $/
    #!/usr/bin/env python
    import os
    import pandas as pd
    import glob
    
    cov_threshold = !{params.mean_cov_threshold}
    nr_samps_threshold = !{params.nr_samps_threshold}
    outdir = "samples"
    
    """
    Input: 
        cov = dataframe of samtools coverage output
        nr_fqs = nr of read files to the sample
        tot_reads = total count of reads from the sample read files
    Returns the weighted mean coverage per read, and the expected average
    coverage for all the reads of that sample to the pangenome.
    """
    def get_weighted_mean(cov, nr_fqs, tot_reads):
        #meandepth is the mean depth of coverage over that contig
        #endpos-startpos=contig length
        cov["contig_length"] = cov.apply(lambda row : (row.endpos - row.startpos)+1, axis=1)
        cov["totaldepth"] = cov.apply(lambda row : row.meandepth*row.contig_length, axis=1) #rounds to nearest int
        weigh_mean=(cov["totaldepth"].sum()/cov["contig_length"].sum())
        #expected average coverage: ((average/million reads)/million*) total number of reads per sample
        weigh_mean_per_read = (weigh_mean/(1000000*nr_fqs))
        exp_cov = weigh_mean_per_read*tot_reads
        return weigh_mean_per_read, exp_cov
    
    count_files = glob.glob("*_readcount.txt")
    readcount = pd.DataFrame()
    for file in count_files:
        print(f"Appending {file}")
        readcount = pd.concat((readcount, pd.read_csv(file, sep='\t')), ignore_index=True)

    cpm_dic = {}
    cov_dic = {}
    cov_files = glob.glob("*_coverage.tsv")
    for file in cov_files:
        print(f"Reading {file}")
        pang_id = file.split("_sub_",1)[0]
        samp_name = file.split("sub_",1)[1].split("_")[0]
        cov = pd.read_csv(file, sep="\t")
        nr_fqs = readcount[readcount["Sample"]==samp_name]["Nr_fastqs"].item()
        tot_reads = readcount[readcount["Sample"]==samp_name]["Total_reads"].item()
        cpm, exp_cov = get_weighted_mean(cov, nr_fqs, tot_reads)
        cpm_dic.setdefault(pang_id, {})[samp_name] = cpm
        cov_dic.setdefault(pang_id, {})[samp_name] = exp_cov
        
        
    #convert dicts to dfs and save to file
    print("Saving coverage and CPM to files.")
    os.makedirs(outdir)
    all_cov = pd.DataFrame.from_dict(cov_dic, orient="index")
    all_cov = all_cov.reset_index().rename(columns={"index": "Pangenome"})
    all_cov.to_csv(f"{outdir}/pangenomes_cov.tsv", sep = '\t', index=False)
    
    all_cpm = pd.DataFrame.from_dict(cpm_dic, orient="index")
    all_cpm = all_cpm.reset_index().rename(columns={"index": "Pangenome"})
    all_cov.to_csv(f"{outdir}/pangenomes_cpm.tsv", sep = '\t', index=False)

    #create new .samples files for pangenomes that fit the coverage criteria
    print("Creating samples files for pangenomes that pass the thresholds.")
    for pang_id in cov_dic.keys():
        ovr_thresh = []
        for sample in cov_dic[pang_id].keys():
            #print(f"Checking {sample} coverage on {pang_id}")
            if cov_dic[pang_id][sample].item() >= cov_threshold:
                ovr_thresh.append(sample)
        if len(ovr_thresh) >= nr_samps_threshold:
            samp_df = pd.DataFrame()
            for sample in ovr_thresh:
                samples_file = pd.read_csv(f"{sample}.samples", sep='\t', names=["sample","read", "pair"])
                samp_df = pd.concat((samples_file, samp_df))
            print(f"Creating new samples file for {pang_id}")
            samp_df.to_csv(f"{outdir}/{pang_id}.samples", header=False, index=None, sep='\t')
            
    if len(glob.glob(f"{outdir}/*.samples")) < 1:
        raise Exception("It seems none of your pangenomes fulfill the thresholds for further analysis. Consider lowering --mean_cov_threshold and/or --nr_samps_threshold, or perhaps using more samples.")
   
    /$
}

/*
Runs squeezemeta to map samples to the pangenome consensus assembly or a reference genome.
Input:
      A tuple with pangenome name, the .samples file for the pangenome/reference genome and the the fasta file for the pangenome/ref genome.
      The directory with raw reads.
Output:
      The new SqueezeMeta output directory.
*/
process pang_to_bams {
    publishDir "${params.project}/pangenomes/sqm", mode: "copy"
    input:
    tuple(val(pang_ID), path(samples), path(pang_fasta))
    path(fastq_dir)
    output:
    path("${pang_ID}", type: "dir", emit: pang_sqm)
    shell:
    '''
    echo "Running SqueezeMeta on pangenome/reference genome !{pang_fasta} to map reads."
    #skips binning, assembly and renaming since we already have these things.
    #Mapping reads with a minimum of 95% identity using bowtie2
    SqueezeMeta.pl -m coassembly -p !{pang_ID} -f !{fastq_dir} -s !{samples} -extassembly !{pang_fasta} -t !{params.threads} --nobins --norename -mapping_options "--ignore-quals --mp 1,1 --np 1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.05" 
    '''
}

/*
This process takes the SqueezeMeta output and downsamples the bam files to get even coverage between samples, in preparation for Variant Calling.
Input: the results directory from SqueezeMeta.
Output: Since there is a possibility that no bams fit the minimum coverage and breadth criteria, this process might have no output or it will send
        a tuple with a fasta file of all contigs longer than 1000 bases from the input pangenome and a merged bam-file from all bams that passed the breadth
        and coverage criteria.
The downsampling shell code is modified from POGENOM's Input_pogenom pipeline by Anders Andersson and Luis F. Delgado
See here: https://github.com/EnvGen/POGENOM/blob/master/Input_POGENOM/src/cov_bdrth_in_dataset.sh
*/
//MAYBE REMOVE THE MERGED_BAMS TSV AND THE CODE FOR IT NOW SINCE VCFLIB SAMPLENAMES WORKS
process downsample_bams_merge {
    label "medium_time"
    input:
    path(pang_sqm)
    output:
    tuple(path("${pang_sqm}_long_contigs.fasta"), path("${pang_sqm}_merged.bam"), optional: true, emit: ref_merged)
    shell:
    '''    
    #Get total length of NBPs longer than 1000
    echo "Counting positions"
    positions=$(grep -oP '(?<=length_).*(?=_cov)' !{pang_sqm}/results/01.*.fasta | awk '{ if( $1*1 >= 1000) {SUM += $1} } END { print SUM }' )
    
    #Create tmp bams
    echo "Creating tmp bams"
    mkdir tmp_bams
    for bam in !{pang_sqm}/data/bam/*.bam;
    do
        #Filter to select only paired reads (-f 2) and avoids optical duplicates (-F 1024)
        samtools view -Sbh -F 1024 -q 20 --threads !{params.threads} $bam > tmp_filtered.bam
        #filter for contigs over 1000 bases put reads aligning to them in tmp_bams
        #names of contigs longer than 1000 in first column, and the length of contig in second column
        samtools idxstats tmp_filtered.bam | awk '$2 >= 1000 { print $0 }' > contigs.tsv
        awk ' { print $1, 1, $2} ' contigs.tsv > contigs.bed
        #create tmp bams
        bam_ID=$(basename $bam .bam)
        samtools view -b -L contigs.bed --threads !{params.threads} tmp_filtered.bam > tmp_bams/${bam_ID}.bam
    done
    
    mkdir -p !{pang_sqm}_mergeable
    echo "Creating mpileup files and checking breadth and coverage."    
    #create mpileup files
    for bam in tmp_bams/*.bam;
    do
        samtools mpileup -d 1000000 -Q 15 -a $bam > tmp.mpileup
    
        # ---- arguments
        mpileupfile=tmp.mpileup
        bamfile=$bam
        outbamfile=$(basename $bam bam)subsampled.bam #name of output
        mag=!{pang_sqm} #pangenome name
        mincov=!{params.min_med_cov}
        minbreadth=!{params.min_breadth}
        samplename=$(basename ${bam#"${mag}."} .bam)

        #--- Median coverage
        cov=$(cut -f4 $mpileupfile | grep -vw "0" | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')

        #---breadth
        non_zero=$(cut -f4 $mpileupfile | grep -cvw "0")
        breadth=$(echo $non_zero*100/$positions | bc -l )

        echo "Genome:" $mag "- Sample:" $samplename "Median_coverage:" $cov " breadth %:" $breadth

        #---selection of BAM files and subsample
        if (( $(echo "$breadth >= $minbreadth" | bc -l) )) && (( $(echo "$cov >= $mincov" | bc -l) )); then
            echo "        Downsampling coverage to $mincov - Genome: $mag - Sample: $samplename "
            limite=$(echo "scale=3; $mincov/$cov" | bc )
            samp=$(echo "scale=3; ($limite)+10" | bc)
            samtools view -Sbh --threads !{params.threads} -s $samp $bamfile | samtools sort -o !{pang_sqm}_mergeable/$outbamfile --threads !{params.threads}
        fi
    done
    
    #Merge bam-files that pass the check, if more than one bam in mergeable/ #*/ what to do if no files?
    #thoughts if at least one file in mergeable, create new fasta with only long contigs
    #also, how to name the output?
    echo "Checking mergeable"
    if [ -z "$(ls -A !{pang_sqm}_mergeable)" ]; then
         echo "No sample fit the alignment criteria. Skipping further analysis for !{pang_sqm}"
    else
        echo "Merging subsampled bams. and creating fasta of pangenome with only NBPs over 1000 bases."
        ls !{pang_sqm}_mergeable/*.bam > bamlist.txt
        samtools merge -o !{pang_sqm}_merged.bam -b bamlist.txt --threads !{params.threads}
        samtools idxstats !{pang_sqm}_merged.bam | awk '$2 >= 1000 { print $0 }' > long_contigs.tsv
        awk '{ print $1 }' long_contigs.tsv > contig_names.tsv
        seqtk subseq !{pang_sqm}/results/01.*.fasta contig_names.tsv > !{pang_sqm}_long_contigs.fasta
    fi
    #add cleanup step?
    '''

}

/*
Run freebayes, simple version for now. Check if can control to run more if it only uses 1 thread?
Input: a pangenome/reference fasta file and a (preferably) downsampled and merged bam-file.
Output: A filtered vcf file.
*/
process detect_variants {
    publishDir "${params.project}/pogenom/vcfs", mode: "copy"
    input:
    tuple(path(pangenome), path(bam))
    output:
    path("*_filtered.vcf", emit: filt_vcf)
    path("*_samples.txt", emit: samps_txt)
    shell:
    '''
    pang_ID=$(basename !{pangenome} _long_contigs.fasta)
    #why cant I find an option to specify threads...
    echo "Detecting variants in !{pangenome}"
    echo "Indexing !{bam}"
    samtools index !{bam}
    echo "Indexing !{pangenome}"
    samtools faidx !{pangenome} -o !{pangenome}.fai
    echo "Running freebayes"
    freebayes -f !{pangenome} -C 4 -p 1 --pooled-continuous --read-max-mismatch-fraction 0.05 --min-alternate-fraction 0.01 -q 15 --report-monomorphic !{bam} > ${pang_ID}_unfiltered.vcf
    echo "Filtering vcf"
    vcffilter -f 'QUAL > 20' ${pang_ID}_unfiltered.vcf > ${pang_ID}_filtered.vcf
    vcfsamplenames ${pang_ID}_filtered.vcf > ${pang_ID}_samples.txt
    echo "Done"
    '''
}

/*

*/
process pogenom {
    publishDir "${params.project}/pogenom", mode: "copy"
    input:
    //vcf file
    path(vcf)
    //pang_sqm for gff file
    path(pang_sqm)
    output:
    path(unknown)
    shell:
    '''
    #Fernando will give me a wrapper script for this, that I will call here
    '''
}





