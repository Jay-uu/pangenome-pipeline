/*
A process that subsamples a specififed number of reads from each sample
Input is a tab delimited samples file and path to the directory with the raw reads.
Output is a tuple of the sample name and the two resulting concatenated subsample files.
*/

process subsample_fastqs {
    label "low_cpu"
    tag "low_cpu"
    publishDir "${params.project}/subsamples/fastqs", mode: "copy", pattern: "*.fq.gz"
    input:
    path(sample)
    path(fastq_dir)
    output:
    tuple(val("${sample.baseName}"), path("sub_*.fq.gz"), emit: sub_reads)
    path("*_readcount.txt", emit: readcount)
    path("*.subsampled.samples", emit: sample_file)
    shell:
    $/
    #!/usr/bin/env python
    import os
    from pathlib import Path
    from subprocess import Popen
    from subprocess import PIPE
    from subprocess import run
    from subprocess import check_output
    import gzip
    import shutil

    NR_SUBSAMP = int("!{params.nr_subsamp}")
    FASTQ_FILES = os.listdir("!{fastq_dir}")
    SAMPLE_ID = Path("!{sample}").stem
    
    """
    A function which takes a list of fastq files and subsamples a number of reads from the combined files.
    Result is a compressed fq file.
    """
    def concat_subtk_compress(file_list, direction, nr_subsamp):
        file_list = ["!{fastq_dir}/" + readfile for readfile in file_list]
        with open(f"sub_{SAMPLE_ID}_{direction}.fq.gz", "w") as subout:
            concat = Popen(["cat", (" ").join(file_list)], stdout=PIPE)
            subtk = Popen(["seqtk", "sample", "-s100", "-", f"{nr_subsamp}"], stdin=concat.stdout, stdout=PIPE)
            run(["gzip"], stdin=subtk.stdout, stdout=subout)
        return f"sub_{SAMPLE_ID}_{direction}.fq.gz"  
    
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
            if filename not in FASTQ_FILES:
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
    
    tot_reads = 0
    for fq in (fwds + revs):
        print(f"Counting reads in {fq}")
        reads_bases = check_output(["seqtk", "size",f"!{fastq_dir}/{fq}"], text=True)
        reads = int(reads_bases.split()[0]) #0 is nr reads, 1 is nr nucleotides
        tot_reads = tot_reads + reads
      
    if tot_reads < NR_SUBSAMP:
        print(f"Less reads available than requested subsampling. Subsampling {tot_reads} number of reads instead.")
        NR_SUBSAMP = tot_reads
    
    #subsampling and concatenating
    rev_out = ""
    #checking if we have reverse reads first to set the subsampling from each file to the correct number
    if len(revs) > 0:
        NR_SUBSAMP = NR_SUBSAMP/2
        rev_out = concat_subtk_compress(revs, "R2", NR_SUBSAMP)
    fwd_out = concat_subtk_compress(fwds, "R1", NR_SUBSAMP)
        
    #saving a new .samples file
    with open(f"{SAMPLE_ID}.subsampled.samples", "w") as outfile:
        outfile.write(f"{SAMPLE_ID}\t{fwd_out}\tpair1")
        if rev_out:
            outfile.write(f"\n{SAMPLE_ID}\t{rev_out}\tpair2")
        
    #write file with samp_name, nr fqs and tot_reads
    with open(f"{SAMPLE_ID}_readcount.txt", "w") as out:
        out.write("Sample\tNr_fastqs\tTotal_reads\n")
        out.write("\t".join([f"{SAMPLE_ID}", str(fq_count), str(tot_reads)+"\n"]))    
    /$
}
