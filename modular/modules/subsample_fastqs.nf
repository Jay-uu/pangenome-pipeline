/*
A process that subsamples a million reads from each raw reads file for a sample,
and then concatenates paired subsampled reads.
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
    from subprocess import call
    from subprocess import check_output
    import gzip
    import shutil

    nr_subsamp = "!{params.nr_subsamp}"
    
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
            call(["seqtk", "sample", "-s100", f"concat_{direction}.fq.gz", nr_subsamp], stdout = subout)
        #Compressing
        with open(f"sub_{sample_ID}_{direction}.fq", "rb") as in_f, gzip.open(f"sub_{sample_ID}_{direction}.fq.gz", "wb") as out_f:
            shutil.copyfileobj(in_f, out_f)
        #Removing intermediate files
        os.remove(f"concat_{direction}.fq.gz")
        os.remove(f"sub_{sample_ID}_{direction}.fq")
        #return the compressed subsampled file name
        return f"sub_{sample_ID}_{direction}.fq.gz"

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
    
    #subsampling and concatenating
    fwd_out = concat_subtk_compress(fwds, "R1")
    rev_out = ""
    if len(revs) > 0:
        rev_out = concat_subtk_compress(revs, "R2")
        
    #saving a new .samples file
    with open(f"{sample_ID}.subsampled.samples", "w") as outfile:
        outfile.write(f"{sample_ID}\t{fwd_out}\tpair1")
        if rev_out:
            outfile.write(f"\n{sample_ID}\t{rev_out}\tpair2")
    
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
