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
