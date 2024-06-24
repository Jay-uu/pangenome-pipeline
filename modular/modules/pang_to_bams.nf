/*
Runs squeezemeta to map samples to the pangenome consensus assembly or a reference genome.
Input:
      A tuple with pangenome name, the .samples file for the pangenome/reference genome and the the fasta file for the pangenome/ref genome.
      The directory with raw reads.
Output:
      The new SqueezeMeta output directory.
*/
process pang_to_bams {
    publishDir "${params.project}/mOTUs/results/${pang_ID}/pangenome", mode: "copy", pattern: "${pang_ID}.zip"
    label "no_label"
    tag "${pang_ID}"
    input:
    tuple(val(pang_ID), path(pang_fasta), path(samples))
    path(fastq_dir)
    output:
    path("${pang_ID}", type: "dir", emit: pang_sqm)
    tuple(val("${pang_ID}"), path("${pang_ID}/results/03.*.gff", type: "file"), path("${pang_fasta}"), emit: id_gff_genome)
    path("${pang_ID}.zip", emit: sqm_zip)
    shell:
    '''
    echo "Running SqueezeMeta on pangenome/reference genome !{pang_fasta} to map reads."
    #skips binning, assembly and renaming since we already have these things.
    #making a dir to be able to use -extbins flag, will make squeezemeta do some extra stuff!
    mkdir extassembly
    cp !{pang_fasta} extassembly/.
    #Mapping reads with a minimum of 95% identity using bowtie2
    SqueezeMeta.pl -m coassembly -p !{pang_ID} -f !{fastq_dir} -s !{samples} -extbins extassembly -t !{task.cpus} --nobins --norename -b !{params.block_size} -mapping_options "--ignore-quals --mp 1,1 --np 1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.05"
    sqm2zip.py !{pang_ID} .
    '''
}
