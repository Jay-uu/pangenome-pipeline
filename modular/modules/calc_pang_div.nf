/*
This process runs Pogenom via a python wrapper.
Input: 
	Tuple of:
	pang_ID string with name of the pan-/reference genome
	vcf: VCF that has been downsampled to even coverage
	gff: gff output from SqueezeMeta
	genome: fasta file of pan-/reference genome.
Output:
	FST files from pogenom.
*/
process calc_pang_div {
    publishDir "${params.project}/mOTUs/results/${pang_ID}/pangenome/pogenom", mode: "copy"
    errorStrategy "ignore"
    tag "no_label"
    input:
    tuple(val(pang_ID), path(vcf), path(gff), path(genome))
    output:
    path("results", type: "dir", emit: pog_dir)
    shell:
    '''
    echo "Running pogenom for !{pang_ID}"
    run-pogenom.py !{vcf} -f !{genome} --gff !{gff} --minCount !{params.min_locus_cov}-t !{task.cpus} -p !{pang_ID} -o results
    #removing tmp dir
    rm -r results/temp
    '''
}
