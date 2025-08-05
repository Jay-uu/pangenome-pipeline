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
    publishDir "${project_path}/mOTUs/results/${pang_ID}/pangenome/pogenom", mode: "copy", pattern: "results"
    errorStrategy "ignore"
    label "calc_pang_div"
    tag "${pang_ID}"
    input:
    tuple(val(pang_ID), path(vcf), path(gff), path(genome))
    val(project_path)
    val(min_locus_cov)
    output:
    path("results", type: "dir", optional: true, emit: pog_dir)
    path("SUCCESS.txt", type: "file", emit: success_message)
    script:
    """ 
    #!/usr/bin/env bash
    echo "Running pogenom for ${pang_ID}"
    run-pogenom.py ${vcf} -f ${genome} --gff ${gff} --minCount ${min_locus_cov} -t ${task.cpus} -p ${pang_ID} -o results
    #removing tmp dir
    rm -r results/temp
    
    if [ -z "\$( ls -A 'results' )" ]; then
        echo "Empty"
        echo "POGENOM could not generate results for ${pang_ID}." > SUCCESS.txt
        rm -r results
    else
        echo "Not Empty"
        echo "POGENOM results available for ${pang_ID}." > SUCCESS.txt
    fi
    """
}
