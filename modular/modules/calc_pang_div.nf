/*

*/
process calc_pang_div {
//publishdir looks like this currently: mOTUs/results/c__Actinomycetia_mOTU_0/pangenome/pogenom/results/motu_name_out/*files*
//files could probably be directly in results?
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
    run-pogenom.py !{vcf} -f !{genome} --gff !{gff} -t !{task.cpus} -p !{pang_ID} -o results
    
    #removing tmp dir
    rm -r results/tmp
    '''
}
