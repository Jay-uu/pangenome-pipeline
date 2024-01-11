/*

*/
process calc_pang_div {
    publishDir "${params.project}/pogenom/", mode: "copy"
    errorStrategy "ignore"
    tag "no_label"
    input:
    tuple(val(pang_ID), path(vcf), path(gff), path(genome))
    output:
    path("${pang_ID}_out", type: "dir", emit: pog_dir)
    shell:
    '''
    echo "Running pogenom for !{pang_ID}"
    run-pogenom.py !{vcf} -f !{genome} --gff !{gff} -t !{params.threads} -p !{pang_ID} -o !{pang_ID}_out
    '''
}
