/*
Run freebayes, simple version for now. Check if can control to run more if it only uses 1 thread?
Input: a pangenome/reference fasta file and a (preferably) downsampled and merged bam-file.
Output: A filtered vcf file.
*/
process detect_variants {
    debug true
    publishDir "${params.project}/pogenom/vcfs", mode: "copy"
    input:
    tuple(path(pangenome), path(bam))
    output:
    tuple(env(pang_ID), path("*_filtered.vcf"), optional: true, emit: filt_vcf)
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
