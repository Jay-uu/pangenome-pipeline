/*
Takes a directory with bins and runs them through SqueezeMeta, resulting in taxonomic classification of the bins etc..
Output is the dir with all SqueezeMeta results, the bins, and the combined checkM and GTDB-Tk results.
*/
process classify_bins {
    publishDir "${params.project}/bins/bintables", mode: "copy", pattern: "${sample.baseName}/results/18.*.bintable", saveAs: { filename -> "${sample.baseName}.bintable" }
    publishDir "${params.project}/bins/fastas", mode: "copy", pattern: "${sample.baseName}/results/bins/*.fa", saveAs: { filename -> filename.split("/")[-1] }
    tag "no_label"
    input:
    path(sample)
    path(in_bins)
    path(fastq_dir)
    output:
    path("${sample.baseName}/results/bins/*.fa", emit: bins)
    path("${sample.baseName}/results/18.*.bintable", emit: bintable)
    shell:
    '''
    echo "The sample file is !{sample.baseName}, the fastq dir is !{fastq_dir}, and the bins dir is !{in_bins}"
    SAMPLE_ID="!{sample.baseName}"
    SqueezeMeta.pl -m coassembly -f !{fastq_dir} -s !{sample} -p $SAMPLE_ID --extbins !{in_bins} --gtdbtk -test 1 -t !{task.cpus}
    17.checkM_batch.pl $SAMPLE_ID/
    18.getbins.pl $SAMPLE_ID/
    '''
}
