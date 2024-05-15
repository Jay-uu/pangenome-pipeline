/*
Takes raw reads and runs them through SqueezeMeta, resulting in bins.
Output is the dir with all SqueezeMeta results, the bins, and the combined checkM and GTDB-Tk results.
*/
process fastq_to_bins {
    publishDir "${params.project}/bins/fastas", mode: "copy", pattern: "${sample.baseName}/results/bins/*.fa", saveAs: { filename -> filename.split("/")[-1] }
    publishDir "${params.project}/bins/bintables", mode: "copy", pattern: "${sample.baseName}/results/18.*.bintable", saveAs: { filename -> "${sample.baseName}.bintable" }
    tag "no_label"
    input:
    path(sample)
    path(fastq_dir)
    output:
    path("${sample.baseName}/results/bins/*.fa", emit: bins)
    path("${sample.baseName}/results/18.*.bintable", emit: bintable)
    shell:
    '''
    echo "The sample file is !{sample.baseName} and the fastq dir is !{fastq_dir}"
    SAMPLE_ID="!{sample.baseName}"
    SqueezeMeta.pl -m coassembly -f !{fastq_dir} -s !{sample} -p $SAMPLE_ID -binners !{params.binners} -t !{task.cpus} --onlybins --gtdbtk
    '''
}
