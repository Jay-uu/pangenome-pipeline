/*
Takes a directory with bins and runs them through SqueezeMeta, resulting in taxonomic classification of the bins etc..
Output is the dir with all SqueezeMeta results, the bins, and the combined checkM and GTDB-Tk results.
*/
process classify_bins {
    publishDir "${project_path}/bins/bintables", mode: "copy", pattern: "${sample.baseName}/results/18.*.bintable", saveAs: { file("${project_path}").baseName + ".bintable" }
    publishDir "${project_path}/bins/fastas", mode: "copy", pattern: "${sample.baseName}/results/bins/*.fa", saveAs: { filename -> filename.split("/")[-1] }
    label "classify_bins"
    tag "${sample.baseName}"
    input:
    path(sample)
    path(in_bins)
    path(fastq_dir)
    val(project_path)
    output:
    path("${sample.baseName}/results/bins/*.fa", emit: bins)
    tuple(val("${sample.baseName}"),path("${sample.baseName}"), emit: sqm_dir)
    //path("${sample.baseName}/results/18.*.bintable", emit: bintable)
    script:
    """
    #!/usr/bin/env bash
    echo "The sample file is ${sample.baseName}, the fastq dir is ${fastq_dir}, and the bins dir is ${in_bins}"
    SAMPLE_ID="${sample.baseName}"
    SqueezeMeta.pl -m extbins -f ${fastq_dir} -s ${sample} -p \$SAMPLE_ID -r ${in_bins} --gtdbtk -test 1 -t ${task.cpus} -–nomarkers
    """
}
