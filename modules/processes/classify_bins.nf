/*
This process takes a directory with bins in fasta format, and formats everything to a SqueezeMeta directory.
Might rename contigs, but keeps bin-names.
*/
process classify_bins {
    publishDir "${project_path}/bins/bintables", mode: "copy", pattern: "${sample.baseName}/results/18.*.bintable", saveAs: { file("${project_path}").baseName + ".bintable" }
    publishDir "${project_path}/bins/fastas", mode: "copy", pattern: "${sample.baseName}/results/bins/*.fa", saveAs: { filename -> filename.split("/")[-1] }
    tag "${sample.baseName}"
    input:
    path(sample)
    path(in_bins)
    path(fastq_dir)
    val(project_path)
    output:
    path("${sample.baseName}/results/bins", emit: bins_dir)
    tuple(val("${sample.baseName}"),path("${sample.baseName}"), emit: sqm_dir)
    path("${sample.baseName}/results/bins/*.fa", emit: bins)
    script:
    """
    #!/usr/bin/env bash
    echo "The sample file is ${sample.baseName}, the fastq dir is ${fastq_dir}, and the bins dir is ${in_bins}"
    SAMPLE_ID="${sample.baseName}" 
    SqueezeMeta.pl -m extbins -f ${fastq_dir} -s ${sample} -p \$SAMPLE_ID -r ${in_bins} --gtdbtk -test 1 -t ${task.cpus} --nomarkers
    """
}
