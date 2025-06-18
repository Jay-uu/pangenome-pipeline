/*
Takes a directory with bins and runs them through SqueezeMeta, resulting in taxonomic classification of the bins etc..
Output is the dir with all SqueezeMeta results, the bins, and the combined checkM and GTDB-Tk results.
*/
process classify_bins {
    //projname = file("${params.project}").baseName 
    //used for saving the bintable previously. testing using it directly within the saveAs closure.
    publishDir "${params.project}/bins/bintables", mode: "copy", pattern: "${sample.baseName}/results/18.*.bintable", saveAs: { file("${params.project}").baseName + ".bintable" }
    publishDir "${params.project}/bins/fastas", mode: "copy", pattern: "${sample.baseName}/results/bins/*.fa", saveAs: { filename -> filename.split("/")[-1] }
    label "classify_bins"
    tag "${sample.baseName}"
    input:
    path(sample)
    path(in_bins)
    path(fastq_dir)
    output:
    path("${sample.baseName}/results/bins/*.fa", emit: bins)
    path("${sample.baseName}/results/18.*.bintable", emit: bintable)
    script:
    """
    #!/usr/bin/env bash
    echo "The sample file is ${sample.baseName}, the fastq dir is ${fastq_dir}, and the bins dir is ${in_bins}"
    SAMPLE_ID="${sample.baseName}"
    SqueezeMeta.pl -m coassembly -f ${fastq_dir} -s ${sample} -p \$SAMPLE_ID --extbins ${in_bins} --gtdbtk -test 1 -t ${task.cpus}
    #17.checkM_batch.pl \$SAMPLE_ID/
    17.checkbins.pl \$SAMPLE_ID/
    18.getbins.pl \$SAMPLE_ID/
    """
}
