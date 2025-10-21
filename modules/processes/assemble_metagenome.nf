/*
Takes raw reads and runs them through SqueezeMeta, resulting in per-sample assemblies.
*/
process assemble_metagenome {
    //should sample sqm.zip be output also? 
    label "assemble_metagenome"
    tag "${sample.baseName}"
    input:
    path(sample)
    path(fastq_dir)
    val(binners)
    val(contig_len)
    output:
    tuple(val("${sample.baseName}"),path("${sample.baseName}"), emit: sqm_dir)
    script:
    """
    echo "The sample file is ${sample.baseName} and the fastq dir is ${fastq_dir}"
    #-test 1 means stopping after assembly.
    SqueezeMeta.pl -m coassembly -f ${fastq_dir} -samples ${sample} -p ${sample.baseName} -binners ${binners} -t ${task.cpus} -contiglen ${contig_len} --onlybins --gtdbtk --nomarkers -test 1
    """
}
