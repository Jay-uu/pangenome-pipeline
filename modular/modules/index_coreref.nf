/*
Indexing the pangenome or reference genome to use for read mapping.
Input is 
Output is .
This is separate from map_subset, since each pangenome only needs to be indexed once but will likely have multiple samples mapped to it.
*/
process index_coreref {
    label "low_cpu"
    label "index_coreref"
    tag "${fasta.baseName}"
    input:
    path(fasta)
    output:
    tuple(path(fasta),path("index*"), env(fasta.baseName), emit: fasta_index_id) 
    script:
    """
    echo "Building index for ${fasta.baseName}"
    bowtie2-build ${fasta} index
    """
}
