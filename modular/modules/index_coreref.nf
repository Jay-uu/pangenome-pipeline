/*
Indexing the pangenome or reference genome to use for read mapping.
Input is 
Output is .
This is separate from map_subset, since each pangenome only needs to be indexed once but will likely have multiple samples mapped to it.
*/
process index_coreref {
    label "low_cpu"
    tag "${fasta.baseName}"
    input:
    path(fasta)
    output:
    tuple(path(fasta),path("index*"), env(fasta_id), emit: fasta_index_id) 
    shell:
    '''
    fasta_id=!{fasta.baseName}
    echo "Building index for ${fasta_id}"
    bowtie2-build !{fasta} index
    '''
}
