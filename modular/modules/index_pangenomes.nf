/*
Indexing the pangenomes to use for read mapping.
Input is the directory with SuperPang output for a pangenome.
Output is the same directory, all of the index files and the environment variable set to either the base pangenome name,
or the base pangenome name + "_consensus" if there was no core genome.
This is separate from map_subset, since each pangenome only needs to be indexed once but will likely have multiple samples mapped to it.
*/
process index_pangenomes {
    label "index_pangenomes"
    tag "${pangenome_dir}.baseName"
    input:
    path(pangenome_dir)
    output:
    tuple(path(pangenome_dir),path("index*"), env(pang_id), emit: pang_index) 
    shell:
    '''
    pang_file=!{pangenome_dir}/*.core.fasta
    pang_id=!{pangenome_dir.baseName}
    #Checking core genome existence
    if [ -s ${pang_file} ]; then
        echo "Mapping samples to core genome"
    else
        echo "The core pangenome file is empty. Using the consensus assembly instead."
        echo "This should not happen unless you're using mock communities."
        pang_file=!{pangenome_dir}/*.NBPs.fasta
        pang_id=$(basename $pang_file .NBPS.fasta)
        pang_id=${pang_id}_consensus
    fi
    
    echo "Building index"
    bowtie2-build $pang_file index --threads !{task.cpus}
    '''
}
