/*
Running checkm on pangenomes
*/
process checkm_pangenomes {
    publishDir "${params.project}/mOTUs/${pangenome_dir.baseName}/checkm", mode: "copy"
    publishDir "${params.project}/checkm_pangenomes", mode: "copy" //remove when finalizing pipeline
    tag "checkm_pangenomes"
    input:
    path(pangenome_dir)
    output:
    path("*_cM1_summary.txt")
    shell:
    '''
    pang_id=!{pangenome_dir.baseName}
    #using the SqueezeMeta installation of checkm
    installpath=$CONDA_PREFIX/SqueezeMeta/bin

    echo "Running checkM on all fastas in $pang_id"
    PATH=$installpath:$installpath/pplacer:$installpath/hmmer:$PATH $installpath/checkm lineage_wf -t !{task.cpus} -x fasta !{pangenome_dir} ${pang_id}_cM1
    PATH=$installpath:$installpath/pplacer:$installpath/hmmer:$PATH $installpath/checkm qa -t !{task.cpus} ${pang_id}_cM1/lineage.ms ${pang_id}_cM1 > ${pang_id}_cM1_summary.txt
    
    '''
}
