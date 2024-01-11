/*
Running checkm on pangenomes
*/
process checkm_pangenomes {
    publishDir "${params.project}/checkm_pangenomes", mode: "copy" //change after deciding whether to use checkm or checkm2
    tag "no_label"
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
    PATH=$installpath:$installpath/pplacer:$installpath/hmmer:$PATH $installpath/checkm lineage_wf -t !{params.threads} -x fasta !{pangenome_dir} ${pang_id}_cM1
    PATH=$installpath:$installpath/pplacer:$installpath/hmmer:$PATH $installpath/checkm qa -t !{params.threads} ${pang_id}_cM1/lineage.ms ${pang_id}_cM1 > ${pang_id}_cM1_summary.txt
    
    '''
}
