/*
This process aligns raw reads to previously indexed genomes.
It uses the names of the read files to determine if it should run in paired-end mode or not.
Input is somewhat complicated. It takes a tuple of a directory with SuperPang output, the index files for the pangenome,
a string with the pangenome name, a string with the name of the sample and the raw reads for that sample.
Output is the coverage information of how well the reads mapped to the genome.
*/
process map_subset {
    label "low_cpu"
    tag "low_cpu"
    input:
    tuple(path(pangenome_dir),path(index), val(pang_id), val(sample_ID), path(sub_reads)) //use the combine operator on the channels in the workflow.
    output:
    path("*_coverage.tsv", emit: coverage)
    shell:
    '''
    #run bowtie2
    reads_id=$(basename sub_*_R1.fq.gz _R1.fq.gz)
    #check if there are sub*R2 reads, if yes:
    if stat --printf='' sub_*_R2.fq.gz 2>/dev/null; then
        echo "Running paired-end mode"
        bowtie2 -x index -1 sub_*_R1.fq.gz -2 sub_*_R2.fq.gz --threads !{task.cpus} | samtools view -bS --threads !{task.cpus} > tmp_alignment.bam
    else
        echo "Running unpaired reads mode"
        bowtie2 -x index -U sub_*_R1.fq.gz --threads !{task.cpus} | samtools view -bS --threads !{task.cpus} > tmp_alignment.bam
    fi
    
    echo "Sorting bam files"
    samtools sort tmp_alignment.bam -O BAM -o !{pang_id}_${reads_id}_alignment.bam --threads !{task.cpus}
    
    echo "Computing coverage"
    #samtools coverage !{pang_id}_${reads_id}_alignment.bam -o !{pang_id}_${reads_id}_coverage.tsv
    #include orphans, max depth 1M, min base qual 15, output all positions
    samtools mpileup -A -d 1000000 -Q 15 -a !{pang_id}_${reads_id}_alignment.bam -o !{pang_id}_${reads_id}_coverage.tsv
    
    echo "Removing .bam files" #to save space
    rm *.bam #the sorted bams are named despite being deleted in case we decide a downstream process needs them.
    
    '''
}
