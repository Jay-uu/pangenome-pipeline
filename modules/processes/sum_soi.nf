/* Summarise species of interest, aka species that have enough coverage in enough samples to be used for variant calling-
*/
process sum_soi {
    label "low_cpu"
    publishDir "${project_path}", mode: "copy"
    input:
    val(project_path)
    path(passed_files)
    val(min_cov)
    val(nr_samps_threshold)
    output:
    path("species_of_interest.txt"), emit: soi_file
    script:
    """
    #concat all files
    echo "#These species/pangenomes passed the thresholds of median coverage over ${min_cov} in ${nr_samps_threshold} different samples and are candidates for variant calling:" > species_of_interest.txt
    cat *_PASSED.txt >> species_of_interest.txt
    """
}