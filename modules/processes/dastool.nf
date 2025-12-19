process dastool {
        tag "${samp_name}"
        input:
        tuple(val(samp_name), path(sqm_dir))
        output:
        path("${samp_name}/results/bins/*.fa", emit: bins)
        tuple(val("${samp_name}"),path("${sqm_dir}"), emit: sqm_dir)
        script:
        """
        echo "Running dastool with {$task.cpus} nr of cpus"
        15.dastool.pl $sqm_dir $task.cpus
        """

}