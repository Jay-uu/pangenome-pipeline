process map_samples {
        label "high_mem"
        tag "${samp_name}"
        input:
        tuple(val(samp_name), path(sqm_dir))
        output:
        tuple(val("${samp_name}"),path("${sqm_dir}"), emit: sqm_dir)
        script:
        """
        echo "Mapping samples to assembly $samp_name, using $task.cpus threads"
        10.mapsamples.pl $sqm_dir 0 $task.cpus
        """
}