process dastool {
        label "dastool"
        tag "${samp_name}"
        input:
        tuple(val(samp_name), path(sqm_dir))
        output:
        tuple(val("${samp_name}"),path("${sqm_dir}"), emit: sqm_dir)
        script:
        """
        echo "Running dastool with {$task.cpus} nr of cpus"
        #sed -i 's#\\(\$numthreads *= *\\).*;#\\1'"$task.cpus;"'#' $sqm_dir/SqueezeMeta_conf.pl
        15.dastool.pl $sqm_dir $task.cpus
        """

}