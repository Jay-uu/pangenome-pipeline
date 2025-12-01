process dastool {
        publishDir "${project_path}/bins/fastas", mode: "copy", pattern: "${samp_name}/results/bins/*.fa", saveAs: { filename -> filename.split("/")[-1] }
        label "dastool"
        tag "${samp_name}"
        input:
        tuple(val(samp_name), path(sqm_dir))
        val(project_path)
        output:
        tuple(val("${samp_name}"),path("${sqm_dir}"), emit: sqm_dir)
        path("${sqm_dir}/results/bins/*.fa", emit: bins)
        script:
        """
        echo "Running dastool with {$task.cpus} nr of cpus"
        #sed -i 's#\\(\$numthreads *= *\\).*;#\\1'"$task.cpus;"'#' $sqm_dir/SqueezeMeta_conf.pl
        15.dastool.pl $sqm_dir $task.cpus
        """

}