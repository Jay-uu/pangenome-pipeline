process checkbins {
        publishDir "${project_path}/bins/fastas", mode: "copy", pattern: "${samp_name}/results/bins/*.fa", saveAs: { filename -> filename.split("/")[-1] }
        publishDir "${project_path}/bins/bintables", mode: "copy", pattern: "${samp_name}/results/18.*.bintable", saveAs: { "${samp_name}.bintable" }
        label "checkbins"
        label "high_mem"
        tag "${samp_name}"
        input:
        tuple(val(samp_name), path(sqm_dir))
        val(project_path)
        output:
        path("${sqm_dir}/results/bins/*.fa", emit: bins)
        path("${sqm_dir}/results/18.*.bintable", emit: bintable)
        script:
        """
        echo "Running checkbins and getbins with {$task.cpus} nr of cpus"
        #sed -i 's#\\(\$numthreads *= *\\).*;#\\1'"$task.cpus;"'#' $sqm_dir/SqueezeMeta_conf.pl
        17.checkbins.pl $sqm_dir $task.cpus
        18.getbins.pl $sqm_dir
        """
}