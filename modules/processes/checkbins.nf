process checkbins {
        publishDir "${project_path}/bins/fastas", mode: "copy", pattern: "${sample.baseName}/results/bins/*.fa", saveAs: { filename -> filename.split("/")[-1] }
        publishDir "${project_path}/bins/bintables", mode: "copy", pattern: "${sample.baseName}/results/18.*.bintable", saveAs: { "${sample.baseName}.bintable" }
        label "checkbins"
        label "high_mem"
        tag "${samp_name}"
        input:
        tuple(val(samp_name), path(sqm_dir))
        output:
        path("${sqm_dir}/results/bins/*.fa", emit: bins)
        path("${sqm_dir}/results/18.*.bintable", emit: bintable)
        script:
        """
        echo "Running checkbins with {$task.cpus} nr of cpus"
        sed -i 's/\($numthreads *= *\).*;/\1$task.cpus;/' ${sqm_dir}/SqueezeMeta_conf.pl
        17.checkbins.pl $sqm_dir
        """
}