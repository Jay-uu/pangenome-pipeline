process checkbins {
        label "binning"
        tag "${samp_name}"
        input:
        tuple(val(samp_name), path(sqm_dir))
        output:
        tuple(val("${samp_name}"),path("${sqm_dir}"), emit: sqm_dir)
        script:
        """
        echo "Running checkbins with {$task.cpus} nr of cpus"
        sed -i 's/\($numthreads *= *\).*;/\1$task.cpus;/' ${sqm_dir}/SqueezeMeta_conf.pl
        17.checkbins.pl $sqm_dir
        """
}