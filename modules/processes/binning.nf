process binning {
        label "binning"
        tag "${samp_name}"
        input:
        tuple(val(samp_name), path(sqm_dir))
        output:
        tuple(val("${samp_name}"),path("${sqm_dir}"), emit: sqm_dir)
        script:
        """
        sed -i 's/\($numthreads *= *\).*;/\1$task.cpus;/' ${sqm_dir}/SqueezeMeta_conf.pl
        14.runbinning.pl $sqm_dir
        """
}