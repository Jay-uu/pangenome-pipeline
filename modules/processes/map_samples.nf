process map_samples {
        label "map_to_asm"
        tag "${samp_name}"
        input:
        tuple(val(samp_name), path(sqm_dir))
        output:
        tuple(val("${samp_name}"),path("${sqm_dir}"), emit: sqm_dir)
        script:
        """
        #sqm inherits the settings from the previous process
        #this sed line changes those to the settings for this process
        sed -i 's/\($numthreads *= *\).*;/\1$task.cpus;/' ${sqm_dir}/SqueezeMeta_conf.pl 
        10.mapsamples.pl $sqm_dir
        """
}