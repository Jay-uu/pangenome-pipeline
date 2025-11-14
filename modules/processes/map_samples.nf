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
        # NB! Since the SqueezeMeta_conf.pl is actually symlinked to the process where it was generated first (assemble_metagenome)
        # it actually modifies the original one. So the last run process settings are the ones that will be in the conf.pl
        echo "Adjusting config settings for $samp_name to use $task.cpus threads"
        #sed -i 's#\\(\$numthreads *= *\\).*;#\\1'"$task.cpus;"'#' $sqm_dir/SqueezeMeta_conf.pl
        echo "Mapping samples to assembly."
        10.mapsamples.pl $sqm_dir 0 $task.cpus
        """
}