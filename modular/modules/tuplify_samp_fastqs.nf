/*
*/

process tuplify_samp_fastqs {
    label "low_cpu"
    tag "${sample.baseName}"
    input:
    path(sample)
    path(fastq_dir)
    output:
    tuple(val("${sample.baseName}"), path("*.*q.gz"), emit: reads)
    shell:
    $/
    #!/usr/bin/env python
    import os
    from pathlib import Path
    from shutil import copy

    FASTQ_FILES = os.listdir("!{fastq_dir}")
    SAMPLE_ID = Path("!{sample}").stem
    
    with open("!{sample}") as infile:
        for line in open("!{sample}"):
            fields = line.strip().split("\t")
            if len(fields) < 3:
                raise Exception(f"Missing columns or wrong delimiter on line: {line}")
            sample, filename, pair, *_ = fields
            if filename not in FASTQ_FILES:
                raise Exception(f"{filename} not found in !{fastq_dir}")
            #copy the link to the fastq into current workdir, so it can be output. They can remain symlinks to save space.
            copy(f"!{fastq_dir}/{filename}", ".", follow_symlinks=False)
            
    
    /$
}
