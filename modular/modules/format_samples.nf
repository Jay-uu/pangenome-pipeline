/*
This process takes a tab-delimited samples file and converts it to individual files per sample.
It also checks that the provided fastq-dir has the files specified in the samples file.
Output is a list of all the individual samples files.
*/
process format_samples {
    label "low_cpu"
    tag "low_cpu"
    input:
    path(samples_file)
    path(fastq_dir)
    output:
    path("singles/*")
    shell:
    $/
    #!/usr/bin/env python
    import os
    from collections import defaultdict
    from collections import Counter
    
    os.makedirs("singles", exist_ok=True)
    fastq_files = os.listdir("!{fastq_dir}")
    
    samples2lines = defaultdict(list)
    with open("!{samples_file}") as infile:
        for line in open("!{samples_file}"):
            fields = line.strip().split("\t")
            if len (fields) < 3:
                raise Exception(f"Missing columns or wrong delimiter on line: {line}") # this halts execution with exit code -1 and emits a message
            sample, filename, pair, *_ = fields
            if pair not in ("pair1", "pair2"):
                raise Exception("Error with provided samples file. Make sure there's no header and that the third column says 'pair1' or 'pair2'") #a header, or spelled wrong
            if filename not in fastq_files:
                 raise Exception(f"{filename} not found in !{fastq_dir}")
            if not filename.endswith(".fastq.gz") and not filename.endswith(".fq.gz"):
                 raise Exception(f"Wrong file format for {filename}. Only .fastq.gz and .fq.gz are accepted.")
            samples2lines[sample].append(line)
            
    paired_flag = set()
    for sample, lines in samples2lines.items():
        pairs = []
        for row in lines:
            pairs.append(row.strip().split("\t")[2])
        count = Counter(pairs)
        if count["pair1"] == count["pair2"]:
            paired_flag.add("paired")
        elif "pair2" not in count:
            paired_flag.add("unpaired")
        else:
            raise Exception(f"Sample {sample} has a paired file mismatch in the samples file.\n Make sure that there are only single reads or only paired end reads.")
        with open(f"singles/{sample}.samples", "w") as outfile:
            for line in lines:
                outfile.write(line)
    /$       
}
