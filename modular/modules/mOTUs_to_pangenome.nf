/*
Creates pangenomes from a directory with bins and a tsv with bin completeness and contamination by running SuperPang.
If there's not enough genomes to a mOTU, it will just copy the genome to results.
Input is a tuple containing a directory with bins (preferably belonging to the same mOTU) and a bintable with completeness and contamination.
Output is a directory with a subdirectory for the mOTU, containing the pangenome fasta file and the NBPs.fasta in an individual channel.
This could possibly be changed to have different script parts, which would mean that I can have the python code for changing the name of the file and the fasta headers (the else part) directly in the process instead of in a separate script.
*/
process mOTUs_to_pangenome {
    publishDir "${params.project}/mOTUs/results/${mOTU_dir}/pangenome", mode: "copyNoFollow",  pattern: "pangenomes/${mOTU_dir}", saveAs: { filename -> "superpang/" }
    publishDir "${params.project}/mOTUs/results/${mOTU_dir}/pangenome/superpang/", mode: "copy", pattern: "*contigs.tsv"
    publishDir "${params.project}/mOTUs/results/${mOTU_dir}/", mode: "copy", pattern: "input_bins.txt", saveAs: { filename -> "pang_bins.txt" }
    label "mOTUs_to_pangenome"
    tag "${mOTU_dir.baseName}"
    input:
    tuple(path(mOTU_dir), path(bintable))
    output:
    path("pangenomes/${mOTU_dir}", type: "dir", emit: pangenome_dir)
    path("pangenomes/${mOTU_dir}/*.NBPs.fasta", emit: NBPs_fasta)
    path("pangenomes/${mOTU_dir}/*.core.fasta", emit: core_fasta)
    path("${mOTU_dir}.core.contigs.tsv"), emit: contigs_tsv
    path("input_bins.txt"), emit: input_bins
    shell:
    $/
    
    #!/usr/bin/env python
    from Bio import SeqIO
    from subprocess import call
    from pathlib import Path
    import os
    import glob
    
    genomes = glob.glob("!{mOTU_dir}/*.fa")
    nr_genomes = len(genomes)
    pg_dir_name = "pangenomes"
    os.makedirs(pg_dir_name)
    core_name = f"{pg_dir_name}/!{mOTU_dir}/" + "!{mOTU_dir}" + ".NBPs.core.fasta"
    
    #Making a file with the bins used for creating the pangenome
    with open("input_bins.txt", "w") as fastas:
        fastas.write("\n".join(genomes))
        
    if nr_genomes > 1:
        print("Enough genomes to run pangenome computation")
        #with open("input_bins.txt", "w") as fastas:
            #fastas.write("\n".join(genomes))
        call(["SuperPang.py", "--fasta", "input_bins.txt","--checkm", "!{bintable}", "--output-dir", f"{pg_dir_name}/!{mOTU_dir}", "--header-prefix", f"!{mOTU_dir}",
        "--output-as-file-prefix", "--nice-headers", "--threads", f"!{task.cpus}"])
        #Check if core file is empty
        if os.stat(core_name).st_size == 0:
            print("Core genome file empty. Will use the consensus assembly for read mapping.")
            print("This shouldn't happen unless you're using mock communities.")
            #remove the old core fasta, since I can't output python variables and this makes the output flexible
            Path.unlink(Path(core_name))
            #make consensus.core.fasta as symlink
            core_name = f"{pg_dir_name}/!{mOTU_dir}/" + "!{mOTU_dir}" + ".consensus.core.fasta"
            Path(core_name).symlink_to("!{mOTU_dir}" + ".NBPs.fasta")
    
    elif nr_genomes == 1:
        print("Only one genome in mOTU. Renaming headers and copying to pangenome dir.")
        outfile= f"{pg_dir_name}/!{mOTU_dir}/" + "!{mOTU_dir}" + ".singlemOTU.core.fasta"
        core_name = outfile
        symfile = f"{pg_dir_name}/!{mOTU_dir}/" + "!{mOTU_dir}" + ".singlemOTU.NBPs.fasta"
        pangenome_file=f"{genomes[0]}"
        os.makedirs(pg_dir_name+"/!{mOTU_dir}")
        with open(pangenome_file, "rt") as pg:
            all_contigs = list(SeqIO.parse(pg, "fasta")) 
            pg_name = pangenome_file.split('.',1)[0]
            for i,s in enumerate(all_contigs):
                s.id =f"{pg_name}_NODE_Sc{i}-0-noinfo_length_{len(s.seq)}_cov_1.tag_noinfo" #renames contigs
                s.description="" #remove the old name
            with open(outfile, "w") as output_handle:
                SeqIO.write(all_contigs, output_handle, "fasta")
        #symlinking so that it can be sent as NBPs_fasta output
        Path(symfile).symlink_to("!{mOTU_dir}" + ".singlemOTU.core.fasta")
    else:
        raise Exception(f"No fastas found in !{mOTU_dir}")
        
    #Getting the core contig names over a certain length in a tsv for the VC workflow
    print("Creating contigs.tsv")
    with open(core_name, "rt") as pg:
        all_contigs = list(SeqIO.parse(pg, "fasta"))
        with open(f"!{mOTU_dir}.core.contigs.tsv", "w") as outfile:
            for i,s in enumerate(all_contigs):
                c_len = len(s.seq)
                if c_len >= !{params.min_contig_len}:
                    outfile.write(f"{s.id}\t{c_len}\n")       
    
         
    /$
}
