/*
This process calculates which pangenomes have enough samples that pass the expected average coverage threshold to create new .samples files
for the pangenomes that will then be used to align the reads of those samples to the pangenomes for further analysis. Input should be the 
collected output from map_subset coverage, single_samps and subsample_fastqs readcount, so that all samples are processed with one process run.
Input: 
       coverage: A file with samtools coverage output, which has the depth of how reads from a sample mapped to a pangenome/reference genome.
       singles: The previously generated .samples files with which reads belong to which sample.
       readcounts: A file named {sample}_readcount.txt which conatains how many fastq files with reads belong to the sample, and the total number of reads 
       in those fastqs.
Output:
       pang_cpm_cov: Two files, one with the calculated results for CPM and the other with coverage.
       pang_samples: The new .samples files for the pangenomes.
*/
process cov_to_pang_samples {
    label "low_cpu"
    tag "low_cpu"
    publishDir "${params.project}/pangenomes", mode: "copy"
    input:
    path(coverage)
    path(samples_file)
    path(readcounts)
    output:
    path("samples/*.tsv", emit: pang_cpm_cov)
    path("samples/*.samples", emit: pang_samples)
    shell:
    $/
    #!/usr/bin/env python
    import os
    import pandas as pd
    import glob
    
    samps_file = "!{samples_file}"
    cov_threshold = !{params.mean_cov_threshold}
    nr_samps_threshold = !{params.nr_samps_threshold}
    outdir = "samples"
    
    """
    Input: 
        cov = dataframe of samtools coverage output
        nr_fqs = nr of read files to the sample
        tot_reads = total count of reads from the sample read files
    Returns the weighted mean coverage per read, and the expected average
    coverage for all the reads of that sample to the pangenome.
    """
    def get_weighted_mean(cov, nr_fqs, tot_reads):
        #meandepth is the mean depth of coverage over that contig
        #endpos-startpos=contig length
        cov["contig_length"] = cov.apply(lambda row : (row.endpos - row.startpos)+1, axis=1)
        cov["totaldepth"] = cov.apply(lambda row : row.meandepth*row.contig_length, axis=1) #rounds to nearest int
        weigh_mean=(cov["totaldepth"].sum()/cov["contig_length"].sum())
        #expected average coverage: ((average/million reads)/million*) total number of reads per sample
        weigh_mean_per_read = (weigh_mean/(1000000*nr_fqs))
        exp_cov = weigh_mean_per_read*tot_reads
        return weigh_mean_per_read, exp_cov
    
    count_files = glob.glob("*_readcount.txt")
    readcount = pd.DataFrame()
    for file in count_files:
        print(f"Appending {file}")
        readcount = pd.concat((readcount, pd.read_csv(file, sep='\t')), ignore_index=True)

    cpm_dic = {}
    cov_dic = {}
    cov_files = glob.glob("*_coverage.tsv")
    for file in cov_files:
        print(f"Reading {file}")
        pang_id = file.split("_sub_",1)[0]
        samp_name = file.split("sub_",1)[1].split("_")[0]
        cov = pd.read_csv(file, sep="\t")
        nr_fqs = readcount[readcount["Sample"]==samp_name]["Nr_fastqs"].item()
        tot_reads = readcount[readcount["Sample"]==samp_name]["Total_reads"].item()
        cpm, exp_cov = get_weighted_mean(cov, nr_fqs, tot_reads)
        cpm_dic.setdefault(pang_id, {})[samp_name] = cpm
        cov_dic.setdefault(pang_id, {})[samp_name] = exp_cov
        
        
    #convert dicts to dfs and save to file
    print("Saving coverage and CPM to files.")
    os.makedirs(outdir)
    all_cov = pd.DataFrame.from_dict(cov_dic, orient="index")
    all_cov = all_cov.reset_index().rename(columns={"index": "Pangenome"})
    all_cov.to_csv(f"{outdir}/pangenomes_cov.tsv", sep = '\t', index=False)
    
    all_cpm = pd.DataFrame.from_dict(cpm_dic, orient="index")
    all_cpm = all_cpm.reset_index().rename(columns={"index": "Pangenome"})
    all_cov.to_csv(f"{outdir}/pangenomes_cpm.tsv", sep = '\t', index=False)

    #create new .samples files for pangenomes that fit the coverage criteria
    print("Creating samples files for pangenomes that pass the thresholds.")
    samples_df = pd.read_csv(samps_file, sep='\t', names=["sample","read", "pair"])
    for pang_id in cov_dic.keys():
        ovr_thresh = []
        for sample in cov_dic[pang_id].keys():
            #print(f"Checking {sample} coverage on {pang_id}")
            if cov_dic[pang_id][sample].item() >= cov_threshold:
                ovr_thresh.append(sample)
        if len(ovr_thresh) >= nr_samps_threshold:
            print(f"Creating new samples file for {pang_id}")
            new_samp_df = samples_df.query('sample in @ovr_thresh')
            new_samp_df.to_csv(f"{outdir}/{pang_id}.samples", header=False, index=None, sep='\t')
            
    if len(glob.glob(f"{outdir}/*.samples")) < 1:
        raise Exception("It seems none of your pangenomes fulfill the thresholds for further analysis. Consider lowering --mean_cov_threshold and/or --nr_samps_threshold, or perhaps using more samples.")
   
    /$
}
