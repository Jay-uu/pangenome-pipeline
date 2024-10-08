/*
This process calculates which pangenomes have enough samples that pass the expected average coverage threshold to create new .samples files
for the pangenomes that will then be used to align the reads of those samples to the pangenomes for further analysis. Input should be the 
collected output from map_subset coverage, single_samps and subsample_fastqs readcount, so that all samples are processed with one process run.
Input: 
       coverage: A file with samtools coverage output, which has the depth of how reads from a sample mapped to a pangenome/reference genome.
       singles: The previously generated .samples files with which reads belong to which sample.
       readcounts: A file named {sample}_readcount.txt which contains the total number of reads per sample.
       in those fastqs.
Output:
       pang_cpm_cov: Two files, one with the calculated results for CPM and the other with coverage.
       pang_samples: The new .samples files for the pangenomes.
       individual_pang_cov: CPM and coverage for individual mOTUs/pangenomes published within their own dirs.
*/
process cov_to_pang_samples {
    label "low_cpu"
    label "cov_to_pang_samples"
    tag "All_mOTUs"
    publishDir "${params.project}/mOTUs", mode: "copy", pattern: "*.*.tsv", failOnError: false
    publishDir "${params.project}/mOTUs/results", mode: "copy", pattern: "pangenome/*.tsv", failOnError: false, saveAs: {"${file(it).getSimpleName()}/pangenome/${file(it).getBaseName()}"}
    publishDir "${params.project}/mOTUs/results", mode: "copy", pattern: "samples/*.samples",failOnError: false, saveAs: {"${file(it).getSimpleName()}/pangenome/${file(it).getSimpleName()}.samples"}
    input:
    path(coverage)
    path(samples_file)
    path(readcounts)
    output:
    path("*.*.tsv", emit: pang_cpm_cov)
    path("samples/*.samples", optional: true, emit: pang_samples)
    path("pangenome/*.tsv"), emit: individual_pang_cov
    path("NONE_PASSED.txt"), optional: true, emit: not_passed_message
    
    shell:
    $/
    #!/usr/bin/env python
    import os
    import pandas as pd
    import glob
    
    SAMPS_FILE = "!{samples_file}"
    COV_THRESHOLD = !{params.min_cov}
    NR_SAMPS_THRESHOLD = !{params.nr_samps_threshold}
    NR_SUBSAMP = !{params.nr_subsamp}
    PROJECT = os.path.basename("!{params.project}".rstrip("/"))
    
    print(f"Project is !{params.project}")
    """
    Input: 
        cov = dataframe of samtools coverage output
        tot_reads = total count of reads from the sample read files
    Returns the weighted mean coverage per read, and the expected average
    coverage for all the reads of that sample to the pangenome.
    """
    def get_expected_cov_cpm(cov, tot_reads, nr_subsamp):
        #median = cov["Depth"].median()
        #CovPM = median*1000000/nr_subsamp
        mean = cov["Depth"].mean() #using mean to avoid getting 0
        CovPM = mean*1000000/nr_subsamp
        exp_cov = CovPM*tot_reads/1000000
        return CovPM, exp_cov
    
    count_files = glob.glob("*_readcounts.tsv")
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
        samp_name = file.split("sub_",1)[1].split("_coverage.tsv")[0]
        header = ["Name", "Position", "Depth"] #samtools depth format
        cov = pd.read_csv(file, sep="\t", names=header)
        tot_reads = readcount[readcount["Sample"]==samp_name]["Total_reads"].item()
        if NR_SUBSAMP > tot_reads:
            print(f"Requested subsampling more than available reads for {samp_name}. Total reads: {tot_reads} used for subsampling and cov estimation instead.")
            nr_subsamp = tot_reads
        else:
            nr_subsamp = NR_SUBSAMP
        cpm, exp_cov = get_expected_cov_cpm(cov, tot_reads, nr_subsamp)
        print(f"Samp: {samp_name}, Median coverage: {exp_cov}, CPM: {cpm}")
        cpm_dic.setdefault(pang_id, {})[samp_name] = cpm
        cov_dic.setdefault(pang_id, {})[samp_name] = exp_cov
    
    #convert dicts to dfs and save to file
    print("Saving coverage and CPM to files.")
    all_cov = pd.DataFrame.from_dict(cov_dic, orient="index")
    all_cov = all_cov.reset_index().rename(columns={"index": "Pangenome"})
    all_cov.to_csv(f"{PROJECT}.cov.tsv", sep = '\t', index=False)
    
    all_cpm = pd.DataFrame.from_dict(cpm_dic, orient="index")
    all_cpm = all_cpm.reset_index().rename(columns={"index": "Pangenome"})
    all_cpm.to_csv(f"{PROJECT}.cpm.tsv", sep = '\t', index=False)
    
    print("Making individual cov and cpm files.")
    os.makedirs("pangenome")
    for pang in all_cov["Pangenome"].unique():
        motu = pang.split(".", 1)[0]
        all_cov[all_cov["Pangenome"] == pang].to_csv(f"pangenome/{motu}.cov.tsv", sep="\t", index = False)
        all_cpm[all_cpm["Pangenome"] == pang].to_csv(f"pangenome/{motu}.cpm.tsv", sep="\t", index = False)

    os.makedirs("samples")
    #create new .samples files for pangenomes that fit the coverage criteria
    print("Creating samples files for pangenomes that pass the thresholds.")
    samples_df = pd.read_csv(SAMPS_FILE, sep='\t', names=["sample","read", "pair"])
    for pang_id in cov_dic.keys():
        ovr_thresh = []
        for sample in cov_dic[pang_id].keys():
            #print(f"Checking {sample} coverage on {pang_id}")
            if cov_dic[pang_id][sample].item() >= COV_THRESHOLD:
                ovr_thresh.append(sample)
        if len(ovr_thresh) >= NR_SAMPS_THRESHOLD:
            print(f"Creating new samples file for {pang_id}")
            new_samp_df = samples_df.query('sample in @ovr_thresh')
            new_samp_df.to_csv(f"samples/{pang_id}.samples", header=False, index=None, sep='\t')

    #This allows the process to finish and publish results, but still printing why the pipeline stops if no pangenomes pass the thresholds
    if len(glob.glob(f"samples/*.samples")) < 1:
        with open("NONE_PASSED.txt", "w") as outfile:
            outfile.write("WARNING: It seems none of your pangenomes fulfill the thresholds for further analysis. Consider lowering --min_cov and/or --nr_samps_threshold, increasing how many reads are subsampled or perhaps using more samples.")

   
    /$
}
