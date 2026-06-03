/*
This progress takes a collection of the output coverages and CPMs from cov_to_pang_samples,
and combines the results from each pangenome into one coverage file and one cpm file
*/
process concat_subsamp_cov {
    label "low_cpu"
    tag "All_pangs"
    publishDir "${project_path}/mOTUs", mode: "copy"
    input:
    val(project_path)
    path(individual_pang_cpm_cov)
    output:
    path("*.expected.*.tsv", emit: aggregated_cpm_cov)
    script:
    """
    #!/usr/bin/env python3
    import os
    import pandas as pd
    import glob

    PROJECT = os.path.basename("${project_path}".rstrip("/"))

    print("=====Process name: concat_subsamp_cov=====")
    #need to concat by column names since order might differ
    print("Appending input files to lists.")
    cov_files = glob.glob("*.cov.tsv")
    cpm_files = glob.glob("*.cpm.tsv")

    print("Sorting lists.")
    cov_files.sort()
    cpm_files.sort()

    assert len(cov_files) == len(cpm_files), "Different number of cov and CPM files."

    nr_files = len(cov_files)
    coverages = pd.DataFrame()
    cpms = pd.DataFrame()
    for progress in range(0,nr_files): #this is done differently. Fix.
        cov_i = cov_files[progress]
        cpm_i = cpm_files[progress]
        print(f"Reading files: {cov_i} and {cpm_i}")
        coverages = pd.concat((coverages, pd.read_csv(cov_i, sep='\t')), ignore_index=True)
        cpms = pd.concat((cpms, pd.read_csv(cpm_i, sep='\t')), ignore_index=True)
        progress = progress + 1
        print(f"Read {progress}/{nr_files} files.")

    print("Saving aggregated results to files.")
    coverages.to_csv(f"{PROJECT}.expected.cov.tsv", sep = '\t', index=False)
    cpms.to_csv(f"{PROJECT}.expected.cpm.tsv", sep = '\t', index=False)

    """

}