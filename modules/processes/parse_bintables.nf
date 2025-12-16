/*
Checks for files with the .bintable extension, reads them and checks if there's checkM2 results for all bins provided.
}
*/
process parse_bintables {
    label "low_cpu"
    tag "All_bins"
    input:
    tuple(val(samp_name), path(sqm_dir))
    path(bintables_dir)
    output:
    tuple(val("${samp_name}"),path("to_check"), optional: true, emit: to_checkbins)
    script:
    """
    #!/usr/bin/env python3
    #go through each bintable and see which bins have results, then compare. 
    #bins w/o bintable gets sent to checkbins
    import glob as glob
    import os
    import shutil

    SQMDIR = "${sqm_dir}"
    BINTABLESDIR = "${bintables_dir}"
    BINSDIR = SQMDIR+"/results/bins/"

    #get bin names
    all_bins = glob.glob(BINSDIR+"*.fa")
    all_bintables = glob.glob(BINTABLESDIR+"/*.bintable")
    for bintable in all_bintables:
        print("Reading: ", bintable)
        if len(all_bins) > 0:
            bins_list = [s.split('\t', 1)[0] for s in open(bintable,"r")]
            [all_bins.remove(BINSDIR+bin_id+".fa") for bin_id in bins_list if BINSDIR+bin_id+".fa" in all_bins]
        else:
            print("All bins found in bintables. No need to run checkM2")
    
    if len(all_bins) > 0: #if bins still in list they need to be checked
        print("Found bins that lack checkM2 results. Will send to checkbins.")
        os.makedirs("to_check/results/bins")
        #copy necessary files to simulate a proper sqm_dir
        os.makedirs("to_check/intermediate")
        shutil.copy(SQMDIR+"/SqueezeMeta_conf.pl", "to_check/SqueezeMeta_conf.pl")
        shutil.copy(SQMDIR+"/parameters.pl", "to_check/parameters.pl")
        for bin_id in all_bins:
            bin_name = os.path.split(bin_id)[1]
            print(f"Copying {bin_name} to output dir for checkbins.")
            shutil.copy(bin_id, "to_check/results/bins/"+bin_name)

    """
}