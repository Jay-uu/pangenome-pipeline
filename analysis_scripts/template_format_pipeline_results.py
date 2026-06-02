# Needed input: path to project directory and where to store results.
#Pipeline needs to have ran step3 so the file completed_pogs.txt exist
PROJPATH=""
POGS=PROJPATH + "/completed_pogs.txt" #file where each line is the path from PROJPATH to a pangenome.samples file - aka only works for pangenomes that passed the coverage and nr_samples thresholds for VC.
#example: mOTUs/results/<pang_name>/pangenome/<pang_name>.samples
#for now it can be created like so in PROJPATH:
# ls mOTUs/results/*/pangenome/*.samples > good_cov_mOTUs.txt
OUTDIR=PROJPATH + "/analysis"
ALL_POGS = True #choose wether to filter based on completeness and contam or not. True means no filtering.
COMPLETENESS=90
CONTAMINATION=0

#1 Create a dict with genome size, completeness, contamination, and cpm per sample, for each POG
#2 For each pog * sample combo, calculate the percentage of reads recruited
#3 Add full GTDB-Tk taxonomy for each POG
#4 write to outdir: one file with all POGs, one with only high quality ones.

# import libraries
import glob
import pandas as pd
import matplotlib
import pathlib

#create outdir if it doesn't exist
pathlib.Path(OUTDIR).mkdir(parents=True, exist_ok=True) 

# Get a list with all the paths to the POGs.
#open the POGS file
print("Reading list of completed POGs")
pog_list = []
with open (POGS, "r") as pogfile:
    for line in pogfile:
        motus,results,pog,*_ = line.strip().split("/") #want to not include linebreaks
        pog_list.append(pog)
        #print(line.split("/"))
        print(pog)
    #format: mOTUs/results/g__BOG-931_mOTU_4/pangenome/pogenom/results\n
    #let's only keep the pog name here.

print("Number of pangenomes for further analysis: ", len(pog_list))

# Then go trhough the CheckM2 results
print("Accessing completeness, contamination and genome size for each POG.\nAlso checking the CPM from subsampled reads.")
pog_dict = {}
high_qual_pogs = []
for pog in pog_list:
    #genome name is pog
    #get genome size
    core_len = "0"
    pog_dict.update({pog: {}})
    with open(PROJPATH + "/mOTUs/results/" + pog + "/pangenome/" + pog + "_cM2_summary.txt", "r") as file:
    # Iterate over each line
        for line in file:
            Name, Completeness, Contamination, Completeness_Model_Used, Translation_Table_Used, Coding_Density, Contig_N50, Average_Gene_Length, Genome_Size, *_ = line.split("\t")
            #print(Name, ": ", Genome_Size)
            if Name.rsplit(".")[-1] == "core":
                pog_dict[pog].update({"Genome_Size":int(Genome_Size)})
                pog_dict[pog].update({"Completeness":float(Completeness)})
                pog_dict[pog].update({"Contamination":float(Contamination)})
                if float(Completeness)>70 and float(Contamination)<10:
                    high_qual_pogs.append(pog)
                #print(f"Core length for {pog} found: {Genome_Size}")
    #now get sample cpms
    with open(PROJPATH + "/mOTUs/results/" + pog + "/pangenome/" + pog + ".cpm", "r") as file:
        allLines = file.readlines()
    samples = allLines[0].split(sep="\t")[1:]
    cpm_values = allLines[1].split(sep="\t")[1:]
    if (len(allLines[0].split(sep="\t"))-1 == len(samples)):
        for i in range(0,len(samples)):
            samp = samples[i].rstrip()
            #ADD A CHECK THAT CPM IS A NUMBER HERE
            if (cpm.replace('.','',1).isdigit()):
                cpm = float(cpm_values[i].rstrip())
            else:
                cpm = 0
            #print(f"Sample is: {samp} and cpm is {cpm}")
            pog_dict[pog].update({samp:cpm})
    else:
        print(f"Warning! Double check input for {pog} cpm values.")

print("Comparing dict keys to number of pogs (should be equal): ", len(pog_dict), len(pog_list))

print("Calculating read recruitment based on CPM for subsampled reads.")
#define samples from .samples file instead?

readsprcnt = {}
for mOTU in pog_list:
    readsprcnt.update({mOTU: {}})
    genome_size = pog_dict[mOTU]["Genome_Size"]
    for sample in samples:
         #transform sample values according to formula
         sample = sample.rstrip()
         cpm = pog_dict[mOTU][sample]
         rpc = cpm*genome_size/150/10000
         readsprcnt[mOTU].update({sample:rpc})
         #percentage of reads: cpm*genome_size/150/10000

print("Comparing dict keys to number of pogs: ", len(pog_dict), len(pog_list))
#Get percentage of reads mapped per genome for all samps, and for all average per sample for all genomes?
nr_genomes = 0
tot_avg_bcp_all_genomes = 0
agg_bcp = 0
for key in pog_dict.keys():
    nr_genomes = nr_genomes + 1
    gsize = pog_dict[key]["Genome_Size"]
    #print(key," genome size is: ", gsize)
    tot_bcp = 0
    tot_samps = 0
    for samp in pog_dict[key].keys():
        #check it's not genome size
        if samp not in ["Genome_Size", "Completeness", "Contamination"]:
            #(cpm*genome_size)/150aka read size = bases_cov/tot_bases_mapped
            cpm = pog_dict[key][samp]
            tot_samps = tot_samps + 1
            #samp_bcp = (cpm*gsize/150)/1000000*100 #percentage of bases covered
            samp_bcp = cpm*gsize/150/1000000 #percentage of reads mapped per million reads
            tot_bcp = tot_bcp + samp_bcp
            #print(f"For genome {key} and sample {samp} the percentage of reads that mapped are: {samp_bcp}")
    #print(f"The total times samples were mapped: {tot_samps} to {key}")
    #print(f"The total bcp for {key} is {tot_bcp}")
    #print(f"Average bcp: {tot_bcp/tot_samps} percentage of bases were mapped to {key}")
    tot_avg_bcp_all_genomes = tot_bcp/tot_samps
    agg_bcp = agg_bcp + tot_bcp
print(f"The average percentage of mapped reads for all genomes with all samples is: {tot_avg_bcp_all_genomes/nr_genomes}")
print(f"While the aggregated percentages of mapped reads to the genomes is: {agg_bcp}")

print("Checking the high quality POGs.")
print("Number of high quality POGs: ",len(high_qual_pogs))
#verify some inputs
for pog in high_qual_pogs:
    comp = pog_dict[pog]["Completeness"]
    cont = pog_dict[pog]["Contamination"]
    print(f"Pog: {pog}, Completeness: {comp}, Contamination: {cont}")

if (ALL_POGS):
    print("Using all POGs with pogenom results, regardless of genome quality.")
    select_pogs = pog_list
elif len(high_qual_pogs) > 0:
    print("Continue analysis with only high quality POGs.")
    select_pogs = high_qual_pogs
else:
    print("Warning! No high quality POGs. Continuing analysis with all POGs.")

#Add taxonomy lineage
sum_bintable = pd.read_csv(PROJPATH+"/bins/summarized_bins.bintable", sep="\t")
#col of interest "Tax GTDB-Tk"
tax = sum_bintable["Tax GTDB-Tk"].unique()
for pog in pog_dict.keys():
    species = pog.split("_")[2]
    lineage = next((s for s in tax if species in s), None)
    #print("pog", pog)
    #print("species:", species)
    #print("lineage", lineage)
    pog_dict[pog].update({"Tax GTDB-Tk":lineage})


print("Saving results in long format for further processing.")

#write a file with the selected pogs?
with open(OUTDIR+"/select_pogs_read_recruitment.tsv", "w") as outfile:
    outfile.write("Sample\tmOTU\tprcnt_mapped\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n")
    for pog in select_pogs:
        for samp in pog_dict[pog].keys():
            if samp not in ['Genome_Size', 'Completeness', 'Contamination', 'Tax GTDB-Tk']: #only want sample names
                val = pog_dict[pog][samp]
                lineage = pog_dict[pog]["Tax GTDB-Tk"]
                Domain, Phylum, Class, Order, Family, Genus, Species = lineage.split(";")
                tax_levels = [Domain, Phylum, Class, Order, Family, Genus, Species]
                #print(tax_levels)
                for i in range(0,len(tax_levels)):
                    if len(tax_levels[i]) == 3: #this means it just has the level signature
                        tax_levels[i] = tax_levels[i]+"Unclassified"
                outfile.write(f"{samp}\t{pog}\t{val}\t"+"\t".join(tax_levels)+"\n")

print(f"File for further analysis created at {OUTDIR}/select_pogs_read_recruitment.tsv")
