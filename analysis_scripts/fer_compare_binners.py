#!/usr/bin/env python3

from os import listdir
from glob import glob
from fastani import fastani

# These variables need to be changed for your purpose:
METHODS = {'metabat2': '/home/jay/data/mock_20240409_metabat', 'all_binners': '/home/jay/data/mock_20240409_all_binners'}
WORKDIR = '/home/jay/data/tmp_binners_eval' #output ends up here
ANI_THRESHOLD = 95

results = {} # This will be a nested dictionary
# LVL1: method
# LVL2: pangenomes within a method
# LVL3: info about each pangenome within a method
# eg {'metabat2' : {'mOTU_00': {'core_completeness': 90, 'core_contamination': 5, 'pangenome_size': 2000000}}}

### Load data
for method, path in METHODS.items():

    results[method] = {}
    mOTUs_path = f'{path}/mOTUs'
    pang_path =  f'{mOTUs_path}/pangenomes'

    # Iterate over each mOTU for this method
    for mOTU in listdir(pang_path):

        results[method][mOTU] = {}
        
        # Get pangenome size for this mOTU, also store the path to the core.NBPs file for calculating ANI later
        NBPs_all_path  = glob(f'{pang_path}/{mOTU}/*NBPs.fasta')[0]
        NBPs_core_path = glob(f'{pang_path}/{mOTU}/*core.fasta')[0]
        results[method][mOTU]['NBPs_core_path'] = NBPs_core_path
        results[method][mOTU]['pangenome_size'] = sum(len(line.strip()) for line in open(NBPs_all_path) if not line.startswith('>'))

        checkm_path = f'{mOTUs_path}/{mOTU}/checkm/{mOTU}_cM2_summary.txt'
        with open(checkm_path) as infile:
            infile.readline() # burn headers
            for line in infile:
                name, completeness, contamination, *_ = line.strip().split('\t')
                if name.endswith('core'):
                    results[method][mOTU]['core_completeness' ] = float(completeness )
                    results[method][mOTU]['core_contamination'] = float(contamination)


### Calculate ANI between pangenomes from both methods

# Initialize ANIs to 0
for method in results:
    for mOTU in results[method]: 
        results[method][mOTU]['best_match']       = ''
        results[method][mOTU]['best_match_ANI']   = 0

# Calculate ANIs
for mOTU_all, info_all in results['all_binners'].items():
    fasta_all = info_all['NBPs_core_path']
    tax_all = mOTU_all.split('_')[2]
    for mOTU_mb2, info_mb2 in results['metabat2'].items():
        tax_mb2 = mOTU_mb2.split('_')[2]
        if tax_all != tax_mb2:
            continue
        fasta_mb2 = info_mb2['NBPs_core_path']
        # Use the all_binners mOTU as the reference, assume ANIs are symmetrical
        fastANI_result = fastani(query=fasta_mb2, reference=fasta_all).as_dict()[fasta_mb2][fasta_all]
        if not fastANI_result: # No hit was found
            continue
        ANI = fastANI_result.ani
        # If this ANI is better than the one recorded previously, update
        if ANI < ANI_THRESHOLD:
            continue
        if ANI > results['all_binners'][mOTU_all]['best_match_ANI']:
            results['all_binners'][mOTU_all]['best_match'    ] = mOTU_mb2
            results['all_binners'][mOTU_all]['best_match_ANI'] = ANI
        if ANI > results['metabat2']   [mOTU_mb2]['best_match_ANI']:
            results['metabat2']   [mOTU_mb2]['best_match'    ] = mOTU_all
            results['metabat2']   [mOTU_mb2]['best_match_ANI'] = ANI


### Report results
with open(f'{WORKDIR}/compare_pangenomes_output.tsv', 'w') as outfile:
    outfile.write('\t'.join(['mOTU_all_binners', 'core_completeness_all_binners', 'core_contamination_all_binners', 'pangenome_length_all_binners'\
                             'mOTU_metabat2', 'core_completeness_metabat2', 'core_contamination_metabat2', 'pangenome_length_metabat2' \
                             'ANI_all_binners_vs_metabat2']) + '\n')
    added_mb2 = set()
    for mOTU_all, info_all in results['all_binners'].items():
        core_comp_all, core_cont_all, pang_length_all = info_all['core_completeness'], info_all['core_contamination'], info_all['pangenome_size']
        mOTU_mb2, ANI = info_all['best_match'], info_all['best_match_ANI']
        if mOTU_mb2: # We got a match for this mOTU
            added_mb2.add(mOTU_mb2) # make a note that this metabat2 mOTU was paired to an all_binners_mOTU
            info_mb2 = results['metabat2'][mOTU_mb2]
            core_comp_mb2, core_cont_mb2, pang_length_mb2 = info_mb2['core_completeness'], info_mb2['core_contamination'], info_mb2['pangenome_size']
        else:
            mOTU_mb2, core_comp_mb2, core_cont_mb2, pang_length_mb2, ANI = ['NA'] * 5
        outfile.write('\t'.join(map(str, [mOTU_all, core_comp_all, core_cont_all, pang_length_all, \
                                          mOTU_mb2, core_comp_mb2, core_cont_mb2, pang_length_mb2, ANI])) + '\n')
    # Add the metabat2 mOTUs that were not paired with ANI
    for mOTU_mb2, info_mb2 in results['metabat2'].items():
        if mOTU_mb2 in added_mb2:
            continue
        core_comp_mb2, core_cont_mb2, pang_length_mb2 = info_mb2['core_completeness'], info_mb2['core_contamination'], info_mb2['pangenome_size']
        assert not info_mb2['best_match'] # There should not be a match here so we double check
        mOTU_all, core_comp_all, core_cont_all, pang_length_all, ANI = ['NA'] * 5
        outfile.write('\t'.join(map(str, [mOTU_all, core_comp_all, core_cont_all, pang_length_all, \
                                          mOTU_mb2, core_comp_mb2, core_cont_mb2, pang_length_mb2, ANI])) + '\n')
