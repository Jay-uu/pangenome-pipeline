#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 13:58:49 2021

@author: chelsea

Modified by Jay. This script requires pandas v.1.9 to work.
"""

import pandas as pd
import numpy as np

#Path to SRA_data and creating dataframe
sra_data = '/data/jay/sra_res/sra_data.csv'
df_sra_data = pd.read_csv(sra_data, sep=',', index_col=1, low_memory=False)

#List of things to be removed from sra_data
if True:
    to_filter_out = ['sediment','virus', 'viral', 'virome', 'transcript', 'marine', 'filtered water', 'aquarium', '^genomic', 'seawater']
else:
    to_filter_out = []

#First filtering of unwanted information, such as non-metagenomes, non-bacteria, non-freshwater
for search_word in to_filter_out:
    df_sra_data = df_sra_data[~df_sra_data.stack().str.contains(search_word, case=False, regex = True).any(level=0)]

#Merge together all coordinate columns to one
coordinate_columns = ['lat_lon', 'latitude', 'longitude', 'geographic location (latitude)', 
                      'geographic location (longitude)', 'latitude start', 'Latitude End', 
                      'Longitude End', 'Latitude Start', 'Longitude Start', 
                      'Latitude and longitude', 'latitude and longitude' ] #create list of coordinate columns
coordinate_columns = [col for col in coordinate_columns if col in df_sra_data.columns]
df_sra_data[coordinate_columns] = df_sra_data[coordinate_columns].fillna('') #Fill NaN with empty string
df_sra_data['coordinates'] = df_sra_data[coordinate_columns].apply(lambda row: ' '.join(row.values.astype(str)), axis=1) #Merge coordinate columns
df_sra_data['coordinates'] = df_sra_data['coordinates'].str.strip() #remove space in the beginning of string

#Filter coordinate column 
df_sra_data = df_sra_data[~df_sra_data.stack().str.contains('°').any(level=0)] #Remove rows with degree sign
df_sra_data = df_sra_data[~df_sra_data.stack().str.contains('º').any(level=0)] #remove rows with ⁰
df_sra_data['coordinates'].replace('', np.nan, inplace = True) #replace empty strings with NaN
df_sra_data = df_sra_data[df_sra_data['coordinates'].notna()] #keep only rows without Nan
df_sra_data = df_sra_data[df_sra_data['coordinates'].str.contains('^[0-9]', regex = True)] #keep only coordinates starting with number

#Filter for freshwater accessions
freshwater_columns = ['taxon', 'Isolation source', 
                      'env_broad_scale', 'env_local_scale', 
                      'env_medium', 'environment (biome)', 
                      'environment (feature)', 'environment (material)', 
                      'env_biome', 'env_feature', 'env_material', 
                      'Environment (Feature)', 'Isolation source', 
                      'environment material','metagenomic-source', 
                      'organism', 'biome', 'environment material'] #list of columns relating to freshwater
if False:
    freshwater = df_sra_data[freshwater_columns].copy() #create df with only freshwater columns
    freshwater = freshwater.fillna('') #replace NaN with empty strings
    freshwater_filtered = freshwater[freshwater.stack().str.contains('fresh').any(level=0)] #keep freshwater rows
else:
    freshwater_filtered = df_sra_data.copy()

#Create dataframes of only accessions
accessions_freshwater = freshwater_filtered.index #create index object of accessions
accessions_freshwater = accessions_freshwater.to_frame(index = False) #create df of accessions
accessions_freshwater = accessions_freshwater.set_index('SRA_ID') #set index to SRA_ID
# accessions = df_sra_data.index 
# accessions = accessions.to_frame(index=False)
# accessions = accessions.set_index('SRA_ID')

#Merge freshwater accessions with sra_id to find the intersection
merge = pd.merge(accessions_freshwater, df_sra_data, how="inner", left_index=True, right_index=True) #Merge triplicates with original data

#Create dataframe with only accessions and coordinates from merged dataframe
id_coord = merge[['coordinates']].copy()

# Convert coordinates
lat = []
long = []
coordinates = []

for index, row in id_coord.iterrows():
    coordinates.append(str(row[0]).split(" "))

# Go through each of the coordinates and check which cardinal direction (letters) which will 
# determine if the value is positive or negative. 
for list in coordinates:
    if len(list) == 2:
        lat_value = str(list[0])
        lat.append(lat_value)
        long_value = str(list[1])
        long.append(long_value)
    else:
        if list[1] == 'S':
            lat_value = '-' + str(list[0])
            lat.append(lat_value)
        else: 
            lat.append(list[0])
        if list[3] == 'W':
            lon_value = '-' + str(list[2])
            long.append(lon_value)
        else: 
            long.append(list[2])

# Add the converted values to the dataset
id_coord['Latitude'] = lat 
id_coord['Longitude'] = long

#Replace commas with dots
id_coord['Latitude'] = id_coord['Latitude'].str.replace(',','.')
id_coord['Longitude'] = id_coord['Longitude'].str.replace(',', '.')

#Remove coordinates column
id_coord = id_coord.drop('coordinates', axis=1)

#Add special cases filtering.
#The lake Most, Czechia coordinates have are formatted as xx.xxx.xxx which causes issues. Need to remove the second dot.
lake_most = ['ERR3719164', 'ERR3719163', 'ERR3719162', 'ERR3719161', 'ERR3719160', 'ERR3719159', 'ERR3719158', 'ERR3719157', 'ERR3719156', 'ERR3719155',
        'ERR3719154', 'ERR3719153', 'ERR3719152', 'ERR3719151', 'ERR3719150', 'ERR3719149', 'ERR3719148', 'ERR3719147', 'ERR3719146', 'ERR3719145',
        'ERR3719144', 'ERR3719143', 'ERR3719142', 'ERR3719141', 'ERR3719140', 'ERR3719139', 'ERR3719138', 'ERR3719137', 'ERR9585962']
for lm in lake_most:
    id_coord.loc[lm]["Latitude"] = "".join(id_coord.loc[lm]["Latitude"].rsplit(".", 1))
    id_coord.loc[lm]["Longitude"] = "".join(id_coord.loc[lm]["Longitude"].rsplit(".", 1))

sowe_river = ['ERR10466844', 'ERR10466843', 'ERR10466842']
for sv in sowe_river:
    id_coord.loc[sv]["Latitude"] = id_coord.loc[sv]["Latitude"].split("DD")[0]
    id_coord.loc[sv]["Longitude"] = id_coord.loc[sv]["Longitude"].split("DD")[0]
"""
These have strange locations. Checking manually.
[0] seems to have been a mixup with coordinates, change longitude to East (positive value)
[1] marine
[2] is fine.
[3] mobile bay, varying salinity: https://www.usgs.gov/publications/mobile-bay, but south of the parkway is saltwater: which this sample is from according to gps
[4] same as [3]
[5] coordinates are likely completely off, and also likely seawater
[6] "sendiment" and ocean coordinates. Dont trust it.
[7] supposedly drinking water (from China?) but location is sea outside singapore, and I can't find more info. Disregarding.
[8] same as [0]. Oh, because they're the same ID....
[9] supposed to be river in San Diego, but coordinates are in the pacific ocean right outside it... Might just be approximate. 
Tijuana River kinda close. But same study also has seawater samples. Hard one. Leave as is? Remove? Change lat to 32.55 and long to -117.1 to point at approx. Tijuana river?
[10] dont trust it, completely off. Including misspeling of metagenomics. "metagenimics"
[11] Drinking water from Hong kong, cant find more info. Cooridnates point to ocean in singapore. Removing.
[12] skin? brain? not water at least. Removing.
[13] estuary, can keep for now. Column: taxon.
[14] estuary, can keep
[15] estuary
[16] south china sea, but classified as freshwater. Removing. Supposedly drinking water.
[17] classified as freshwater, but description is duck feces and environment, and location is the south china sea.  Removing.
[18] ballast water. Removing.
[19] says reservoir, but coords is in sea/ocean. It is in right are of world though... Removing. 
[20]
"""
mystery_loc = ["DRR095146", "SRR6370751", "SRR10066355", "SRR3820960", "SRR8003412", "SRR24075716", "SRR26197971", "SRR6986811", "DRR095146", "ERR4702269",
               "SRR24075715", "SRR6797150", "SRR24991626", "SRR20245410", "SRR20245412", "SRR15213102", "SRR6797136", "SRR19631146", "SRR17478312",
              "SRR19503657", "SRR11267126", "SRR6797128", "SRR17405562", "DRR095147", "SRR3534995", "SRR3820958", "SRR8040758", "SRR24075714", "SRR2138582", ]

#how to check example:
merge.loc[mystery_loc[15]]
id_coord.loc[mystery_loc[15]]
merge.loc[mystery_loc[15]]["taxon"]

#How check specific studies example:
#merge[merge.apply(lambda r: r.str.contains('Amazon Continuum Metagenomes').any(), axis=1)] #any column
merge.loc[merge["study"]=="Amazon Continuum Metagenomes"]
id_coord.loc[merge.loc[merge["study"]=="Amazon Continuum Metagenomes"].index]

#Microbial metagenome of suspended particulate matter in the Pearl River Estuary
#"the freshwater extending as far as 55 km offshore" https://www.nature.com/articles/s41467-023-39507-0. Not sure where shore starts, but from
#the inner location Xiahuipai to the bridge it's less than 55, and the coordinates are somewhere in the middle. Keeping these samples.
#Metagenomic sequencing at PRE are also Pearl River estuary. Some are within the 55km and some outside. I can list the ones that are outside to remove.
pearl_river_ovr55 = ["SRR15213101", "SRR15213102", "SRR15213103", "SRR15213104"]

#amazon continuum are seawater samples
amazon_co = id_coord.loc[merge.loc[merge["study"]=="Amazon Continuum Metagenomes"].index].index.to_list()

#Rongjiang River and South China Sea Raw sequence reads
#maybe using the outer bridge as a freshwater border? https://link.springer.com/article/10.1007/s12237-021-00981-8 according to this there's saltwater intrusion happening but oh well, I need to choose something.
rongjiang = id_coord.loc[merge.loc[merge["study"]=="Rongjiang River and South China Sea Raw sequence reads"].index]
#rongjiang[["Longitude", "Latitude"]] = rongjiang[["Longitude", "Latitude"]].apply(pd.to_numeric)
#rongjiang.sort_values(by=["Latitude", "Longitude"])
rongjiang_not_fresh = ["SRR21201134", "SRR21201132", "SRR21201130", "SRR21201110", "SRR21201129", "SRR21201128", "SRR21201119", "SRR21201127",
                      "SRR21201126", "SRR21201117", "SRR21201121", "SRR21201125", "SRR21201118", "SRR21201124", "SRR21201122", "SRR21201123"]

#Metagenomic exploration of antibiotic resistance genes and their hosts in aquaculture waters of Dongshan Bay (China)
#I think it is technically marine, but there is a mix with freshwater and both brackish and freshwater species of phytoplankton has been found
#so maybe other interesting bacteria too? altough all of these samples are kinda close to the outer part of the bay, but close to land. So I'll just remove them.
dongshan = id_coord.loc[merge.loc[merge["study"]=="Metagenomic exploration of antibiotic resistance genes and their hosts in aquaculture waters of Dongshan Bay (China)"].index].index.to_list()

#oil metagenome Raw sequence reads
#seems to be related to some industrial thing
oil = id_coord.loc[merge.loc[merge["study"]=="oil metagenome Raw sequence reads"].index].index.to_list()

#wetland microbial community structure and function diversity
#classified as aquatic, but is wetland. Removing
wmc = id_coord.loc[merge.loc[merge["study"]=="wetland microbial community structure and function diversity"].index].index.to_list()

#study = Mine wastewater metagenomics
#study = 	mining impacted wastewater
#study 	Black Sea metagenomic datasets from 5, 30, 150 and 750 m depths Raw sequence reads
#	BLACK-OMICS: Black Sea metagenomics
#

not_fresh = [mystery_loc[1], mystery_loc[3], mystery_loc[4], mystery_loc[5], mystery_loc[6], mystery_loc[7], mystery_loc[10], mystery_loc[11],
             mystery_loc[12], mystery_loc[16], mystery_loc[17], mystery_loc[18], mystery_loc[19]] + amazon_co + pearl_river_ovr55 + \
            rongjiang_not_fresh + dongshan + oil + wmc
id_coord = id_coord.drop(index=not_fresh)

#fix mystery_loc 0 and jiulong river, the coordinates were put as western longitude when it should be eastern longitude. Convert to positive.
jiulong_id = merge.loc[merge["study"]=="Jiulong River 201209-201306 Metagenome"].index
rev_lon = [mystery_loc[0]] + jiulong_id.to_list()
for lon in rev_lon:
    id_coord.loc[lon]["Longitude"] = id_coord.loc[lon]["Longitude"].split("-")[1]

#save accessions and latitudes and longitudes to csv
id_coord.to_csv('/home/jay/pangenome-pipeline/SRA_scripts/results/accessions_coordinates.csv', encoding='utf-8')


#metagenome_columns = ['LIBRARY_SOURCE','sample-type', 'metagenome-source', 'sample_type', 'Omics']
#metagenome = df_sra_data[metagenome_columns]
#metagenome = metagenome.fillna('')
#metagenome_filtered = metagenome[metagenome['LIBRARY_SOURCE'].str.contains('META')]
