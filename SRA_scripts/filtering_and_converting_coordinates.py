#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 13:58:49 2021

@author: chelsea
"""
import pandas as pd
import numpy as np

#Path to SRA_data and creating dataframe
sra_data = '/home/fer/.local/opt/SRAnsack/workflow/scripts/output/sra_data.csv'
df_sra_data = pd.read_csv(sra_data, sep=',', index_col=1)

#List of things to be removed from sra_data
if True:
    to_filter_out = ['sediment','virus', 'viral', 'virome', 'transcript', 'marine', 'tap water', 'filtered water', 'aquarium', '^genomic', 'seawater', 'river']
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
    freshwater = df_sra_data[freshwater_columns] #create df with only freshwater columns
    freshwater = freshwater.fillna('') #replace NaN with empty strings
    freshwater_filtered = freshwater[freshwater.stack().str.contains('fresh').any(level=0)] #keep freshwater rows
else:
    freshwater_filtered = df_sra_data

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
id_coord = merge[['coordinates']]

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

#save accessions and latitudes and longitudes to csv
id_coord.to_csv('/home/fer/.local/opt/SRAnsack/workflow/scripts/output/accessions_coordinates.csv', encoding='utf-8')


#metagenome_columns = ['LIBRARY_SOURCE','sample-type', 'metagenome-source', 'sample_type', 'Omics']
#metagenome = df_sra_data[metagenome_columns]
#metagenome = metagenome.fillna('')
#metagenome_filtered = metagenome[metagenome['LIBRARY_SOURCE'].str.contains('META')]
