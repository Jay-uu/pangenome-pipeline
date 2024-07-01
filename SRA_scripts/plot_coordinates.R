# RUN
# python3 ~/.local/opt/Population-dynamics-of-freshwater-microorganisms/scripts/Data_collection/filtering_and_converting_coordinates.py
# THIS NEEDS MANUAL EDITING TO SELECT INPUT OUTPUT AND ENABLE/DISABLE FILTERS!!!!
# ALSO IT CURRENTLY CONSIDERS CERTAIN COLUMNS RELATED TO GPS COORDINATES, BUT THERE MIGHT BE MORE. WE HAVE TO CHECK MANUALLY AFTER WE CREATE OUR sra_data.csv

library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(mapview)
library(ggplot2)
ta = read.table('accessions_coordinates.csv', sep=',', header=T, row.names=1)
ta2 = read.table("sra_data.csv", header=T, sep=',', row.names=2, quote='"', comment.char='', as.is=T)
all(rownames(ta) %in% rownames(ta2)) # TRUE

#count how many have geotags
ta2[ta2==''] = NA
gps=colnames(ta2)[grepl('gps|ongitude|atitude|other_gps_coordinates|lat_lon', colnames(ta2))]
ta2 = ta2[,gps]
table(rowSums(!is.na(ta2))>0) # 9281 F 603 T



ta$bases = ta2[rownames(ta),'nb_bases']
ta$contact_name = ta2[rownames(ta),'contact_name']
ta$sample_name = ta2[rownames(ta),'sample_name']
ta$study = ta2[rownames(ta),'study']
ta$title = ta2[rownames(ta),'Title']
ta$collection_date = ta2[rownames(ta),'collection_date']

samples_sf = st_as_sf(ta, coords = c("Longitude", "Latitude"),  crs = 4326)
options(browser = 'firefox')
mapviewOptions(fgb = FALSE)
mapview(samples_sf, zcol='bases')

theme_set(theme_bw())
world = ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world) + geom_sf(fill="floralwhite", alpha=0.5, col="lightgray")+geom_sf(data=samples_sf, shape=18, col='deepskyblue')+coord_sf(xlim=c(-180, 180), ylim=c(-90,90))+scale_size_continuous(range = c(2, 10))
