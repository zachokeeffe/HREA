##########################################################################################
# Title:       download_HREA_COGs.R
# Description: Example script to view STAC information and download HREA data.
# Author:      Zach O'Keeffe <zokeeffe@umich.edu>
# Date:        2024-05-26
##########################################################################################


# List of packages used in this example
package_list <- c('jsonlite', 'curl', 'data.table', 'terra')
# Install packages
install.packages(package_list)
# Load packages
sapply(package_list, require, character.only=TRUE)

# Set printing options to view wide tables
options(width=200)

# Read lookup table
(LUDT <- fread('https://globalnightlight.s3.amazonaws.com/HREA_aux_data/HREA_place_lookup.csv', header=TRUE, sep=','))

# Find available countries
print(unique(LUDT[,list(country_iso3, country_name)]), nrows=nrow(LUDT))

# Select a country of interest (here we'll use Fiji as an example)
(focal_country_iso3 <- unique(LUDT[country_name=='Fiji'][['country_iso3']]))

# Find cog_ids for a particular country
(cog_ids <- LUDT[country_iso3==focal_country_iso3][['cog_id']])


###################
# STAC navigation #
###################

# Read top-level HREA COG catalogue
(HREA_catalog <- fromJSON(txt='https://globalnightlight.s3.amazonaws.com/HREAv1.1_COGs/catalog.json'))

# Read top-level catalogue for country of interest
(country_catalog <- fromJSON(txt=paste0('https://globalnightlight.s3.amazonaws.com/HREAv1.1_COGs/', focal_country_iso3, '/catalog.json')))

# Read catalogues for each cog_id
(cog_id_catalog <- lapply(cog_ids, function(cogid) fromJSON(txt=paste0('https://globalnightlight.s3.amazonaws.com/HREAv1.1_COGs/', focal_country_iso3, '/', cogid, '/catalog.json'))))

# Select a particular cog_id of interest
(focal_cog_id <- cog_ids[1L])

# Read the settlement population json for that cog_id
(set_pop_item <- fromJSON(txt=paste0('https://globalnightlight.s3.amazonaws.com/HREAv1.1_COGs/', focal_country_iso3, '/', focal_cog_id, '/', focal_cog_id, '_set_pop.json')))

# Get the filepath for the settlement layer:
(set_pop_URL <- set_pop_item$assets$population$href)

# Alternatively, one can find the settlement layer by building the URL:
(set_pop_URL_alt <- paste0('https://globalnightlight.s3.amazonaws.com/HREAv1.1_COGs/', focal_country_iso3, '/', focal_cog_id, '/', focal_cog_id, '_set_pop.tif'))

identical(set_pop_URL, set_pop_URL_alt)

# One can also view the collection of COGs of some file_type (rade9lnmu, set_zscore, set_lightscore, set_prplit) for a cog_id:
(file_type <- c('rade9lnmu', 'set_zscore', 'set_lightscore', 'set_prplit')[1L])

(file_type_collection <- fromJSON(txt=paste0('https://globalnightlight.s3.amazonaws.com/HREAv1.1_COGs/', focal_country_iso3, '/', focal_cog_id, '/', file_type, '/collection.json')))


#####################
# Downloading files #
#####################

# Set path for downloading files (the following is an example that follows the directory structure on AWS)
(COG_local_dir <- paste0('~/HREAv1.1_COGs/', focal_country_iso3, '/', focal_cog_id, '/'))
dir.create(COG_local_dir, FALSE, TRUE)

# Download the settlement population layer
download.file(set_pop_URL, destfile=paste0(COG_local_dir, focal_cog_id, '_set_pop.tif'))
# Or do a system call with wget
system(paste('wget -P', COG_local_dir, set_pop_URL))
# Or use curl
curl_download(set_pop_URL, destfile=paste0(COG_local_dir, focal_cog_id, '_set_pop.tif'))

# To download a particular file_type for a given year, select the year (2013-2020):
(year <- (2013:2020)[1L])

# Build the URL:
(COG_URL <- paste0('https://globalnightlight.s3.amazonaws.com/HREAv1.1_COGs/', focal_country_iso3, '/', focal_cog_id, '/', file_type, '/', focal_cog_id, '_', file_type, '_', year, '.tif'))

# One can directly read the COG from the URL in terra:
(HREA_COG <- rast(COG_URL))
