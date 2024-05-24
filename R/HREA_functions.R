##########################################################################################
# Title:       HREA_functions.R
# Description: Core functions to replicate HREA.
# Author:      Zach O'Keeffe <zokeeffe@umich.edu>
# Date:        2024-05-24
##########################################################################################


# ****************************************************************************************
# Name:    make_country_boundary
# Purpose: Creates a simplified country boundary vector. One can supply a boundary, or it
#          will automatically pull data from GADM. If a HRSL file is provided, it will
#          crop the boundary to the layer.
# Inputs:  HRSL_path (character, optional) - path to settlement layer
#          country_boundary_path (character, optional) - path to vector boundary file
#          focal_country_iso (character, required) - 3 letter ISO code for country
#          cog_id (character, required) - COG ID of unit
# Output:  NULL on success
# Notes:   Run at the cog_id level.
# ****************************************************************************************

make_country_boundary <- function(HRSL_path=NA_character_, 
                         country_boundary_path=NA_character_,
                         focal_country_iso=NA_character_, cog_id=NA_character_) {
  # Load necessary packages
  require(geodata)
  
  # Stop function if missing necessary inputs
  if(is.na(focal_country_iso3)) stop('Requires the country ISO3.')
  if(is.na(cog_id)) stop('Requires the country COG ID.')
  # Make directories for files
  out_dir_dat <- paste0('data', focal_country_iso3, sep = '/')
  dir.create(out_dir_dat, FALSE, TRUE)
  
  # Read in country boundary
  if(file.exists(country_boundary_path)){
    country_boundary <- vect(country_boundary_path)
  } else {
    country_boundary <- gadm(focal_country_iso, level=0, path=tempdir(), version="3.6")
  }
  
  # If HRSL is present, crop
  if(file.exists(HRSL_path)){
    # Read in HRSL settlement raster
    set_rast <- rast(HRSL_path)
    
    # Crop to settlement layer
    country_boundary <- crop(country_boundary, set_rast)
  }
  
  # Simplify
  country_boundary <- disagg(country_boundary)
  country_boundary <- simplifyGeom(country_boundary, tolerance=0.1, preserveTopology=TRUE,
                      makeValid=TRUE)

  # Save vector to disk
  saveRDS(country_boundary, paste0(out_dir_dat, '/', cog_id, '_boundary.rds'))
  
  # Return NULL
  return(NULL)
}


# ****************************************************************************************
# Name:    make_grid_pop_rast
# Purpose: Creates constant geospatial grids for a location (usually a country) and 
#          assigns them population data from the settlement layer, as well as land cover
#          values from MODIS. One must specify either HRSL_path or GHSL_path to point to
#          the settlement layer to use.
# Inputs:  HRSL_path (character, optional) - path to HRSL file
#          GHSL_path (character, optional) - path to GHSL file
#          focal_country_iso (character, required) - 3 letter ISO code for country
#          cog_id (character, required) - COG ID of unit
# Output:  NULL on success
# Notes:   Run at the cog_id level. Requires that the cog_id boundary has been generated.
# ****************************************************************************************

make_grid_pop_rast <- function(HRSL_path=NA_character_, GHSL_path=NA_character_,
                      focal_country_iso=NA_character_, cog_id=NA_character_) {
  # Load necessary packages
  require(MODISTools)
  require(data.table)
  
  # Stop function if missing necessary inputs
  if(is.na(focal_country_iso3)) stop('Requires the country ISO3.')
  if(is.na(cog_id)) stop('Requires the country COG ID.')
  # Make directories for files
  out_dir_dat <- paste0('data', focal_country_iso3, sep = '/')
  dir.create(out_dir_dat, FALSE, TRUE)
  out_dir_cog <- paste0('data', focal_country_iso3, cog_id, sep = '/')
  dir.create(out_dir_cog, FALSE, TRUE)
  
  # Read in country boundary (fails if doesn't exist)
  country_boundary <- readRDS(paste0(out_dir_dat, '/', cog_id, '_boundary.rds')))
  
  if(file.exists(HRSL_path)){
    # Read in HRSL settlement raster
    set_rast <- rast(HRSL_path)
  
    # Make 1as raster with same extent as settlement layer
    rast_hr <- rast(resolution=c(0.0002777778, 0.0002777778), extent=ext(set_rast)
               crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
    # Save empty raster to populate later
    writeRaster(rast_hr, paste0(out_dir_dat, '/', cog_id, '_empty_grid.tif'), 
      format="GTiff", overwrite=TRUE, dataType='FLT4S', setStatistics=FALSE, 
      NAflag=-999999, options=c("COMPRESS=LZW"))
    
    # Fill with 0s
    rast_hr[] <- 0
    # Mask with shape
    rast_hr <- mask(rast_hr, country_boundary)
    
    # Get non-missing cell numbers
    cells_hr <- cells(rast_hr, na.rm=TRUE)
    # Get raster x-y values
    xys_hr <- crds(rast_hr, na.rm=TRUE)
    # Get HRSL cell numbers
    set_cells <- cellFromXY(set_rast, xys_hr)
    # Get population values
    pop_hr <- set_rast[set_cells]
    # Make data.table of values
    cellpop_hr <- data.table(cellnum_hr = cells_hr, pop = pop_hr)
    # Set missing population to 0
    cellpop_hr[is.na(pop), pop:=0]
  } else if(file.exists(GHSL_path)){
    # Read in GHSL settlement raster
    set_rast <- rast(GHSL_path)
    
    # Reproject country boundary
    country_boundary_proj <- project(shp, set_rast)
    
    # Crop and mask
    rast_hr <- crop(set_rast, country_boundary_proj)
    rast_hr <- mask(rast_hr, country_boundary_proj)
    
    # Reproject raster
    rast_hr <- project(rast_hr, 
               '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0',
               res=c(0.0020833333333333333,0.0020833333333333333), method='bilinear')
    
    # Get non-missing cell numbers
    cells_hr <- cells(rast_hr, na.rm=TRUE)
    # Get raster x-y values
    xys_hr <- crds(rast_hr, na.rm=TRUE)
    # Get HRSL cell numbers
    set_cells <- cellFromXY(set_rast, xys_hr)
    # Get population values
    pop_hr <- set_rast[set_cells]
    # Make data.table of values
    cellpop_hr <- data.table(cellnum_hr = cells_hr, pop = pop_hr)
    # Set missing population to 0
    cellpop_hr[is.na(pop)|pop<2, pop:=0]
    # Set everything to missing
    rast_hr[] <- NA
    # Save empty raster to populate later
    writeRaster(rast_hr, paste0(out_dir_dat, '/', cog_id, '_empty_grid.tif'), 
      format="GTiff", overwrite=TRUE, dataType='FLT4S', setStatistics=FALSE, 
      NAflag=-999999, options=c("COMPRESS=LZW"))
  } else {
    stop('Requires either HRSL or GHSL layer.')
  }
  
  # Make 15as grid
  rast_15as <- rast(resolution=c(0.0041666667, 0.0041666667), extent=ext(rast_hr)
               crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
  
  # Save empty raster to populate later
  writeRaster(rast_15as, paste0(out_dir_dat, '/', cog_id, '_empty_grid_15as.tif'), 
    format="GTiff", overwrite=TRUE, dataType='FLT4S', setStatistics=FALSE, 
    NAflag=-999999, options=c("COMPRESS=LZW"))
  
  # Get 15as cell values from 1as coordinates
  set(cellpop_hr, NULL, 'cellnum_15as', cellFromXY(rast_15as, xys_hr))
  
  # Set data.table keys
  setkey(cellpop_hr, cellnum_15as, cellnum_hr)
  
  # Make 15as data.table
  cellpop_15as <- cellpop_hr[, list(pop=sum(pop)), by='cellnum_15as']
  
  # Get 15as cell centroid coordinates and add to data.table
  cell_15as_xy <- xyFromCell(rast_15as, cellpop_15as[['cellnum_15as']])
  set(cellpop_15as, NULL, 'x', cell_15as_xy[,1L])
  set(cellpop_15as, NULL, 'y', cell_15as_xy[,2L])
  
  # Make bounding box polygon of grid
  rast_bbox <- vect(ext(rast_15as))
  # Get centroid of country
  rast_cent <- centroids(rast_bbox), inside=FALSE)
  # Get distance from points
  country_km_lr <- ceiling(distance(x=crds(rast_bbox)[1,,drop=F], 
                    y=crds(rast_bbox)[4,,drop=F], lonlat=TRUE, pairwise=TRUE)/1000)
  country_km_ab <- ceiling(distance(x=crds(rast_bbox)[1,,drop=F], 
                    y=crds(rast_bbox)[2,,drop=F], lonlat=TRUE, pairwise=TRUE)/1000)
  
  # Get MODIS land cover data
  LC <- mt_subset(product = "MCD12Q1", lat = crds(rast_cent)[2], lon = crds(rast_cent)[1], 
        band = "LC_Type1", start = "2012-01-01", end = "2012-12-31",
        km_lr = country_km_lr, km_ab = country_km_ab, site_name = "tempsite", 
        internal = TRUE, progress = FALSE)
  # Convert to raster
  LC_r <- mt_to_terra(df = LC, reproject = TRUE, method = "near")
  # Get land cover cells
  cells_LC <- cellFromXY(LC_r, cell_15as_xy)
  # Get land cover values
  LC_vals <- LC_r[cells_LC]
  set(cellpop_15as, NULL, 'LC',LC_vals)
  set(cellpop_15as, NULL, 'LC_fac', factor(cellpop_15as[['LC']], levels=1:17, 
    labels=c('Evergreen needleleaf forest', 'Evergreen broadleaf forest',
    'Deciduous needleleaf forest', 'Deciduous broadleaf forest', 'Mixed forest',
    'Closed shrubland', 'Open shrubland', 'Woody savanna', 'Savanna', 'Grassland',
    'Permanent wetland', 'Cropland', 'Urban and built-up',
    'Cropland/natural vegetation mosaic', 'Snow and ice', 'Barren', 'Water')))
  
  # Save data.tables
  saveRDS(cellpop_hr, paste0(out_dir_dat, '/', cog_id, '_cell_pop.rds'))
  saveRDS(cellpop_15as, paste0(out_dir_dat, '/', cog_id, '_cell_pop_xy_15as.rds'))
  
  # Discard cells with zero population
  cellpop_hr <- cellpop_hr[pop!=0]
  
  # Assign values to raster
  rast_hr[cellpop_hr[['cellnum_hr']]] <- cellpop_hr[['pop']]
  
  # Save raster  
  writeRaster(rast_hr, paste0(out_dir_cog, '/', cog_id, '_set_pop.tif'), format="GTiff", 
    overwrite=TRUE, dataType='FLT4S', setStatistics=TRUE,  NAflag=-999999, 
    options=c("COMPRESS=LZW"))
  
  # Return NULL
  return(NULL)
}


# ****************************************************************************************
# Name:    find_iso_nsets
# Purpose: Identified non-settlements that are at least one 15as pixel away from 
#          settlement cells.
# Inputs:  focal_country_iso (character, required) - 3 letter ISO code for country
#          cog_id (character, required) - COG ID of unit
# Output:  NULL on success
# Notes:   Run at the cog_id level. Requires that the grid and settlement data have been
#          generated.
# ****************************************************************************************

find_iso_nsets <- function(focal_country_iso=NA_character_, cog_id=NA_character_) {
  # Load necessary packages
  require(terra)
  require(data.table)
  
  # Stop function if missing necessary inputs
  if(is.na(focal_country_iso3)) stop('Requires the country ISO3.')
  if(is.na(cog_id)) stop('Requires the country COG ID.')
  # Make directories for files
  out_dir_dat <- paste0('data', focal_country_iso3, sep = '/')
  dir.create(out_dir_dat, FALSE, TRUE)
  out_dir_cog <- paste0('data', focal_country_iso3, cog_id, sep = '/')
  dir.create(out_dir_cog, FALSE, TRUE)
  
  # Read in country cell info (fails if doesn't exist)
  cellpop_15as <- readRDS(paste0(out_dir_dat, '/', cog_id, '_cell_pop_xy_15as.rds'))
  # Specify logical whether cell is a settlement
  cellpop_15as[,set:=pop>0]
  # Read in 15as raster (fails if doesn't exist)
  rast_15as <- rast(paste0(out_dir_dat, '/', cog_id, '_empty_grid_15as.tif'))
  
  # Add column and row numbers to data.table
  set(cellpop_15as, NULL, 'col', colFromCell(rast_15as, cellpop_15as[['cellnum_15as']]))
  set(cellpop_15as, NULL, 'row', rowFromCell(rast_15as, cellpop_15as[['cellnum_15as']]))
  # Get total number of rows and columns
  rastNrow <- nrow(rast_15as)
  rastNcol <- ncol(rast_15as)
  
  # Specify edge rows and columns to ignore
  edgerows <- c(1, rastNrow)
  edgecols <- c(1, rastNcol)
  
  # Make settlement and initial isolated non-settlement data.tables
  setDT <- cellDTf[set==TRUE, list(set,col,row)])
  isoDT <- cellDTf[set==FALSE&!(row%in%edgerows|col%in%edgecols), list(cellnum_15as,
           col, row)]
  
  # Specify columns to walk through
  cs2do <- -1:1
  # Walk through columns
  for(tc in cs2do){
    # Specify rows to walk through (ignore origin 0,0)
    if(tc==0L){
      (rs2do <- setdiff(cs2do,0))
    } else {
      (rs2do <- cs2do)
    }
    # Walk through rows
    for(r in rs2do){
      # Add row and column to settlement table
      setcellDT <- setDT[,list(set,col=col+tc,row=row+r)]
      # Merge with isolated non-settlement table
      isoDT <- merge(isoDT, setcellDT, by=c('col','row'), all.x=T)
      # Remove non-settlement cells next to settlements
      isoDT <- isoDT[is.na(set), !'set', with=F]
      # Stop if nothing left
      if(nrow(isoDT)==0) stop('No isolated non-settlements.')
    }
  }
  
  # Extract isolated non-settlement cells
  iso_cells <- isoDT[['cellnum_15as']]
  
  # Save isolated non-settlement cells
  saveRDS(iso_cells, paste0(out_dir_dat, '/', cog_id, '_iso_cells.rds'))
  
  # Return NULL
  return(NULL)
}


# ****************************************************************************************
# Name:    get_VIIRS_dat
# Purpose: Downloads and processes VIIRS-DNB data from the LEN repository for a given
#          cog_id and month.
# Inputs:  focal_country_iso (character, required) - 3 letter ISO code for country
#          cog_id (character, required) - COG ID of unit
#          year (integer, required) - year in YYYY format
#          month (integer, required) - month in numeric MM format
#          good_vflag_path (character, required) - path to good vflag values
#          viirs_dat_dir (character, required) - path to downloaded VIIRS data
#          tiff_list_path (character, required) - path to overlapping VIIRS TIFFs
# Output:  NULL on success
# Notes:   Run at the cog_id-month level. Requires that grid has been generated,
#          overlapping LEN TIFFs have been identified, and the orbital strip data has
#          been downloaded.
# ****************************************************************************************

get_VIIRS_dat <- function(focal_country_iso=NA_character_, cog_id=NA_character_,
                 year=NA_integer_, month=NA_integer_, good_vflag_path=NA_character_,
                 viirs_dat_dir=NA_character_, tiff_list_path=NA_character_) {
  # Load necessary packages
  require(data.table)
  require(terra)
  
  # Stop function if missing necessary inputs
  if(is.na(focal_country_iso3)) stop('Requires the country ISO3.')
  if(is.na(cog_id)) stop('Requires the country COG ID.')
  if(is.na(yearmo)) stop('Requires the month of observation.')
  
  # Make directories for files
  out_dir_dat <- paste0('data', focal_country_iso3, sep = '/')
  dir.create(out_dir_dat, FALSE, TRUE)
  
  # Make satellite year-month string
  if (year>=2013L & year<=2017L) {
    YYYYMM <- paste0(year, formatC(month, width = 2, format = "d", flag = "0"))
  } else {
    YYYYMM <- paste0('npp_', year, formatC(month, width = 2, format = "d", flag = "0"))
  }
  # Read in cell info (fails if doesn't exist)
  cellpop_15as <- readRDS(paste0(out_dir_dat, '/', cog_id, '_cell_pop_xy_15as.rds'))
  # Make second offset
  cellpop_15as[,secOffset:=x*240]
  # Make x-y point matrices with an offset for longitude values above 180 and extents
  xyMat <- as.matrix(cellpop_15as[,list(x, y)])
  xyMat2 <- as.matrix(cellpop_15as[,list(x=x%%360, y)])
  xyExt <- ext(min(xyDT[['x']])-.1, max(xyDT[['x']])+.1, min(xyDT[['y']])-.1,
           max(xyDT[['y']])+.1)
  xyExt2 <- ext(sort(xyExt[1:2]%%360),xyExt[3],xyExt[4])
  cellpop_15as <- cellpop_15as[, list(cellnum_15as, secOffset)]
  
  # Read in good-quality vflags (fails if doesn't exist)
  good_vflags <- readRDS(good_vflag_path)
  
  # Read in table of overlapping TIFFs
  tiff_list <- fread(tiff_list_path, header=TRUE, sep=',')
  nstrips <- nrow(tiff_list)
  
  # Make empty list to populate
  finDat <- vector('list', nstrips)
  # For each row in the table of TIFFs:
  for(r in 1:nstrips){
    # Make empty list to populate
    DT <- lapply(c("li","vflag","rade9","samples"), function(ft) {
      trast <- rast(paste(tiff_list_path, tiff_list[r][[ft]], sep='/'))
      rast_ext <- ext(trast)
      if(rast_ext[2]>180&xyExt[2]<0){
        trast <- tryCatch(crop(trast,xyExt2),error=function(e) return(NULL))
        if(is.null(trast)) return(NULL)
        return(extract(trast, xyMat2, cells=FALSE, method="simple"))
      } else {
        trast <- tryCatch(crop(trast,xyExt),error=function(e) return(NULL))
        if(is.null(trast)) return(NULL)
        return(extract(trast, xyMat, cells=FALSE, method="simple"))
      }
    })
    if(any(sapply(DT,is.null))) next
    
    # Combine and subset
    DT <- do.call(cbind, DT)
    DT <- as.data.table(DT)
    setnames(DT, filetypes)
    DT <- cbind(DT, cellpop_15as)
    DT <- DT[vflag%in%good_vflags]
    DT[,vflag:=NULL]
    DT <- DT[!is.na(rade9)]
    DT <- DT[!is.na(li)]
    DT <- DT[rade9>(-1.5)]
    DT <- DT[li<=.005]
    if(nrow(DT)==0) next
    
    # Get time
    year <- as.integer(substr(vfn,6L,9L)))
    month <- as.integer(substr(vfn,10L,11L)))
    day <- as.integer(substr(vfn,12L,13L)))
    stime <- substr(vfn,16L,22L))
    shr <- as.integer(substr(vfn,16L,17L)))
    smin <- as.integer(substr(vfn,18L,19L)))
    ssec <- as.integer(substr(vfn,20L,21L)))
    ehr <- as.integer(substr(vfn,25L,26L)))
    emin <- as.integer(substr(vfn,27L,28L)))
    esec <- as.integer(substr(vfn,29L,30L)))
    sdatetime <- as.POSIXct(paste(paste(year,month,day,sep='-'), paste(shr, smin, ssec, 
                 sep=':')), format='%Y-%m-%d %H:%M:%S', tz="UTC"))
    edatetime <- as.POSIXct(paste(paste(year,month,day,sep='-'), paste(ehr, emin, esec,
                 sep=':')), format='%Y-%m-%d %H:%M:%S', tz="UTC"))
    if(edatetime<sdatetime) edatetime <- edatetime+86400]
    setimediff <- edatetime-sdatetime
    meantime <- sdatetime + setimediff/2
    set(DT, NULL, 'mtimeloc', meantime+DT[['secOffset']])
    DT[,secOffset:=NULL]
    
    # Return strip data
    finDat[[r]] <- DT
    rm(DT);gc()
  }
  # Combine data
  finDat <- rbindlist(finDat, fill=TRUE)
  if(nrow(finDat)==0) stop('No data extracted.')
  # Save each column as a vector
  for(tcol in c('cellnum_15as','li','rade9','samples','mtimeloc')){
    saveRDS(finDat[[tcol]], paste0(out_dir_dat, cog_id, '_', YYYYMM, '_', tcol, '.rds'))
  }
  
  # Return NULL
  return(NULL)
}


# ****************************************************************************************
# Name:    make_HREA_mod
# Purpose: Identifies isolated non-settlements, removes outliers, and randomly samples
#          them for use in the regression. Runs a linear mixed model predicting radiance
#          for these points using the HREA model specification. Saves the model to disk.
# Inputs:  focal_country_iso (character, required) - 3 letter ISO code for country
#          year (integer, required) - year for HREA regression
#          cog_lookup_path (character, required) - location of country lookup table
# Output:  NULL on success
# Notes:   Run at the country_iso3-year level. Requires that all VIIRS data have been
#          downloaded for the year, and that the isolated non-settlement points have been
#          identified.
# ****************************************************************************************

make_HREA_mod <- function(focal_country_iso=NA_character_, year=NA_integer_,
                 cog_lookup_path=NA_character_) {
  # Load necessary packages
  require(data.table)
  require(stackoverflow)
  require(lme4)
  RNGkind("L'Ecuyer-CMRG")
  
  # Stop function if missing necessary inputs
  if(is.na(focal_country_iso3)) stop('Requires the country ISO3.')
  if(is.na(cog_id)) stop('Requires the country COG ID.')
  if(is.na(year)) stop('Requires the year of observation.')
  
  # Make directories for files
  out_dir_dat <- paste0('data', focal_country_iso3, sep = '/')
  dir.create(out_dir_dat, FALSE, TRUE)
  
  # Specify months
  if (year>=2013L & year<=2017L) {
    YYYYMMs <- paste0(year, formatC(1:12, width = 2, format = "d", flag = "0"))
  } else {
    YYYYMMs <- paste0('npp_', year ,formatC(1:12, width = 2, format = "d", flag = "0"))
  }
  nYYYYMMs <- length(YYYYMMs)
  
  # Check if country has multiple parts
  cog_lookup <- fread(cog_lookup_path, header=TRUE, sep=',')
  cog_ids <- cog_lookup[country_iso3==focal_country_iso][['cog_id']]
  ncog_ids <- length(cog_ids)
  
  if(ncog_ids==1){
    cog_id <- cog_ids
    # Read in cell info (fails if doesn't exist)
    cellpop_15as <- readRDS(paste0(out_dir_dat, '/', cog_id, '_cell_pop_xy_15as.rds'))
    # Read isolated non-settlement cells
    iso_cells <- readRDS(paste0(out_dir_dat, '/', cog_id, '_iso_cells.rds'))
    # Keep only isolated non-settlements
    cellpop_15as <- cellpop_15as[cellnum_15as%in%iso_cells]
    # Get unique land cover values
    lcs2keep <- unique(cellpop_15as[['LC']])    
    # Make empty vector for data
    nsmsDT <- vector('list', nYYYYMMs)
    # Read in and manipulate data
    for(i in 1:nYYYYMMs){
      m <- YYYYMMs[i]
      if(!file.exists(paste0(out_dir_dat, cog_id, '_', m, "_cellnum_15as.rds"))){
        nsmsDT[[i]] <- data.table(NULL)
        next
      }
      tdt <- data.table(cellnum_15as=readRDS(paste0(out_dir_dat, cog_id, '_', m,
             "_cellnum_15as.rds")))
      set(tdt, NULL, 'r9', readRDS(paste0(out_dir_dat, cog_id, '_', m, "_rade9.rds")))
      set(tdt, NULL, 'li', readRDS(paste0(out_dir_dat, cog_id, '_', m, "_li.rds")))
      set(tdt, NULL, 'mtimeloc', readRDS(paste0(out_dir_dat, cog_id, '_', m, 
        "_mtimeloc.rds")))
      tdt <- merge(tdt, cellpop_15as, by='cellnum_15as')
      if(nrow(tdt)==0) next
      set(tdt, NULL, 'lirescale', sqrt(tdt[['li']]/.005))
      set(tdt, NULL, 'li',NULL)
      set(tdt, NULL, 'locdatechar', format.POSIXct(tdt[['mtimeloc']], format='%Y-%m-%d',
        tz='UTC'))
      set(tdt, NULL, 'timehour', hour(tdt[['mtimeloc']]) + minute(tdt[['mtimeloc']])/60 +
        second(tdt[['mtimeloc']])/3600)
      set(tdt, NULL, 'timefix', ((tdt[['timehour']]-23)%%24)/5)
      tdt <- tdt[timefix<=1]
      set(tdt, NULL, 'mtimeloc', NULL)
      set(tdt, NULL, 'timehour', NULL)
      nsmsDT[[i]] <- tdt
      rm(tdt);gc()
    }
    matchvars <- 'cellnum_15as'
  } else {
    lcs2keep <- cellpop_15as <- vector('list', ncog_ids)
    for(j in 1:ncog_ids){
      cog_id <- cog_ids[j]
      tdt <- readRDS(paste0(out_dir_dat, '/', cog_id, '_cell_pop_xy_15as.rds'))
      tcells <- readRDS(paste0(out_dir_dat, '/', cog_id, '_iso_cells.rds'))
      tdt <- tcells[cellnum_15as%in%tcells]
      set(tdt, NULL, 'COGID', cog_id)
      setkey(tdt,NULL,COGID,cellnum_15as)
      cellpop_15as[[j]] <- tdt
      lcs2keep[[j]] <- unique(tdt[['LC']])
    }
    lcs2keep <- unique(unlist(lcs2keep))
    # Make empty vector for data
    nsmsDT <- vector('list', nYYYYMMs)
    # Read in and manipulate data
    for(i in 1:nYYYYMMs){
      nsmsDT[[i]] <- vector('list', ncog_ids)
      for(j in 1:ncog_ids){
        m <- YYYYMMs[i]
        cog_id <- cog_ids[j]
        if(!file.exists(paste0(out_dir_dat, cog_id, '_', m, "_cellnum_15as.rds"))){
          nsmsDT[[i]][[j]] <- data.table(NULL)
          next
        }
        tdt <- data.table(cellnum_15as=readRDS(paste0(out_dir_dat, cog_id, '_', m,
               "_cellnum_15as.rds")))
        set(tdt, NULL, 'r9', readRDS(paste0(out_dir_dat, cog_id, '_', m, "_rade9.rds")))
        set(tdt, NULL, 'li', readRDS(paste0(out_dir_dat, cog_id, '_', m, "_li.rds")))
        set(tdt, NULL, 'mtimeloc', readRDS(paste0(out_dir_dat, cog_id, '_', m,
          "_mtimeloc.rds")))
        tdt <- merge(tdt, cellpop_15as[[j]], by='cellnum_15as')
        tdt <- tdt[li<=.005]
        if(nrow(tdt)==0) next
        set(tdt, NULL, 'lirescale', sqrt(tdt[['li']]/.005))
        set(tdt, NULL, 'li',NULL)
        set(tdt, NULL, 'locdatechar', format.POSIXct(tdt[['mtimeloc']], format='%Y-%m-%d',
          tz='UTC'))
        set(tdt, NULL, 'timehour', hour(tdt[['mtimeloc']]) + minute(tdt[['mtimeloc']])/60+
          second(tdt[['mtimeloc']])/3600)
        set(tdt, NULL, 'timefix', ((tdt[['timehour']]-23)%%24)/5)
        tdt <- tdt[timefix<=1]
        set(tdt, NULL, 'mtimeloc', NULL)
        set(tdt, NULL, 'timehour', NULL)
        nsmsDT[[i]][[j]] <- tdt
        rm(tdt);gc()
      }
    }
    matchvars <- c('COGID', 'cellnum_15as')
  }
  
  # Get means and standard deviations of radiance by cell
  if(ncog_ids==1){
    Nrows <- sum(as.numeric(unlist(lapply(nsmsDT,nrow))))
    if(Nrows>.Machine$integer.max){
      Nsplits <- ceiling(Nrows/.Machine$integer.max)
      idchunks <- chunk2(iso_cells, Nsplits)
      nsmsDT2 <- vector('list', Nsplits)
      for(chun in 1:Nsplits){
        tmpDT <- rbindlist(lapply(nsmsDT, function(xx) 
                 tryCatch(xx[cellnum_15as%in%idchunks[[chun]], list(cellnum_15as, r9)],
                 error=function(e) NULL)), fill=T)
        tmpDT <- tmpDT[,list(r9m=mean(r9), r9s=sd(r9), N=.N), by='cellnum_15as']
        nsmsDT2[[chun]] <- tmpDT
        rm(tmpDT);gc()
      }
      nsmsDT2 <- rbindlist(nsmsDT2, fill=TRUE)
    } else {
      nsmsDT2 <- rbindlist(nsmsDT, fill=TRUE)
      nsmsDT2 <- nsmsDT2[,list(r9m=mean(r9), r9s=sd(r9), N=.N), by='cellnum_15as']
    }
  } else {
    nsmsDT2 <- vector('list', ncog_ids)
    for(j in 1:ncog_ids){
      tdt <- lapply(nsmsDT, function(xx) tryCatch(xx[[j]][,c(matchvars,'r9'), with=FALSE],
             error=function(e) NULL))
      tdt <- rbindlist(tdt, fill=TRUE)
      tdt <- tdt[,list(r9m=mean(r9), r9s=sd(r9), N=.N), by=matchvars]
      nsmsDT2[[j]] <- tdt
    }
    nsmsDT2 <- rbindlist(nsmsDT2, fill=TRUE)
    cellpop_15as <- rbindlist(cellpop_15as, fill=TRUE)
  }
  nsmsDT2 <- na.omit(nsmsDT2[N>=(nYYYYMMs*3)])
  stopifnot(nrow(nsmsDT2)>0)
  
  # Randomly sample cells to keep after outlier removal
  set.seed(48109)
  ids2keep <- lapply(1:length(lcs2keep), function(i) {
    lctt <- lcs2keep[i]
    tnsmsdt <- merge(nsmsDT2, cellpop_15as[LC==lctt,c(matchvars),with=F], by=matchvars)
    tmquants <- quantile(tnsmsdt[['r9m']],c(.01,.5),na.rm=T)
    tmquantdif <- tmquants[[2]]-tmquants[[1]]
    tmthreshes <- c(tmquants[1L],tmquants[2L]+tmquantdif)
    tsquants <- quantile(tnsmsdt[['r9s']],c(.01,.5),na.rm=T)
    tsquantdif <- tsquants[[2]]-tsquants[[1]]
    tsthreshes <- c(tsquants[1L],tsquants[2L]+tsquantdif)
    tids2keep <- tnsmsdt[r9m>tmthreshes[1L]&r9m<tmthreshes[2L]&r9s>tsthreshes[1]&
                 r9s<tsthreshes[2L]][,c(matchvars),with=F]
    tids2keep <- tids2keep[sample(.N, min(.N,1000), replace=F)]
    return(tids2keep)
  })
  ids2keep <- rbindlist(ids2keep, fill=TRUE)
  setkeyv(ids2keep, matchvars)
  
  # Subset data
  for(i in 1:nYYYYMMs){
    if(ncog_ids==1){
      if(nrow(nsmsDT)==0) next
      nsmsDT[[i]] <- merge(nsmsDT[[i]], ids2keep, by=matchvars)
    } else {
      for(j in 1:ncog_ids){
        if(nrow(nsmsDT[[i]][[j]])==0) next
        nsmsDT[[i]][[j]] <- merge(nsmsDT[[i]][[j]], ids2keep, by=matchvars)
      }
      nsmsDT[[i]] <- rbindlist(nsmsDT[[i]], fill=TRUE)
    }
  }
  nsmsDT <- rbindlist(nsmsDT, fill=TRUE)
  
  # Additional outlier removal
  set(nsmsDT,NULL,'r9l',log(nsmsDT[['r9']]+2.5))
  r9med <- median(nsmsDT[['r9l']])
  r9sd <- sd(nsmsDT[['r9l']])
  nsmsDT <- nsmsDT[r9l<=(r9med+4*r9sd),!'r9l',with=F]
  lcdates<-unique(nsmsDT[,list(LC,locdatechar)])
  setorder(lcdates,LC,locdatechar)
  regDT <- lapply(1:nrow(lcdates), function(r) {
    tdt <- nsmsDT[LC==lcdates[r][['LC']]&locdatechar==lcdates[r][['locdatechar']]]
    if(nrow(tdt)>=5L){
      r9m<-mean(tdt[['r9']])
      r9s<-sd(tdt[['r9']])
      tdt<-tdt[r9<=(r9m+4*r9s)]
    }
    return(tdt)
  })
  regDT <- rbindlist(regDT, fill=TRUE)
  
  # Create factor variables
  set(regDT, NULL, 'monthfac', factor(as.integer(substr(regDT[['locdatechar']],6L,7L)),
    levels=1:12, labels=month.abb))
  
  # Run model
  mod <- lmer(r9 ~ lirescale + timefix + monthfac + LC_fac + LC_fac:lirescale +
         (1|locdatechar), data=regDT, REML=TRUE))
  # Save model
  saveRDS(mod, paste0(out_dir_dat, '/', focal_country_iso3, '_HREA_LMM_', year, '.rds'))
  
  return(NULL)
}


# ****************************************************************************************
# Name:    make_HREA_data
# Purpose: Creates three HREA metrics (mean z-score, lightscore, and proportion of nights
#          lit) by reading in settlement data, generating predicted values from the HREA
#          model, removing extreme outliers, and generating statistics. Saves the HREA
#          metrics to disk.
# Inputs:  focal_country_iso (character, required) - 3 letter ISO code for country
#          cog_id (character, required) - COG ID of unit
#          year (integer, required) - year of HREA regression
#          cog_lookup_path (character, required) - location of country lookup table
#          LC_sub_path (character, required) - location of land cover substitution table
# Output:  NULL on success
# Notes:   Run at the cog_id-year level. Requires that the appropriate country model has
#          been generated, and all VIIRS data have been downloaded. Land cover
#          substitution file must also be available.
# ****************************************************************************************

make_HREA_data <- function(focal_country_iso=NA_character_, cog_id=NA_character_,
                  year=NA_integer_, cog_lookup_path=NA_character_,
                  LC_sub_path=NA_character_) {
  # Load necessary packages
  require(data.table)
  require(lme4)
  
  # Make directories for files
  out_dir_dat <- paste0('data', focal_country_iso3, sep = '/')
  dir.create(out_dir_dat, FALSE, TRUE)
  
  # Specify months
  if (year>=2013L & year<=2017L) {
    YYYYMMs <- paste0(year, formatC(1:12, width = 2, format = "d", flag = "0"))
  } else {
    YYYYMMs <- paste0('npp_', year ,formatC(1:12, width = 2, format = "d", flag = "0"))
  }
  nYYYYMMs <- length(YYYYMMs)
  
  # Specify country model (small island countries require other models)
  if(focal_country_iso %in% c('NRU',"PLW",'FSM')) {
    mod_iso <- 'SLB'
  } else {
    mod_iso <- focal_country_iso
  }
  
  # Read in model
  mod <- readRDS(paste0('data/', mod_iso, '/', mod_iso, '_HREA_LMM_', year, '.rds'))
  mod_sig <- sigma(mod)
  
  # Get months of data
  monthLabels <- levels(mod@frame$monthfac)
  nmonthLabels <- length(monthLabels)
  
  # Get land cover types
  lcLabels <- levels(mod@frame$LC_fac)
  nlcLabels <- length(lcLabels)
    
  # Read in country cell info (fails if doesn't exist)
  cellpop_15as <- readRDS(paste0(out_dir_dat, '/', cog_id, '_cell_pop_xy_15as.rds'))
  # Keep only settlements
  cellpop_15as <- cellpop_15as[pop>0]
  # Get unique land cover types
  set_LCs <- unique(cellpop_15as[['LC_fac']])
  # Use substitution if necessary
  if(!all(set_LCs%in%lcLabels)){
    # Read in land cover substitution table
    LClookup <- readRDS(LC_sub_path)
    setnames(cellpop_15as,'LC_fac','LC_fac_orig')
    LCmatch1 <- LClookup[LC_fac_orig%in%intersect(set_LCs,lcLabels), 
                list(LC=origLC, LC_fac=LC_fac_orig)]
    LCmatch2 <- LClookup[LC_fac_orig%in%setdiff(set_LCs,lcLabels) & 
                LC_fac_alt%in%lcLabels, list(LC=origLC, LC_fac=LC_fac_alt)]
    setorder(LCmatch2, LC, ind)
    LCmatch2[,indfin:=1:.N,by='LC']
    LCmatch2 <- LCmatch2[indfin==1, list(LC, LC_fac)]
    LCmatch <- rbindlist(list(LCmatch1, LCmatch2), fill=TRUE)
    cellpop_15as <- merge(cellpop_15as, LCmatch, by='LC')
  }
 
  # Specify months to read
  YYYYMMs <- as.character(match(monthLabels,month.abb)+year*100)
  if(YYYY>=2018L) YYYYMMs <- paste0('npp_',YYYYMMs)
  nYYYYMMs <- length(YYYYMMs)
  
  # Make empty vector for data
  predDT <- vector('list', nYYYYMMs)
  # Read in and manipulate data
  for(i in 1:nYYYYMMs){
    m <- YYYYMMs[i]
    if(!file.exists(paste0(out_dir_dat, cog_id, '_', m, "_cellnum_15as.rds"))){
      predDT[[i]] <- data.table(NULL)
      next
    }
    tdt <- data.table(cellnum_15as=readRDS(paste0(out_dir_dat, cog_id, '_', m,
           "_cellnum_15as.rds")))
    set(tdt, NULL, 'r9', readRDS(paste0(out_dir_dat, cog_id, '_', m, "_rade9.rds")))
    set(tdt, NULL, 'li', readRDS(paste0(out_dir_dat, cog_id, '_', m, "_li.rds")))
    set(tdt, NULL, 'mtimeloc', readRDS(paste0(out_dir_dat, cog_id, '_', m,
      "_mtimeloc.rds")))
    tdt <- merge(tdt, cellpop_15as, by='cellnum_15as')
    if(nrow(tdt)==0) next
    set(tdt, NULL, 'lirescale', sqrt(tdt[['li']]/.005))
    set(tdt, NULL, 'li',NULL)
    set(tdt, NULL, 'locdatechar', format.POSIXct(tdt[['mtimeloc']], format='%Y-%m-%d',
      tz='UTC'))
    set(tdt, NULL, 'timehour', hour(tdt[['mtimeloc']]) + minute(tdt[['mtimeloc']])/60 +
      second(tdt[['mtimeloc']])/3600)
    set(tdt, NULL, 'timefix', ((tdt[['timehour']]-23)%%24)/5)
    tdt <- tdt[timefix<=1]
    set(tdt, NULL, 'mtimeloc', NULL)
    set(tdt, NULL, 'timehour', NULL)
    set(tdt, NULL, 'monthfac', factor(as.integer(substr(tdt[['locdatechar']],6L,7L)),
      levels=1:12, labels=month.abb))
    set(tdt,NULL,'r9pred',predict(mod, newdata=tdt, re.form = NULL, random.only=FALSE,
      type='response', allow.new.levels=TRUE))
    set(tdt,NULL,'r9rs',(tdt[['r9']]-tdt[['r9pred']])/mod_sig)
    tdt <- tdt[,list(cellnum_15as,r9,r9pred,r9rs)]
    predDT[[i]] <- tdt
    rm(tdt);gc()
  }
  predDT <- rbindlist(predDT, fill=TRUE)
  
  # Remove outliers
  setkey(predDT,cellnum_15as)
  predDT[,r9med:=as.numeric(median(r9)), by='cellnum_15as']
  predDT[,r9sd:=sd(r9), by='cellnum_15as']
  predDT <- predDT[r9<r9med+6*r9sd]
  
  # Create HREA metrics
  predDT[, zmu:=mean(r9rs), by='cellnum_15as']
  set(predDT, NULL, 'lit', predDT[['r9rs']] > qnorm(.9))
  predDT[, prplit:=mean(lit), by='cellnum_15as']
  set(predDT, NULL, 'lightscore', (pnorm(predDT[['zscore']])-.5)/.5)
  set(predDT, predDT[,.I[lightscore<0]], which(names(predDT)=='lightscore'), 0)
  
  # Subset and save file
  predDT <- predDT[, list(cellnum_15as, zmu, prplit, lightscore)]
  predDT <- unique(predDT, by='cellnum_15as')
  saveRDS(predDT, paste0(out_dir_dat, '/', cog_id, '_HREA_', year, '.rds'))
  
  return(NULL)
}


# ****************************************************************************************
# Name:    make_HREA_cogs
# Purpose: Creates Cloud-Optimized GeoTIFFs of HREA data by cog_id and year.
# Inputs:  focal_country_iso (character, required) - 3 letter ISO code for country
#          cog_id (character, required) - COG ID of unit
#          year (integer, required) - year of HREA regression
# Output:  NULL on success
# Notes:   Run at the cog_id-year level. Requires that HREA statistics have been
#          generated.
# ****************************************************************************************

make_HREA_cogs <- function(focal_country_iso=NA_character_, cog_id=NA_character_,
                  year=NA_integer_) {
  # Load necessary packages
  require(terra)
  require(data.table)
  
  # Stop function if missing necessary inputs
  if(is.na(focal_country_iso3)) stop('Requires the country ISO3.')
  if(is.na(cog_id)) stop('Requires the country COG ID.')
  # Make directories for files
  out_dir_dat <- paste0('data', focal_country_iso3, sep = '/')
  dir.create(out_dir_dat, FALSE, TRUE)
  out_dir_cog <- paste0('data', focal_country_iso3, cog_id, sep = '/')
  dir.create(out_dir_cog, FALSE, TRUE)
  
  # Read empty raster grids
  rast_hr <- rast(paste0(out_dir_dat, '/', cog_id, '_empty_grid.tif'))
  rast_15as <- rast(paste0(out_dir_dat, '/', cog_id, '_empty_grid_15as.tif'))
  # Read cell info
  cellpop_hr <- readRDS(paste0(out_dir_dat, '/', cog_id, '_cell_pop.rds'))
  
  # Read HREA data
  predDT <- readRDS(paste0(out_dir_dat, '/', cog_id, '_HREA_', year, '.rds'))
  predDT[zmu<=-999999, zmu:=-999998.9]
  # Merge with settlement cell info
  predDT <- merge(predDT, cellpop_hr, by='cellnum_15as')
  rm(cellpop_hr);gc()
  
  # Make mean z-score COG
  rast_hr[predDT[['cellnum_hr']]] <- predDT[['zmu']]
  writeRaster(rast_hr, paste0(out_dir_cog, '/', cog_id, '_set_zscore_', year, '.tif'),
    format="GTiff", overwrite=TRUE, dataType='FLT4S', setStatistics=TRUE,  NAflag=-999999, 
    options=c("COMPRESS=LZW"))
  
  # Make lightscore COG
  rast_hr[predDT[['cellnum_hr']]] <- predDT[['lightscore']]
  writeRaster(rast_hr, paste0(out_dir_cog, '/', cog_id, '_set_lightscore_', year, '.tif'),
    format="GTiff", overwrite=TRUE, dataType='FLT4S', setStatistics=TRUE,  NAflag=-999999, 
    options=c("COMPRESS=LZW"))
  
  # Make proportion of nights lit COG
  rast_hr[predDT[['cellnum_hr']]] <- predDT[['zmu']]
  writeRaster(rast_hr, paste0(out_dir_cog, '/', cog_id, '_set_prplit_', year, '.tif'),
    format="GTiff", overwrite=TRUE, dataType='FLT4S', setStatistics=TRUE,  NAflag=-999999, 
    options=c("COMPRESS=LZW"))
  rm(rast_hr, predDT);gc()
  
  # Specify months
  if (year>=2013L & year<=2017L) {
    YYYYMMs <- paste0(year, formatC(1:12, width = 2, format = "d", flag = "0"))
  } else {
    YYYYMMs <- paste0('npp_', year ,formatC(1:12, width = 2, format = "d", flag = "0"))
  }
  nYYYYMMs <- length(YYYYMMs)
  
  # Collect annual data
  avgDT <- vector('list', nYYYYMMs)
  for(i in 1:nYYYYMMs){
    m <- YYYYMMs[i]
    if(!file.exists(paste0(out_dir_dat, cog_id, m, "_cellnum_15as.rds"))){
      predDT[[i]] <- data.table(NULL)
      next
    }
    tdt <- data.table(cellnum_15as=readRDS(paste0(out_dir_dat, cog_id, '_', m,
           "_cellnum_15as.rds")))
    set(tdt, NULL, 'r9', readRDS(paste0(out_dir_dat, cog_id, '_', m, "_rade9.rds")))
    set(tdt, NULL, 'li', readRDS(paste0(out_dir_dat, cog_id, '_', m, "_li.rds")))
    set(tdt, NULL, 'sam', readRDS(paste0(out_dir_dat, cog_id, '_', m, "_samples.rds")))
    tdt <- tdt[li<=.0005]
    tdt <- tdt[sam>=225&sam<=3840]
    set(tdt, NULL, 'r9l', log(tdt[['r9']]+2.5))
    set(tdt, NULL, 'li', NULL)
    set(tdt, NULL, 'sam', NULL)
    set(tdt, NULL, 'r9', NULL)
    predDT[[i]] <- tdt
    rm(tdt);gc()
  }
  avgDT <- rbindlist(avgDT, fill=TRUE)
  setkey(avgDT, cellnum_15as)
  avgDT[, r9lmu:=mean(r9l), by='cellnum_15as']
  set(avgDT, NULL, 'r9l', NULL)
  avgDT <- unique(avgDT, by='cellnum_15as')
  
  # Make average radiance raster
  rast_15as[avgDT[['cellnum_15as']]] <- avgDT[['r9lmu']]
  writeRaster(rast_15as, paste0(out_dir_cog, '/', cog_id, '_rade9lnmu_', year, '.tif'),
    format="GTiff", overwrite=TRUE, dataType='FLT4S', setStatistics=TRUE,  NAflag=-999999, 
    options=c("COMPRESS=LZW"))
  rm(rast_15as, avgDT);gc()
  
  # Return NULL
  return(NULL)
}
