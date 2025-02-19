country <- "Guinea"
gadm_abbrev <- "GIN"
# Load libraries and info ----------------------------------------------------------
options(gsubfn.engine = "R")
library(rgdal)
options(warn=0)
library(spdep)
library(SUMMER)
library(geosphere)
library(stringr)
library(tidyverse)
#devtools::install_github("ropensci/rdhs")
library(rdhs)
library(terra)

#### GADM DATA #################################################################
poly_layer_adm0 <- paste('gadm41', gadm_abbrev,
                         '0', sep = "_")
poly_layer_adm1 <- paste('gadm41', gadm_abbrev,
                         '1', sep = "_")
poly_layer_adm2 <- paste('gadm41', gadm_abbrev,
                         '2', sep = "_")
poly_path <- "data/Guinea/gadm41_GIN_shp/"

poly_adm0 <- st_read(dsn = poly_path,
                     layer = as.character(poly_layer_adm0)) 
# use encoding to read special characters
poly_adm1 <-  st_read(dsn = poly_path,
                      layer = as.character(poly_layer_adm1)) 
if(sum(grepl(paste('gadm41', gadm_abbrev,
                   '2', sep = "_"), list.files(poly_path))) != 0){
  poly_adm2 <-  st_read(dsn = poly_path,
                        layer = as.character(poly_layer_adm2)) 
}
if (exists("poly_adm2")) {
  st_crs(poly_adm0) <- st_crs(poly_adm1)  <- st_crs(poly_adm2)
}else {
  st_crs(poly_adm0) <- st_crs(poly_adm1)
}
poly_adm2$NAME_1 <- as.character(poly_adm2$NAME_1)
poly_adm1$NAME_1 <- as.character(poly_adm1$NAME_1)

#### WORLDPOP POPULATION DATA ##################################################
frame_year <- 2017 # sampling frame based on 2017 census
svy_year <- 2018
pop_abbrev <- "gin"

for (year in c(frame_year, svy_year)) {
  options(timeout=1000)
  file <- paste0("data/", country, "/Population/", pop_abbrev, '_ppp_', year, '_1km_Aggregated_UNadj.tif')
  
  if(!file.exists(file)){
    url <- paste0("https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj/", 
                  year, "/", toupper(pop_abbrev),"/",      
                  pop_abbrev,'_ppp_', year,'_1km_Aggregated_UNadj.tif')
    
    download.file(url, file, method = "libcurl", mode="wb")
  }
  for (s in c("f", "m")) {
    file <- paste0("data/", country, "/Population/", pop_abbrev, '_', s, '_1_', year,'.tif')
    if(!file.exists(file)){
      url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020/",
                    year, "/", toupper(pop_abbrev),"/",      
                    pop_abbrev, '_', s, '_1_', year,'.tif')
      download.file(url, file, method = "libcurl", mode="wb")
      
    }
  }
}
file <- paste0("data/", country, "/Covariates/", pop_abbrev,"_viirs_100m_2016.tif")
if(!file.exists(file)){
  url <- "https://data.worldpop.org/GIS/Covariates/Global_2000_2020/GIN/VIIRS//gin_viirs_100m_2016.tif"
  download.file(url, file, method = "libcurl", mode="wb")
}
#### COVARIATE DATA ############################################################
# access raster
access <- rast("data/Nigeria/2015_accessibility_to_cities_v1.0.tif")


# nighttime lights raster
ntl <- rast(paste0("data/", country, "/Covariates/",
                      pop_abbrev,"_viirs_100m_2016.tif"))

pop_data_folder <- "data/Guinea/Population/"
cen_f1_pop <- rast(paste0(pop_data_folder, pop_abbrev, "_f_1_2017.tif"))
cen_m1_pop <- rast(paste0(pop_data_folder, pop_abbrev, "_m_1_2017.tif"))
cen_1_5_pop <- cen_f1_pop + cen_m1_pop

svy_f1_pop <- rast(paste0(pop_data_folder, pop_abbrev, "_f_1_2018.tif"))
svy_m1_pop <- rast(paste0(pop_data_folder, pop_abbrev, "_m_1_2018.tif"))
svy_1_5_pop <- svy_f1_pop + svy_m1_pop

# Total population rasters
# worldpop estimates for year of census
cen_pop <- rast(paste0(pop_data_folder, pop_abbrev, '_ppp_2017_1km_Aggregated_UNadj.tif'))

# worldpop estimates for year of survey
svy_pop <- rast(paste0(pop_data_folder, pop_abbrev, '_ppp_2018_1km_Aggregated_UNadj.tif'))
cen_1_5_pop <- aggregate(cen_1_5_pop, 10)
cen_1_5_pop <- project(cen_1_5_pop, cen_pop)
svy_1_5_pop <- aggregate(svy_1_5_pop, 10)
svy_1_5_pop <- project(svy_1_5_pop, cen_pop)
pop_stk <- c(cen_1_5_pop, svy_1_5_pop, cen_pop, svy_pop)
pop_stk <- c(pop_stk, project(access, pop_stk))
pop_stk <- c(pop_stk, project(ntl, pop_stk))
names(pop_stk) <- c("cen_1_5_pop", "svy_1_5_pop", "cen_pop", "svy_pop", "access", "ntl")

writeRaster(pop_stk, "data/Guinea/pop_stk_1km.tif", overwrite = TRUE)


#### DHS DATA ##################################################################
# set API to get DHS data -- you will need to change this to your information!
set_rdhs_config(email = "petergao@uw.edu",
                project = "Small Area Estimation with Covariates")

# get country ID
surveys <- dhs_datasets(countryIds = "GN", surveyYearStart = 2000) %>%
  dplyr::filter(
    (FileType == 'Births Recode' & 
       FileFormat=='Stata dataset (.dta)') | 
      (FileType == 'Geographic Data' & 
         FileFormat =='Flat ASCII data (.dat)')
    )
dhs_paths <- get_datasets(surveys[surveys$SurveyYear==2018,]$FileName, clear_cache = T)
in_dat <- readRDS(paste0(dhs_paths[1]))
geo_dat <- readRDS(paste0(dhs_paths[2]))
svy_dat = data.frame(
  cluster = in_dat$v001, 
  hshold = in_dat$v002,
  stratum = in_dat$v022, 
  mcv_yes = 1 * (in_dat$h9 == 1 | # "vaccination date on card"
                   in_dat$h9 == 2 | # "reported by mother" 
                   in_dat$h9 == 3), # "vaccination marked on card"
  alive = in_dat$b5 == 1,
  age = in_dat$b19,
  wt = in_dat$v005/1000000
)
svy_dat <- subset(svy_dat, age <= 23 & age >= 12)
svy_dat <- subset(svy_dat, !is.na(mcv_yes))
svy_dat = subset(svy_dat, alive)

svy_dat %>% group_by(stratum) %>% summarize(mcv = mean(mcv_yes), n = n())
  
ea_locs_path <- "data/Guinea/GNGE71FL/"
ea_locs <- st_read(ea_locs_path)
ea_dat <- data.frame(cluster = ea_locs$DHSCLUST,
                     urban = ea_locs$URBAN_RURA,
                     lon = ea_locs$LONGNUM,
                     lat = ea_locs$LATNUM)
svy_dat = merge(svy_dat, ea_dat, by = "cluster")



svy_dat <- st_as_sf(svy_dat, coords = c("lon", "lat"))
st_crs(svy_dat) <- st_crs(poly_adm0)
svy_dat <- st_join(svy_dat,
                   select(poly_adm1, NAME_1))
svy_dat <- st_join(svy_dat,
                   select(poly_adm2, NAME_2))


svy_dat$admin1 <- match(svy_dat$NAME_1, poly_adm1$NAME_1)
svy_dat$admin1_char <- paste0("admin1_", svy_dat$admin1)
svy_dat$admin1_name <- poly_adm1$NAME_1
svy_dat$admin2 <- match(svy_dat$NAME_2, poly_adm2$NAME_2)
svy_dat$admin2_char <- paste0("admin2_", svy_dat$admin2)
svy_dat$admin2_name <- poly_adm1$NAME_2
saveRDS(svy_dat,
        file = paste0("data/Guinea/Guinea_DHS_2018_MCV.rds"))




#### CONVERT FROM RASTER TO POINT-REFERENCED DATA ####
#### CALCULATE URBAN/RURAL FRACTION ####
# this function computes the population density threshold t such that 
# pixels with population density > t are classified as urban
get_thresh <- function(this_area_pop, adm1, poppa) {
  # based on John Paige's getRegionThresh from SUMMER package
  prop_urb <- poppa$propUrb[poppa[["GADMarea"]] == adm1] 
  prop_rural <- 1 - prop_urb
  if (round(prop_rural, 10) == 1) {
    return(-Inf)
  } else if (round(prop_rural, 10) == 0) {
    return(Inf)
  }
  area_pop <- this_area_pop
  total_area_pop <- sum(area_pop)
  
  area_pop_sorted = sort(area_pop)
  cumsum_area_pop = cumsum(area_pop_sorted)
  thresh_idx = match(1, cumsum_area_pop >= total_area_pop * prop_rural)
  
  if ((thresh_idx != 1) & (thresh_idx != length(area_pop))) {
    thresh = area_pop_sorted[thresh_idx]
  } else {
    # make sure not all pixels are urban or all are rural
    if(thresh_idx == 1) {
      thresh = mean(c(area_pop_sorted[1], area_pop_sorted [2]))
    } else {
      thresh = mean(c(area_pop_sorted [length(area_pop)], 
                      area_pop_sorted[length(area_pop) - 1]))
    }
  }
  return(thresh)
}

# return both population summary as well as urban/rural raster 
stk_to_pix_dat <- function(pop_stk, poly_adm1, poly_adm2, svy_poppa, path) {

  # make pixel table
  adm1_names <- as.character(poly_adm1$NAME_1)
  
  adm1_pix_list <- list()
  adm1_pop_totals_list <- list()
  for (i in 1:length(adm1_names)) {
    print(i)
    adm1 <- adm1_names[i]
    print(adm1)
    if (!file.exists(paste0(path, "gin_pix_tbl_admin1_", i, ".rds"))) {
      adm1_pix <-
        mask(crop(pop_stk, subset(poly_adm2, NAME_1 == adm1), snap='out'),
             vect(subset(poly_adm2, NAME_1 == adm1)))
      adm1_pix <- st_join(st_as_sf(as.points(adm1_pix)),
                          st_as_sf(poly_adm2)[, c("NAME_1", "NAME_2")])
      adm1_pix <- rename(adm1_pix, "admin1_name" = "NAME_1")
      adm1_pix <- rename(adm1_pix, "admin2_name" = "NAME_2")
      adm1_pix <- adm1_pix[adm1_pix$admin1_name == adm1,]
      adm1_pix$cen_pop[is.na(adm1_pix$cen_pop)] <- 0
      # USE SURVEY POPULATION TO GET URBAN/RURAL THRESHOLD
      thresh <- get_thresh(adm1_pix$cen_pop, adm1, svy_poppa)
      adm1_pix$urban <- adm1_pix$cen_pop >= thresh

      adm1_pix_coords <- st_coordinates(st_as_sf(adm1_pix))
      adm1_pix$lon <- adm1_pix_coords[, 1] 
      adm1_pix$lat <- adm1_pix_coords[, 2]
      
      adm1_pix$access <- ifelse(is.na(adm1_pix$access), mean(adm1_pix$access, na.rm = T), adm1_pix$access)
      adm1_pix$l1a <- log(1 + adm1_pix$access)
      adm1_pix$ntl <- ifelse(is.na(adm1_pix$ntl), mean(adm1_pix$ntl, na.rm = T), adm1_pix$ntl)
      adm1_pix <- st_set_geometry(adm1_pix, NULL)
      adm1_pix <- adm1_pix[complete.cases(adm1_pix), ]
      saveRDS(adm1_pix, 
              paste0(path, "gin_pix_tbl_admin1_", i, ".rds"))
    } else {
      adm1_pix <- readRDS(paste0(path, "gin_pix_tbl_admin1_", i, ".rds"))
    }
    adm1_pix_list[[i]] <- adm1_pix
    adm1_pop_totals <-
      aggregate(.~admin1_name + admin2_name,
                dplyr::select(adm1_pix, -lon, -lat, -urban, -access, -ntl), sum)
    adm1_urban_pop_totals <-
      aggregate(.~admin1_name + admin2_name + urban,
                dplyr::select(adm1_pix, -lon, -lat, -access, -ntl), sum)
    adm1_urban_pop_totals <- subset(adm1_urban_pop_totals, urban == T)
    adm1_urban_pop_totals <- dplyr::select(adm1_urban_pop_totals, -urban)
    colnames(adm1_urban_pop_totals)[3:6] <- 
      paste0(colnames(adm1_urban_pop_totals)[3:6], "_urb")
    adm1_pop_totals <- merge(adm1_pop_totals, 
                             adm1_urban_pop_totals[, ], all.x = T)
    adm1_pop_totals[is.na(adm1_pop_totals)] <- 0
    adm1_pop_totals_list[[i]] <- adm1_pop_totals
    
    
  }
  saveRDS(do.call(rbind, adm1_pix_list), 
          paste0(path, "gin_pix_tbl.rds"))
  adm1_pop_totals <- do.call(rbind, adm1_pop_totals_list)
  adm1_pop_totals <- aggregate(.~admin1_name + admin2_name, data = adm1_pop_totals, sum)
  saveRDS(adm1_pop_totals, 
          paste0(path, "gin_pop_totals.rds"))
}



svy_poppa <- read.csv("data/Guinea/gin_2018_poppa.csv")
stk_to_pix_dat(pop_stk, poly_adm1, poly_adm2, svy_poppa, "data/Guinea/pix_tbl_1km/")

