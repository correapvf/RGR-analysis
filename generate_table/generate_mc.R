
# Script to load raw annotations, track, multibeam and CTD data and generate mc and envs files used in the analysis

#### load raw data ####
# (currently not available in this repository)
library(readxl)
library(data.table)

tabela <- paste0(Sys.getenv("GoogleDrive"), "\\Doutorado\\Videos HyBIS\\tabela.xlsm")

habitat_raw <- setDT(read_excel(tabela, sheet = "habitat"))
habitatType_raw <- setDT(read_excel(tabela, sheet = "habitatType"))
annotations_raw <- setDT(read_excel(tabela, sheet = "main"))
morphotypes_raw <- setDT(read_excel(tabela, sheet = "morphotypes"))
locations_raw <- setDT(read_excel(tabela, sheet = "locations"))

track <- fread(paste0(Sys.getenv("GoogleDrive"), "\\Doutorado\\Videos HyBIS\\track clean\\python\\track_clean.csv"))
track[, Time := as.character(Time, format = "%H:%M:%S")]

# set parameters to define segments - as selected in the choose-th.R script
min_ninds <- 20
thr <- 120
min_thr <- 0.8
max_thr <- 1.0


# get habitat and annotations table
source('generate_table/wrap_functions.R')
habitat <- merge_habitat_track(habitat_raw, habitatType_raw, locations_raw, track)
annot <- merge_annotation(annotations_raw, morphotypes_raw, locations_raw)


samples <- create_sample(habitat, annot, thr, min_ninds)
annotations <- samples[1][[1]]
tmp <- samples[2][[1]]
habitat <- samples[3][[1]]


# remove all segmentes with distance lower than thr and small Ninds
tmp2 <- tmp[keep == TRUE]

# get median lat and long of all files
track_median <- habitat[tmp2, get_mids(lon, lat), by = .(sample), on = .(sample)]


# get cumulative distance for the annotation table
habitat2 <- copy(habitat)
habitat2 <- habitat2[, dist_total := cumsum(distance)]
habitat2 <- habitat2[, dist_habitat := cumsum(distance), by = .(short_name)]

# SAVE LONG DATA
annotations[habitat2, c('dist_total', 'dist_habitat') := .(dist_total, dist_habitat), on = c('station', 'time')]
annotations$morphotype <- NULL
write.csv(annotations, 'long_benthic.csv', quote = FALSE, row.names = FALSE)


#### create a community matrix ####
# convert annotation table to community matrix
long <- dcast(annotations, sample ~ code, fun.aggregate = sum, value.var = 'Ninds')

# remove all 'samples' with distance lower than thr
long <- long[tmp2, .SD, on = .(sample)]

# create the community matrix
mc <- long[, -1]

# remove spp that has no observations
mc <- mc[, apply(mc, 2, sum) > 0, with = FALSE]

# write csv
write.csv(mc, 'mc.csv', quote = FALSE, row.names = FALSE)



#### get data from rasters ####
library(raster)
rasters_name <- c("depth", "slope", "broadbpi_20_200", "finebpi_3_30",
            "cos_thr2", "sin_thr2", "curvature", "vrm_011", "backscatter")

rasters <- paste0(Sys.getenv("GoogleDrive"), "\\Doutorado\\images and data\\BTM_DY94\\", rasters_name, ".tif")
r <- lapply(rasters, raster)

# extract rasters values
envs <- sapply(r, raster::extract, y = habitat[, c('lon', 'lat')])
colnames(envs) <- rasters_name

habitat <- cbind(habitat, envs)



# get environment values by sample
habitat_group <- habitat[, lapply(.SD, mean), by = .(sample), .SDcols = 14:22]

## get the environments
envs <- habitat_group[tmp2, .SD, on = .(sample)]

# join medianlat and long
envs[track_median, c("lat", "lon") := .(lat, lon), on = .(sample)]
envs$sample <- NULL

# factors environments
envs_f <- annotations[, lapply(.SD, first), by = .(sample), .SDcols = c(1, 2, 12, 14, 15, 4, 5)]
envs_f <- envs_f[tmp2, .SD, on = .(sample)]
envs <- cbind(envs_f, envs)

setcolorder(envs, c("sample", "station", "time", "lon", "lat"))




#### add CTD measurements ####
files <- list.files(path = paste0(Sys.getenv("GoogleDrive"), "\\Doutorado\\cruzeiro RGR1\\CTD"),
                   pattern = "\\.csv$")

files <- files[as.numeric(substr(files, 3, 5)) >= 468] # select CTDs only near HyBIS dives
files <- files[-12] # imcomplete CTD 493

min.depth <- 500
max.depth <- 1500

x <- NULL
for (i in files) {
  tmp <- data.table(read.delim(paste0(Sys.getenv("GoogleDrive"), "\\Doutorado\\cruzeiro RGR1\\CTD\\", i), stringsAsFactors = F))
  colnames(tmp) <- sapply(strsplit(colnames(tmp), '\\.'), `[`, 1)
  tmp <- tmp[Depth >= min.depth]
  tmp <- tmp[Depth <= max.depth]

  tmp$station <- substr(i, 1, 5)
  x <- rbind(x, tmp)
}

x$Depth_round <- round(x$Depth / 10) * 10
x <- x[, lapply(.SD, mean), by = .(station, Depth_round), .SDcols = c(6:8, 10, 11, 13)]

envs$depth_round <- round(envs$depth / 10) * 10


# interpolate CTD data every 10 m depth and store results in env table
library(gstat)
library(sp)
for (d in unique(envs$depth_round)) {
  for (variable in c('Temperature', 'Salinity', 'Oxygen', 'Chlorophyll')) {
    sample <- x[Depth_round == d * -1, c('Latitude', 'Longitude', variable), with = FALSE]
    colnames(sample) <- c('lat', 'lon', 'variable')
    coordinates(sample) <- ~ lon + lat
    proj4string(sample) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
    sample <- spTransform(sample, CRS("+proj=utm +zone=25 +south +datum=WGS84 +units=m +no_defs"))

    index <- envs$depth_round == d
    loc <- envs[index, c('lon', 'lat')]
    coordinates(loc)  <- ~ lon + lat
    proj4string(loc) <- CRS("+proj=utm +zone=25 +south +datum=WGS84 +units=m +no_defs")

    oo <- idw(formula = variable ~ 1, locations = sample, newdata = loc, idp = 2.0)
    envs[index, variable] <- oo@data$var1.pred
  }
}

envs$depth_round <- NULL


# WRITE ENVS MATRIX
write.csv(envs, 'envs.csv', quote = FALSE, row.names = FALSE)


#### Generate tables for paper ####

# Table 1
sec2time <- function(t) {
  h <- t %/% (60 * 60)
  m <- t %/% 60 %% 60
  s <- t %% 60
  return(sprintf("%02d:%02d:%02d", h, m, s))
}

substrRight <- function(x, n) {
  nt <- nchar(x)
  substr(x, nt - n + 1, nt)
}

annot2 <- annot[, .(Ninds = sum(Ninds)), by = .(station, time)]

habitat2 <- copy(habitat)
habitat2$short_name <- factor(habitat2$short_name, levels =
                        c('RiftDeb', 'RiftSed', 'RiftHole', 'RiftRock', 'RiftCr', 'Sed', 'CalPav', 'CrSed', 'CrPav', 'Cr', 'CrCob'))
habitat2[annot2, Ninds := Ninds, on = .(station, time)]

x <- habitat2[, .(habitat = paste(first(name), " (", first(short_name), ")", sep = ""),
           dives = paste(substr(unique(station), 3, 4), collapse = ", "),
           segments = 0,
           time = length(time),
           distance = round(sum(distance), 1),
           Ninds = sum(Ninds, na.rm = TRUE)), by = .(short_name, Region)]

tmp3 <- tmp2[, .(segments2 = .N), by = .(habitat, region)]
x[tmp3, segments := segments2, on = c(short_name = 'habitat', Region = 'region')]

setorder(x, short_name, Region)
setcolorder(x, c("habitat", "Region"))
x$short_name <- NULL

total <- colSums(x[, 4:7])
x <- rbind(x, data.frame('Total', NA, NA, t(total)), use.names = FALSE)

x[, time := sec2time(time)]

write.csv(x, "Table 1.csv", na = "", row.names = FALSE)


# statistics for Table 2
summary_mean <- function(x) {
  return(sprintf("%.01f; %.f-%.f", mean(x), min(x), max(x)))
}

x2 <- habitat2[, .(depth = summary_mean(depth * -1), slope = summary_mean(slope)), by = .(short_name)]
setorder(x2, short_name)
write.csv(x2, "Table 2.csv", row.names = FALSE)
