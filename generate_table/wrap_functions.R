library(data.table)
library(readxl)

shift <- data.table::shift
## load and merge habitat table with track by time
merge_habitat_track <- function(habitat, habitatType, locations, track) {
  # calculate distance
  track[, distance := sqrt((shift(Lon) - Lon)^2 + (shift(Lat) - Lat)^2 + (shift(Depth) - Depth)^2), by = .(Dive)]

  # merge habitat with habitatType
  habitat[habitatType, c('name', 'short_name', 'Subsystem', 'Substrate') := .(name, sigla, localidade, substrato), on = "habitat"]

  # remove unwanted columns
  habitat[, c('filename', 'comments', 'time')] <- NULL
  names(habitat)[2] <- "start_time"

  # create a column with end time of each habitat
  habitat[, duration := shift(start_time, -1) - start_time, by = .(station)]
  habitat <- habitat[, end_time := start_time + duration]
  habitat <- habitat[!is.na(duration)]

  # expand rows based on start and end times, jumping one second each row
  habitat <- habitat[, .SD[rep(1:.N, duration)]][,
                                                Date := seq(start_time, end_time - 1, by = 'sec'),
                                                by = .(station, start_time, end_time)][]

  # remove columns
  habitat[, time := as.character(Date, format = '%H:%M:%S')]
  habitat <- habitat[!is.na(name)]
  habitat[, c('start_time', 'habitat', 'duration', 'end_time', 'Date')] <- NULL

  # merge habitat and track
  habitat[track, c('distance', 'lon', 'lat') := .(distance, Lon, Lat), on = .('station' = Dive, 'time' = Time)]
  habitat$distance[is.na(habitat$distance)] <- 0

  habitat[locations, Region := location, on = "station"]

  return(habitat)
}


## load and merge the annotation table
merge_annotation <- function(annotations, morphotypes, locations) {
  # remove unwanted columns
  annotations[, c('filename', 'substrato', 'coments', 'imagelink')] <- NULL
  annotations[, time := as.character(time, format = '%H:%M:%S')]

  # merge habitat with annotations
  annotations[habitat, c('habitat', 'habitat_code') := .(name, short_name), on = .(station, time)]


  # merge morphotypes with annotations
  annotations[morphotypes, c('phylum', 'class', 'order', 'group', 'code', 'OTU', 'zone') :=
                .(phylum, class, order, group, code, taxonName, zone), on = .(morphotype)]

  # remove unwanted annotations
  annotations <- annotations[group != 'ignorar' & zone == 'b�ntico']
  annotations$zone <- NULL


  # merge locations with annotations
  annotations[locations, region := location, on = 'station']

  return(annotations)
}


ceiling_dec <- function(x, level=3) round(x + 5 * 10^(-level - 1), level)

divide_sample_helper <- function(distance, Ninds, thr_cor, dir, thr) {
  if (dir == "F") {
    sample_id <- cumsum(distance) %/% thr_cor
  } else {
    sample_id <- rev(cumsum(rev(distance)) %/% thr_cor)
  }
  df <- data.table(distance, Ninds, sample_id)
  df <- df[, lapply(.SD, sum), by = .(sample_id), .SDcols = c('distance', 'Ninds')]
  df <- df[distance >= thr * min_thr]
  ninds <- sum(df$Ninds)
  return(ninds)
}


divide_sample <- function(distance, thr, Ninds) {
  total_disntace <- sum(distance)
  nsamples <- total_disntace %/% thr

  if (nsamples == 0) { # if distance is lower than thr
    return("00")
  }

  # calculate number of samples
  thr_cor_max <- ceiling_dec(min(total_disntace / nsamples, thr * max_thr + .5))
  thr_cor_min <- ceiling_dec(max(total_disntace / (nsamples + 1), thr * min_thr + .5))
  Ninds[is.na(Ninds)] <- 0

  opcoes <- expand.grid(thr_cor = c(thr_cor_max, thr_cor_min), dir = c('F', 'R'), ninds = 0)

  for (i in 1:nrow(opcoes)) {
    opcoes$ninds[i] <- divide_sample_helper(distance, Ninds, opcoes$thr_cor[i], opcoes$dir[i], thr)
  }

  # check which orientation will get more abundance
  i <- which.max(opcoes$ninds)
  if (opcoes$dir[i] == "F") {
    sample_id <- cumsum(distance) %/% opcoes$thr_cor[i]
  } else {
    sample_id <- rev(cumsum(rev(distance)) %/% opcoes$thr_cor[i])
  }


  # return sample id
  return(sprintf("%02d", sample_id))

}

## create samples along tracks with confirmed distances
create_sample <- function(habitat, annot, thr, min_ninds) {

  # get a sample id by each station and each habitat
  habitat2 <- copy(habitat)
  habitat2[, sample_id := sprintf("%03d", rleid(station, short_name))]

  annot2 <- annot[, .(Ninds = sum(Ninds)), by = .(station, time)]
  habitat2[annot2, Ninds_tmp := Ninds, on = .(station, time)]

  # divide samples id by the distance
  habitat2[, sample := paste(sample_id, divide_sample(distance, thr, Ninds_tmp), sep = '-'), by = .(sample_id)]

  dists <- habitat2[, .(dists = sum(distance)), by = .(sample)]

  annot2 <- habitat2[annot, .(sample, code, Ninds, lon, lat), on = .(station, time)]

  tmp <- annot2[, .(Nspp = length(unique(code)), Ninds = sum(Ninds), lon = mean(lon), lat = mean(lat)), by = .(sample)]
  tmp[dists, dists := dists, on = .(sample)]

  tmp[habitat2, c('station', 'habitat', 'region') := .(station, short_name, Region), on = .(sample)]
  tmp[, keep := dists >= thr * min_thr & Ninds >= min_ninds]

  annot2 <- copy(annot)
  annot2[habitat2, c('sample', 'Subsystem', 'Substrate') := .(sample, Subsystem, Substrate), on = .(station, time)]

  return(list(annot2, tmp, habitat2))
}


## get midpoint of the line in one segment
get_mids <- function(lon, lat) {
  distance <- sqrt((diff(lon)^2 + (diff(lat))^2))
  coords <- data.frame(lon = lon, lat = lat)
  dist_mid <- sum(distance) / 2
  dist_cum <- cumsum(distance)
  end_index <- which(dist_cum > dist_mid)[1]
  end <- coords[end_index, ]
  start_index <- end_index - 1L

  if (distance[start_index] == 0) {
    return(list(lon = end[1, 1], lat = end[1, 2]))
  } else {
    start <- coords[start_index, ]
    dist_remaining <- dist_mid - dist_cum[start_index]
    mid <- start + (end - start) * (dist_remaining / distance[start_index])
    return(list(lon = mid[1, 1], lat = mid[1, 2]))
  }
}
