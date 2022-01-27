library(data.table)
library(readxl)
library(ggplot2)

# Script to load raw annotations and track and generate segments using multiple distance thresholds

# load raw data (currently not available in this repository)
tabela <- paste0(Sys.getenv("GoogleDrive"), "\\Doutorado\\Videos HyBIS\\tabela.xlsm")

habitat_raw <- setDT(read_excel(tabela, sheet = "habitat"))
habitatType_raw <- setDT(read_excel(tabela, sheet = "habitatType"))
annotations_raw <- setDT(read_excel(tabela, sheet = "main"))
morphotypes_raw <- setDT(read_excel(tabela, sheet = "morphotypes"))
locations_raw <- setDT(read_excel(tabela, sheet = "locations"))

track <- fread(paste0(Sys.getenv("GoogleDrive"), "\\Doutorado\\Videos HyBIS\\track clean\\python\\track_clean.csv"))
track[, Time := as.character(Time, format = "%H:%M:%S")]

# load functions
source('generate_table/wrap_functions.R')

theme_artigo <- function() {
  theme_classic(base_size = 8, base_line_size = 0.25) %+replace%    #replace elements we want to change
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = rel(0.9))
    )
}


# get habitat and annotations table
habitat <- merge_habitat_track(habitat_raw, habitatType_raw, locations_raw, track)
annot <- merge_annotation(annotations_raw, morphotypes_raw, locations_raw)


# set parameters to define segments
min_ninds <- 20 # segments below this number of individuals will be rejected
distances <- seq(50, 300, 10) # distance thresholds
min_thr <- 0.8 # minimal distance tolerance for the segments
max_thr <- 1.0 # maximum distance tolerance



# iterate for each distance threshold
results <- data.frame(w = distances, # the thr used
                     max_richness = 0, # max richness of the samples
                     min_richness = 0, # min richness of the samples
                     n_samples = 0, # number of samples
                     n_samples_min_inds = 0, # number of samples above min_ninds
                     n_inds = 0, # number of inds
                     distance = 0) # distance that was considered in the samples

for (i in seq_along(distances)) { # for each thr
  # get samples ids
  samples <- create_sample(habitat, annot, distances[i], min_ninds)
  tmp <- samples[2][[1]]

  # select samples above thr
  tmp3 <- tmp[keep == TRUE]

  # store results
  results$max_richness[i] <- max(tmp3$Nspp)
  results$median_richness[i] <- median(tmp3$Nspp)
  results$min_richness[i] <- min(tmp3$Nspp)
  results$n_samples[i] <- nrow(tmp)
  results$n_samples_min_inds[i] <- nrow(tmp3)
  results$n_inds[i] <- sum(tmp3$Ninds)
  results$distance[i] <- sum(tmp3$dists)
  results$n_habitat[i] <- sum(table(tmp3$habitat) > 0)
  results$n_station[i] <- sum(table(tmp3$station) > 0)

}


thr <- 120 # threshold to be displayed in the plot

# PLOT - Figure 2 of the paper

g <- ggplot(results, aes(w, n_inds / 1000)) + geom_line(size = 0.3) +
  geom_vline(xintercept = thr, color = 'gray', linetype = "dashed", size = 0.5) +
  scale_x_continuous(breaks = seq(50, 300, 50)) + xlab("Segment Size (m)") +
  geom_point(size = 0.9) + theme_artigo() + ylab("Total abundance (x1000)") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))

ggsave("Figure 2.svg", g, height = 6, width = 9, units = 'cm')



# correlation of Ninds and Nspp by distance - should not be significant
samples <- create_sample(habitat, annot, thr, min_ninds)
tmp <- samples[2][[1]]
tmp3 <- tmp[keep == TRUE]


summary(lm(Ninds ~ dists, tmp3))
summary(lm(Nspp ~ dists, tmp3))
