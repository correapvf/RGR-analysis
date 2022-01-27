library(data.table)

# Figure 4 - richness vs distance
random_acc <- function(code, nperm = 999) {
  results <- replicate(nperm, cumsum(!duplicated(sample(code))))
  return(apply(results, 1, mean))
}

annotations <- fread('long_benthic.csv')
annotations$habitat_code <- factor(annotations$habitat_code, levels =
                                    c('RiftDeb', 'RiftSed', 'RiftHole', 'RiftRock', 'RiftCr', 'Sed', 'CalPav', 'CrSed', 'CrPav', 'Cr', 'CrCob'))
annotations$region <- factor(annotations$region, levels = c('Northeast', 'Center', 'South', 'Southwest'))
setorder(annotations, Subsystem, habitat_code, dist_habitat)


annotations[, Nspp_acum := random_acc(code), by = .(habitat_code)]
annotations2 <- annotations[, .(Subsystem, dist_habitat, Nspp_acum, habitat_code)]

fwrite(annotations2[Subsystem == 'Plateau', !"Subsystem"], "richness_by_distance_plateau.csv")
fwrite(annotations2[Subsystem != 'Plateau', !"Subsystem"], "richness_by_distance_rifte.csv")



# make morphotype table
morphs <- annotations[, .(.N, taxonName = first(OTU)), by = .(phylum, code)]
write.csv(morphs, "morphotypes.csv", row.names = FALSE, quote = FALSE)
