
# ggplot theme
theme_artigo <- function() {
  theme_classic(base_size = 8, base_line_size = 0.25) %+replace%    #replace elements we want to change
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = rel(0.9)),
      plot.tag = element_text(size = 9.5, face = "bold", hjust = -0.1),
      plot.tag.position = c(0, 1)
    )
}

font_size <- 7 / ggplot2:::.pt # font size of the geom_text

# Dendogram
dendro_data_k <- function(hc, k) {

  hcdata    <-  ggdendro::dendro_data(hc, type = "rectangle")
  seg       <-  hcdata$segments
  labori    <-  cutree(hc, k)
  labclust  <-  labori[hc$order]
  segclust  <-  rep(0L, nrow(seg))
  heights   <-  sort(hc$height, decreasing = TRUE)
  height    <-  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)

  for (i in 1:k) {
    xi      <-  hcdata$labels$x[labclust == i]
    idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3    <-  seg$yend < height
    idx     <-  idx1 & idx2 & idx3
    segclust[idx] <- i
  }

  idx                    <-  which(segclust == 0L)
  segclust[idx]          <-  segclust[idx + 1L]
  hcdata$segments$clust  <-  segclust
  hcdata$segments$line   <-  as.integer(segclust < 1L)
  hcdata$labels$clust    <-  labclust
  hcdata$labels$original <-  labori

  return(hcdata)
}


plot_ggdendro <- function(hcdata,
                          scale.color = NULL,
                          branch.size = 1,
                          expand.y = c(0, 0)) {

  ybreaks   <- pretty(segment(hcdata)$y, n = 5)
  ymax      <- max(segment(hcdata)$y)

  ## branches
  p <- ggplot() +
    geom_segment(data         =  segment(hcdata),
                 aes(x        =  x,
                     y        =  y,
                     xend     =  xend,
                     yend     =  yend,
                     linetype =  "solid",
                     colour   =  factor(clust)),
                 lineend      =  "round",
                 show.legend  =  FALSE,
                 size         =  branch.size)

  ## orientation
  p <- p + scale_x_continuous(breaks = NULL)

  p <- p + scale_y_continuous(breaks = ybreaks, expand = expand.y)


  # colors and limits
  if (!is.null(scale.color)) {
    p <- p + scale_color_manual(values = scale.color)
  }


  return(p)
}


grpdist <- function(X) {
  require(cluster)
  gr <- as.data.frame(as.factor(X))
  distgr <- daisy(gr, "gower")
  return(distgr)
}


# Indval
permute_custom <- function(nr, control, selected_stations) {
  x <- shuffle(nr, control)
  to <- sample(1:nr, length(selected_stations))
  x[to] <- selected_stations
  x[selected_stations] <- to
  return(x)
}


create_other <- function(x2) {
  x2 <- x2[order(-x2)]
  index <- sum(cumsum(x2) <= 80) + 1

  if ((length(x2) - index) > 5) {
    x2 <- head(x2, index)
    x2 <- c(x2, 100 - sum(x2))
    names(x2)[length(x2)] <- "Others"
  }
  return(x2)
}


# Richness plots
richness_plot <- function(mc, perm = 100) {
  pool <- estaccumR(mc, permutations = perm, parallel = 2)
  pool2 <- data.table(Samples = pool$N,
                      Sobs = rowMeans(pool$S),
                      Chao1 = rowMeans(pool$chao),
                      ACE = rowMeans(pool$ace, na.rm = TRUE))

  # correct for a few NAs in the ACE
  rows_NA <- which(apply(pool$ace, 1, function(x) sum(is.na(x))) > 0)
  for (i in rows_NA) {
    pool$ace[i, is.na(pool$ace[i, ])] <- mean(pool$ace[i, ], na.rm = TRUE)
  }

  pool_incidence <- poolaccum(mc, permutations = perm)
  pool_incidence2 <- data.table(Samples = pool_incidence$N,
                                Chao2 = rowMeans(pool_incidence$chao),
                                Jack1 = rowMeans(pool_incidence$jack1),
                                Jack2 = rowMeans(pool_incidence$jack2),
                                Boot = rowMeans(pool_incidence$boot))

  pool_spp <- pool_incidence2[pool2, .(Samples, Sobs, Chao1, Chao2, Jack1, Jack2, ACE, Boot), on = .(Samples)]
  pool_spp <- melt(pool_spp, id.vars = c('Samples'), na.rm = TRUE)
  return(list('pool' = pool, 'pool_incidence' = pool_incidence, 'pool_spp' = pool_spp))
}



theme_artigo2 <- function(font_size = 6.5) {
  theme_artigo() %+replace%    #replace elements we want to change
    theme(
      legend.title = element_blank(),
      panel.grid.major.y = element_line(color = 'gray80'),
      legend.key.height = unit(3, 'mm'),
      legend.justification = c(1, 0),
      legend.position = c(1.0, 0.02)
    )
}


# Table 3
mean_sd1 <- function(x) {
  if (length(x) == 1) {
    return(as.character(x))
  } else {
    return(sprintf("%.2f � %.1f", mean(x), sd(x)))
  }
}

mean_sd2 <- function(x) {
  if (length(x) == 1) {
    return(sprintf("%.2f", x))
  } else {
    return(sprintf("%.2f � %.2f", mean(x), sd(x)))
  }
}

mean_sd3 <- function(x) {
  if (length(x) == 1) {
    return(sprintf("%.1f", x))
  } else {
    return(sprintf("%.1f � %.1f", mean(x), sd(x)))
  }
}


# box-plot with group letters
tri.to.squ <- function(x) {
  rn <- row.names(x)
  cn <- colnames(x)
  an <- unique(c(cn, rn))
  myval <- x[!is.na(x)]
  mymat <- matrix(1, nrow = length(an), ncol = length(an), dimnames = list(an, an))
  for (ext in 1:length(cn)) {
    for (int in 1:length(rn)) {
      if (is.na(x[row.names(x) == rn[int], colnames(x) == cn[ext]])) next
      mymat[row.names(mymat) == rn[int], colnames(mymat) == cn[ext]] <- x[row.names(x) == rn[int], colnames(x) == cn[ext]]
      mymat[row.names(mymat) == cn[ext], colnames(mymat) == rn[int]] <- x[row.names(x) == rn[int], colnames(x) == cn[ext]]
    }

  }
  return(mymat)
}


# MVPART
dendro_data_mvpart <- function(res.part, format_labels = function(x) x) {
  res <- list()
  parms <- list(uniform = FALSE, branch = 1, nspace = -1, minbranch = 0.3)

  # coordinates for segments
  xy <- mvpart:::rpartco(res.part, parms)
  node <- as.numeric(row.names(res.part$frame))
  res$segments <- ggdendro:::rpart_segments(ggdendro:::rpart.branch(xy$x, xy$y, node, parms$branch))
  res$segments$n <- NULL
  res$segments <- rbind(res$segments, data.frame(x = xy$x[1], y = 1.03, xend = xy$x[1], yend = 1))

  # labels for tree
  left.child <- match(2L * node, node)
  left.child <- left.child[!is.na(left.child)]
  right.child <- match(node * 2L + 1L, node)
  right.child <- right.child[!is.na(right.child)]
  node.left <- node[node %% 2L == 0L]
  parent <- match(node.left / 2L, node)

  # get labels and format text
  t <- labels(res.part,  pretty = 1)
  t <- format_labels(t)

  res$labels.left <- cbind(data.frame(xy)[parent, ], label = t[left.child])
  res$labels.right <- cbind(data.frame(xy)[parent, ], label = t[right.child])


  # labels for leaves
  leaves <- res.part$frame$var == "<leaf>"
  res$leaves <- cbind(data.frame(xy)[leaves, ],
                      label = paste0(round(res.part$frame$dev[leaves], 1), " : n=", res.part$frame$n[leaves]),
                      where = (1:nrow(res.part$frame))[leaves])


  # multivariate means
  yval <- as.data.frame(res.part$frame$yval2[res$leaves$where, ])
  colnames(yval) <- colnames(res.part$y)
  yval$where <- res$leaves$where
  res$bars <- reshape2::melt(yval, id.vars = "where")

  return(res)
}


# format labels to plot decision tree
format_labels <- function(t) {
  t <- gsub("depth", "depth ", t)
  t <- gsub("slope", "slope ", t)
  ind <- grep("^habitat_code=", t)
  t <- gsub("^region=", "", t)
  t <- gsub("^habitat_code=", "", t)
  t <- gsub("=", "= ", t)

  # format habitat
  for (i in ind) {
    temp <- strsplit(t[i], ",")[[1]]
    temp <- paste0('*', temp, '*')
    if (length(temp) >= 4) temp[3] <- paste0("<br>", temp[3])
    t[i] <- paste(temp, collapse = ",")
  }
  return(t)
}
