
library(geomorph)

# indicates the working directory

# carapace ----

# imputing lm data
lm_carapace <- readland.tps("carapace_all.TPS", specID = "imageID")

# applying gpa
gpa_carapace <- gpagen(lm_carapace)

# looking do landmark coord variation
plotAllSpecimens(gpa_carapace$coords)

# reading classifier object
classifier_carapace <- read.csv("carapace_classifier.csv", sep=";")

# applying procD.lm (equivalent to MANCOVA)
mancova_carapace <- procD.lm(
  gpa_carapace$coords ~ log(gpa_carapace$Csize) * island, 
  data = classifier_carapace, 
  random = ~id, 
  iter = 999
)

# summary results of mancova
summary(mancova_carapace)

# posthoc pairwise comparison
pairwise_carapace <- pairwise(mancova_carapace,
  groups = interaction(classifier_carapace$island))

# summary results of pairwise comparison
summary(pairwise_carapace)

# abdomen ----

# imputing lm data
lm_abdomen <- readland.tps("abdomen_all.TPS", specID = "imageID")

# applying gpa
gpa_abdomen <- gpagen(lm_abdomen)

# looking do landmark coord variation
plotAllSpecimens(gpa_abdomen$coords)

# reading classifier object
classifier_abdomen <- read.csv("abdomen_classifier.csv", sep=";")

# applying procD.lm (equivalent to MANCOVA)
mancova_abdomen <- procD.lm(
  gpa_abdomen$coords ~ log(gpa_abdomen$Csize) * island, 
  data = classifier_abdomen, 
  random = ~id, 
  iter = 999
)

# summary results of mancova
summary(mancova_abdomen)

# posthoc pairwise comparison
pairwise_abdomen <- pairwise(mancova_abdomen,
                              groups = interaction(classifier_abdomen$island))

# summary results of pairwise comparison
summary(pairwise_abdomen)

# chela ----

# imputing lm data
lm_chela <- readland.tps("chela_all.TPS", specID = "imageID")

# applying gpa
gpa_chela <- gpagen(lm_chela)

# looking do landmark coord variation
plotAllSpecimens(gpa_chela$coords)

# reading classifier object
classifier_chela <- read.csv("chela_classifier.csv", sep=";")

# applying procD.lm (equivalent to MANCOVA)
mancova_chela <- procD.lm(
  gpa_chela$coords ~ log(gpa_chela$Csize) * chela_id * morph * island, 
  data = classifier_chela, 
  random = ~id, 
  iter = 999
)

# summary results of mancova
summary(mancova_chela)

# posthoc pairwise comparison
pairwise_chela <- pairwise(mancova_chela,
  groups = interaction(classifier_chela$island,
                       classifier_chela$morph))

# summary results of pairwise comparison
summary(pairwise_chela)

# ploting mean shape by island ----

# indicating island colors
island_colors <- c("RA" = '#3677b2', "FN" = '#3f642f', "TR" = '#b96848')

# indicating link between landmarks to draw connection between dots

# carapace 
link_carapace <- matrix(c(1, 2,
                  2, 3,
                  3, 4,
                  4, 5,
                  5, 6,
                  6, 7,
                  7, 8,
                  8, 9,
                  9, 10,
                  10, 11,
                  11, 12,
                  12, 13,
                  13, 1,
                  14, 15,
                  15, 16,
                  16, 17, 
                  17, 18,
                  18, 19,
                  19, 14
), ncol = 2, byrow = TRUE)

# calculating average shape for each island
mean_shapes_by_area <- list()

for (current_area in unique(classifier_carapace$island)) {
  inds_area <- which(classifier_carapace$island == current_area)
  mean_shape_area <- apply(gpa_carapace$coords[,,inds_area], 
                           c(1, 2), mean)
  mean_shapes_by_area[[as.character(current_area)]] <- mean_shape_area
}

# ploting
plot(NA, 
     xlim = range(gpa_carapace$coords[, 1, ]), 
     ylim = range(gpa_carapace$coords[, 2, ]), 
     asp = 1, xlab = "", ylab = "", axes = FALSE)

# drawing shape for each island
for (current_area in names(mean_shapes_by_area)) {
  shape <- mean_shapes_by_area[[current_area]]
  for (i in 1:nrow(link_carapace)) {
    segments(
      shape[link_carapace[i, 1], 1], shape[link_carapace[i, 1], 2],
      shape[link_carapace[i, 2], 1], shape[link_carapace[i, 2], 2],
      col = adjustcolor(island_colors[current_area], alpha.f = 0.6),
      lwd = 6
    )
  }
  points(shape, pch = 16, cex = 2, 
         col = island_colors[current_area])
}

# abdomen
link_abdomen <- matrix(c(1, 2,
                  2, 3,
                  3, 4,
                  4, 5,
                  5, 6,
                  6, 7,
                  7, 8,
                  8, 9,
                  9, 10,
                  10, 11,
                  11, 12,
                  12, 1), 
                ncol = 2, byrow = TRUE)

# calculating average shape for each island
mean_shapes_by_area <- list()

for (current_area in unique(classifier_abdomen$island)) {
  inds_area <- which(classifier_abdomen$island == current_area)
  mean_shape_area <- apply(gpa_abdomen$coords[,,inds_area], 
                           c(1, 2), mean)
  mean_shapes_by_area[[as.character(current_area)]] <- mean_shape_area
}

# ploting
plot(NA, 
     xlim = range(gpa_abdomen$coords[, 1, ]), 
     ylim = range(gpa_abdomen$coords[, 2, ]), 
     asp = 1, xlab = "", ylab = "", axes = FALSE)

# drawing shape for each island
for (current_area in names(mean_shapes_by_area)) {
  shape <- mean_shapes_by_area[[current_area]]
  for (i in 1:nrow(link_abdomen)) {
    segments(
      shape[link_abdomen[i, 1], 1], shape[link_abdomen[i, 1], 2],
      shape[link_abdomen[i, 2], 1], shape[link_abdomen[i, 2], 2],
      col = adjustcolor(island_colors[current_area], alpha.f = 0.6),
      lwd = 6
    )
  }
  points(shape, pch = 16, cex = 2, 
         col = island_colors[current_area])
}

# chela
link_chela <- matrix(c(1, 2,
                  2, 3,
                  3, 4,
                  4, 5,
                  5, 6,
                  6, 7,
                  7, 8,
                  8, 9,
                  9, 10,
                  10, 11,
                  11, 12,
                  12, 1), 
                ncol = 2, byrow = TRUE)

# calculating average shape for each island
mean_shapes_by_area <- list()

for (current_area in unique(classifier_chela$island)) {
  inds_area <- which(classifier_chela$island == current_area)
  mean_shape_area <- apply(gpa_chela$coords[,,inds_area], 
                           c(1, 2), mean)
  mean_shapes_by_area[[as.character(current_area)]] <- mean_shape_area
}

# ploting
plot(NA, 
     xlim = range(gpa_chela$coords[, 1, ]), 
     ylim = range(gpa_chela$coords[, 2, ]), 
     asp = 1, xlab = "", ylab = "", axes = FALSE)

# drawing shape for each island
for (current_area in names(mean_shapes_by_area)) {
  shape <- mean_shapes_by_area[[current_area]]
  for (i in 1:nrow(link_chela)) {
    segments(
      shape[link_chela[i, 1], 1], shape[link_chela[i, 1], 2],
      shape[link_chela[i, 2], 1], shape[link_chela[i, 2], 2],
      col = adjustcolor(island_colors[current_area], alpha.f = 0.6),
      lwd = 6
    )
  }
  points(shape, pch = 16, cex = 2, 
         col = island_colors[current_area])
}
