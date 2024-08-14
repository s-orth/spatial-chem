# TITLE: Create Figures for Paper
# AUTHOR: Sarah Orth, Ostling Lab
# DATE: August 13, 2024

library(matrixStats) # load necessary packages
library(Rfast)
library(ggplot2)
library(ggpubr)

setwd("/Users/Sarah/Library/CloudStorage/Box-Box/spatial-chem/") # set working directory

grid.dir <- "spatial-nulls-2010/genus-psychotria/10m-nullPAs/"
list.all <- list.files(path = grid.dir) # list all files names in grid.dir folder
print(list.all) # check file names

# Read in presence-absence matrices. Cell id columns will have turned into "X1", "X2"... and so on but this is OK. R doesn't like variable names to be a number. We use the slower read.table() instead of fread() because fread can't understand a row name.
null.mat <- lapply(X = paste(grid.dir, list.all, sep = ""), 
                   FUN = read.table, sep = "\t", header = T) # tab-delimited file with a header

# census data
bci.census <- read.csv(file = "clean-BCI-by-quadrat-2010/clean_BCI_census_10m.txt",
                       sep = "\t", header = T) # tab-delimited file with a header

p.1 <- bci.census[grep("psycli|psycho", bci.census$Mnemonic), ] # things with micone or micoaf
p.1 <- p.1[, -c(3:10)] # remove stem information
p.1 <- unique(p.1) # make it so we just have one cell obs per species per cell
p.1 <- p.1[, -c(1:2)] # isolate just cell information
p.1 <- p.1[duplicated(p.1), ] # removing unique observations so we know which cells BOTH WERE IN

all.psyc <- bci.census[grep("Psychotria", bci.census$Latin), ]
all.psyc <- all.psyc[, -c(3:10)] # remove unique information
all.psyc <- unique(all.psyc) # keep just one cell per species
all.psyc <- all.psyc[, -c(1:2)] # isolate just cell information
all.psyc <- all.psyc[duplicated(all.psyc), ] # we only want cells that had at least 2 psyc
all.psyc <- unique(all.psyc) # now, remove dupicates (would occur when 3+ psyc present)

# number of obs in all.psyc should match number of cells recorded in results table.

# remove rows of all.psyc that are in p.1
inds <- which(all.psyc$quadrat_id %in% p.1$quadrat_id)
all.psyc <- all.psyc[-inds, ]

# assign level for shading
p.1$shade <- 1
all.psyc$shade <- 0.5

all.shading.dat <- rbind(p.1, all.psyc)

p1 <- ggplot(all.shading.dat, aes(x = quadrat_locX,
                                  y = quadrat_locY,
                                  fill = shade)) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "lightgray", high = "black") +
  theme_classic() +
  theme(aspect.ratio = 1/2) +
  scale_x_continuous(name = "X", expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(name = "Y", expand = c(0, 0), limits = c(0, NA)) +
  ggtitle("a) P. horizontalis & P. limonensis") +
  labs(subtitle = "less co-occupancy than expected, 4/270")

# getting quadrat locations for null PA matrix
quad.locs <- bci.census[, c(17:19)]
quad.locs <- unique(quad.locs) # just quadrat IDs and location! I think can just cbind this to the PA matrix

null.1 <- null.mat[[1]]
p.1 <- null.1[, c(8, 9)] # psycli and psycho

null.1$sums <- rowSums(null.1)
p.1$sums <- rowSums(p.1)

inds <- which(p.1$sums == 2)
p.1 <- p.1[inds, ]
p.1 <- quad.locs[as.numeric(rownames(p.1)), ]
p.1$shade <- 1

inds <- which(null.1$sums >= 2)
all.psyc <- null.1[inds, ]
all.psyc <- quad.locs[as.numeric(rownames(all.psyc)), ]
all.psyc$shade <- 0.5

inds <- which(all.psyc$quadrat_id %in% p.1$quadrat_id) # remove duplicates
all.psyc <- all.psyc[-inds, ]

all.shading.dat <- rbind(p.1, all.psyc)

p2 <- ggplot(all.shading.dat, aes(x = quadrat_locX,
                                  y = quadrat_locY,
                                  fill = shade)) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "lightgray", high = "black") +
  theme_classic() +
  theme(aspect.ratio = 1/2) +
  scale_x_continuous(name = "X", expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(name = "Y", expand = c(0, 0), limits = c(0, NA)) +
  ggtitle("b) Null P. horizontalis & P. limonensis") +
  labs(subtitle = "14/270 quadrats")

# another
null.1 <- null.mat[[2]]
p.1 <- null.1[, c(8, 9)] # psycli and psycho

null.1$sums <- rowSums(null.1)
p.1$sums <- rowSums(p.1)

inds <- which(p.1$sums == 2)
p.1 <- p.1[inds, ]
p.1 <- quad.locs[as.numeric(rownames(p.1)), ]
p.1$shade <- 1

inds <- which(null.1$sums >= 2)
all.psyc <- null.1[inds, ]
all.psyc <- quad.locs[as.numeric(rownames(all.psyc)), ]
all.psyc$shade <- 0.5

inds <- which(all.psyc$quadrat_id %in% p.1$quadrat_id) # remove duplicates
all.psyc <- all.psyc[-inds, ]

all.shading.dat <- rbind(p.1, all.psyc)

p3 <- ggplot(all.shading.dat, aes(x = quadrat_locX,
                                  y = quadrat_locY,
                                  fill = shade)) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "lightgray", high = "black") +
  theme_classic() +
  theme(aspect.ratio = 1/2) +
  scale_x_continuous(name = "X", expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(name = "Y", expand = c(0, 0), limits = c(0, NA)) +
  ggtitle("c) Null P. horizontalis & P. limonensis") +
  labs(subtitle = "10/270 quadrats")

# another
null.1 <- null.mat[[3]]
p.1 <- null.1[, c(8, 9)] # psycli and psycho

null.1$sums <- rowSums(null.1)
p.1$sums <- rowSums(p.1)

inds <- which(p.1$sums == 2)
p.1 <- p.1[inds, ]
p.1 <- quad.locs[as.numeric(rownames(p.1)), ]
p.1$shade <- 1

inds <- which(null.1$sums >= 2)
all.psyc <- null.1[inds, ]
all.psyc <- quad.locs[as.numeric(rownames(all.psyc)), ]
all.psyc$shade <- 0.5

inds <- which(all.psyc$quadrat_id %in% p.1$quadrat_id) # remove duplicates
all.psyc <- all.psyc[-inds, ]

all.shading.dat <- rbind(p.1, all.psyc)

p4 <- ggplot(all.shading.dat, aes(x = quadrat_locX,
                                  y = quadrat_locY,
                                  fill = shade)) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "lightgray", high = "black") +
  theme_classic() +
  theme(aspect.ratio = 1/2) +
  scale_x_continuous(name = "X", expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(name = "Y", expand = c(0, 0), limits = c(0, NA)) +
  ggtitle("d) Null P. horizontalis & P. limonensis") +
  labs(subtitle = "7/270 quadrats")

# repeat it all again for the other pair... P. marginata and P. acuminata
p.1 <- bci.census[grep("psycma|psycac", bci.census$Mnemonic), ] # things with psycma or psycac
p.1 <- p.1[, -c(3:10)] # remove stem information
p.1 <- unique(p.1) # make it so we just have one cell obs per species per cell
p.1 <- p.1[, -c(1:2)] # isolate just cell information
p.1 <- p.1[duplicated(p.1), ] # removing unique observations so we know which cells BOTH WERE IN

all.psyc <- bci.census[grep("Psychotria", bci.census$Latin), ]
all.psyc <- all.psyc[, -c(3:10)] # remove unique information
all.psyc <- unique(all.psyc) # keep just one cell per species
all.psyc <- all.psyc[, -c(1:2)] # isolate just cell information
all.psyc <- all.psyc[duplicated(all.psyc), ] # we only want cells that had at least 2 psyc
all.psyc <- unique(all.psyc) # now, remove dupicates (would occur when 3+ psyc present)

# number of obs in all.psyc should match number of cells recorded in results table.

# remove rows of all.psyc that are in p.1
inds <- which(all.psyc$quadrat_id %in% p.1$quadrat_id)
all.psyc <- all.psyc[-inds, ]

# assign level for shading
p.1$shade <- 1
all.psyc$shade <- 0.5

all.shading.dat <- rbind(p.1, all.psyc)

p5 <- ggplot(all.shading.dat, aes(x = quadrat_locX,
                                  y = quadrat_locY,
                                  fill = shade)) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "lightgray", high = "black") +
  theme_classic() +
  theme(aspect.ratio = 1/2) +
  scale_x_continuous(name = "X", expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(name = "Y", expand = c(0, 0), limits = c(0, NA)) +
  ggtitle("e) P. acuminata & P. marginata") +
  labs(subtitle = "more co-occupancy than expected, 6/270")

# last pair, psycma and psycac
null.1 <- null.mat[[1]]
p.1 <- null.1[, c(1, 10)] # psycac and psycma

null.1$sums <- rowSums(null.1)
p.1$sums <- rowSums(p.1)

inds <- which(p.1$sums == 2)
p.1 <- p.1[inds, ]
p.1 <- quad.locs[as.numeric(rownames(p.1)), ]
p.1$shade <- 1

inds <- which(null.1$sums >= 2)
all.psyc <- null.1[inds, ]
all.psyc <- quad.locs[as.numeric(rownames(all.psyc)), ]
all.psyc$shade <- 0.5

inds <- which(all.psyc$quadrat_id %in% p.1$quadrat_id) # remove duplicates
all.psyc <- all.psyc[-inds, ]

all.shading.dat <- rbind(p.1, all.psyc)

p6 <- ggplot(all.shading.dat, aes(x = quadrat_locX,
                                  y = quadrat_locY,
                                  fill = shade)) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "lightgray", high = "black") +
  theme_classic() +
  theme(aspect.ratio = 1/2) +
  scale_x_continuous(name = "X", expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(name = "Y", expand = c(0, 0), limits = c(0, NA)) +
  ggtitle("f) Null P. acuminata & P. marginata") +
  labs(subtitle = "4/270 quadrats")

# another
null.1 <- null.mat[[2]]
p.1 <- null.1[, c(1, 10)] # psycac and psycma

null.1$sums <- rowSums(null.1)
p.1$sums <- rowSums(p.1)

inds <- which(p.1$sums == 2) 
p.1 <- p.1[inds, ]
p.1 <- quad.locs[as.numeric(rownames(p.1)), ]
p.1$shade <- 1

inds <- which(null.1$sums >= 2)
all.psyc <- null.1[inds, ]
all.psyc <- quad.locs[as.numeric(rownames(all.psyc)), ]
all.psyc$shade <- 0.5

inds <- which(all.psyc$quadrat_id %in% p.1$quadrat_id) # remove duplicates
all.psyc <- all.psyc[-inds, ]

all.shading.dat <- rbind(p.1, all.psyc)

p7 <- ggplot(all.shading.dat, aes(x = quadrat_locX,
                                   y = quadrat_locY,
                                   fill = shade)) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "lightgray", high = "lightgray") +
  theme_classic() +
  theme(aspect.ratio = 1/2) +
  scale_x_continuous(name = "X", expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(name = "Y", expand = c(0, 0), limits = c(0, NA)) +
  ggtitle("g) Null P. acuminata & P. marginata") +
  labs(subtitle = "2/270 quadrats")

# another one
null.1 <- null.mat[[3]]
p.1 <- null.1[, c(1, 10)] # psycac and psycma

null.1$sums <- rowSums(null.1)
p.1$sums <- rowSums(p.1)

inds <- which(p.1$sums == 2) # no co-occurrences
p.1 <- p.1[inds, ]
p.1 <- quad.locs[as.numeric(rownames(p.1)), ]
p.1$shade <- 1

inds <- which(null.1$sums >= 2)
all.psyc <- null.1[inds, ]
all.psyc <- quad.locs[as.numeric(rownames(all.psyc)), ]
all.psyc$shade <- 0.5

inds <- which(all.psyc$quadrat_id %in% p.1$quadrat_id) # remove duplicates
all.psyc <- all.psyc[-inds, ]

all.shading.dat <- rbind(p.1, all.psyc)

p8 <- ggplot(all.shading.dat, aes(x = quadrat_locX,
                                   y = quadrat_locY,
                                   fill = shade)) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "lightgray", high = "lightgray") +
  theme_classic() +
  theme(aspect.ratio = 1/2) +
  scale_x_continuous(name = "X", expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(name = "Y", expand = c(0, 0), limits = c(0, NA)) +
  ggtitle("h) Null P. acuminata & P. marginata") +
  labs(subtitle = "1/270 quadrats")

ggarrange(p1, p5,
          p2, p6, 
          p3, p7,
          p4, p8,
          ncol = 2, nrow = 4)

