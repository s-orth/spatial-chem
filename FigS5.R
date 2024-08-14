# TITLE: Create Figures for Paper
# AUTHOR: Sarah Orth, Ostling Lab
# DATE: August 13, 2024

library(matrixStats) # load necessary packages
library(Rfast)
library(ggplot2)
library(ggpubr)

setwd("/Users/Sarah/Library/CloudStorage/Box-Box/spatial-chem/") # set working directory

grid.dir <- "spatial-nulls-2010/genus-miconia/10m-nullPAs/"
list.all <- list.files(path = grid.dir) # list all files names in grid.dir folder
print(list.all) # check file names

# Read in presence-absence matrices. Cell id columns will have turned into "X1", "X2"... and so on but this is OK. R doesn't like variable names to be a number. We use the slower read.table() instead of fread() because fread can't understand a row name.
null.mat <- lapply(X = paste(grid.dir, list.all, sep = ""), 
                  FUN = read.table, sep = "\t", header = T) # tab-delimited file with a header

# census data
bci.census <- read.csv(file = "clean-BCI-by-quadrat-2010/clean_BCI_census_10m.txt",
                       sep = "\t", header = T) # tab-delimited file with a header

# repeat it all again for the other pair... M. affinis and M. argentea
p.1 <- bci.census[grep("micoar|micoaf", bci.census$Mnemonic), ] # things with micoar or micoaf
p.1 <- p.1[, -c(3:10)] # remove stem information
p.1 <- unique(p.1) # make it so we just have one cell obs per species per cell
p.1 <- p.1[, -c(1:2)] # isolate just cell information
p.1 <- p.1[duplicated(p.1), ] # removing unique observations so we know which cells BOTH WERE IN

all.miconia <- bci.census[grep("Miconia", bci.census$Latin), ]
all.miconia <- all.miconia[, -c(3:10)] # remove unique information
all.miconia <- unique(all.miconia) # keep just one cell per species
all.miconia <- all.miconia[, -c(1:2)] # isolate just cell information
all.miconia <- all.miconia[duplicated(all.miconia), ] # we only want cells that had at least 2 miconia
all.miconia <- unique(all.miconia) # now, remove dupicates (would occur when 3+ miconia present)

# number of obs in all.miconia should match number of cells recorded in results table.

# remove rows of all.miconia that are in p.1
inds <- which(all.miconia$quadrat_id %in% p.1$quadrat_id)
all.miconia <- all.miconia[-inds, ]

# assign level for shading
p.1$shade <- 1
all.miconia$shade <- 0.5

all.shading.dat <- rbind(p.1, all.miconia)

p1 <- ggplot(all.shading.dat, aes(x = quadrat_locX,
                                  y = quadrat_locY,
                                  fill = shade)) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "lightgray", high = "black") +
  theme_classic() +
  theme(aspect.ratio = 1/2) +
  scale_x_continuous(name = "X", expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(name = "Y", expand = c(0, 0), limits = c(0, NA)) +
  ggtitle("a) M. affinis and M. argentea") +
  labs(subtitle = "less co-occupancy than expected, 47/169")

# NULL FOR MICOAF AND MICOAR
# getting quadrat locations for null PA matrix
quad.locs <- bci.census[, c(17:19)]
quad.locs <- unique(quad.locs) # just quadrat IDs and location! I think can just cbind this to the PA matrix

null.1 <- null.mat[[1]]
p.1 <- null.1[, c(1:2)] # micoaf and micoar

null.1$sums <- rowSums(null.1)
p.1$sums <- rowSums(p.1)

inds <- which(p.1$sums == 2)
p.1 <- p.1[inds, ]
p.1 <- quad.locs[as.numeric(rownames(p.1)), ]
p.1$shade <- 1

inds <- which(null.1$sums >= 2)
all.miconia <- null.1[inds, ]
all.miconia <- quad.locs[as.numeric(rownames(all.miconia)), ]
all.miconia$shade <- 0.5

inds <- which(all.miconia$quadrat_id %in% p.1$quadrat_id) # remove duplicates
all.miconia <- all.miconia[-inds, ]

all.shading.dat <- rbind(p.1, all.miconia)

p2 <- ggplot(all.shading.dat, aes(x = quadrat_locX,
                                  y = quadrat_locY,
                                  fill = shade)) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "lightgray", high = "black") +
  theme_classic() +
  theme(aspect.ratio = 1/2) +
  scale_x_continuous(name = "X", expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(name = "Y", expand = c(0, 0), limits = c(0, NA)) +
  ggtitle("b) Null M. affinis and M. argentea") +
  labs(subtitle = "70/169 quadrats")

# NULL FOR MICOAF AND MICOAR
null.1 <- null.mat[[2]]
p.1 <- null.1[, c(1:2)] # micoaf and micoar

null.1$sums <- rowSums(null.1)
p.1$sums <- rowSums(p.1)

inds <- which(p.1$sums == 2)
p.1 <- p.1[inds, ]
p.1 <- quad.locs[as.numeric(rownames(p.1)), ]
p.1$shade <- 1

inds <- which(null.1$sums >= 2)
all.miconia <- null.1[inds, ]
all.miconia <- quad.locs[as.numeric(rownames(all.miconia)), ]
all.miconia$shade <- 0.5

inds <- which(all.miconia$quadrat_id %in% p.1$quadrat_id) # remove duplicates
all.miconia <- all.miconia[-inds, ]

all.shading.dat <- rbind(p.1, all.miconia)

p3 <- ggplot(all.shading.dat, aes(x = quadrat_locX,
                                  y = quadrat_locY,
                                  fill = shade)) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "lightgray", high = "black") +
  theme_classic() +
  theme(aspect.ratio = 1/2) +
  scale_x_continuous(name = "X", expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(name = "Y", expand = c(0, 0), limits = c(0, NA)) +
  ggtitle("c) Null M. affinis and M. argentea") +
  labs(subtitle = "64/169 quadrats")

null.1 <- null.mat[[3]]
p.1 <- null.1[, c(1:2)] # micoaf and micoar

null.1$sums <- rowSums(null.1)
p.1$sums <- rowSums(p.1)

inds <- which(p.1$sums == 2)
p.1 <- p.1[inds, ]
p.1 <- quad.locs[as.numeric(rownames(p.1)), ]
p.1$shade <- 1

inds <- which(null.1$sums >= 2)
all.miconia <- null.1[inds, ]
all.miconia <- quad.locs[as.numeric(rownames(all.miconia)), ]
all.miconia$shade <- 0.5

inds <- which(all.miconia$quadrat_id %in% p.1$quadrat_id) # remove duplicates
all.miconia <- all.miconia[-inds, ]

all.shading.dat <- rbind(p.1, all.miconia)

p4 <- ggplot(all.shading.dat, aes(x = quadrat_locX,
                                  y = quadrat_locY,
                                  fill = shade)) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "lightgray", high = "black") +
  theme_classic() +
  theme(aspect.ratio = 1/2) +
  scale_x_continuous(name = "X", expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(name = "Y", expand = c(0, 0), limits = c(0, NA)) +
  ggtitle("d) Null M. affinis and M. argentea") +
  labs(subtitle = "78/169 quadrats")

# NEXT GROUP
p.1 <- bci.census[grep("micone|micoaf", bci.census$Mnemonic), ] # things with micone or micoaf
p.1 <- p.1[, -c(3:10)] # remove stem information
p.1 <- unique(p.1) # make it so we just have one cell obs per species per cell
p.1 <- p.1[, -c(1:2)] # isolate just cell information
p.1 <- p.1[duplicated(p.1), ] # removing unique observations so we know which cells BOTH WERE IN

all.miconia <- bci.census[grep("Miconia", bci.census$Latin), ]
all.miconia <- all.miconia[, -c(3:10)] # remove unique information
all.miconia <- unique(all.miconia) # keep just one cell per species
all.miconia <- all.miconia[, -c(1:2)] # isolate just cell information
all.miconia <- all.miconia[duplicated(all.miconia), ] # we only want cells that had at least 2 miconia
all.miconia <- unique(all.miconia) # now, remove dupicates (would occur when 3+ miconia present)

# number of obs in all.miconia should match number of cells recorded in results table.

# remove rows of all.miconia that are in p.1
inds <- which(all.miconia$quadrat_id %in% p.1$quadrat_id)
all.miconia <- all.miconia[-inds, ]

# assign level for shading
p.1$shade <- 1
all.miconia$shade <- 0.5

all.shading.dat <- rbind(p.1, all.miconia)

p5 <- ggplot(all.shading.dat, aes(x = quadrat_locX,
                                  y = quadrat_locY,
                                  fill = shade)) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "lightgray", high = "black") +
  theme_classic() +
  theme(aspect.ratio = 1/2) +
  scale_x_continuous(name = "X", expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(name = "Y", expand = c(0, 0), limits = c(0, NA)) +
  ggtitle("e) M. affinis and M. nervosa") +
  labs(subtitle = "more co-occupancy than expected, 69/169")

# NULL PLOTS
null.1 <- null.mat[[1]]
p.1 <- null.1[, c(1, 7)] # micoaf and micone

null.1$sums <- rowSums(null.1)
p.1$sums <- rowSums(p.1)

inds <- which(p.1$sums == 2)
p.1 <- p.1[inds, ]
p.1 <- quad.locs[as.numeric(rownames(p.1)), ]
p.1$shade <- 1

inds <- which(null.1$sums >= 2)
all.miconia <- null.1[inds, ]
all.miconia <- quad.locs[as.numeric(rownames(all.miconia)), ]
all.miconia$shade <- 0.5

inds <- which(all.miconia$quadrat_id %in% p.1$quadrat_id) # remove duplicates
all.miconia <- all.miconia[-inds, ]

all.shading.dat <- rbind(p.1, all.miconia)

p6 <- ggplot(all.shading.dat, aes(x = quadrat_locX,
                                  y = quadrat_locY,
                                  fill = shade)) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "lightgray", high = "black") +
  theme_classic() +
  theme(aspect.ratio = 1/2) +
  scale_x_continuous(name = "X", expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(name = "Y", expand = c(0, 0), limits = c(0, NA)) +
  ggtitle("f) Null M. affinis and M. nervosa") +
  labs(subtitle = "40/169 quadrats")

# ANOTHER NULL FOR MICOAF AND MICONE
null.1 <- null.mat[[2]]
p.1 <- null.1[, c(1, 7)] # micoaf and micone

null.1$sums <- rowSums(null.1)
p.1$sums <- rowSums(p.1)

inds <- which(p.1$sums == 2)
p.1 <- p.1[inds, ]
p.1 <- quad.locs[as.numeric(rownames(p.1)), ]
p.1$shade <- 1

inds <- which(null.1$sums >= 2)
all.miconia <- null.1[inds, ]
all.miconia <- quad.locs[as.numeric(rownames(all.miconia)), ]
all.miconia$shade <- 0.5

inds <- which(all.miconia$quadrat_id %in% p.1$quadrat_id) # remove duplicates
all.miconia <- all.miconia[-inds, ]

all.shading.dat <- rbind(p.1, all.miconia)

p7 <- ggplot(all.shading.dat, aes(x = quadrat_locX,
                                  y = quadrat_locY,
                                  fill = shade)) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "lightgray", high = "black") +
  theme_classic() +
  theme(aspect.ratio = 1/2) +
  scale_x_continuous(name = "X", expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(name = "Y", expand = c(0, 0), limits = c(0, NA)) +
  ggtitle("g) Null M. affinis and M. nervosa") +
  labs(subtitle = "50/169 quadrats")

# ANOTHER NULL FOR MICOAF AND MICONE
null.1 <- null.mat[[3]]
p.1 <- null.1[, c(1, 7)] # micoaf and micone

null.1$sums <- rowSums(null.1)
p.1$sums <- rowSums(p.1)

inds <- which(p.1$sums == 2)
p.1 <- p.1[inds, ]
p.1 <- quad.locs[as.numeric(rownames(p.1)), ]
p.1$shade <- 1

inds <- which(null.1$sums >= 2)
all.miconia <- null.1[inds, ]
all.miconia <- quad.locs[as.numeric(rownames(all.miconia)), ]
all.miconia$shade <- 0.5

inds <- which(all.miconia$quadrat_id %in% p.1$quadrat_id) # remove duplicates
all.miconia <- all.miconia[-inds, ]

all.shading.dat <- rbind(p.1, all.miconia)

p8 <- ggplot(all.shading.dat, aes(x = quadrat_locX,
                                  y = quadrat_locY,
                                  fill = shade)) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradient(low = "lightgray", high = "black") +
  theme_classic() +
  theme(aspect.ratio = 1/2) +
  scale_x_continuous(name = "X", expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(name = "Y", expand = c(0, 0), limits = c(0, NA)) +
  ggtitle("h) Null M. affinis and M. nervosa") +
  labs(subtitle = "42/169 quadrats")

ggarrange(p1, p5,
          p2, p6, 
          p3, p7,
          p4, p8,
          ncol = 2, nrow = 4)
