# TITLE: Create Figures for Paper - Occupancy vs. Disparity
# AUTHOR: Sarah Orth, Ostling Lab
# LAST MODIFIED: August 12, 2024

library(matrixStats) # load necessary packages
library(tidyverse)
library(car)
library(lmtest)
library(gtools)

setwd("/Users/Sarah/Library/CloudStorage/Box-Box/spatial-chem/") # set working directory

grid.dir <- "BCI-PA-matrices-2010/" # set folder containing presence-absence matrices
list.all <- list.files(path = grid.dir) # list all files names in grid.dir folder
print(list.all) # check order of file names
list.all <- list.all[c(2:4, 1)] # put files in order from smallest to largest spatial scale
print(list.all) # confirm correct order of files

bci.mat <- lapply(X = paste(grid.dir, list.all, sep = ""), 
                  FUN = read.table, sep = "\t", header = T) # tab-delimited file with a header

bci.cd <- read.table(
  file = "clean-BCI-by-quadrat-2010/disparity_matrix_defense_only.txt",
  sep = "\t", # tab-delimited
  row.names = 1, # row names are in column 1
  header = T) # header present

diag(bci.cd) <- NA # ignore conspecific relationships, as this may dilute metrics.

Get.Sub <- function(d, g = NULL) {
  # d is a data frame containing the presence-absence matrix
  # g is a string, the name of a genus of interest, such as "psychotria". Default value is NULL
  # this function will be used in conjunction with lapply() to be iterated over a list of data frames.
  
  # filter data frame for a genus of interest, set by providing a string for g
  # if there is a string for "g"
  if(length(g) > 0) { 
    
    # indicate rows including genus of interest
    inds <- grepl(g, d$Latin, ignore.case = TRUE)
    
    # filter
    d <- d[inds, ]}
  
  # remove extraneous identification information  
  d <- d[ , -c(1:2)] 
  
  # convert d to matrix format
  d <- as.matrix(d) 
  
  # transform so that species ID is column and quadrat is row.
  d <- t(d) 
}

Sub.Mat <- lapply(X = bci.mat, FUN = Get.Sub, 
                  g = "inga") # change for different genera (e.g. set g = "miconia" or g = "psychotria")

Rinds <- which(rownames(bci.cd) %in% colnames(Sub.Mat[[4]]))
Cinds <- which(colnames(bci.cd) %in% colnames(Sub.Mat[[4]]))
sub.cd <- data.matrix(bci.cd[Rinds, Cinds])

occ.rate <- list()
cell.sizes <- c("10m", "20m", "50m", "100m")
cell_length <- c(10, 20, 50, 100)

for(i in 1:4){
  d <- Sub.Mat[[i]] # extract one spatial scale
  tot <- matrixStats::count(matrixStats::rowSums2(d, na.rm = TRUE) > 0) # count total number of occupied cells
  d <- colSums(d, na.rm = TRUE) # calculate number of occupations for each species
  d <- as.data.frame(d/tot) # calculate relative occupancy
  names(d)[1] <- cell.sizes[i] # add a variable for spatial scale (cell size)
  occ.rate[[i]] <- d # add data to occ.rate list
}

occ <- occ.rate[[1]] # create a dataframe
for(i in 2:4){
  occ <- merge(occ, occ.rate[[i]], by = 'row.names') # merge occupancy rates at all spatial scales
  row.names(occ) <- occ$Row.names
  occ <- occ[, -c(1)]
}

avg.disp <- as.data.frame(matrixStats::rowMeans2(sub.cd, na.rm = TRUE))
row.names(avg.disp) <- row.names(sub.cd)
names(avg.disp)[1] <- "avg_disp"

all.dat <- list()
for(i in 1:4){
  d <- merge(occ.rate[[i]], avg.disp, by = 'row.names')
  rownames(d) <- d$Row.names
  d <- d[, -c(1)]
  
  colnames(d) <- c("occ_rate", "avg_disp")
  d$cell_length <- cell_length[i]
  
  all.dat[[i]] <- d
}

# fit the models for plotting

# function to extract overall p value
overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# set up plotting layout
m <- rbind(c(1, 4, 5, 6, 7),
           c(2, 3, 3, 3, 3))
layout(m,
       heights = c(0.6, 0.1),
       widths = c(0.1, 0.5, 0.5, 0.5, 0.5))
layout.show(7)
par(mar = c(0, 0, 0, 0)) # set margins to zero when doing titles

# plot y axis
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "") 
text(1, 1, "occupancy rate", cex = 1.5, font = 2, srt = 90)

# plot empty square
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "") 

# plot shared x axis label
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "") 
text(1, 1, "mean defense chemical disparity", cex = 1.5, font = 2)

par(mar = c(1., 1.5, 1.5, 1), cex.main = 1, cex.axis = 1) # Make the margins non-zero for plots

# text plotting stuff

for(i in 1:4){
  
  ad <- all.dat[[i]]$avg_disp
  or <- all.dat[[i]]$occ_rate
  fit <- lm(or ~ ad)
  cf <- round(coef(fit), 2)
  
  eq.dat <- paste0("y = ", cf[1],
               ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x ")
  rs <- bquote(italic(R)^2 == .(format(summary(fit)$adj.r.squared, digits = 2)))
  
  plot(x = ad,
       y = or,
       pch = 16, cex = 0.75,
       ylim = c(0, 1),
       xlim = c(0, 0.8),
       xlab = "",
       ylab = "",
       xaxs = "i",
       yaxs = "i",
       las = 1,
       main = cell.sizes[i], cex.main = 1.5,
       cex.axis = 0.9)
  abline(fit, col = "gray") # comment out for miconia and psychotria, where there is not a good fit.
  mtext(eq.dat, side = 3, line = -2, cex = 0.75, adj = 0.1)
  mtext(rs, side = 3, line = -3, cex = 0.75, adj = 0.1)
}



