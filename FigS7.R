# TITLE: Create Figures for Paper
# AUTHOR: Sarah Orth, Ostling Lab
# LAST MODIFIED: August 12, 2024

library(matrixStats) # load necessary packages
library(Rfast)
library(ggplot2)
library(ggpubr)

setwd("/Users/Sarah/Library/CloudStorage/Box-Box/spatial-chem/") # set working directory

grid.dir <- "BCI-PA-matrices-2010/" # set folder containing presence-absence matrices
list.all <- list.files(path = grid.dir) # list all files names in grid.dir folder
print(list.all) # check order of file names
list.all <- list.all[c(2:4, 1)] # put files in order from smallest to largest spatial scale
print(list.all) # confirm correct order of files

# Read in presence-absence matrices. Cell id columns will have turned into "X1", "X2"... and so on but this is OK. R doesn't like variable names to be a number. We use the slower read.table() instead of fread() because fread can't understand a row name.
bci.mat <- lapply(X = paste(grid.dir, list.all, sep = ""), 
                  FUN = read.table, sep = "\t", header = T) # tab-delimited file with a header

bci.cd <- read.table(
  file = "clean-BCI-by-quadrat-2010/disparity_matrix_defense_only.txt",
  sep = "\t", # tab-delimited
  row.names = 1, # row names are in column 1
  header = T) # header present

diag(bci.cd) <- NA # ignore conspecific relationships, as this may dilute metrics.

# read in null co-occupancy data
my.dir <- paste0("spatial-nulls-2010/genus-inga/") # subset path
list.all <- list.files(path = my.dir) # list files
print(list.all)
null.stats <- list.all[c(15:17, 14)] # reorder
null.zscores <- list.all[c(23:25, 22)]

null.stats <- lapply(X = paste(my.dir, null.stats, sep = ""), # read in data
                     FUN = read.table, 
                     sep = "\t", header = T) # tab-delim file w/ header

null.zscores <- lapply(X = paste(my.dir, null.zscores, sep = ""), # read in data
                       FUN = read.table, 
                       sep = "\t", header = T) # tab-delim file w/ header

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

Sub.Mat <- lapply(X = bci.mat, FUN = Get.Sub, g = "inga")

Rinds <- which(rownames(bci.cd) %in% colnames(Sub.Mat[[4]])) 
Cinds <- which(colnames(bci.cd) %in% colnames(Sub.Mat[[4]]))
sub.cd <- data.matrix(bci.cd[Rinds, Cinds])

cell.sizes <- c("10m", "20m", "50m", "100m")
cell_length <- c(10, 20, 50, 100)

# calculate co-occupancy rates
pair.occ <- matrix(data = NA, nrow = nrow(sub.cd), ncol = ncol(sub.cd), 
                   dimnames = list(rownames(sub.cd), colnames(sub.cd)))
all.pair.occ <- list()

for(i in 1:4){
  d <- Sub.Mat[[i]]
  for(xi in 1:nrow(sub.cd)){
    for(yi in 1:ncol(sub.cd)){
      pair.cd <- sub.cd[xi, yi]
      sp.x <- rownames(sub.cd)[xi]
      sp.y <- colnames(sub.cd)[yi]
      tot.cell <- matrixStats::count(rowSums(d) > 0)
      sp.xy <- d[, c(sp.x, sp.y)]
      sp.xy <- matrixStats::count(rowSums(sp.xy) == 2)
      pair.occ[xi, yi] <- sp.xy/tot.cell
    }
  }
  
  pair.occ[upper.tri(x = pair.occ, diag = TRUE)] <- NA
  all.pair.occ[[i]] <- pair.occ
  
}

occ.rate <- matrix(data = NA, nrow = nrow(sub.cd), ncol = ncol(sub.cd), 
                   dimnames = list(rownames(sub.cd), colnames(sub.cd)))

occ.rate.prods <- list() # PRODUCTS of occupancy rates

for(i in 1:4){
  d <- Sub.Mat[[i]]
  for(xi in 1:nrow(sub.cd)){
    for(yi in 1:ncol(sub.cd)){
      sp.x <- rownames(sub.cd)[xi]
      sp.y <- colnames(sub.cd)[yi]
      tot.cell <- matrixStats::count(rowSums(d) > 0)
      sp.x <- d[ , sp.x]
      sp.x <- sum(sp.x) / tot.cell
      sp.y <- d[ , sp.y]
      sp.y <- sum(sp.y) / tot.cell
      occ.rate[xi, yi] <- sp.x * sp.y
    }
  }
  
  occ.rate[upper.tri(x = occ.rate, diag = TRUE)] <- NA
  occ.rate.prods[[i]] <- occ.rate
  
}

cd <- as.data.frame(as.table(sub.cd))
cd <- na.omit(cd)
names(cd)[1] <- "sp_i"
names(cd)[2] <- "sp_j"
names(cd)[3] <- "disp_ij"

all.dat <- list()

for(i in 1:4){
  paired <- as.data.frame(as.table(all.pair.occ[[i]]))
  paired <- na.omit(paired)
  names(paired)[1] <- "sp_i"
  names(paired)[2] <- "sp_j"
  names(paired)[3] <- "occ_ij"
  
  ij <- as.data.frame(as.table(occ.rate.prods[[i]]))
  ij <- na.omit(ij)
  names(ij)[1] <- "sp_i"
  names(ij)[2] <- "sp_j"
  names(ij)[3] <- "i_times_j"
  
  dat <- merge(paired, ij, by = c("sp_i", "sp_j"))
  dat <- merge(dat, cd, by = c("sp_i", "sp_j"))
  
  all.dat[[i]] <- dat
  
}

# organize it all for plotting
plot.dat <- list()

for(i in 1:4){
  obs <- all.dat[[i]]
  stats <- null.stats[[i]]
  dists <- null.zscores[[i]]
  
  # give null data col names to help match up with the observed data
  names(stats)[names(stats) == "Var1"] <- "sp_i"
  names(stats)[names(stats) == "Var2"] <- "sp_j"
  names(dists)[names(dists) == "Var1"] <- "sp_i"
  names(dists)[names(dists) == "Var2"] <- "sp_j"
  
  # ensure data all sorted in same order
  obs <- obs[with(obs, order(sp_i, sp_j)), ]
  stats <- stats[with(stats, order(sp_i, sp_j)), ]
  dists <- dists[with(dists, order(sp_i, sp_j)), ]
  
  nums <- as.matrix(dists[, -c(1:2)])
  
  zscore <- (obs$occ_ij - stats$null_occRate_mean) / stats$null_occRate_sd
  null.z.low <- rownth(nums, elems = as.integer(rep(25, dim(nums)[1])))
  null.z.high <- rownth(nums, elems = as.integer(rep(25, dim(nums)[1])), descending = TRUE)
  
  obs$zscore <- zscore
  obs$null.z.low <- null.z.low
  obs$null.z.high <- null.z.high
  obs$null.occ.ij <- stats$null_occRate_mean
  
  plot.dat[[i]] <- obs
}

dat <- plot.dat[[1]] # select one spatial scale

p1 <- ggplot(dat, aes(x = disp_ij, y = zscore, col = null.occ.ij)) +
  geom_point(size = 2) +
  scale_color_gradient(low = "gray87", high = "black") +
  theme_classic() +
  labs(title = "Co-occupancy rates of Inga, 10m scale",
       color = "mean null co-occupancy") +
  xlab("pairwise chemical disparity") +
  ylab("co-occupancy rate z-score") +
  geom_line(aes(x = disp_ij, y = null.z.low), color = "gray") +
  geom_line(aes(x = disp_ij, y = null.z.high), color = "gray")

p1

