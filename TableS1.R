# TITLE: Compare Observed to Null and create table for paper.
# AUTHOR: Sarah Orth, Ostling Lab
# LAST MODIFIED: August 14, 2024

# load necessary packages
library(tidyr)
library(ggpubr)
library(Rfast)
library(dplyr)
library(gtools)

# set working directory
setwd("/Users/Sarah/Library/CloudStorage/Box-Box/spatial-chem/")

# list of folders to read in
s <- c("2005", "kursar-2005") # list of subsets

obs.dat <- list() # list to store data lists (lists within a list)

trait.stats <- list() # list to store null summary data lists
spatial.stats <- list()

trait.dists <- list() # list to store null distribution lists
spatial.dists <- list()

# for loop to read in relevant data
for(i in 1:2){
  my.dir <- paste0("quadrat-values-", s[i], "/genus-inga/") # subset path
  list.all <- list.files(path = my.dir) # list files
  list.all <- list.all[c(2:4, 1)] # reorder
  obs.dat[[i]] <- lapply(X = paste(my.dir, list.all, sep = ""), # read in data
                         FUN = read.table, 
                         sep = "\t", header = T) # tab-delim file w/ header
  
  my.dir <- paste0("trait-nulls-", s[i], "/genus-inga/") # subset path for null
  list.all <- list.files(path = my.dir) # list files
  list.trait.stats <- list.all[c(26:28, 25)] # select null stats
  list.trait.dists <- list.all[c(6:8, 5)] # select null distributions
  trait.stats[[i]] <- lapply(X = paste(my.dir, list.trait.stats, sep = ""), # read in data
                            FUN = read.table, 
                            sep = "\t", header = T)
  trait.dists[[i]] <- lapply(X = paste(my.dir, list.trait.dists, sep = ""), # read in data
                            FUN = read.table, 
                            sep = "\t", header = T)
  
  my.dir <- paste0("spatial-nulls-", s[i], "/genus-inga/") # subset path for null
  list.all <- list.files(path = my.dir) # list files
  list.spatial.stats <- list.all[c(27:29, 26)] # select null stats
  list.spatial.dists <- list.all[c(7:9, 6)] # select null distributions
  spatial.stats[[i]] <- lapply(X = paste(my.dir, list.spatial.stats, sep = ""), # read in data
                            FUN = read.table, 
                            sep = "\t", header = T)
  spatial.dists[[i]] <- lapply(X = paste(my.dir, list.spatial.dists, sep = ""), # read in data
                            FUN = read.table, 
                            sep = "\t", header = T)
}

# add names to help keep organized
names(obs.dat) <- s
names(trait.stats) <- s
names(trait.dists) <- s
names(spatial.stats) <- s
names(spatial.dists) <- s

scl <- c("10m", "20m", "50m", "100m") # scales for forming table

# big place to save everything
all.table.dat <- list()

for(i in 1:2){ # going through each subset...
  
  # pull out one subset for the observed and each null. Now we have plain old lists, not nested lists.
  obs <- obs.dat[[i]]
  tst <- trait.stats[[i]]
  tdi <- trait.dists[[i]]
  sst <- spatial.stats[[i]]
  sdi <- spatial.dists[[i]]
  
  table.dat <- list() # to store organized data
  
  for(j in 1:4){ # now going through each spatial scale at the subset...
    
    # select one spatial scale for each thing. now we have plain old vectors
    obs.x <- obs[[j]]
    tst.x <- tst[[j]]
    tdi.x <- tdi[[j]]
    sst.x <- sst[[j]]
    sdi.x <- sdi[[j]]
    
    # compute observed NDIs and NNIs for each null
    t.obs.NDI <- (obs.x$mean_disp - tst.x$null_meanDisp_mean) / tst.x$null_meanDisp_sd
    t.obs.NDI <- na.omit(t.obs.NDI)
    
    t.obs.NNI <- (obs.x$mean_NN - tst.x$null_meanNN_mean) / tst.x$null_meanNN_sd
    t.obs.NNI <- na.omit(t.obs.NNI)
    
    s.obs.NDI <- (obs.x$mean_disp - sst.x$null_meanDisp_mean) / sst.x$null_meanDisp_sd
    s.obs.NDI <- na.omit(s.obs.NDI)
    
    s.obs.NNI <- (obs.x$mean_NN - sst.x$null_meanNN_mean) / sst.x$null_meanNN_sd
    s.obs.NNI <- na.omit(s.obs.NNI)
    
    # conduct t-tests
    t.NDI.ttest <- t.test(t.obs.NDI, mu = 0, alternative = "two.sided")
    t.NNI.ttest <- t.test(t.obs.NNI, mu = 0, alternative = "two.sided")
    
    s.NDI.ttest <- t.test(s.obs.NDI, mu = 0, alternative = "two.sided")
    s.NNI.ttest <- t.test(s.obs.NNI, mu = 0, alternative = "two.sided")
    
    # compare observed NDIs to null NDIs
    t.NDI.mean <- mean(t.obs.NDI, na.rm = T)
    t.NDI.phi <- sum(tdi.x$null_NDI_mean > t.NDI.mean) / 1000
    t.NDI.plo <- sum(tdi.x$null_NDI_mean < t.NDI.mean) / 1000
    
    t.NNI.mean <- mean(t.obs.NNI, na.rm = T)
    t.NNI.phi <- sum(tdi.x$null_NNI_mean > t.NNI.mean) / 1000
    t.NNI.plo <- sum(tdi.x$null_NNI_mean < t.NNI.mean) / 1000
    
    s.NDI.mean <- mean(s.obs.NDI, na.rm = T)
    s.NDI.phi <- sum(sdi.x$null_NDI_mean > s.NDI.mean) / 1000
    s.NDI.plo <- sum(sdi.x$null_NDI_mean < s.NDI.mean) / 1000
    
    s.NNI.mean <- mean(s.obs.NNI, na.rm = T)
    s.NNI.phi <- sum(sdi.x$null_NNI_mean > s.NNI.mean) / 1000
    s.NNI.plo <- sum(sdi.x$null_NNI_mean < s.NNI.mean) / 1000
    
    # number of cells at this scale
    n.cells <- length(s.obs.NDI)
    
    table.dat[[j]] <- c(scl[j], n.cells, 
                        round(t.NDI.mean, 3), 
                        stars.pval(t.NDI.ttest$p.value), stars.pval(min(t.NDI.plo, t.NDI.phi)), 
                        round(t.NNI.mean, 3), 
                        stars.pval(t.NNI.ttest$p.value), stars.pval(min(t.NNI.plo, t.NNI.phi)),
                        round(s.NDI.mean, 3), 
                        stars.pval(s.NDI.ttest$p.value), stars.pval(min(s.NDI.plo, s.NDI.phi)),
                        round(s.NNI.mean, 3), 
                        stars.pval(s.NNI.ttest$p.value), stars.pval(min(s.NNI.plo, s.NNI.phi)))
  }
  # form all spatial scale data for one subset into one item
  table.dat <- do.call(rbind, table.dat)
  
  colnames(table.dat) <- c("scale", "n.cells", 
                           "t.meanNDI", "t.NDI.t-test", "t.NDI.null-test",
                           "t.meanNNI", "t.NNI.t-test", "t.NNI.null-test",
                           "s.meanNDI", "s.NDI.t-test", "s.NDI.null-test",
                           "s.meanNNI", "s.NNI.t-test", "s.NNI.null-test")
  
  rownames(table.dat) <- rep(s[i], 4)
  
  # store it into the BIG TABLE
  all.table.dat[[i]] <- table.dat
} 

# put it all together
all.table.dat <- do.call(rbind, all.table.dat)

all.table.dat

