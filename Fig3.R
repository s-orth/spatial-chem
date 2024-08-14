# TITLE: Create Figures for Paper
# AUTHOR: Sarah Orth, Ostling Lab
# LAST MODIFIED: August 12, 2024

library(gtools) # stars.pvals function
library(Rfast)
library(ggplot2)
library(tidyr)
library(forcats)
library(dplyr)
library(gridExtra)

setwd("/Users/Sarah/Library/CloudStorage/Box-Box/spatial-chem/") # set working directory

s <- c("all-species", "genus-inga", "genus-miconia", "genus-psychotria") # list of subsets
obs.dat <- list() # list to store data lists (lists within a list)
null.stats <- list() # list to store null summary data lists
null.dists <- list() # list to store null distribution lists

# for loop to read in relevant data
for(i in 1:4){
  my.dir <- paste0("quadrat-values-2010/", s[i], "/") # subset path
  list.all <- list.files(path = my.dir) # list files
  list.all <- list.all[c(2:4, 1)] # reorder
  obs.dat[[i]] <- lapply(X = paste(my.dir, list.all, sep = ""), # read in data
                    FUN = read.table, 
                    sep = "\t", header = T) # tab-delim file w/ header
  
  my.dir <- paste0("spatial-nulls-2010/", s[i], "/") # subset path for null
  list.all <- list.files(path = my.dir) # list files
  list.null.stats <- list.all[c(27:29, 26)] # select null stats
  list.null.dists <- list.all[c(7:9, 6)] # select null distributions
  null.stats[[i]] <- lapply(X = paste(my.dir, list.null.stats, sep = ""), # read in data
                      FUN = read.table, 
                      sep = "\t", header = T)
  null.dists[[i]] <- lapply(X = paste(my.dir, list.null.dists, sep = ""), # read in data
                     FUN = read.table, 
                     sep = "\t", header = T)
}

# name list items, just to help keep things straight
names(obs.dat) <- s
names(null.stats) <- s
names(null.dists) <- s
plot.titles <- c("a) All Species", "b) Inga", "c) Miconia", "d) Psychotria")

all.box.dat <- list()
all.plot.dat <- list()

for(i in 1:4){ # going through each subset...
  
  # pull out one subset for the observed and each null. Now we have plain old lists, not nested lists.
  obs <- obs.dat[[i]]
  sta <- null.stats[[i]]
  dis <- null.dists[[i]]
  
  plot.dat <- list()
  
  for(j in 1:4){ # now going through each spatial scale at the subset...
    
    # select one spatial scale for each thing. now we have plain old vectors
    obs.x <- obs[[j]]
    sta.x <- sta[[j]]
    dis.x <- dis[[j]]
    
    # compute observed NDIs - FOR MAIN TEXT FIGURE
    obs.NDI <- (obs.x$mean_disp - sta.x$null_meanDisp_mean) / sta.x$null_meanDisp_sd
    obs.NDI <- na.omit(obs.NDI)

    # compare observed NDIs to null NDIs
    obs.NDI.mean <- mean(obs.NDI, na.rm = T)
    p_hi <- sum(dis.x$null_NDI_mean > obs.NDI.mean) / 1000
    p_lo <- sum(dis.x$null_NDI_mean < obs.NDI.mean) / 1000

    plot.dat[[j]] <- c(obs.NDI.mean, stars.pval(p_lo), stars.pval(p_hi))
  }
  
  # back to our per subset situation... lets prepare plotting data.
  box.dat <- cbind(dis[[1]]$null_NDI_mean, dis[[2]]$null_NDI_mean, dis[[3]]$null_NDI_mean, dis[[4]]$null_NDI_mean)
  box.dat <- as.data.frame(box.dat)
  colnames(box.dat) <- c("10m", "20m", "50m", "100m")
  box.dat <- gather(box.dat, "scale", "value", 1:4)
  box.dat <- box.dat %>%
    mutate(scale = fct_relevel(scale, "10m", "20m", "50m", "100m"))
  
  all.box.dat[[i]] <- box.dat

  plot.dat <- do.call(rbind, plot.dat)
  plot.dat <- as.data.frame(plot.dat)
  colnames(plot.dat) <- c("obs.NDI", "pval.lo", "pval.hi")
  plot.dat$scale <- c("10m", "20m", "50m", "100m")
  plot.dat$obs.NDI <- as.numeric(plot.dat$obs.NDI)
  plot.dat <- plot.dat %>%
    mutate(scale = fct_relevel(scale, "10m", "20m", "50m", "100m"))
  
  all.plot.dat[[i]] <- plot.dat
  
}

names(all.box.dat) <- s
names(all.plot.dat) <- s

p1 <- ggplot(all.box.dat[[1]], 
             aes(x = scale, 
                 y = value)) +
  geom_violin(fill = "lightgray") +
  stat_summary(fun = "mean", 
               geom = "crossbar",
               width = 0.4, 
               fatten = 0.7, 
               color = "black") +
  geom_point(data = all.plot.dat[[1]], 
             aes(x = scale, 
                 y = obs.NDI), 
             color = "red", 
             size = 1) +
  ylab("mean NDI") +
  ylim(-0.05, 0.05) +
  geom_text(data = all.plot.dat[[1]], 
            aes(label = pval.hi, 
                y = 0.045), 
            size = 8) +
  geom_text(data = all.plot.dat[[1]], 
            aes(label = pval.lo, 
                y = 0.045),
            size = 8) +
  ggtitle(plot.titles[1]) +
  theme_bw()

p2 <- ggplot(all.box.dat[[2]], 
             aes(x = scale, 
                 y = value)) +
  geom_violin(fill = "lightgray") +
  stat_summary(fun = "mean", 
               geom = "crossbar",
               width = 0.4, 
               fatten = 0.7, 
               color = "black") +
  geom_point(data = all.plot.dat[[2]], 
             aes(x = scale, 
                 y = obs.NDI), 
             color = "red", 
             size = 1) +
  ylab("mean NDI") +
  ylim(-0.4, 0.4) +
  geom_text(data = all.plot.dat[[2]], 
            aes(label = pval.hi, 
                y = 0.35), 
            size = 8) +
  geom_text(data = all.plot.dat[[2]], 
            aes(label = pval.lo, 
                y = 0.35),
            size = 8) +
  ggtitle(plot.titles[2]) +
  theme_bw()

p3 <- ggplot(all.box.dat[[3]], 
             aes(x = scale, 
                 y = value)) +
  geom_violin(fill = "lightgray") +
  stat_summary(fun = "mean", 
               geom = "crossbar",
               width = 0.4, 
               fatten = 0.7, 
               color = "black") +
  geom_point(data = all.plot.dat[[3]], 
             aes(x = scale, 
                 y = obs.NDI), 
             color = "red", 
             size = 1) +
  ylab("mean NDI") +
  ylim(-0.4, 0.4) +
  geom_text(data = all.plot.dat[[3]], 
            aes(label = pval.hi, 
                y = 0.35), 
            size = 8) +
  geom_text(data = all.plot.dat[[3]], 
            aes(label = pval.lo, 
                y = 0.35),
            size = 8) +
  ggtitle(plot.titles[3]) +
  theme_bw()

p4 <- ggplot(all.box.dat[[4]], 
             aes(x = scale, 
                 y = value)) +
  geom_violin(fill = "lightgray") +
  stat_summary(fun = "mean", 
               geom = "crossbar",
               width = 0.4, 
               fatten = 0.7, 
               color = "black") +
  geom_point(data = all.plot.dat[[4]], 
             aes(x = scale, 
                 y = obs.NDI), 
             color = "red", 
             size = 1) +
  ylab("mean NDI") +
  ylim(-0.4, 0.4) +
  geom_text(data = all.plot.dat[[4]], 
            aes(label = pval.hi, 
                y = 0.35), 
            size = 8) +
  geom_text(data = all.plot.dat[[4]], 
            aes(label = pval.lo, 
                y = 0.35),
            size = 8) +
  ggtitle(plot.titles[4]) +
  theme_bw()


# plot.outputs
grid.arrange(p1, p2, p3, p4,
             nrow = 2)
