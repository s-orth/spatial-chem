# TITLE: Create Figure 1 for Paper
# AUTHOR: Sarah Orth, Ostling Lab
# LAST MODIFIED: July 31, 2024
library(rlist)

setwd("/Users/Sarah/Library/CloudStorage/Box-Box/spatial-chem/simulation-results") # set working directory

get.dir <- "scenario-I/" # set folder containing scenario results
list.all <- list.files(path = get.dir) # list all files names in get.dir folder
print(list.all)
I <- list.all[c(1:10)]
I.slopes <- list.all[c(21:30)]

# Read in simulation results. tab-delimited file with a header
I <- lapply(X = paste(get.dir, I, sep = ""),
            FUN = read.csv)
I.slopes <- lapply(X = paste(get.dir, I.slopes, sep = ""),
                   FUN = read.csv)

get.dir <- "scenario-II/" # set folder containing scenario results
list.all <- list.files(path = get.dir) # list all files names in get.dir folder
print(list.all)
II <- list.all[c(1:10)]
II.slopes <- list.all[c(21:30)]

II <- lapply(X = paste(get.dir, II, sep = ""), 
             FUN = read.csv)
II.slopes <- lapply(X = paste(get.dir, II.slopes, sep = ""),
                    FUN = read.csv)

get.dir <- "scenario-III/" # set folder containing scenario results
list.all <- list.files(path = get.dir) # list all files names in get.dir folder
print(list.all)
III <- list.all[c(1:10)]
III.slopes <- list.all[c(21:30)]

III <- lapply(X = paste(get.dir, III, sep = ""), 
              FUN = read.csv)
III.slopes <- lapply(X = paste(get.dir, III.slopes, sep = ""),
                     FUN = read.csv)

get.dir <- "scenario-IV/" # set folder containing scenario results
list.all <- list.files(path = get.dir) # list all files names in get.dir folder
print(list.all)
IV <- list.all[c(1:10)]
IV.slopes <- list.all[c(21:30)]

IV <- lapply(X = paste(get.dir, IV, sep = ""), 
             FUN = read.csv)
IV.slopes <- lapply(X = paste(get.dir, IV.slopes, sep = ""),
                      FUN = read.csv)

I <- list.rbind(I) # Scenario I
II <- list.rbind(II) # Scenario II
III <- list.rbind(III) # Scenario III
IV <- list.rbind(IV) # Scenario IV

all.scenarios <- list(I, II, III, IV)
names(all.scenarios) <- c("I", "II", "III", "IV")

II.slopes <- list.rbind(II.slopes)
III.slopes <- list.rbind(III.slopes)
IV.slopes <- list.rbind(IV.slopes)
I.slopes <- list.rbind(I.slopes) # bind together slopes

# summing up simulations finding over-dispersion.
plot.dat <- list()

for(i in 1:4){
  
  dat <- all.scenarios[[i]] # select current scenario
  
  inds <- which(dat$UN.NDI.means > 0 & dat$UN.NDI.t.pvals < 0.05) # trait and t test
  tt <- dat[inds, ]

  inds <- which(dat$UN.NDI.means > 0 & dat$UN.NDI.dist.pvals.hi < 0.05) # trait and dist test
  td <- dat[inds, ]

  inds <- which(dat$CN.NDI.means > 0 & dat$CN.NDI.t.pvals < 0.05) # spatial and t test
  st <- dat[inds, ]

  inds <- which(dat$CN.NDI.means > 0 & dat$CN.NDI.dist.pvals.hi < 0.05) # trait and dist test
  sd <- dat[inds, ]

  # to look at detected levels of underdispersion
  inds <- which(dat$UN.NDI.means < 0 & dat$UN.NDI.t.pvals < 0.05) # trait and t test
  # tt <- dat[inds, ]
  # 
  # inds <- which(dat$UN.NDI.means < 0 & dat$UN.NDI.dist.pvals.lo < 0.05) # trait and dist test
  # td <- dat[inds, ]
  # 
  # inds <- which(dat$CN.NDI.means < 0 & dat$CN.NDI.t.pvals < 0.05) # spatial and t test
  # st <- dat[inds, ]
  # 
  # inds <- which(dat$CN.NDI.means < 0 & dat$CN.NDI.dist.pvals.lo < 0.05) # trait and dist test
  # sd <- dat[inds, ]
  
  tots <- c(nrow(tt), nrow(td), nrow(st), nrow(sd))
  plot.dat[[i]] <- tots

}
plot.dat <- list.rbind(plot.dat)
plot.dat <- as.data.frame(plot.dat)
rownames(plot.dat) <- names(all.scenarios)
colnames(plot.dat) <- c("tt", "td", "st", "sd")

# set colors for plotting
#pal <- palette.colors(palette = "Set1", n = 4)
pal = c("lightgreen", "purple3", "cornflowerblue", "red")

par(mfrow = c(3, 2),
    mar = c(5, 4, 1.5, 1.5),
    cex.axis = 1,
    cex.lab = 1.2,
    cex.main = 1.4)

# conceptual figure, inspired by actual saved slopes
I.slopes <- seq(from = min(I.slopes$x), to = max(I.slopes$x), by = 0.1)
II.slopes <- seq(from = min(II.slopes$x), to = max(II.slopes$x), by = 0.001)
III.slopes <- seq(from = min(III.slopes$x), to = max(III.slopes$x), by = 0.1)
IV.slopes <- seq(from = min(IV.slopes$x), to = max(IV.slopes$x), by = 0.5)

plot(x = 1, 
     type = "n",
     xlim = c(0, 1),
     ylim = c(-5, 6),
     xlab = "species mean chemical disparity",
     ylab = "species occupancy rate",
     main = "a) Simulation scenarios",
     xaxs = "i",
     yaxs = "i", bty = "n")
for(i in 1:length(I.slopes)){
  abline(a = 0.5, b = I.slopes[i], col = pal[1])
}
for(i in 1:length(II.slopes)){
  abline(a = 0.5, b = II.slopes[i], col = pal[2])
}
for(i in 1:length(III.slopes)){
  abline(a = 0.5, b = III.slopes[i], col = pal[3])
}
for(i in 1:length(IV.slopes)){
  abline(a = 0.5, b = IV.slopes[i], col = pal[4])
}

plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1)) 
legend("bottomleft",
       legend = c("I: no constraints",
                  "II: -0.001 < m < 0.001",
                  "III: 0.5 < m < 1.5",
                  "IV: m > 1.5"),
       pch = 16,
       pt.cex = 2.5,
       cex = 1.5,
       bty = "n",
       col = pal)
text(x = 0.2, y = 0.75, cex = 1.8, labels = "Scenarios")

# trait null and t test figure
b <- barplot(plot.dat$tt, 
             names = rownames(plot.dat),
             col = pal,
             xlab = "Scenario", 
             ylab = "# sims sig. overdispersed",
             # ylab = "# sims sig. underdispersed",
             main = "b) trait null, t-test",
             las = 1,
             ylim = c(0, 1000))
text(x = b, y = plot.dat$tt,
     labels = plot.dat$tt,
     pos = 3,
     # offset = 1.5) # for underdisperion fig.
     offset = -1.5)

# spatial null and t test figure
b <- barplot(plot.dat$st, 
             names = rownames(plot.dat),
             col = pal,
             xlab = "Scenario", 
             ylab = "# sims sig. overdispersed",
             # ylab = "# sims sig. underdispersed",
             main = "c) spatial null, t-test",
             las = 1,
             ylim = c(0, 1000))
text(x = b, y = plot.dat$st,
     labels = plot.dat$st,
     pos = 3,
     offset = 1.5)

# trait null and distribution test figure
b <- barplot(plot.dat$td, 
             names = rownames(plot.dat),
             col = pal,
             xlab = "Scenario", 
             ylab = "# sims sig. overdispersed",
             # ylab = "# sims sig. underdispersed",
             main = "d) trait null, distribution test",
             las = 1,
             ylim = c(0, 1000))
text(x = b, y = plot.dat$td,
     labels = plot.dat$td,
     pos = 3,
     offset = 1.5)

# spatial and dist test figure
b <- barplot(plot.dat$sd, 
             names = rownames(plot.dat),
             col = pal,
             xlab = "Scenario", 
             ylab = "# sims sig. overdispersed",
             # ylab = "# sims sig. underdispersed",
             main = "e) spatial null, distribution test",
             las = 1,
             ylim = c(0, 1000))
text(x = b, y = plot.dat$sd,
     labels = plot.dat$sd,
     pos = 3,
     offset = 1.5)
