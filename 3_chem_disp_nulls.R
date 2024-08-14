# TITLE: BCI Chemical Disparity Nulls
# AUTHOR: Sarah Orth, Ostling Lab
# LAST MODIFIED: August 12, 2024

# DESCRIPTION: Generate null distributions by randomizing which species are in which quadrats. In the trait null, only quadrat species richness will be preserved. In the spatial null, quadrat species richness will be preserved, as well as species occupancy rates. Presence-Absence matrices will be randomized using a very slightly modified version of the curveball algorithm from the 'curveball' package in R for the spatial null. See Strona et al. 2014. 

# SET-UP ----

# load necessary packages
library(matrixStats) 
library(Rfast)
library(cooccur)
library(reshape2)

# set working directory
setwd("/Users/Sarah/Library/CloudStorage/Box-Box/spatial-chem/") 

# set folder containing presence-absence matrices
pres.abs.dir <- "BCI-PA-matrices-2010/" 

# list all files names in pres.abs.dir folder
list.all <- list.files(path = pres.abs.dir) 

# check order of file names
print(list.all)

# put files in order from smallest to largest spatial scale
list.all <- list.all[c(2:4, 1)]

# confirm correct order of files
print(list.all) 

# Read in presence-absence matrices.
bci.mat <- lapply(X = paste(pres.abs.dir, list.all, sep = ""), 
                  FUN = read.table, # function to apply to each list item
                  sep = "\t", # tab-delimited
                  header = T) # header present

# read in chemical disparity matrix
bci.cd <- read.table(
   file = "clean-BCI-by-quadrat-2010/disparity_matrix_defense_only.txt",
  # file = "raw-data/kursar_disparity_matrix_all_defense.txt",
  sep = "\t", # tab-delimited
  row.names = 1, # row names are in column 1
  header = T) # header present

# ignore conspecific relationships, as this may dilute metrics.
diag(bci.cd) <- NA

# number of times to run the nulls
num.nulls <- 1000

# creating directories for saving null results
dir.create(paste0(getwd(), "/spatial-nulls-2010"))
dir.create(paste0(getwd(), "/trait-nulls-2010"))

# FUNCTION FOR SUBSETTING THE PRESENCE-ABSENCE DATA ----
Get.Sub <- function(d, g = NULL) {
  # d is a data frame containing the presence-absence matrix
  # g is a string, the name of a genus of interest, such as "inga". Default value is NULL
  
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

# FUNCTION FOR CALCULATING CELL-LEVEL CHEMICAL DISPARITY ----

# This function will obtain mean nearest neighbor and mean chemical disparity scores for each quadrat, at each spatial scale (when used in conjunction with lapply()). 
# The required input is a list of presence absence matrices that have species IDs as the column names and cell IDs as the rows.

# we will use lapply() to iterate this function over each PA matrix in the list
Get.Null.Vals <- function(d) { 
  
  # create an empty list to hold species IDs present in each quadrat
  sp.list <- list() 
  
  # for each row (quadrat) in the matrix
  # create an entry in sp.list containing a list of the species present in the quadrat
  for(row in 1:dim(d)[1]) { sp.list[[row]] <- (which(d[row, ] == 1)) } 
  
  # create an empty list to hold chemical data fpr the quadrat
  chem.list <- list() 
  
  # for each quadrat stored in sp.list...
  for(q in 1:length(sp.list)) { 
    
    # if there are 2 or more species in the quadrat
    if(length(sp.list[[q]]) > 1) { 
      
      # filter chem matrix have only species in quadrat "q" and convert to data matrix
      quad.chm <- data.matrix(sub.cd[names(sp.list[[q]]), # indicate rows to keep
                                     names(sp.list[[q]])]) # indicate cols to keep
      
      # extract row minimum - this is each species nearest neighbor
      nn <- matrixStats::rowMins(quad.chm, na.rm = TRUE)
      
      # calculate the mean of the nearest neighbor scores
      mean.nn <- mean(nn) 
      
      # make the entire upper triangle of the quad.chm matrix NA to remove duplicates
      quad.chm[upper.tri(x = quad.chm, diag = TRUE)] <- NA 
      
      # calculate the mean chemical disparity among species present in the quadrat
      mean.cd <- mean(as.matrix(quad.chm), 
                      na.rm = T) # ignore NAs
      
      # store mean NN and mean chemical disparity in quad.chm
      quad.chm <- c(mean.nn, mean.cd)
      
      # if there are 0 or 1 species in the quadrat, store NA values for chemistry in quad.chm  
    } else { quad.chm <- c(NA, NA) } 
    
    # store finished chem information for quadrat "q" in chem.list 
    chem.list[[q]] <- quad.chm 
  } 
  
  # turn list of quadrat information into a data frame for all quadrats at a spatial scale
  d <- as.data.frame(matrix(unlist(chem.list), 
                            ncol = 2, # number of columns in new data frame
                            nrow = length(chem.list), # number of rows in new data frame
                            byrow = T, # list by row
                            dimnames = list(NULL, c("mean_NN", "mean_disp")))) # name columns
  
  # create a quadrat ID variable
  d$quadrat_id <- c(1:nrow(d)) 
  
  # re-arrange columns
  d <- d[ , c(3, 1:2)] 
  
}

# FUNCTION FOR SHUFFLING PRESENCE-ABSENCE MATRIX ----
# modification of the curveball function in the R package "backbone"
Curveball.Mod <- function(m){
  
  # extract the dimensions (rows and columns) of the matrix m. returns vector of 2
  RC <- dim(m)
  
  # R is the number of rows in matrix m, first element of RC. returns integer.
  R <- RC[1]
  
  # C is the number of columns in matrix m, second element of RC. returns integer.
  C <- RC[2]
  
  # create empty list called hp
  hp <- list()
  
  # save matrix names
  mat.nms <- colnames(m) 
  
  # Next part loops down all rows in matrix m
  # for each row in data matrix m, creates a list indicating which species are in the quadrat.
  # This is a vector of integers. 
  # In their paper, Strona et. al. say that this should be inverted if there are more quadrats than species - such that there would be a list of which quadratss a species is present in rather than which species are present in a quadrat. 
  # To do this, simply iterate through columns instead of rows. 
  # See alternative code below. SO automated this with an if/else statement that asks if number of rows (quadrats/sites) is greater than number of columns (species). Mostly there will be more rows than columns (more quadrats than species)
  
  # if there are more species than quadrats, create lists of which species are in each quadrat
  if(C >= R) { 
    for (row in 1:dim(m)[1]) { hp[[row]] = (which(m[row, ] == 1)) }
    
  # if there are more quadrats than species, create lists of which quadrats each species is in
  } else { 
    for (col in 1:dim(m)[2]) { hp[[col]] = (which(m[, col] == 1)) }
  }
  
  # returns the length of list hp as integer.
  l_hp <- length(hp) 
  
  # notes refer to example where cells > species. 
  # Invert language for instance of cells < species.
  
  # for how many repetitions... do what is below
  for (rep in 1:(5*l_hp)) { 
    
    # pick two random species presences lists from hp and store integers in AB
    AB <- sample(1:l_hp, 2)
    
    # for the first random species (A), list the quadrats that species is present in
    a <- hp[[AB[1]]] 
    
    # for the second random species (B), list the quadrats that species is present in 
    b <- hp[[AB[2]]] 
    
    # find the quadrats that both species A and B are present in and return vector of integers
    ab <- intersect(a, b)
    
    # number of quadrats that species A and B have in common (length of ab) - integer
    l_ab <- length(ab)
    
    # number of quadrats that species A is present in - integer
    l_a <- length(a) 
    
    # number of quadrats that species B is present in - integer
    l_b = length(b) 
    
    # As long as the number of quadrats that A and B have in common is not equal to the number of quadrats A is in or the number of quadrats B is in, do the following:
    if ((l_ab %in% c(l_a, l_b)) == F){
      
      # vector of the quadrats that species A and B are NOT both present in
      tot <- setdiff(c(a,b), ab)
      
      # number of quadrats that species A and B are NOT both present in
      l_tot <- length(tot)
      
      # shuffle the quadrats that are in tot
      tot <- sample(tot, l_tot, replace = FALSE, prob = NULL)
      
      # how many quadrats for species A do we need to fill?
      L <- l_a - l_ab 
      
      # for the first randomly selected presences list, replace in hp a new shuffled list in which the quadrats that both species were in are all there, and then from the shuffled set of quadrats that the species did not have in common put in enough to make sure the richness of that quadrat/occupancy rate of the species is maintained.
      hp[[AB[1]]] <- c(ab, tot[1:L])
      
      # do the same thing for the second randomly selected presences list but use the things from tot that did not get used in the first list.
      hp[[AB[2]]] <- c(ab,tot[(L+1):l_tot])} 
  }
  
  # create an empty matrix rm (randomized matrix) with the same dimensions as original matrix m
  rm <- matrix(0, R, C)
  
  # for all of the rows in rm, fill in with the new randomized rows that are stored as vectors in list hp. This also needs to be changed if you had more cells than species. 
  if( C >= R ) { 
    for (row in 1:R) { rm[row, hp[[row]]] = 1 }
  } else { 
    for (col in 1:C) { rm[hp[[col]], col] = 1 }}
  
  # add names back to matrix
  colnames(rm) <- mat.nms 
  m <-rm
}

# APPLYING FUNCTION Get.Sub ----

# NOTE: This code is not set up to do all nulls for an entire group of subsets, because the nulls are more time consuming. User will need to return to the code to specify the next subset.

# Leave g blank if you want to do all species
Sub.Mat <- lapply(X = bci.mat, 
                  FUN = Get.Sub) # function to apply
                  #g = "inga") # g can be set to "inga", "miconia", or "psychotria". g must be commented out for all species.

# unique identifier for subset to be used in writing files
s <- "all-species" # suggested conventions: "all-species", "genus-inga", "genus-miconia", "genus-psychotria"

# create the folders where outputs will be saved
dir.create(paste0(getwd(), "/spatial-nulls-2010/", s))
dir.create(paste0(getwd(), "/trait-nulls-2010/", s))

# grab row and column indices of species in subset in chemistry matrix
Rinds <- which(rownames(bci.cd) %in% colnames(Sub.Mat[[4]]))
Cinds <- which(colnames(bci.cd) %in% colnames(Sub.Mat[[4]]))

# filter chemistry matrix to have only species in subset
sub.cd <- data.matrix(bci.cd[Rinds, Cinds])

# filter Sub.Mat to have only species in the disparity matrix (may need these lines if using a different chemical disparity matrix or census year)
# inds <- which(colnames(Sub.Mat[[4]]) %in% colnames(bci.cd))
# Sub.Mat <- Map(function(x){x[ , inds]}, Sub.Mat)


# GENERATE NULL PRESENCE ABSENCE MATRICES ----

# Before we start running nulls, let's generate some null presence-absence matrices. These will be used to compare observed spatial patterns to null spatial patterns, namely occupancy rates. They come into play for some figures that are in the supplement. 

dir.create(paste0(getwd(), "/spatial-nulls-2010/", s, "/", "10m-nullPAs")) # create a folder for them

spatial.Null.PAs <- list() # create a list to store outputs

d <- Sub.Mat[[1]] # 1 = 10m spatial scale. This is the only scale we are interested in

for(i in 1:10){ # just saving 10
  d <- Curveball.Mod(d) # shuffle the matrix 
  spatial.Null.PAs[[i]] <- d # store it
}

nms <- c("nullPA_10m_1", "nullPA_10m_2", "nullPA_10m_3", "nullPA_10m_4", "nullPA_10m_5",
         "nullPA_10m_6", "nullPA_10m_7", "nullPA_10m_8", "nullPA_10m_9", "nullPA_10m_10")

names(spatial.Null.PAs) <- nms # assign the names

lapply(names(spatial.Null.PAs), function(x) { # actual chunk to save code
  f <- spatial.Null.PAs[[x]]
  write.table(f, file = paste0(getwd(),'/spatial-nulls-2010/',
                               s, '/', "10m-nullPAs/", x, '.txt'),
              sep = "\t", # tab delimited
              row.names = F, # no row names
              col.names = T) # yes column names
})


# RUN spatial NULLS ----

# create empty listS to store null results for each spatial scale
spatial.Null.Results <- list() # summary statistics of mean NN and mean Disparity by quadrat
spatial.Null.Distributions <- list() # null distribution of NNI and NDI for statistical tests
spatial.Null.NNIs <- list() # full output of null NNIs by quadrat
spatial.Null.NDIs <- list() # full output of null NDIs by quadrat
spatial.Null.Occupancies <- list() # summary statistics of null occupancy rates
spatial.Null.Occupancy.Full <- list() # full output of each pair's occupancy rate in null
spatial.Null.OccZscores <- list() # null distribution of z-scored occupancy rates

# loop to calculate all spatial scales #### WILL TAKE SEVERAL MINUTES
for(scale in 1:length(Sub.Mat)) {

  # return which spatial scale is being worked on
  print(paste("Working on spatial null spatial scale ", scale, "...", sep = ""))

  # pull out one spatial scale at a time
  d <- Sub.Mat[[scale]]

  # create empty list to store null mean NN values and mean disparity values
  Null.meanNN.list <- list()
  Null.meanDisp.list <- list()
  Null.occRate.list <- list()

  # for each run of the null
  for(run.num in 1:num.nulls){

    # Shuffle the presence-absence matrix
    d <- Curveball.Mod(d)

    # apply Get.Null.Vals function to return quadrat-level metrics
    Null.Vals <- Get.Null.Vals(d)

    # store mean NN values and mean disparity values in their respective lists
    Null.meanNN.list[[run.num]] <- Null.Vals$mean_NN
    Null.meanDisp.list[[run.num]] <- Null.Vals$mean_disp

    # also, calculate null co-occupancy rates and store for later use
    occ.rates <- create.N.matrix(t(d)) # returns a pairwise matrix of number of co-occupancies
    colnames(occ.rates) <- colnames(d) # add column names
    rownames(occ.rates) <- colnames(d) # add row names
    occ.rates <- occ.rates/(matrixStats::count(rowSums(d) > 0)) # convert to rate, with # cells occupied by any member of the group as denom

    occ.rates[upper.tri(x = occ.rates, diag = TRUE)] <- NA # convert upper triangle to NA to avoid duplicates
    occ.rates <- reshape2::melt(occ.rates) # transform into a data frame
    occ.rates <- na.omit(occ.rates) # remove NA values

    Null.occRate.list[[run.num]] <- occ.rates$value # only save values in list, will add names later
  }

  # bind mean NN and mean Disp lists into one (does this for a whole spatial scale)
  # rows will be cell IDs, columns are the null run
  Null.meanNN <- matrix(unlist(Null.meanNN.list), ncol = num.nulls)
  Null.meanDisp <- matrix(unlist(Null.meanDisp.list), ncol = num.nulls)
  Null.occRate <- matrix(unlist(Null.occRate.list), ncol = num.nulls)

  # summarize spatial null data and organize into a savable format. Metrics summarize the nulls of a given quadrat

  # calculate mean, standard deviation, and percentiles of null mean NN values
  null_meanNN_mean <- rowMeans2(Null.meanNN)
  null_meanNN_sd <- rowSds(Null.meanNN)
  null_meanNN_5pct <- rownth(Null.meanNN,
                             elems = as.integer(rep(50, dim(Null.meanNN)[1])))
  null_meanNN_95pct <- rownth(Null.meanNN,
                              elems = as.integer(rep(50, dim(Null.meanNN)[1])),
                              descending = TRUE)
  null_meanNN_2.5pct <- rownth(Null.meanNN,
                               elems = as.integer(rep(25, dim(Null.meanNN)[1])))
  null_meanNN_97.5pct <- rownth(Null.meanNN,
                                elems = as.integer(rep(25, dim(Null.meanNN)[1])),
                                descending = TRUE)

  # calculate mean, standard deviation, and percentiles of null mean disparity values
  null_meanDisp_mean <- rowMeans2(Null.meanDisp)
  null_meanDisp_sd <- rowSds(Null.meanDisp)
  null_meanDisp_5pct <- rownth(Null.meanDisp,
                               elems = as.integer(rep(50, dim(Null.meanDisp)[1])))
  null_meanDisp_95pct <- rownth(Null.meanDisp,
                                elems = as.integer(rep(50, dim(Null.meanDisp)[1])),
                                descending = TRUE)
  null_meanDisp_2.5pct <- rownth(Null.meanDisp,
                                 elems = as.integer(rep(25, dim(Null.meanDisp)[1])))
  null_meanDisp_97.5pct <- rownth(Null.meanDisp,
                                  elems = as.integer(rep(25, dim(Null.meanDisp)[1])),
                                  descending = TRUE)

  # calculate mean, standard deviation, and percentiles of null co-ocurrences
  null_occRate_mean <- rowMeans2(Null.occRate)
  null_occRate_sd <- rowSds(Null.occRate)
  null_occRate_5pct <- rownth(Null.occRate,
                              elems = as.integer(rep(50, dim(Null.occRate)[1])))
  null_occRate_95pct <- rownth(Null.occRate,
                               elems = as.integer(rep(50, dim(Null.occRate)[1])),
                               descending = TRUE)
  null_occRate_2.5pct <- rownth(Null.occRate,
                                elems = as.integer(rep(25, dim(Null.occRate)[1])))
  null_occRate_97.5pct <- rownth(Null.occRate,
                                 elems = as.integer(rep(25, dim(Null.occRate)[1])),
                                 descending = TRUE)

  # compare nulls to themselves to generate a distribution of expected NNI and NDI scores as well as z-scored occupancy rates
  
  null_NNI_median <- vector() # spaces to store data
  null_NDI_median <- vector()

  null_NNI_mean <- vector()
  null_NDI_mean <- vector()

  null_NNIs <- list()
  null_NDIs <- list()
  null_occ_zscores <- list()

  for(i in 1:num.nulls){
    # select one column (one run of the null) and compare to the null distribution as a whole, creating a set of 1000 "null" NNI and NDI values, and z-scored occupancy rates. For each run of the null, the list will be the number of cells at that scale. Take the median to find a distribution of median NNI and NDI values.

    null_NNIs[[i]] <- (Null.meanNN[, i] - null_meanNN_mean) / null_meanNN_sd
    null_NDIs[[i]] <- (Null.meanDisp[, i] - null_meanDisp_mean) / null_meanDisp_sd
    null_NNI_median[i] <- median((Null.meanNN[, i] - null_meanNN_mean) / null_meanNN_sd,
                             na.rm = T)
    null_NDI_median[i] <- median((Null.meanDisp[, i] - null_meanDisp_mean) / null_meanDisp_sd,
                             na.rm = T)
    null_NNI_mean[i] <- mean((Null.meanNN[, i] - null_meanNN_mean) / null_meanNN_sd,
                                 na.rm = T)
    null_NDI_mean[i] <- mean((Null.meanDisp[, i] - null_meanDisp_mean) / null_meanDisp_sd,
                                 na.rm = T)
    
    null_occ_zscores[[i]] <- (Null.occRate[, i] - null_occRate_mean) / null_occRate_sd
  }

  # combine all results into one data frame
  results <- cbind(null_meanNN_mean, null_meanNN_sd,
                   null_meanNN_5pct, null_meanNN_95pct,
                   null_meanNN_2.5pct, null_meanNN_97.5pct,
                   null_meanDisp_mean, null_meanDisp_sd,
                   null_meanDisp_5pct, null_meanDisp_95pct,
                   null_meanDisp_2.5pct, null_meanDisp_97.5pct)

  null_NNI_NDI_distributions <- cbind(null_NNI_median, null_NDI_median,
                                      null_NNI_mean, null_NDI_mean)

  null_NNIs <- do.call(cbind, null_NNIs)
  null_NDIs <- do.call(cbind, null_NDIs)

  null_occRates_save <- cbind(occ.rates[, c(1:2)], null_occRate_mean, null_occRate_sd, # also save pair ids for matching
                              null_occRate_5pct, null_occRate_95pct,
                              null_occRate_2.5pct, null_occRate_97.5pct)
  
  null_occRate_distributions <- cbind(occ.rates[, c(1:2)], Null.occRate)
  null_occRate_zscores <- do.call(cbind, null_occ_zscores)

  spatial.Null.Results[[scale]] <- results
  spatial.Null.Distributions[[scale]] <- as.data.frame(null_NNI_NDI_distributions)
  spatial.Null.NNIs[[scale]] <- as.data.frame(null_NNIs)
  spatial.Null.NDIs[[scale]] <- as.data.frame(null_NDIs)
  spatial.Null.Occupancies[[scale]] <- null_occRates_save
  spatial.Null.Occupancy.Full[[scale]] <- null_occRate_distributions
  spatial.Null.OccZscores[[scale]] <- cbind(occ.rates[, c(1:2)], null_occRate_zscores)

}

# SAVE spatial NULL RESULTS!

# names for each item in spatial.Null.NNIs
nms <- c("null_NNIs_10m", "null_NNIs_20m", "null_NNIs_50m", "null_NNIs_100m")

names(spatial.Null.NNIs) <- nms

# names for each item in spatial.Null.NDIs
nms <- c("null_NDIs_10m", "null_NDIs_20m", "null_NDIs_50m", "null_NDIs_100m")

names(spatial.Null.NDIs) <- nms

# names for each item in spatial.Null.Results
nms <- c("null_stats_10m", "null_stats_20m", "null_stats_50m", "null_stats_100m")

names(spatial.Null.Results) <- nms

# names of reach item in spatial.Null.Distributions
nms <- c("null_NNI_NDI_dist_10m", "null_NNI_NDI_dist_20m", "null_NNI_NDI_dist_50m",
         "null_NNI_NDI_dist_100m")

names(spatial.Null.Distributions) <- nms

# names for each item in spatial.Null.Occupancies
nms <- c("null_occ_stats_10m", "null_occ_stats_20m", "null_occ_stats_50m", 
         "null_occ_stats_100m")

names(spatial.Null.Occupancies) <- nms

# names for each item in spatial.Null.Occupancy.Full
nms <- c("null_occs_10m", "null_occs_20m", "null_occs_50m", "null_occs_100m")

names(spatial.Null.Occupancy.Full) <- nms

# names for each item in spatial.Null.OccZscores
nms <- c("null_occs_null_Zscores_10m", "null_occs_null_Zscores_20m", 
         "null_occs_null_Zscores_50m", "null_occs_null_Zscores_100m")

names(spatial.Null.OccZscores) <- nms

# save spatial null results as separate files for each spatial scale

# for each item in spatial.Null.NNIs
lapply(names(spatial.Null.NNIs), function(x) {
  f <- spatial.Null.NNIs[[x]]
  write.table(f, file = paste0(getwd(),'/spatial-nulls-2010/',
                               s, '/', x, '.txt'),
              sep = "\t", # tab delimited
              row.names = F, # no row names
              col.names = T) # yes column names
})

# for each item in spatial.Null.NDIs
lapply(names(spatial.Null.NDIs), function(x) {
  f <- spatial.Null.NDIs[[x]]
  write.table(f, file = paste0(getwd(),'/spatial-nulls-2010/',
                               s, '/', x, '.txt'),
              sep = "\t", # tab delimited
              row.names = F, # no row names
              col.names = T) # yes column names
})

# for each item in spatial.Null.Results
lapply(names(spatial.Null.Results), function(x) {
  f <- spatial.Null.Results[[x]]
  write.table(f, file = paste0(getwd(),'/spatial-nulls-2010/',
                               s, '/', x, '.txt'),
              sep = "\t", # tab delimited
              row.names = F, # no row names
              col.names = T) # yes column names
})

# for each item in spatial.Null.Distributions
lapply(names(spatial.Null.Distributions), function(x) {
  f <- spatial.Null.Distributions[[x]]
  write.table(f, file = paste0(getwd(),'/spatial-nulls-2010/',
                               s, '/', x, '.txt'),
              sep = "\t", # tab delimited
              row.names = F, # no row names
              col.names = T) # yes column names
})

# for each item in spatial.Null.Occupancies
lapply(names(spatial.Null.Occupancies), function(x) {
  f <- spatial.Null.Occupancies[[x]]
  write.table(f, file = paste0(getwd(),'/spatial-nulls-2010/',
                               s, '/', x, '.txt'),
              sep = "\t", # tab delimited
              row.names = F, # no row names
              col.names = T) # yes column names
})

# for each item in spatial.Null.Occupancy.Full
lapply(names(spatial.Null.Occupancy.Full), function(x) {
  f <- spatial.Null.Occupancy.Full[[x]]
  write.table(f, file = paste0(getwd(),'/spatial-nulls-2010/',
                               s, '/', x, '.txt'),
              sep = "\t", # tab delimited
              row.names = F, # no row names
              col.names = T) # yes column names
})

# for each item in spatial.Null.Occ.Zscores
lapply(names(spatial.Null.OccZscores), function(x) {
  f <- spatial.Null.OccZscores[[x]]
  write.table(f, file = paste0(getwd(),'/spatial-nulls-2010/',
                               s, '/', x, '.txt'),
              sep = "\t", # tab delimited
              row.names = F, # no row names
              col.names = T) # yes column names
})

# RUN trait NULLS ----

# create empty lists to store trait null results for each spatial scale. Naming conventions similar to that in the spatial null.
trait.Null.Results <- list() 
trait.Null.Distributions <- list()
trait.Null.NNIs <- list()
trait.Null.NDIs <- list()
trait.Null.Occupancies <- list()
trait.Null.Occupancy.Full <- list()
trait.Null.OccZscores <- list()

# loop over each item in Sub.Mat to calculate all spatial scales
for(scale in 1:length(Sub.Mat)) { 
  
  # return which spatial scale is being worked on
  print(paste("Working on trait null spatial scale ", scale, "...", sep = "")) 
  
  # pull out one spatial scale at a time
  d <- Sub.Mat[[scale]] 
  
  # create empty list to store null mean NN values and mean disparity values
  Null.meanNN.list <- list() 
  Null.meanDisp.list <- list()
  Null.occRate.list <- list()
  
  # for each run of the null
  for(run.num in 1:num.nulls){
    
    # Shuffle the presence/absence matrix
    colnames(d) <- sample(colnames(d))
    
    # sort into alphabetical order
    new_order <- sort(colnames(d))
    d <- d[ , new_order]
    
    # apply Get.Null.Vals function to return quadrat-level metrics
    Null.Vals <- Get.Null.Vals(d) 
    
    # store mean NN values and mean disparity values in their lists
    Null.meanNN.list[[run.num]] <- Null.Vals$mean_NN
    Null.meanDisp.list[[run.num]] <- Null.Vals$mean_disp
    
    # also, calculate null co-occupancy rates and store for later use
    occ.rates <- create.N.matrix(t(d)) # returns a pairwise matrix of number of co-occupancies
    colnames(occ.rates) <- colnames(d) # add column names
    rownames(occ.rates) <- colnames(d) # add row names
    occ.rates <- occ.rates/(matrixStats::count(rowSums(d) > 0)) # convert to rate, with # cells occupied by any member of the group as denom
    
    occ.rates[upper.tri(x = occ.rates, diag = TRUE)] <- NA # convert upper triangle to NA to avoid duplicates
    occ.rates <- reshape2::melt(occ.rates) # transform into a data frame
    occ.rates <- na.omit(occ.rates) # remove NAs
    
    Null.occRate.list[[run.num]] <- occ.rates$value # only save values in list, will add names later (alphabetical!)
  }
  
  # bind mean NN and mean Disp lists into one (does this for a whole spatial scale)
  # rows will be cell IDs, columns are the null run
  Null.meanNN <- matrix(unlist(Null.meanNN.list), ncol = num.nulls)
  Null.meanDisp <- matrix(unlist(Null.meanDisp.list), ncol = num.nulls)
  Null.occRate <- matrix(unlist(Null.occRate.list), ncol = num.nulls)
  
  # summarize spatial null data and organize into a savable format. Metrics summarize all the nulls of a given quadrat
  
  # calculate mean, standard deviation, and percentiles of null mean NN values
  null_meanNN_mean <- rowMeans2(Null.meanNN)
  null_meanNN_sd <- rowSds(Null.meanNN)
  null_meanNN_5pct <- rownth(Null.meanNN, 
                             elems = as.integer(rep(50, dim(Null.meanNN)[1])))
  null_meanNN_95pct <- rownth(Null.meanNN, 
                              elems = as.integer(rep(50, dim(Null.meanNN)[1])), 
                              descending = TRUE)
  null_meanNN_2.5pct <- rownth(Null.meanNN, 
                               elems = as.integer(rep(25, dim(Null.meanNN)[1])))
  null_meanNN_97.5pct <- rownth(Null.meanNN, 
                                elems = as.integer(rep(25, dim(Null.meanNN)[1])), 
                                descending = TRUE) 
  
  # calculate mean, standard deviation, and percentiles of null mean disparity values
  null_meanDisp_mean <- rowMeans2(Null.meanDisp)
  null_meanDisp_sd <- rowSds(Null.meanDisp)
  null_meanDisp_5pct <- rownth(Null.meanDisp, 
                               elems = as.integer(rep(50, dim(Null.meanDisp)[1])))
  null_meanDisp_95pct <- rownth(Null.meanDisp, 
                                elems = as.integer(rep(50, dim(Null.meanDisp)[1])), 
                                descending = TRUE)
  null_meanDisp_2.5pct <- rownth(Null.meanDisp, 
                                 elems = as.integer(rep(25, dim(Null.meanDisp)[1])))
  null_meanDisp_97.5pct <- rownth(Null.meanDisp,
                                  elems = as.integer(rep(25, dim(Null.meanDisp)[1])), 
                                  descending = TRUE)
  
  # calculate mean, standard deviation, and percentiles of null co-ocurrences
  null_occRate_mean <- rowMeans2(Null.occRate) 
  null_occRate_sd <- rowSds(Null.occRate)
  null_occRate_5pct <- rownth(Null.occRate, 
                              elems = as.integer(rep(50, dim(Null.occRate)[1]))) 
  null_occRate_95pct <- rownth(Null.occRate, 
                               elems = as.integer(rep(50, dim(Null.occRate)[1])), 
                               descending = TRUE)
  null_occRate_2.5pct <- rownth(Null.occRate, 
                                elems = as.integer(rep(25, dim(Null.occRate)[1]))) 
  null_occRate_97.5pct <- rownth(Null.occRate,
                                 elems = as.integer(rep(25, dim(Null.occRate)[1])), 
                                 descending = TRUE)
  
  # compare nulls to themselves
  null_NNI_median <- vector()
  null_NDI_median <- vector()
  
  null_NNI_mean <- vector()
  null_NDI_mean <- vector()
  
  null_NNIs <- list()
  null_NDIs <- list()
  
  null_occ_zscores <- list()
  
  for(i in 1:1000){
    # select one column (one run of the null) and compare to the null distribution as a whole, creating a set of 1000 "null" NNI and NDI values and occupancy rate z-scores. For each run of the null, the list will be the number of cells at that scale. Take the median to find a distribution of medial NNI and NDI values and occupancy rate z-scores.
    
    null_NNIs[[i]] <- (Null.meanNN[, i] - null_meanNN_mean) / null_meanNN_sd
    null_NDIs[[i]] <- (Null.meanDisp[, i] - null_meanDisp_mean) / null_meanDisp_sd
    null_NNI_median[i] <- median((Null.meanNN[, i] - null_meanNN_mean) / null_meanNN_sd, 
                                 na.rm = T)
    null_NDI_median[i] <- median((Null.meanDisp[, i] - null_meanDisp_mean) / null_meanDisp_sd, 
                                 na.rm = T)
    null_NNI_mean[i] <- mean((Null.meanNN[, i] - null_meanNN_mean) / null_meanNN_sd, 
                             na.rm = T)
    null_NDI_mean[i] <- mean((Null.meanDisp[, i] - null_meanDisp_mean) / null_meanDisp_sd, 
                             na.rm = T)
    
    null_occ_zscores[[i]] <- (Null.occRate[, i] - null_occRate_mean) / null_occRate_sd
    
  }
  
  # combine all results into one data frame
  results <- cbind(null_meanNN_mean, null_meanNN_sd, 
                   null_meanNN_5pct, null_meanNN_95pct,
                   null_meanNN_2.5pct, null_meanNN_97.5pct, 
                   null_meanDisp_mean, null_meanDisp_sd,
                   null_meanDisp_5pct, null_meanDisp_95pct, 
                   null_meanDisp_2.5pct, null_meanDisp_97.5pct)
  
  null_NNI_NDI_distributions <- cbind(null_NNI_median, null_NDI_median,
                                      null_NNI_mean, null_NDI_mean)
  
  null_NNIs <- do.call(cbind, null_NNIs)
  null_NDIs <- do.call(cbind, null_NDIs)
  
  null_occRates_save <- cbind(occ.rates[, c(1:2)], null_occRate_mean, null_occRate_sd, # also save pair ids for matching
                              null_occRate_5pct, null_occRate_95pct,
                              null_occRate_2.5pct, null_occRate_97.5pct)
  
  null_occRate_distributions <- cbind(occ.rates[, c(1:2)], Null.occRate)
  
  null_occRate_zscores <- do.call(cbind, null_occ_zscores)
  
  trait.Null.Results[[scale]] <- results
  trait.Null.Distributions [[scale]] <- as.data.frame(null_NNI_NDI_distributions)
  trait.Null.NNIs[[scale]] <- as.data.frame(null_NNIs)
  trait.Null.NDIs[[scale]] <- as.data.frame(null_NDIs)
  trait.Null.Occupancies[[scale]] <- null_occRates_save
  trait.Null.Occupancy.Full[[scale]] <- null_occRate_distributions
  trait.Null.OccZscores[[scale]] <- cbind(occ.rates[, c(1:2)], null_occRate_zscores)
}

 # SAVE RESULTS FOR THE trait NULL!!
# names for each item in trait.Null.NNIs
nms <- c("null_NNIs_10m", "null_NNIs_20m", "null_NNIs_50m", "null_NNIs_100m")

names(trait.Null.NNIs) <- nms

# names for each item in trait.Null.NDIs
nms <- c("null_NDIs_10m", "null_NDIs_20m", "null_NDIs_50m", "null_NDIs_100m")

names(trait.Null.NDIs) <- nms

# names for each item in trait.Null.Results
nms <- c("null_stats_10m", "null_stats_20m", "null_stats_50m", "null_stats_100m")

names(trait.Null.Results) <- nms

# names of reach item in trait.Null.Distributions
nms <- c("null_NNI_NDI_dist_10m", "null_NNI_NDI_dist_20m", "null_NNI_NDI_dist_50m",
         "null_NNI_NDI_dist_100m")

names(trait.Null.Distributions) <- nms

# names for each item in trait.Null.Occupancies
nms <- c("null_occ_stats_10m", "null_occ_stats_20m", "null_occ_stats_50m", 
         "null_occ_stats_100m")

names(trait.Null.Occupancies) <- nms

# names for each item in trait.Null.Occupancy.Full
nms <- c("null_occs_10m", "null_occs_20m", "null_occs_50m", "null_occs_100m")

names(trait.Null.Occupancy.Full) <- nms

# names for each item in trait.Null.OccZscores
nms <- c("null_occs_null_Zscores_10m", "null_occs_null_Zscores_20m", 
         "null_occs_null_Zscores_50m", "null_occs_null_Zscores_100m")

names(trait.Null.OccZscores) <- nms

# save trait null results as separate files for each trait scale

# for each item in trait.Null.NNIs
lapply(names(trait.Null.NNIs), function(x) {
  f <- trait.Null.NNIs[[x]]
  write.table(f, file = paste0(getwd(),'/trait-nulls-2010/',
                               s, '/', x, '.txt'),
              sep = "\t", # tab delimited
              row.names = F, # no row names
              col.names = T) # yes column names
})

# for each item in trait.Null.NDIs
lapply(names(trait.Null.NDIs), function(x) {
  f <- trait.Null.NDIs[[x]]
  write.table(f, file = paste0(getwd(),'/trait-nulls-2010/',
                               s, '/', x, '.txt'),
              sep = "\t", # tab delimited
              row.names = F, # no row names
              col.names = T) # yes column names
})

# for each item in trait.Null.Results
lapply(names(trait.Null.Results), function(x) {
  f <- trait.Null.Results[[x]]
  write.table(f, file = paste0(getwd(),'/trait-nulls-2010/',
                               s, '/', x, '.txt'),
              sep = "\t", # tab delimited
              row.names = F, # no row names
              col.names = T) # yes column names
})

# for each item in trait.Null.Distributions
lapply(names(trait.Null.Distributions), function(x) {
  f <- trait.Null.Distributions[[x]]
  write.table(f, file = paste0(getwd(),'/trait-nulls-2010/',
                               s, '/', x, '.txt'),
              sep = "\t", # tab delimited
              row.names = F, # no row names
              col.names = T) # yes column names
})

# for each item in trait.Null.Occupancies
lapply(names(trait.Null.Occupancies), function(x) {
  f <- trait.Null.Occupancies[[x]]
  write.table(f, file = paste0(getwd(),'/trait-nulls-2010/',
                               s, '/', x, '.txt'),
              sep = "\t", # tab delimited
              row.names = F, # no row names
              col.names = T) # yes column names
})

# for each item in trait.Null.Occupancy.Full
lapply(names(trait.Null.Occupancy.Full), function(x) {
  f <- trait.Null.Occupancy.Full[[x]]
  write.table(f, file = paste0(getwd(),'/trait-nulls-2010/',
                               s, '/', x, '.txt'),
              sep = "\t", # tab delimited
              row.names = F, # no row names
              col.names = T) # yes column names
})

# for each item in trait.Null.Occ.Zscores
lapply(names(trait.Null.OccZscores), function(x) {
  f <- trait.Null.OccZscores[[x]]
  write.table(f, file = paste0(getwd(),'/trait-nulls-2010/',
                               s, '/', x, '.txt'),
              sep = "\t", # tab delimited
              row.names = F, # no row names
              col.names = T) # yes column names
})

