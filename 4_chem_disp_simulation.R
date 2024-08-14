# TITLE: Simulation
# AUTHOR: Sarah Orth, Ostling Lab
# LAST MODIFIED: August 1, 2024
# DESCRIPTION: We suspect that one reason we find differences in our results depending on the null model used is that the trait null is reflecting some kind of pattern in individual species occupancy rates and their disparity with the other species in question. When we see overdispersion in the trait null that goes away with the spatial null, we expect that this is actually reflecting an underlying pattern in high occupancy species also being more different than the other species in their group on average.

# to test this, we will generate data in which there is a known pattern with occupancy and disparity to see what happens if we do the analysis.

library(matrixStats) # load necessary packages. You may need to install these on your system.
library(Rfast)

setwd("/Users/Sarah/Library/CloudStorage/Box-Box/spatial-chem/")

grid.dir <- "BCI-PA-matrices-2010/" # set folder containing presence-absence matrices
list.all <- list.files(path = grid.dir) # list all files names in grid.dir folder
print(list.all)
list.all <- list.all[c(2:4, 1)] # put in order from smallest to largest scale

# Read in presence-absence matrices. tab-delimited file with a header
bci.mat <- lapply(X = paste(grid.dir, list.all, sep = ""), 
                  FUN = read.table, sep = "\t", header = T) 

# read in chemical disparity matrix. tab-delimited file w/ row and column names
bci.cd <- read.table(file = "clean-BCI-by-quadrat-2010/disparity_matrix_defense_only.txt",
                     sep = "\t", row.names = 1, header = T) 

diag(bci.cd) <- NA # ignore conspecific relationships, as this may dilute metrics.

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

# FUNCTION FOR SHUFFLING PRESENCE-ABSENCE MATRIX ----
# modification of the curveball function in the R package "backbone"
# in this simulation, we will use this to simulate species locations on the plot. 
Curveball.Mod <- function(m){
  
  RC <- dim(m) # extract the dimensions (rows and columns) of the matrix m. returns vector of 2
  R <- RC[1] # R is the number of rows in matrix m, first element of RC. returns integer.
  C <- RC[2] # C is the number of columns in matrix m, second element of RC. returns integer.
  hp <- list() # create empty list called hp
  mat.nms <- colnames(m) # save matrix names
  
  # Next part of the code loops down all rows in matrix m (total number of rows comes from dim(m)[1]) for each row (site) in data matrix m, creates a list item in hp indicating which species are in the cell. This is a vector of integers. In their paper, Strona et. al. say that this should be inverted if there are more cells than species - such that there would be a list of which cells a species is present in rather than which species are present in a cell. To change this I believe one would simply iterate through columns instead of rows. See alternative code below. I automated this with an if/else statement that asks if number of rows (cells/sites) is greater than number of columns (species). Mostly there will be more rows than columns (more sites than species)
  
  if(C >= R) { # if there are more species than cells, create lists of which species are in each cell
    for (row in 1:dim(m)[1]) { hp[[row]] = (which(m[row, ] == 1)) }
  } else { # if there are more cells than species, create lists of which cells each species is in
    for (col in 1:dim(m)[2]) { hp[[col]] = (which(m[, col] == 1)) }
  }
  
  l_hp <- length(hp) # returns the length of list hp as integer.
  
  # notes refer to example where cells > species. Invert language for instance of cells < species.
  for (rep in 1:(5*l_hp)) { # for how many repetitions... do what is below
    
    AB <- sample(1:l_hp, 2) # pick two random species presences lists from hp and store integers in AB
    a <- hp[[AB[1]]] # for the first random species (A), list the cells that species is present in
    b <- hp[[AB[2]]] # for the second random species (B), list the cells that species is present in 
    ab <- intersect(a, b) # find the cells that both species A and B are present in and return vector of integers
    l_ab <- length(ab) # number of cells that species A and B have in common (length of ab) - integer
    l_a <- length(a) # number of cells that species A is present in - integer
    l_b = length(b) # number of cells that species B is present in - integer
    
    # As long as the number of cells that A and B have in common is not equal to the number of cells A is in or the number of cells B is in, do the following:
    if ((l_ab %in% c(l_a, l_b)) == F){
      
      tot <- setdiff(c(a,b), ab) # vector of the cells that species A and B are NOT both present in
      l_tot <- length(tot) # number of cells that species A and B are NOT both present in
      tot <- sample(tot, l_tot, replace = FALSE, prob = NULL) # shuffle the cells that are in tot
      L <- l_a - l_ab # how many cells for species A do we need to fill?
      
      # for the first randomly selected presences list, replace in hp a new shuffled list in which the cells that both species were in are all there, and then from the shuffled set of cells that the species did not have in common put in enough to make sure the richness of that cell/occurrence frequency of the species is maintained.
      hp[[AB[1]]] <- c(ab, tot[1:L])
      hp[[AB[2]]] <- c(ab,tot[(L+1):l_tot])} # do the same thing for the second randomly selected presences list but use the things from tot that did not get used in the first list.
  }
  rm <- matrix(0, R, C) # create an empty matrix rm (randomized matrix) with the same dimensions as original matrix m
  
  # for all of the rows in rm, fill in with the new now randomized rows that are stored as vectors in list hp. This also needs to be changed if you had more cells than species. I think I wrote it right!
  if( C >= R ) { for (row in 1:R) { rm[row, hp[[row]]] = 1 }
  } else { for (col in 1:C) { rm[hp[[col]], col] = 1 }}
  
  colnames(rm) <- mat.nms # add names back to matrix
  m <-rm
}

# FUNCTION FOR CALCULATING CELL-LEVEL CHEMICAL DISPARITY ----

# This function will obtain mean nearest neighbor and mean chemical disparity scores for each grid cell. The required input is a list of presence absence matrices that have species IDs as the column names and cell IDs as the rows. 1s indicate that a species is present in a cell, a 0 indicates its absence. For THIS simulation analysis, be sure to use the simulated chemical disparity matrix. The code is set up to use a disparity data table named "simulated.cd"

Get.Cell.Vals <- function(d) { # d is the presence-absence matrix
  
  sp.list <- list() # create an empty list to hold species IDs present in each cell
  
  for(row in 1:dim(d)[1]) { # for each row (cell) in the matrix
    sp.list[[row]] <- (which(d[row, ] == 1))} # create an entry in sp.list containing a list of the species present in the cell
  
  chem.list <- list() # create an empty list to hold chemical data
  
  # ********SLOW********  
  for(cel in 1:length(sp.list)) { # for each cell stored in sp.list
    
    if(length(sp.list[[cel]]) > 1) { # if there are 2 or more species in the cell
      
      cell.chm <- data.matrix(simulated.cd[names(sp.list[[cel]]), names(sp.list[[cel]])]) # filter chem matrix have only species in cell i
      
      nn <- matrixStats::rowMins(cell.chm, na.rm = TRUE) # extract row minimum - this is each species nearest neighbor
      
      mean.nn <- mean(nn) # calculate the mean of the nearest neighbor scores
      
      cell.chm[upper.tri(x = cell.chm, diag = TRUE)] <- NA # make the entire upper triangle of the matrix NA to remove duplicates
      
      mean.cd <- mean(as.matrix(cell.chm), na.rm = T) # calculate the mean chemical disparity in a cell
      
      cell.chm <- c(mean.nn, mean.cd) # store mean NN and mean chemical disparity in cell.chm
      
    } else { cell.chm <- c(NA, NA) } # if there are 0 or 1 species in the cell, store NA values for chemistry in cell.chm
    
    chem.list[[cel]] <- cell.chm # store chem information for cell in chem.list
  } 
  
  d <- as.data.frame(matrix(unlist(chem.list), # turn list of cell information into a data frame for all cells at a spatial scale
                            ncol = 2, nrow = length(chem.list), byrow = T,
                            dimnames = list(NULL, c("mean_NN", "mean_disp")))) # name columns
  
  d$cell_id <- c(1:nrow(d)) # create a cell ID variable
  d <- d[ , c(3, 1:2)] # re-arrange columns
  
} # set up to use simulated.cd

# APPLYING FUNCTION Get.Sub ----
folder <- "/simulation-results" # name of folder to save files
dir.create(paste0(getwd(), folder)) # create folder. Will return a warning if folder already exists.

# we will be using the inga data as a starting point
Sub.Mat <- lapply(X = bci.mat, FUN = Get.Sub, g = "inga") # apply Get.Sub function to each spatial scale in bci.mat

# filtering chemistry matrix to match chosen subset
Rinds <- which(rownames(bci.cd) %in% colnames(Sub.Mat[[4]])) # row indices of species in subset 
Cinds <- which(colnames(bci.cd) %in% colnames(Sub.Mat[[4]])) # column indices of species within subset
sub.cd <- data.matrix(bci.cd[Rinds, Cinds]) # filter chemistry matrix to have only species in subset

# calculate occupancy rates for each species using the presence-absence matrix
occ.rate <- list() # generate an empty list to store values
cell.sizes <- c("m10", "m20", "m50", "m100") # cell sizes to name columns
cell_length <- c(10, 20, 50, 100) # numerical values for cell lengths

for(i in 1:4){ # for each spatial scale...
  d <- Sub.Mat[[i]] # extract one spatial scale
  tot <- matrixStats::count(matrixStats::rowSums2(d, na.rm = TRUE) > 0) # count number of occupied cells
  d <- colSums(d, na.rm = TRUE) # calculate number of occupations for each species
  d <- as.data.frame(d/tot) # calculate relative occupancy
  names(d)[1] <- cell.sizes[i] # add a variable for spatial scale (cell size)
  occ.rate[[i]] <- d # add data to occ.rate list
}

occ <- occ.rate[[1]] # create a dataframe with the first spatial scale to build off of 
for(i in 2:4){ # add the rest of the spatial scales, matching by species identity
  occ <- merge(occ, occ.rate[[i]], by = 'row.names') # merge occupancy rates at all spatial scales
  row.names(occ) <- occ$Row.names # assign row names to full table
  occ <- occ[, -c(1)] # remove extra column
}

# GENERATING RANDOM INGA DATA FOR SIMULATION ####
# For this simulation, I will not be setting up to do multiple spatial scales in one go. The simulations are more time consuming and thus it makes sense to just run the spatial scales of particular interest. Right now, everything is set up to run at the 20m spatial scale

# First, we will be randomizing the chemical relationships by shuffling the internal structure of the pairwise chemical disparity matrix. We will keep diagonals as NA and the resulting matrix will be symmetric across the diagonal. This will vary the chemical relationships of the high occupancy species.

n <- nrow(sub.cd) # number of rows in matrix
num.sims <- 100 # number of times to run the simulation
num.nulls <- 1000 # number of times to run the nulls

sim.mean.disp <- list() # save species mean pairwise disparity to ensure unique sims.
sim.slopes <- vector()

UN.NDI.means <- vector() # initialize vectors to store test results for each simulation
UN.NNI.means <- vector()
UN.NDI.t.pvals <- vector()
UN.NNI.t.pvals <- vector()
UN.NNI.dist.pvals.hi <- vector()
UN.NDI.dist.pvals.hi <- vector()
UN.NNI.dist.pvals.lo <- vector()
UN.NDI.dist.pvals.lo <- vector()

CN.NDI.means <- vector()
CN.NNI.means <- vector()
CN.NDI.t.pvals <- vector()
CN.NNI.t.pvals <- vector()
CN.NDI.dist.pvals.hi <- vector()
CN.NNI.dist.pvals.hi <- vector()
CN.NDI.dist.pvals.lo <- vector()
CN.NNI.dist.pvals.lo <- vector()

for(sim.num in 1:num.sims){ # for each repetition of the simulation
  
  # I - no constraints (no initialized or while loop required)
  # II - -0.001 to 0.001
  # III - 0.5 to 1.5 # these constraints chosen based on original inga slope of 1.12
  # IV - greater than 1.5
  
  print(paste("Working on simulation", sim.num, "..."))
  # initialized value for the while loop
  m <- 0 # CHANGE ME TO m <- 1 for Scenario I and II, m <- 0 for III and IV
  
  # Edit the while loop based on which scenario.
  # For Scenario I: No while loop. Comment out the closing bracket as well.
  # while(m < -0.001 || m > 0.001){ # Scenario II: only accept slope between -0.001 and 0.001
  # while(m < 0.5 || m > 1.5){ # Scenario III: only accept when slope between 0.5 and 1.5
  while(m < 1.5){ # Scenario IV: only accept when slope greater than 1.5
  
  rand.cd <- sub.cd # re-assign chemistry matrix to start fresh
  rand.cd[lower.tri(rand.cd)] <- NA # convert lower triangle to all NA (maintain symmetry)
  
  rand.cd[upper.tri(rand.cd)] <- sample(rand.cd[upper.tri(rand.cd)], # shuffle values
                                        (n*(n-1))/2, replace = FALSE) # w/o replacement
  
  for(i in 1:nrow(rand.cd)){ # copy upper triangle of matrix to the lower triangle (symmetry)
    for(j in 1:i){
      rand.cd[i, j] <- rand.cd[j, i]
    }
  }
  
  avg.disp <- as.data.frame(rowMeans2(rand.cd, na.rm = T)) # calc sp disparity averages
  lm1 <- lm(occ[, 2] ~ avg.disp[, 1]) # basic linear model occupancy rate (at 20m scale) vs avg disparity
  m <- lm1$coefficients[2] # extract slope
  sim.slopes[sim.num] <- m # place to be saved
   } # comment out this bracket if not using the while loop (Scenario I)
  
  simulated.cd <- rand.cd # save the version of the randomized cd that worked here.
  sim.mean.disp[[sim.num]] <- avg.disp # save the species average disparity scores that were simulated
  
  Sim.Mat <- Curveball.Mod(Sub.Mat[[2]]) # shuffle the PA matrix (20m scale) to get simulated combinations of species in their cells.
  
  Sim.Cell.Vals <- Get.Cell.Vals(Sim.Mat) # apply Get.Cell.Vals to presence-absence matrix (uses simulated.cd)
  
  # GENERATING NULL DISTRIBUTIONS FOR THE SIMULATED RESULTS
  # TRAIT (UNCONSTRAINED) NULL
  
  UN.Null.meanNN <- list() # create empty list to store null mean NN values
  UN.Null.meanDisp <- list() # create empty list to store null mean Disparity values
  
  for(run.num in 1:num.nulls){
    
    # Shuffle the chemistry matrix, starting with the simulated one
    rownames(simulated.cd) <- sample(rownames(simulated.cd))
    colnames(simulated.cd) <- rownames(simulated.cd) # simulated.cd is now randomized
    
    UN.Null.Vals <- Get.Cell.Vals(Sim.Mat)
    
    UN.Null.meanNN[[run.num]] <- UN.Null.Vals$mean_NN # store mean NN values in one list
    UN.Null.meanDisp[[run.num]] <- UN.Null.Vals$mean_disp # store mean disparity values in another list
  }
  
  # bind meanNN and meanDisp lists into one. rows will be cell IDs, columns are the null run
  UN.Null.meanNN <- matrix(unlist(UN.Null.meanNN), ncol = num.nulls)
  UN.Null.meanDisp <- matrix(unlist(UN.Null.meanDisp), ncol = num.nulls)
  
  # summarize null data by cell and organize into a save-able format
  # mean NN metric
  UN_meanNN_mean <- rowMeans2(UN.Null.meanNN) # calculate cell mean of null mean NN values
  UN_meanNN_sd <- rowSds(UN.Null.meanNN) # calculate cell standard deviation of null mean NN values
  
  # mean disparity metric
  UN_meanDisp_mean <- rowMeans2(UN.Null.meanDisp) # calculate cell mean of null mean disparity values
  UN_meanDisp_sd <- rowSds(UN.Null.meanDisp) # calculate cell standard deviation
  
  # compare the nulls to themselves to create distribution of null NNI, NDI scores
  null_UN_NNI_means <- vector()
  null_UN_NDI_means <- vector()
  
  for(i in 1:num.nulls){ # calculate each cell's NNI and NDI, then take the mean across all cells

    null_UN_NNI_means[i] <- mean((UN.Null.meanNN[, i] - UN_meanNN_mean) / UN_meanNN_sd,
                             na.rm = T)
    null_UN_NDI_means[i] <- mean((UN.Null.meanDisp[, i] - UN_meanDisp_mean) / UN_meanDisp_sd,
                             na.rm = T)
                        
  } # at the end we will have 1000 null mean NNIs and NDIs to compare our "observed" mean NNI and NDI to

  
  # SPATIAL (CONSTRAINED) NULL  
  CN.Null.meanNN <- list() # create empty list to store null mean NN values
  CN.Null.meanDisp <- list() # create empty list to store null mean Disparity values
  simulated.cd <- rand.cd # re-assign original simulated matrix values
  
  for(run.num in 1:num.nulls){
    
    d <- Curveball.Mod(Sim.Mat) # shuffle presence absence matrix
    CN.Null.Vals <- Get.Cell.Vals(d) # apply Get.Null.Vals function to return cell-level metrics
    
    CN.Null.meanNN[[run.num]] <- CN.Null.Vals$mean_NN # store mean NN values in one list
    CN.Null.meanDisp[[run.num]] <- CN.Null.Vals$mean_disp # store mean disparity values in another list
  }
  
  # bind mean NN and mean Disp lists into one. rows will be cell IDs, columns are the null run
  CN.Null.meanNN <- matrix(unlist(CN.Null.meanNN), ncol = num.nulls)
  CN.Null.meanDisp <- matrix(unlist(CN.Null.meanDisp), ncol = num.nulls)
  
  # summarize null data and organize into a savable format
  CN_meanNN_mean <- rowMeans2(CN.Null.meanNN)
  CN_meanNN_sd <- rowSds(CN.Null.meanNN)
  
  CN_meanDisp_mean <- rowMeans2(CN.Null.meanDisp)
  CN_meanDisp_sd <- rowSds(CN.Null.meanDisp)
  
  # compare the nulls to themselves to create distribution of null NNI, NDI scores
  null_CN_NNI_means <- vector()
  null_CN_NDI_means <- vector()
  
  for(i in 1:num.nulls){
    
    null_CN_NNI_means[i] <- mean((CN.Null.meanNN[, i] - CN_meanNN_mean) / CN_meanNN_sd,
                                 na.rm = T)
    null_CN_NDI_means[i] <- mean((CN.Null.meanDisp[, i] - CN_meanDisp_mean) / CN_meanDisp_sd,
                                 na.rm = T)
    
  }
  
  # calculate comparison metrics - NDI and NNI.
  # Sim.Cell.Vals are the "observed" values for this round of simulation

  UN_NDI <- (Sim.Cell.Vals$mean_disp - UN_meanDisp_mean) / UN_meanDisp_sd
  UN_NDI <- na.omit(UN_NDI)
  UN_NNI <- (Sim.Cell.Vals$mean_NN - UN_meanNN_mean) / UN_meanNN_sd
  UN_NNI <- na.omit(UN_NNI)
  
  CN_NDI <- (Sim.Cell.Vals$mean_disp - CN_meanDisp_mean) / CN_meanDisp_sd
  CN_NDI <- na.omit(CN_NDI)
  CN_NNI <- (Sim.Cell.Vals$mean_NN - CN_meanNN_mean) / CN_meanNN_sd
  CN_NNI <- na.omit(CN_NNI)
  
  # conduct t test on NDI and NNI 
  UN.NDI.test <- t.test(UN_NDI, mu = 0, alternative = "two.sided")
  UN.NNI.test <- t.test(UN_NNI, mu = 0, alternative = "two.sided")
  
  UN.NDI.means[sim.num] <- UN.NDI.test[[5]] # store results 
  UN.NNI.means[sim.num] <- UN.NNI.test[[5]]
  UN.NDI.t.pvals[sim.num] <- UN.NDI.test[[3]]
  UN.NNI.t.pvals[sim.num] <- UN.NNI.test[[3]]
  
  CN.NDI.test <- t.test(CN_NDI, mu = 0, alternative = "two.sided")
  CN.NNI.test <- t.test(CN_NNI, mu = 0, alternative = "two.sided")
  
  CN.NDI.means[sim.num] <- CN.NDI.test[[5]] # store results
  CN.NNI.means[sim.num] <- CN.NNI.test[[5]]
  CN.NDI.t.pvals[sim.num] <- CN.NDI.test[[3]]
  CN.NNI.t.pvals[sim.num] <- CN.NNI.test[[3]]
  
  # compare the "observed" simulated median NNI, NDI to the distribution of null NNI, NDIs

  UN_NNI_mean <- mean(UN_NNI, na.rm = T)
  UN_NNI_p_hi <- sum(null_UN_NNI_means > UN_NNI_mean) / 1000
  UN_NNI_p_lo <- sum(null_UN_NNI_means < UN_NNI_mean) / 1000

  UN_NDI_mean <- mean(UN_NDI, na.rm = T)
  UN_NDI_p_hi <- sum(null_UN_NDI_means > UN_NDI_mean) / 1000
  UN_NDI_p_lo <- sum(null_UN_NDI_means < UN_NDI_mean) / 1000

  CN_NNI_mean <- mean(CN_NNI, na.rm = T)
  CN_NNI_p_hi <- sum(null_CN_NNI_means > CN_NNI_mean) / 1000
  CN_NNI_p_lo <- sum(null_CN_NNI_means < CN_NNI_mean) / 1000

  CN_NDI_mean <- mean(CN_NDI, na.rm = T)
  CN_NDI_p_hi <- sum(null_CN_NDI_means > CN_NDI_mean) / 1000
  CN_NDI_p_lo <- sum(null_CN_NDI_means < CN_NDI_mean) / 1000

  UN.NNI.dist.pvals.hi[sim.num] <- UN_NNI_p_hi
  UN.NDI.dist.pvals.hi[sim.num] <- UN_NDI_p_hi
  CN.NNI.dist.pvals.hi[sim.num] <- CN_NNI_p_hi
  CN.NDI.dist.pvals.hi[sim.num] <- CN_NDI_p_hi
  
  UN.NNI.dist.pvals.lo[sim.num] <- UN_NNI_p_lo
  UN.NDI.dist.pvals.lo[sim.num] <- UN_NDI_p_lo
  CN.NNI.dist.pvals.lo[sim.num] <- CN_NNI_p_lo
  CN.NDI.dist.pvals.lo[sim.num] <- CN_NDI_p_lo
  
}

all.sim.results <- cbind(UN.NNI.means, UN.NNI.t.pvals, UN.NNI.dist.pvals.hi, UN.NNI.dist.pvals.lo,
                         UN.NDI.means, UN.NDI.t.pvals, UN.NDI.dist.pvals.hi, UN.NDI.dist.pvals.lo,
                         CN.NNI.means, CN.NNI.t.pvals, CN.NNI.dist.pvals.hi, CN.NNI.dist.pvals.lo,
                         CN.NDI.means, CN.NDI.t.pvals, CN.NDI.dist.pvals.hi, CN.NDI.dist.pvals.lo)

# SAVING DATA
# suggested naming conventions: "scenario_I", "scenario_II", "scenario_III", "scenario_IV"
folder <- "/scenario-IV" # name of folder to save files, CHANGE FOR SCENARIO
dir.create(paste0(getwd(), "/simulation-results", folder))

ver <- "scenario_IV" # CHANGE FOR SCENARIO. 
set <- "set1.csv" # CHANGE FOR EACH SET OF 100

# save simulation outputs (see what is in all.sim.results)
write.csv(all.sim.results, 
            file = paste0(getwd(), "/simulation-results", folder, "/", 
                          ver, "_results_", set))

sim.mean.disp.save <- as.data.frame(matrix(unlist(sim.mean.disp),
                                           ncol = num.sims, # number of columns in new data frame
                                           nrow = n, # number of rows in new data frame
                                           byrow = F))
# save mean disparities
write.csv(sim.mean.disp.save,
          file = paste0(getwd(), "/simulation-results", folder, "/",
                        ver, "_sim_mean_disp_", set))
# save slopes
write.csv(sim.slopes,
          file = paste0(getwd(), "/simulation-results/", folder, "/",
                        ver, "_slopes_", set))
