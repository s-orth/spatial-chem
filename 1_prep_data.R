# TITLE: Prepare Data
# AUTHOR: Sarah Orth, Ostling Lab
# LAST MODIFIED: August 12, 2024

# DESCRIPTION: 
# This script prepares data for use in analyses of the spatial chemistry project.
# LOAD PACKAGES ----

# load necessary packages
library(stringr)
library(data.table)
library(tidyr)

# CLEAN BCI CENSUS DATA ----

# set working directory to read in raw data
setwd("/Users/Sarah/Library/CloudStorage/Box-Box/spatial-chem/raw-data/") # directory name

# load BCI census data (full measurement table), repeat this code with different census years as needed.
bci <- read.table(file = "bci_census_2010.txt", # file name (change if using a different census year)
                  sep = "\t", # tab-delimited
                  quote = "", # no quoting
                  h = T) # header exists

# remove duplicate entries. 'match' returns a vector of positions
bci <- bci[match(x = unique(bci[ ,8]), # unique TreeID values
                 table = bci$TreeID), ]

# read in taxonomy information to help in matching any species in chemical data that are "missing" from census
taxa <- read.table(file = "bci_taxa.txt",
                   sep = "\t", quote = "", h = T)

# load Wright BCI trait data
traits <- read.csv(file = "wright_BCI_traits_20101220.csv", # file name
                   header = T) # header exists

traits <- traits[ ,c(7, 28)] # keep species ID and DBH_AVG, which is the avg max DBH of 6 largest,
names(traits)[1] <- "Mnemonic" # rename species ID column to match BCI census name
traits$Mnemonic <- str_to_lower(traits$Mnemonic) # make species ID lowercase
traits$adult_DBH <- traits$DBH_AVG / 2 # calculate estimated DBH for an individual to be an adult

# CHEMICAL DATA ----

# load CSCS data from Brian Sedio (contains data on 133 species)
load("sedio_BCI_CSCS_defense_20210630.RData") # defense only data

# remove individuals from census that are not part of metabolomics data
inds <- which(bci$Mnemonic %in% rownames(cscs)) # indicate rows containing species to keep
bci.filtered <- bci[inds, c(2:3, 6:8, 12:13)] # keep only indicated rows and keep only relevant columns

# add a column of the 1/2 * DBH_AVG
bci.filtered <- merge(x = bci.filtered, 
                      y = traits, 
                      by = "Mnemonic", 
                      all.x = T) # keep all obs. from bci even if missing trait values

# identify which species in the CSCS matrix did not show up in the BCI census data and examine
missing.tax <- setdiff(rownames(cscs), bci.filtered$Mnemonic)

# Notes on the missing taxa: 
  # "clidse" falls out of 2010 census
  # "hieral" appears to be under "hyeral" in the census data codes for 2005 but we will ignore it 
  # "hybpar" is a typo for "hybapr" and can be ignored
  # "nectsp" is a different mnemonic for "nects1", but "nects1" occurs only twice and therefore can be ignored
  # "ocotpe", "pipeda", "pipeim", and "psycfu" perhaps were on the island but did not recruit into the census in 2005
  # "pyscg2" is a typo for "psycg2" and can be ignored

# CLEAN CSCS MATRIX ----

Rinds <- which(rownames(cscs) %in% bci.filtered$Mnemonic) # which species are in BCI rows
Cinds <- which(colnames(cscs) %in% bci.filtered$Mnemonic) # which species are in BCI columns
cscs.filtered <- cscs[Rinds,Cinds] # include only species in BCI - remove those that are not

# convert CSCS scores to "chemical disparity" scores by subtracting from 1
bci.cd <- 1 - cscs.filtered

# to look at adults or juveniles only, edit and run this line. Make sure to update names for saving to reflect!!
# bci.filtered <- subset(bci.filtered, bci.filtered$DBH < 100) # flip sign for adults/juveniles
# bci.filtered <- subset(bci.filtered, bci.filtered$DBH > bci.filtered$adult_DBH) # flip sign for adults/juveniles

# ASSIGN GRIDS ----

# READ ME: The code below will assign a quadrat identity for each stem in the census at specified quadrat sizes. Some parts of the code may need to be modified if you are changing which quadrat sizes (length) you are using: the list "quadrats" and the list "nms". Analyses will use this information to compare the chemical diversity among species co-occurring within a quadrat. Additionally, if using a different census year, you may wish to change the naming convention used in "nms".

# specify quadrat sizes
quadrats <- c(10, 20, 50, 100)

# set up an empty list to store data with quadrat assignments
all.list <- list()

# for each scale specified in 'quadrats'...
for(i in 1:length(quadrats)){
  
  # select quadrat size i to work with
  qdr.size <- quadrats[[i]] 
  
  # Set up a grid based on current quadrat size
  # generate locations of quadrat corners on x axis
  x.coord <- seq(0, # starting point of x axis
                 1000, # end point of x axis (length of BCI plot)
                 by = qdr.size) # count by value in qdr.size
  
  # generate locations of quadrat corners on y axis
  y.coord <- seq(0, # starting point of y axis
                 500, # end point of y axis (width of BCI plot)
                 by = qdr.size) # count by value in qdr.size
  
  # initialize counter to be used later for quadrat id
  counter <- 0
  
  # calculate number of quadrats at spatial scale i
  num.quadrats <- (length(x.coord) - 1) * (length(y.coord) - 1)
  
  # create empty list to hold information for individual quadrats
  quadrat.list <- list()
  
  # for each x-coordinate, and for each y-coordinate...
  for(xi in 1:(length(x.coord) - 1)){ 
    for(yi in 1:(length(y.coord) - 1)){
      
      # Get (x, y) limits for the quadrat at (xi, yi)
      xmin <- x.coord[xi]
      xmax <- x.coord[xi + 1]
      ymin <- y.coord[yi]
      ymax <- y.coord[yi + 1]
      
      # Select individuals in the quadrat located at (xi, yi) from the clean census data
      my.quadrat <- subset(bci.filtered,
                           bci.filtered$PX >= xmin & # individuals within xmin 
                           bci.filtered$PX < xmax & # and xmax
                           bci.filtered$PY >= ymin & # individuals within ymin 
                           bci.filtered$PY < ymax) # and ymax
      
      # quadrat identifier
      counter <- counter + 1
      
      # calculate quadrat species richness
      unique.spp <- length(unique(my.quadrat$Mnemonic))
      
      # for occupied quadrats (those with at least one species in them), record the following:
      if(unique.spp >= 1){
        
        # create a table of the # of stems (abundance) for each species in the quadrat
        abund <- as.data.frame(table(my.quadrat$Mnemonic)) 
        
        names(abund)[1] <- "Mnemonic" # rename column with species IDs
        names(abund)[2] <- "spp_abund_in_quadrat" # rename column with quadrat abundances
        
        # add quadrat-level species abudances to my.quadrat by merging by mnemonic
        my.quadrat <- merge(x = my.quadrat, 
                            y = abund, 
                            by = "Mnemonic") 
        
        my.quadrat$xmin <- xmin # store xmin coordinate
        my.quadrat$xmax <- xmax # store xmax coordinate
        my.quadrat$ymin <- ymin # store ymin coordinate     
        my.quadrat$ymax <- ymax # store ymax coordinate
        my.quadrat$quadrat_length <- qdr.size # store length of quadrat
        my.quadrat$quadrat_size <- qdr.size^2 # store quadrat area
        my.quadrat$quadrat_locX <- mean(c(xmin, xmax)) # store quadrat x center point
        my.quadrat$quadrat_locY <- mean(c(ymin, ymax)) # store quadrat y center point
        my.quadrat$quadrat_id <- counter # store quadrat id
        my.quadrat$quadrat_richness <- unique.spp # store quadrat species richness
        
        # store full quadrat information in the list
        quadrat.list[[counter]] <- my.quadrat
        
        
      # if the quadrat is unoccupied, store location but otherwise fill with NAs
      } else { 
        
        my.quadrat[nrow(my.quadrat) + 1, ] <- NA # create a row to store info in
        my.quadrat$spp_abund_in_quadrat <- NA # no species abundance to record
        my.quadrat$xmin <- xmin # store xmin coordinate
        my.quadrat$xmax <- xmax # store xmax coordinate
        my.quadrat$ymin <- ymin # store ymin coordinate     
        my.quadrat$ymax <- ymax # store ymax coordinate
        my.quadrat$quadrat_length <- qdr.size # store length of quadrat
        my.quadrat$quadrat_size <- qdr.size^2 # store quadrat area
        my.quadrat$quadrat_locX <- mean(c(xmin, xmax)) # store quadrat x center point
        my.quadrat$quadrat_locY <- mean(c(ymin, ymax)) # store quadrat y center point
        my.quadrat$quadrat_id <- counter # store quadrat id
        my.quadrat$quadrat_richness <- unique.spp # store richness (0 for empty quadrat)
        
        # store quadrat information in repository for quadrat information
        quadrat.list[[counter]] <- my.quadrat 
        
      }
    }
  }
  
  # bind information together to have all quadrats at a spatial scale
  all <- do.call(rbind, quadrat.list)
  
  # store each spatial scale in the list 'all.list'
  all.list[[i]] <- all 
  
}

# CREATE PRESENCE-ABSENCE MATRICES ----

# Below is a function, Pres.Abs, that will take the census data with quadrat assignments at each spatial scale and convert it into a species presence-absence matrix (0 indicates absence, 1 indicates presence). 
# The matrix will include species ID and other identifiers that can be used in sub-setting the data for various types of analysis later.

# this function will be used in conjunction with lapply() to be iterated over a list of data frames. "d" is a placeholder for a data frame.
Pres.Abs <- function(d){ 
  
  # silently set d to have data frame structure
  setDF(d) 
  
  # remove variables that uniquely identify a stem and quadrat location information
  d <- d[, -c(3:18)] 
  
  # removing duplicates lists each species that is present in a quadrat.
  # this allows us to track occupancy rate rather than abundance
  d <- unique(d) 
  
  # each quadrat id becomes its own variable, with the value being the quadrat richness
  d <- spread(d, quadrat_id, quadrat_richness) 
  
  # isolate the quadrat id variables from the data matrix d
  quadrat.ids <- d[ , c(3:length(d))]
  
  # if the values in the quadrat id columns are greater than 0, convert them to 1s (presence)
  quadrat.ids[quadrat.ids > 0] <- 1 
  
  # if the values in the quadrat id columns are NA, convert them to Os (absence)
  quadrat.ids[is.na(quadrat.ids)] <- 0 
  
  # rejoin the quadrat id columns with the species information from d and store back in d
  d <- cbind(d[, c(1:2)], quadrat.ids)
  
  # identify row of NA species ID (an artifact of storing empty quadrats). 
  # Returns TRUE for non-NA species
  inds <- !is.na(d$Mnemonic) 
  
  # remove the row with NA for species ID
  d <- d[inds, ] 
  
  # make the species ID the row name. Keep other identifiers for now.
  rownames(d) <- d$Mnemonic
  
  # return presence-absence data frame
  d 
}

# apply the Pres.Abs function to each cleaned census in all.list.
PresAbs.Matrix <- lapply(X = all.list, FUN = Pres.Abs) 


# SAVING CLEANED DATA ----

# names for each list item in all.list. Edit to use whatever naming convention you'd like for saving the files.
nms <- c("clean_BCI_census_10m", 
         "clean_BCI_census_20m",
         "clean_BCI_census_50m", 
         "clean_BCI_census_100m")

# add names to list for use in saving
names(all.list) <- nms

# set wd for saving data
setwd("/Users/Sarah/Library/CloudStorage/Box-Box/spatial-chem/")

# name of folder to create to save data in
folder <- "/clean-BCI-by-quadrat-2010"

# create the folder to save data in. Returns a warning if folder already exists.
dir.create(paste0(getwd(), folder))

# save cleaned and quadrat-assigned census data as separate files for each spatial scale

# for each data frame in all.list
lapply(names(all.list), function(x) { 
  f <- all.list[[x]] # select the data frame, then write as a table
  write.table(f, file = paste0(getwd(), # save location
                               folder, "/", x, '.txt'), 
              sep = "\t", # tab-delimited
              row.names = F, # no row names
              col.names = T) # yes column names
})

# save chemical disparity matrix as a table
write.table(bci.cd, 
            file = paste0(getwd(), # save location
                          folder, "/", 'disparity_matrix_defense_only.txt' ),
            sep = "\t", # tab-delimited
            row.names = T, # yes row names
            col.names = T) # yes column names


# names for each item in PresAbs.Matrix
nms <- c("BCI_PA_10m", 
         "BCI_PA_20m",
         "BCI_PA_50m", 
         "BCI_PA_100m")

# apply names to PresAbs.Matrix list for use in saving
names(PresAbs.Matrix) <- nms 

# name of folder to save presence-absence matries in
folder <- "/BCI-PA-matrices-2010" 

# create folder. Will return a warning if folder already exists.
dir.create(paste0(getwd(), folder)) 

# save presence-absence matrices as separate files for each spatial scale
# for each data frame in PresAbs.Matrix
lapply(names(PresAbs.Matrix), function(x) { 
  f <- PresAbs.Matrix[[x]] # select the data frame then write as a table
  write.table(f, file = paste0(getwd(), # save location
                               folder, "/", x, '.txt'),
              sep = "\t", # tab delimited
              row.names = T, # yes row names
              col.names = T) # yes column names
})


