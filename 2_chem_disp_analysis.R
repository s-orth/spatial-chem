# TITLE: BCI Chemical Disparity Analysis
# AUTHOR: Sarah Orth, Ostling Lab
# LAST MODIFIED: August 12, 2024

# DESCRIPTION: This code contains functions to subset the BCI presence/absence matrix and to calculate mean chemical disparity and mean nearest neighbor for each quadrat. Application of functions and code to save data is at the end of this script.

# SET-UP ----

# load necessary packages
library(matrixStats) 

# set working directory
setwd("/Users/Sarah/Library/CloudStorage/Box-Box/spatial-chem/")

# set folder containing presence-absence matrices
pres.abs.dir <- "BCI-PA-matrices-2010/" 

# list all file names present in pres.abs.dir folder
list.all <- list.files(path = pres.abs.dir)

# View order of file names
print(list.all)

# put files in order from smallest to largest spatial scale
list.all <- list.all[c(2:4, 1)]

# confirm correct order of files
print(list.all) 

# Read in presence-absence matrices. 
# Quadrat ID columns will have turned into "X1", "X2"... but this is OK. 
# R doesn't like variable names to be a number.
bci.mat <- lapply(X = paste(pres.abs.dir, list.all, 
                            sep = ""), 
                  FUN = read.table, # function to apply to each item in list
                  sep = "\t", # tab-delimited
                  header = T) # header included

# read in chemical disparity matrix
bci.cd <- read.table(file = "clean-BCI-by-quadrat-2010/disparity_matrix_defense_only.txt",
                     # file = "raw-data/kursar_disparity_matrix_all_defense.txt",
                     sep = "\t", # tab-delimited
                     row.names = 1, # row names are in the first column
                     header = T) # header included

# ignore conspecific relationships by changing the number to NA
diag(bci.cd) <- NA

# FUNCTION FOR SUBSETTING THE PRESENCE-ABSENCE DATA ----

# This part of the code selects certain combinations of species. Filtering options will be based on genus.

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

# FUNCTION FOR CALCULATING quadrat-LEVEL CHEMICAL DISPARITY ----

# This function will obtain mean nearest neighbor and mean chemical disparity scores for each quadrat, at each spatial scale (when used in conjunction with lapply()). The required input is a list of presence absence matrices that have species IDs as the column names and quadrat IDs as the rows.

# we will use lapply() to iterate this function over each PA matrix in the list
Get.Quadrat.Vals <- function(d) { 

  # create an empty list to hold species IDs present in each quadrat
  sp.list <- list() 

  # for each row (quadrat) in the matrix
  # create an entry in sp.list containing a list of the species present in the quadrat
  for(row in 1:dim(d)[1]) { sp.list[[row]] <- (which(d[row, ] == 1)) } 
  
  # create an empty list to hold chemical data for the quadrat
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

# APPLYING FUNCTIONS Get.Sub AND Get.Quadrat.Vals ----
# A series of loops will be used to automate calculating for all subsets in one go.

# name of folder to save values in
folder <- "quadrat-values-2010"

# create folder. Will return a warning if folder already exists.
dir.create(paste0(getwd(),"/", folder)) 

# vector listing genera of interest
genera <- c("inga", "miconia", "psychotria")

# for each genus listed in "genera"
for(genus in 1:length(genera)){
  
  # define g as the current genus we are on
  g.current <- genera[genus]
  
  # apply the Sub.Mat function to get current genus of interest
  Sub.Mat <- lapply(X = bci.mat, # apply Sub.Mat to each item in bci.mat
                    FUN = Get.Sub, # use the Get.Sub function
                    g = g.current) # options for Get.Sub function
  
  # grab row and column indices of species in subset in chemistry matrix
  Rinds <- which(rownames(bci.cd) %in% colnames(Sub.Mat[[4]]))
  Cinds <- which(colnames(bci.cd) %in% colnames(Sub.Mat[[4]]))
  
  # filter chemistry matrix to have only species in subset
  sub.cd <- data.matrix(bci.cd[Rinds, Cinds])
  
  # filter Sub.Mat to have only species in the disparity matrix (may need these lines if using a different chemical disparity matrix or census year)
  # inds <- which(colnames(Sub.Mat[[4]]) %in% colnames(bci.cd))
  # Sub.Mat <- Map(function(x){x[ , inds]}, Sub.Mat)
  
  # apply Get.Quadrat.Vals to list of subsetted presence-absence matrices
  Quadrat.Vals <- lapply(X = Sub.Mat, FUN = Get.Quadrat.Vals)
  
  # create a folder to store quadrat values. "folder" created earlier
  dir.create(paste0(getwd(),"/", folder, "/genus-", g.current)) 
  
  # name each item in Quadrat.Vals
  nms <- c(paste0(g.current, "_quad-vals_10m"), 
           paste0(g.current, "_quad-vals_20m"),
           paste0(g.current, "_quad-vals_50m"), 
           paste0(g.current, "_quad-vals_100m"))
  
  # apply names to list for use in saving
  names(Quadrat.Vals) <- nms 
  
  # save quadrat value information as separate files for each spatial scale
  # for each data frame in Quadrat.Vals
  lapply(names(Quadrat.Vals), function(x) { 
    
    # select the data frame
    f <- Quadrat.Vals[[x]]
    
    # write as a table
    write.table(f, file = paste0(getwd(),'/', folder, '/genus-', g.current, '/', x, '.txt'),
                sep = "\t", # tab delimited
                row.names = F, # no row names
                col.names = T) # yes column names
  })
  
}

# apply the Get.Sub function to all species to format the data.
Sub.Mat <- lapply(X = bci.mat, FUN = Get.Sub) 

# assign chemistry matrix to sub.cd (EDIT ME FOR DESIRED DEFENSE MATRIX)
sub.cd <- data.matrix(bci.cd) 

# apply Get.Quad.Vals to list of presence-absence matrices
Quadrat.Vals <- lapply(X = Sub.Mat, FUN = Get.Quadrat.Vals) 

# create a folder to store quadrat values. "folder" created earlier
dir.create(paste0(getwd(), "/", folder, "/all-species")) 

# name each item in Quadrat.Vals
nms <- c("all_species_quad_vals_10m",
         "all_species_quad_vals_20m",
         "all_species_quad_vals_50m",
         "all_species_quad_vals_100m")

# apply names to list for use in saving
names(Quadrat.Vals) <- nms 

# save presence-absence matrices as separate files for each spatial scale
# for each data frame in Quadrat.Vals
lapply(names(Quadrat.Vals), function(x) { 
  
  # select the data frame
  f <- Quadrat.Vals[[x]] 
  
  # write as a table
  write.table(f, file = paste0(getwd(),'/', folder, '/all-species/', x, '.txt'),
              sep = "\t", # tab delimited
              row.names = F, # no row names
              col.names = T) # yes column names
})
