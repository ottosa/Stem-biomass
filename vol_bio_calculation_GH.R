
#################################################################################################
### This code:                                                                                ### 
### Calculating stemcurves, volumesa and bimass for tree stem measurements generated whit PCT ###
#################################################################################################

# PCT = Point-Cloud-Tools developed by Tuomas Yrttimaa, School of Forest Sciences, University of Eastern Finland
# "Automatic point cloud processing tools for Matlab to characterize trees from terrestrial laser scanning point clouds"
# https://github.com/tuomasyr/Point-Cloud-Tools/tree/main

# Workflow
# 1. Generate stemcurves with smooth spline for each successfully segmented tree
# 2. Generate stemcurves for manually segmented and fixed trees
# 3. Calculate volume with cylinder fitting
# 4. Calculate stem biomass 

# Empty environment 
rm(list = ls())

# libraries
library(stringr)
library(lme4)
library(tidyverse)

### 1. Stemcurves for successfully segmented trees ###

# Directory of the segmented CCF and RF plots
dir_CCF <- "F:/Tree_segments/CCF/segment_plots/"
dir_RF <- "F:/Tree_segments/RF/"
folders_CCF <- list.files(dir_CCF, all.files = F) # list all plot folders in the directory
folders_RF <- list.files(dir_RF, all.files = F) # list all plot folders in the directory
dir_manual <- "F:/Tree_segments/CCF/error_trees/segment_results/" # Directory where manualy segmented trees are 

dir <- "F:/Tree_segments/" #Directory to save the final data

# List of all the CCF and RF trees. Containing plot name, tree name, height and diameter of the trees 
# Note! There is missing diameters for some of the trees which are added later
DataTrees_CCF <- read.csv("F:/Tree_segments/CCF/TreeData/trees_list.csv",sep=";")
densTree_CCF <- read.table('F:/Tree_segments/Article2_RFiles/data/dens_mats/CCF_mat.txt', sep = " ", header = TRUE) # CCF density matrix 

DataTrees_RF_sample <- read.csv("F:/Tree_segments/RFSample.csv",sep=",")
DataTrees_RF <- read.csv("F:/Tree_segments/RF_trees.csv")
DataTrees_RF <- distinct(DataTrees_RF)

densTree_RF <- read.table('F:/Tree_segments/Article2_RFiles/data/dens_mats/RF_mat.txt', sep = " ", header = TRUE) # RF density matrix
RF_h <- read.csv("F:/Tree_segments/RF_trees.csv",sep=",") # Data frame with corrected heights for RF trees 

# Drop out OMT forest types from RF data,as we are interested in MT forest types (plots = 2, 5 and 86)
DataTrees_RF <- DataTrees_RF[!(DataTrees_RF$plot %in% "2_1" |
                                 DataTrees_RF$plot %in% "2_2" | 
                                 DataTrees_RF$plot %in% "2_3" |
                                 DataTrees_RF$plot %in% "2_4" |
                                 DataTrees_RF$plot %in% "2_5" |
                                 DataTrees_RF$plot %in% "86_1" |
                                 DataTrees_RF$plot %in% "86_2" |
                                 DataTrees_RF$plot %in% "5_1" |
                                 DataTrees_RF$plot %in% "5_2" |
                                 DataTrees_RF$plot %in% "5_3" |
                                 DataTrees_RF$plot %in% "5_4" |
                                 DataTrees_RF$plot %in% "5_5"),]

################################
### Stemcurves for CCF trees ###
################################

# Data frame to store the smooth spline stemcurves 
trees_CCF <- data.frame()

# Nested loop, looping through plots and trees inside each plot
# p = plots (first loop) and n = trees (second loop)
for (p in levels(as.factor(folders_CCF))) {
  files <- list.files(paste0(dir_CCF,"/",p,"/stemcurves_raw_olrm/"), pattern = ".txt") # stem measurements from plot p
  treeData <- DataTrees_CCF[DataTrees_CCF$plot == p,] # Subsetting trees from tree list
  treesNum <- as.character(treeData$tree) # Taking tree numbers as a vector for the nested loop
  
  for (n in levels(as.factor(treesNum))) {
    
    # If tree n is in files list smooth spline is fitted to stem measurements
    if("TRUE"%in%str_detect(files, paste0("tree",n,"_")) == TRUE){ # detect if tree n is in stemcurve list (files)
      file <- files[grep(paste0("tree",n,"_"), files, fixed = TRUE)] # tree n detected as file
      tree <- read.delim(paste0(dir_CCF,p,"/stemcurves_raw_olrm/",file)) # read tree n measurements
      
      #treeFilter <- subset(tree, tree$d > low & tree$d < up) # Filter outliers
      treeFilter <- tree
      # Getting tree height from tree list and adding to the tree data
      h <- c(h = treeData$h[treeData$tree==as.numeric(n)], d = 0)
      treeFilter <- rbind(treeFilter, h)
      
      # if tree is more than 10 measurements from the stem is smooth spline fitted 
      if(length(treeFilter$d) > 10){
        # fitting smooth spline to the tree data. Note! lambada set to 0.005, so that the spline is more robust
        spl <- smooth.spline(treeFilter$h, treeFilter$d, lambda = 0.005) 
        # Data frame to store the predicted stem diameter values 
        treeStem <- data.frame(h = seq(from = 0, to = max(treeFilter$h), by = 0.1), 
                               plot = p,
                               tree = as.numeric(n))
        treeStem <- treeStem[, c('plot', 'tree', 'h')]
        # Predicting the stem diameter values at 0.1m intervals
        treeStem$d <- predict(spl, treeStem$h)$y
        treeStem$d[treeStem$h == max(treeStem$h)] <- 0 # diameter at max height set to 0

        # Bind the data to tree data.frame()
        trees_CCF <- rbind(trees_CCF, treeStem)
        
        # If the if statements are not TRUE = iterate the next tree
      }else{
        next
      }
    }else{
      next
    }
  }
}


###################################
### Stemcurves for ALL RF trees ###
###################################

trees_RF_ALL <- data.frame()

for (p in 1:length(DataTrees_RF$plot)) {
  plot <- DataTrees_RF[p,]$plot
  id <- DataTrees_RF[p,]$tree
  plot_id <- paste0(plot, "_", id)
  file_olrm <- paste0(dir_RF,plot,"/stemcurves_raw_olrm/","tree",id,"_stemcurve_raw_olrm.txt") # stem measurements from plot p
  file <- paste0(dir_RF,plot,"/stemcurves/","tree",id,"_stemcurve.txt") # stem measurements from plot p
  if(file.exists(file)){
    tree_olrm <- read.delim(file_olrm)
    tree <- read.delim(file)
  }else{
    next
  }
  
  DataTrees_RF[DataTrees_RF$plot == plot & DataTrees_RF$tree == id,]
  # Removing outliers using robust Standard deviation method
  treeFilter <- tree_olrm
  # Getting tree height from tree list and adding to the tree data
  if (is.na(DataTrees_RF$height[DataTrees_RF$plot == plot & DataTrees_RF$tree == id])){
    h <- c(h = max(tree$h), 
           d = 0)
    DataTrees_RF$height[DataTrees_RF$plot == plot & DataTrees_RF$tree == id] <- h
  }else{
    h <- c(h = DataTrees_RF$height[DataTrees_RF$plot == plot & DataTrees_RF$tree == id],
           d = 0)
  }
  
  treeFilter <- rbind(treeFilter, h)
  
  # if tree has more than 10 measurements from the stem is smooth spline fitted
  if(length(treeFilter$d) > 10){
    # fitting smooth spline to the tree data. Note! lambada set to 0.005, so that the spline is more robust
    spl <- smooth.spline(treeFilter$h, treeFilter$d, lambda = 0.005) 
    # Data frame to store the predicted stem diameter values 
    treeStem <- data.frame(h = seq(from = 0, to = max(treeFilter$h), by = 0.1), 
                           plot = plot,
                           tree = id)
    treeStem <- treeStem[, c('plot', 'tree', 'h')]
    # Predicting the stem diameter values at 0.1m intervals
    treeStem$d <- predict(spl, treeStem$h)$y
    treeStem$d[treeStem$h == max(treeStem$h)] <- 0 # diameter at max height set to 0
    
    # Bind the data to tree data.frame()
    trees_RF_ALL <- rbind(trees_RF_ALL, treeStem)
    
    # If the if statements are not TRUE = iterate the next tree
  }else{
    next
  }
}


#################################################################
### 2. Stemcurves for fixed and manually segmented  CCF trees ###
#################################################################

# Basically same thing than previous loop, but data is taken from other dir 

folders <- list.files(dir_manual, all.files = F)

for (p in levels(as.factor(folders))) {
  files <- list.files(paste0(dir_manual,"/",p,"/stemcurves_raw_olrm/"), pattern = ".txt")
  treeData <- DataTrees_CCF[DataTrees_CCF$plot == p,]
  treesNum <- as.character(treeData$tree)
  
  for (n in levels(as.factor(treesNum))) {
    
    if("TRUE"%in%str_detect(files, paste0("tree",n,"_")) == TRUE){
      file <- files[grep(paste0("tree",n,"_"), files, fixed = TRUE)]
      tree <- read.delim(paste0(dir_manual,p,"/stemcurves_raw_olrm/",file))
      treeFilter <- tree
      
      h <- c(h = treeData$h[treeData$tree==as.numeric(n)], d = 0)
      treeFilter <- rbind(treeFilter, h)
      
      if(length(treeFilter$d) > 10){
        spl <- smooth.spline(treeFilter$h, treeFilter$d, lambda = 0.005)
        
        treeStem <- data.frame(h = seq(from = 0, to = max(treeFilter$h), by = 0.1), 
                               plot = p,
                               tree = as.numeric(n))
        treeStem <- treeStem[, c('plot', 'tree', 'h')]
        
        treeStem$d <- predict(spl, treeStem$h)$y
        treeStem$d[treeStem$h == max(treeStem$h)] <- 0
        
        ### NOTE! if tree data is already in tree data.frame() previously calculated diameters are deleted and new data is stored 
        if(unique(trees_CCF$plot[trees_CCF$plot == p]) == unique(treeStem$plot) & 
           unique(trees_CCF$tree[trees_CCF$tree == unique(treeStem$tree)]) == unique(treeStem$tree)){
          trees_CCF <- trees_CCF[!(trees_CCF$plot == p & trees_CCF$tree == as.numeric(n)),]
        }
        
        trees_CCF <- rbind(trees_CCF, treeStem)
      }else{
        next
      }
    }else{
      next
    }
  }
}

############################################
### 3. Huber formula to calculate volume ###
############################################

# Create column with plot and tree id in same column
trees_CCF$treePlot <- paste0(trees_CCF$plot,"_",trees_CCF$tree)

##################
### Volume CCF ###
##################

treeCyls_CCF <- data.frame()

for (i in levels(factor(trees_CCF$treePlot))){
  # subset tree by id
  stem <- trees_CCF[trees_CCF$treePlot == i,]
  stem$relH <- stem$h/max(stem$h)
  vol <- 0 # value to store the tree volume 
  
  # loop through each height starting from iteration 2 
  for (j in 2:length(stem$h)) {
    cyl_h <- stem$h[j] - stem$h[j-1] # cylinder height (meters)
    cyl_r <- mean(c(stem$d[j], stem$d[j-1]))/2 # cylinder radius (meters)
    cyl_vol <- pi*cyl_r^2*cyl_h # calculate cylinder volume (m^3)
    vol <- vol + cyl_vol # add cylinder volume to tree whole volume
    # Store all the data to data.frame
    iteration <- stem[j,]
    iteration$cyl_h <- cyl_h
    iteration$cyl_r <- cyl_r
    iteration$cyl_vol <- cyl_vol
    treeCyls_CCF <- rbind(treeCyls_CCF,iteration)
  }
  DataTrees_CCF$vol[DataTrees_CCF$plot == unique(stem$plot) 
                & DataTrees_CCF$tree == unique(stem$tree)] <- vol # store the stem volume to DataTrees data frame
  DataTrees_CCF$dbh_splin[DataTrees_CCF$plot == unique(stem$plot) 
                      & DataTrees_CCF$tree == unique(stem$tree)] <- stem$d[stem$h == 1.3] * 100 # store diameter at 1.3m to DataTrees data frame
}

DataTrees_CCF <- DataTrees_CCF[!is.na(DataTrees_CCF$dbh_splin),]
DataTrees_CCF$vol_Laasis <- 0.022927*(DataTrees_CCF$dbh_splin^1.99146)*((0.98600)^DataTrees_CCF$dbh_splin)*(DataTrees_CCF$h^2.82541)*((DataTrees_CCF$h-1.3)^-1.53547)

DataTrees_CCF$vol_Laasis <- DataTrees_CCF$vol_Laasis/1000
plot(DataTrees_CCF$vol, DataTrees_CCF$vol_Laasis)
abline(0,1, col = "red")


###########################
### Volume RF ALL Trees ###
###########################

# Create column with plot and tree id in same column
trees_RF_ALL$treePlot <- paste0(trees_RF_ALL$plot,"_",trees_RF_ALL$tree)
unique(trees_RF_ALL$plot)

treeCyls_RF <- data.frame()

for (i in levels(factor(trees_RF_ALL$treePlot))){
  # subset tree by id
  stem <- trees_RF_ALL[trees_RF_ALL$treePlot == i,]
  stem$relH <- stem$h/max(stem$h)
  vol <- 0 # value to store the tree volume 
  
  # loop through each height starting from iteration 2 
  for (j in 2:length(stem$h)) {
    cyl_h <- stem$h[j] - stem$h[j-1] # cylinder height (meters)
    cyl_r <- mean(c(stem$d[j], stem$d[j-1]))/2 # cylinder radius (meters)
    cyl_vol <- pi*cyl_r^2*cyl_h # calculate cylinder volume (m^3)
    vol <- vol + cyl_vol # add cylinder volume to tree whole volume
    # Store all the data to data.frame
    iteration <- stem[j,]
    iteration$cyl_h <- cyl_h
    iteration$cyl_r <- cyl_r
    iteration$cyl_vol <- cyl_vol
    treeCyls_RF <- rbind(treeCyls_RF,iteration)
  }
  DataTrees_RF$vol[DataTrees_RF$plot == unique(stem$plot) 
                    & DataTrees_RF$tree == unique(stem$tree)] <- vol # store the stem volume to DataTrees data frame
  DataTrees_RF$dbh_splin[DataTrees_RF$plot == unique(stem$plot) 
                          & DataTrees_RF$tree == unique(stem$tree)] <- stem$d[stem$h == 1.3] * 100 # store diameter at 1.3m to DataTrees data frame
}

DataTrees_RF <- DataTrees_RF[!is.na(DataTrees_RF$dbh_splin),]
DataTrees_RF$vol_Laasis <- 0.022927*(DataTrees_RF$dbh_splin^1.99146)*((0.98600)^DataTrees_RF$dbh_splin)*(DataTrees_RF$h^2.82541)*((DataTrees_RF$h-1.3)^-1.53547)

DataTrees_RF$vol_Laasis <- DataTrees_RF$vol_Laasis/1000
plot(DataTrees_RF$vol, DataTrees_RF$vol_Laasis)
abline(0,1, col = "red")

##########################################################################################################################################################

#############################
### 4.Biomass calculation ###
#############################

# Treatment based biomass

colnames(densTree_CCF) <- gsub('X','', colnames(densTree_CCF))
densTree_CCF <- as.matrix(densTree_CCF)

biomass_table_CCF_treat <- data.frame() 

for (i in levels(as.factor(treeCyls_CCF$treePlot))) {
  stem <- treeCyls_CCF[treeCyls_CCF$treePlot == i,] # stem = individual tree containing volume of each cyl
  # Looping each cylinder 
  for (j in levels(as.factor(stem$h))) {
    iteration <- stem[stem$h == j,] # taking in each cyl individualy
    r_cyl <- iteration$cyl_r*100 # calculate cyl dbh (cm) 
    r_split <- seq(from = 1, to = r_cyl+1, by = 1) # split 
    k <- max(r_split)
    bio <- 0
    while(k > 0){
      if (k == max(r_split)){
        d <- r_cyl
      }
      else{
        d <- k
      }
      if(r_split[k] > max(as.numeric(colnames(densTree_CCF)))){
        c <- max(as.numeric(colnames(densTree_CCF)))
      } else{
        c <- r_split[k]
      }
      if(iteration$h > max(as.numeric(row.names(densTree_CCF)))){
        r <- max(as.numeric(row.names(densTree_CCF)))
      } else{
        r <- iteration$h
      } 
      dens <- densTree_CCF[as.character(r),as.character(c)]
      
      vol <- (pi*d^2*(iteration$cyl_h*100)) - (pi*(k-1)^2*(iteration$cyl_h*100))
      bio_part <- dens * vol 
      bio <- bio + bio_part
      
      k <- k-1
    }
    iteration$biomass_treat <- bio/1000
    biomass_table_CCF_treat <- rbind(biomass_table_CCF_treat,iteration)
  }
}

for (i in levels(as.factor(biomass_table_CCF_treat$treePlot))) {
  tree <- biomass_table_CCF_treat[biomass_table_CCF_treat$treePlot == i,]
  bio <- sum(tree$biomass)
  DataTrees_CCF$biomass_treat[DataTrees_CCF$plot == unique(tree$plot) & DataTrees_CCF$tree == unique(tree$tree)] <- bio 
}

##################
### Biomass RF ###
##################

colnames(densTree_RF) <- gsub('X','', colnames(densTree_RF))
densTree_RF <- as.matrix(densTree_RF)

biomass_table_RF_treat <- data.frame() 

treeCyls_RF$d[treeCyls_RF$d < 0] <- 0
treeCyls_RF$cyl_r[treeCyls_RF$cyl_r < 0] <- 0

for (i in levels(as.factor(treeCyls_RF$treePlot))) {
  stem <- treeCyls_RF[treeCyls_RF$treePlot == i,]
  for (j in levels(as.factor(stem$h))) {
    iteration <- stem[stem$h == j,]
    r_cyl <- iteration$cyl_r*100
    r_split <- seq(from = 1, to = r_cyl+1, by = 1)
    k <- max(r_split)
    bio <- 0
    while(k > 0){
      if (k == max(r_split)){
        d <- r_cyl
      }else{
        d <- k
      }
      if(r_split[k] > max(as.numeric(colnames(densTree_RF)))){
        c <- max(as.numeric(colnames(densTree_RF)))
      } else{
        c <- r_split[k]
      }
      if(iteration$h > max(as.numeric(row.names(densTree_RF)))){
        r <- max(as.numeric(row.names(densTree_RF)))
      } else{
        r <- iteration$h
      } 
      dens <- densTree_RF[as.character(r),as.character(c)]
      
      vol <- (pi*d^2*(iteration$cyl_h*100)) - (pi*(k-1)^2*(iteration$cyl_h*100))
      bio_part <- dens * vol 
      bio <- bio + bio_part
      
      k <- k-1
    }
    iteration$biomass <- bio/1000
    biomass_table_RF_treat <- rbind(biomass_table_RF_treat,iteration)
  }
}

for (i in levels(as.factor(biomass_table_RF_treat$treePlot))) {
  tree <- biomass_table_RF_treat[biomass_table_RF_treat$treePlot == i,]
  bio <- sum(tree$biomass)
  DataTrees_RF$biomass_treat[DataTrees_RF$plot == unique(tree$plot) & DataTrees_RF$tree == unique(tree$tree)] <- bio 
}


# Save data 
write.csv(DataTrees_CCF, paste0(dir,"TreeDataCCF.csv"))
write.csv(DataTrees_RF, paste0(dir,"TreeDataRF.csv"))

write.csv(treeCyls_CCF, paste0(dir,"TreeCylsCCF.csv"))
write.csv(treeCyls_RF, paste0(dir,"TreeCylsRF.csv"))






