
# Load packages
rm(list = ls())

library(ggplot2)
library(raster)
library(maps)
library(lme4)
library(nlme)
library(performance)
library(rsample)
library(yardstick)
library(dplyr)
library(Metrics)


# load data
# CCF = Continuous Cover Forestry 
ringD_CCF <- read.delim("F:/Tree_segments/Article2_RFiles/data/rings_ERI.txt") # CCF ring data
TLSTrees_CCF <- read.csv("F:/Tree_segments/CCF/TreeData/treeData.csv", sep = ";") # Tree data CCF
TLSTrees_CCF <- TLSTrees_CCF[!is.na(TLSTrees_CCF$vol),] # drop rows with na values in volume col

# RF
ringD_RF <- read.delim("F:/Tree_segments/Article2_RFiles/data/rings_HON.txt") # RF ring data
treeD <- read.csv("F:/Tree_segments/Article2_RFiles/data/Trees_fieldM.csv", sep = ";") # Fiel data CCF
treeD_RF <- read.delim("F:/Tree_segments/Article2_RFiles/data/field_RF.txt") # Field data RF
Field_RF <- read.delim("F:/Tree_segments/Article2_RFiles/data/field_RF.txt")

#######################################
### CCF ###############################
### Density prediction for CCF data ###
#######################################

# Remove year 2021 and 2020, as there seams to be error in measurement  
ringD_CCF <- ringD_CCF[!ringD_CCF$Year == max(ringD_CCF$Year),]
ringD_CCF <- ringD_CCF[!ringD_CCF$Year == max(ringD_CCF$Year),]

#### Ring sample data preparation ###

# Loop calculates:
# 1. distance from pith to each outer edge of ring sample
# 2. relative height of sample disc 
for(i in levels(as.factor(ringD_CCF$Tree))){
  
  tree <- ringD_CCF[ringD_CCF$Tree == i,]
  
  for (j in levels(as.factor(tree$Sample))) {
    sample <- tree[tree$Sample == j,]
    sample <- sample[order(sample$CA),]
    
    s = 0
    for (k in 1:length(sample$Tree)) {
      # Distance to the outer edge of sample
      r <- sample[k,]
      s <- s + r$RW
      ringD_CCF$R_new[ringD_CCF$Tree == i & 
                        ringD_CCF$Sample == j &
                        ringD_CCF$CA == r$CA] <- s
      # relative height of sample disc
      h <- treeD$H[treeD$ID==i]
      rel_h <- r$sampleHeight/h
      ringD_CCF$rel_sampleH[ringD_CCF$Tree == i & 
                              ringD_CCF$Sample == j &
                              ringD_CCF$CA == r$CA] <- rel_h
    }
  }
}

# Change treeD data.frame column names to match ringD data.frame
colnames(treeD)[colnames(treeD) == "Tree"] <- "PlotTreeID"
colnames(treeD)[colnames(treeD) == "ID"] <- "Tree"

ringD_CCF <- merge(ringD_CCF, treeD[,c("Tree", "Plot")], by = "Tree", all.x = TRUE)

# Remove EVO2 Tree 13 sample 1 and tree 14 sample 1
# As there seams to be some problem with measurements
ringD_CCF <- ringD_CCF[!(ringD_CCF$Stand == "EVO2" & ringD_CCF$Tree == 13 & ringD_CCF$Sample == 1), ]
ringD_CCF <- ringD_CCF[!(ringD_CCF$Stand == "EVO2" & ringD_CCF$Tree == 14 & ringD_CCF$Sample == 1), ]

### Ring density model ### 

m_CCF <- lme(
  log(RD) ~ log(R_new) + log(rel_sampleH),
  random = ~1|Plot/Tree,
  correlation = corCAR1(form = ~CA|Plot/Tree/Sample),
  data = ringD_CCF,
  na.action = na.omit
)

# Model summary
summary(m_CCF)
model_performance(m_CCF)
plot(m_CCF)
icc(m_CCF, by_group = T)
VarCorr(m_CCF)

### Ring density Predicting ### 

# Define ring density classes for measured density
ringD_CCF$RD_cut <- cut(ringD_CCF$RD, breaks = c(0, 0.40, 0.42, 0.44, 0.46, 0.48, 0.50,
                                                 0.52, 0.54, 0.56, 0.58, 0.60, 0.62), labels = F)

# Change plot names in TLS data to correspond with ring data
{
  TLSTrees_CCF$plot[TLSTrees_CCF$plot == "evo2_1" | TLSTrees_CCF$plot == "evo2_2"] <- "EVO2"
  TLSTrees_CCF$plot[TLSTrees_CCF$plot == "evo3_1" | TLSTrees_CCF$plot == "evo3_2"] <- "EVO3"
  TLSTrees_CCF$plot[TLSTrees_CCF$plot == "evo4"] <- "EVO4"
  TLSTrees_CCF$plot[TLSTrees_CCF$plot == "vep2_1" | TLSTrees_CCF$plot == "vep2_2"] <- "VEP2"
  TLSTrees_CCF$plot[TLSTrees_CCF$plot == "vep4_1" | TLSTrees_CCF$plot == "vep4_2"] <- "VEP4"
  TLSTrees_CCF$plot[TLSTrees_CCF$plot == "ves16_1" | TLSTrees_CCF$plot == "ves16_2"] <- "VES16"
  TLSTrees_CCF$plot[TLSTrees_CCF$plot == "ves2_1" | TLSTrees_CCF$plot == "ves2_2"] <- "VES2"
}


### Create CCF density matrix ### 

# Define max height and dbh for density matrix 
rows <- rev(seq(from = 0.1, to = max(TLSTrees_CCF$h), by = 0.1)) # Max height from max height of TLS tree
cols <- seq(from = 1, to = (max(ringD_CCF$R_new)/10)+1, by = 1) # max r from max r of r

# Create the matrix based on rows and column defined above
densTree_CCF <- matrix(, nrow = length(rows), ncol = length(cols)) # empty matrix

colnames(densTree_CCF) <- cols # define col name
rownames(densTree_CCF) <- rows # define row name

# Predict wood density at each defined height and dbh
for (j in levels(as.factor(cols))) {
  for (k in levels(as.factor(rows))) {
    frame <- data.frame(R_new = as.numeric(j), rel_sampleH = as.numeric(k)/max(rows))
    value <- exp(predict(m_CCF,frame, level = 0))
    densTree_CCF[rownames(densTree_CCF) == k, colnames(densTree_CCF) == j] = as.numeric(value)
  }
}

# Save the matrix (Change dir) 
write.table(densTree_CCF, "F:/Tree_segments/Article2_RFiles/data/dens_mats/CCF_mat.txt")




##############################
### RF #######################
### RF = Rotation Forestry ###
##############################

#### Same analysis for RF as was done above for CCF ###

ringD_RF <- ringD_RF[!ringD_RF$Tree == 1235,]

for(i in levels(as.factor(ringD_RF$Tree))){
  
  tree <- ringD_RF[ringD_RF$Tree == i,]
  
  for (j in levels(as.factor(tree$Sample))) {
    sample <- tree[tree$Sample == j,]
    sample <- sample[order(sample$CA),]
    s = 0
    for (k in 1:length(sample$Tree)) {
      r <- sample[k,]
      s <- s + r$RW
      ringD_RF$R_new[ringD_RF$Tree == i & 
                       ringD_RF$Sample == j &
                       ringD_RF$CA == r$CA] <- s
      
      h <- treeD_RF$h_field[treeD_RF$TreeID==i]
      rel_h <- r$sampleHeight/h
      ringD_RF$rel_sampleH[ringD_RF$Tree == i & 
                             ringD_RF$Sample == j &
                             ringD_RF$CA == r$CA] <- rel_h  
      
      
    }
  }
}
ringD_RF$R_class <- cut(ringD_RF$R_new,
                        c(seq(0,265,5),Inf), labels = F)

# remove last measurement
ringD_RF <- ringD_RF[!ringD_RF$Year == max(ringD_RF$Year),]

colnames(treeD_RF)[colnames(treeD_RF) == "TreeID"] <- "Tree"

ringD_RF <- merge(ringD_RF, treeD_RF[,c("Tree", "Plot")], by = "Tree", all.x = TRUE)

ringD_RF$Plot <- gsub("\\_.*", "", ringD_RF$Plot)

# Remove OMT forest stands (We are only interested in MT stands) 
ringD_RF <- ringD_RF[!(ringD_RF$Plot %in% "2" |
                         ringD_RF$Plot %in% "86" |
                         ringD_RF$Plot %in% "5"),]


ringD_RF$rel_sampleH[ringD_RF$rel_sampleH == 0] <- 0.0001 #change rel_sample h 0 to 0.0001, because log(0) = -inf
m_RF <-  lme(
  log(RD) ~ log(R_new) + log(rel_sampleH),
  random = ~1|Plot/Tree,
  correlation = corCAR1(form = ~CA|Plot/Tree/Sample),
  data = ringD_RF,
  na.action = na.omit
)
summary(m_RF)
model_performance(m_RF)
plot(m_RF)

ringD_RF$RD_cut <- cut(ringD_RF$RD, breaks = c(0, 0.40, 0.42, 0.44, 0.46, 0.48, 0.50,
                                               0.52, 0.54, 0.56, 0.58, 0.60, 0.62), labels = F)

# RF treatment matrix
rows <- rev(seq(from = 0.1, to = max(treeD_RF$h_field), by = 0.1)) # Max height from max height of TLS tree
cols <- seq(from = 1, to =( max(ringD_RF$R_new)/10)+1, by = 1) # max r from max r of r

densTree_RF <- matrix(, nrow = length(rows), ncol = length(cols)) # empty matrix

colnames(densTree_RF) <- cols # define col name
rownames(densTree_RF) <- rows # define row name

for (j in levels(as.factor(cols))) {
  for (k in levels(as.factor(rows))) {
    frame <- data.frame(R_new = as.numeric(j), rel_sampleH = as.numeric(k)/max(rows))
    value <- exp(predict(m_RF,frame, level = 0))
    densTree_RF[rownames(densTree_RF) == k, colnames(densTree_RF) == j] = as.numeric(value)
  }
}

### Save the data (Change dir) ###
write.table(densTree_RF, "F:/Tree_segments/Article2_RFiles/data/dens_mats/RF_mat.txt")


###############################################
### Create density matrix for visualization ###
###############################################

###############
### For CCF ###
###############

# Example tree stem curve (Laasasenaho model/param)
laasis <- function(relH) 2.3366*(1-relH)^1 + -3.2684*(1-relH)^2 + 3.6513*(1-relH)^3 +-2.2608*(1-relH)^5 + 0*(1-relH)^8 + 2.1501*(1-relH)^13 + -2.7412*(1-relH)^21 + 1.8876*(1-relH)^34

ex_tree <- data.frame(h = seq(from = 0, to = max(as.numeric(row.names(densTree_CCF))), by = 0.1)) 
ex_tree$relH <- ex_tree$h/max(ex_tree$h)
ex_tree$relD <- laasis(ex_tree$relH)
ex_tree$dbh <- ex_tree$relD*max(as.numeric(colnames(densTree_CCF)))/ex_tree$relD[ex_tree$h==0]

plot(ex_tree$dbh,ex_tree$h)

rows <- rev(seq(from = 0.1, to = max(TLSTrees_CCF$h), by = 0.1)) # Max height from max height of TLS tree
cols <- seq(from = 1, to = (max(ringD_CCF$R_new)/10)+1, by = 1) # max r from max r of r

densTree_CCF_2 <- matrix(, nrow = length(rows), ncol = length(cols)) # empty matrix

colnames(densTree_CCF_2) <- cols # define col name
rownames(densTree_CCF_2) <- rows # define row name

for (i in levels(as.factor(ex_tree$h))) {
  if(i == 0){
    next
  }
  dbh = round(ex_tree$dbh[ex_tree$h == i], digits = 0)
  if(dbh == 0){
    next
  }
  row <- densTree_CCF[as.character(i),]
  row <- unlist(split(row,ceiling(seq_along(row)/dbh))[1],use.names = F)
  for (j in 1:length(row)) {
    val <- row[j]
    densTree_CCF_2[rownames(densTree_CCF_2) == i, colnames(densTree_CCF_2) == j] = as.numeric(val)
  }
}

##############
### For RF ###
##############

# Example tree stem curve (Laasasenaho model/param)
laasis <- function(relH) 2.3366*(1-relH)^1 + -3.2684*(1-relH)^2 + 3.6513*(1-relH)^3 +-2.2608*(1-relH)^5 + 0*(1-relH)^8 + 2.1501*(1-relH)^13 + -2.7412*(1-relH)^21 + 1.8876*(1-relH)^34

ex_tree <- data.frame(h = seq(from = 0, to = max(as.numeric(row.names(densTree_RF))), by = 0.1)) 
ex_tree$relH <- ex_tree$h/max(ex_tree$h)
ex_tree$relD <- laasis(ex_tree$relH)
ex_tree$dbh <- ex_tree$relD*max(as.numeric(colnames(densTree_RF)))/ex_tree$relD[ex_tree$h==0]

plot(ex_tree$dbh,ex_tree$h)

rows <- rev(seq(from = 0.1, to = max(treeD_RF$h_field), by = 0.1)) # Max height from max height of TLS tree
cols <- seq(from = 1, to =( max(ringD_RF$R_new)/10)+1, by = 1) # max r from max r of r

densTree_RF_2 <- matrix(, nrow = length(rows), ncol = length(cols)) # empty matrix

colnames(densTree_RF_2) <- cols # define col name
rownames(densTree_RF_2) <- rows # define row name

for (i in levels(as.factor(ex_tree$h))) {
  if(i == 0){
    next
  }
  dbh = round(ex_tree$dbh[ex_tree$h == i], digits = 0)
  if(dbh == 0){
    next
  }
  row <- densTree_RF[as.character(i),]
  row <- unlist(split(row,ceiling(seq_along(row)/dbh))[1],use.names = F)
  for (j in 1:length(row)) {
    val <- row[j]
    densTree_RF_2[rownames(densTree_RF_2) == i, colnames(densTree_RF_2) == j] = as.numeric(val)
  }
}


### Save CCF and RF visualization matrix (Change dir) ###
write.table(densTree_CCF_2, "F:/Tree_segments/Article2_RFiles/data/dens_mats/CCF_mat_2.txt")
write.table(densTree_RF_2, "F:/Tree_segments/Article2_RFiles/data/dens_mats/RF_mat_2.txt")

# Check mean and weighted mean density for CCF and RF
mean(ringD_CCF$RD)
weighted.mean(ringD_CCF$RD, ringD_CCF$RW)

mean(ringD_RF$RD)
weighted.mean(ringD_RF$RD, ringD_RF$RW)


