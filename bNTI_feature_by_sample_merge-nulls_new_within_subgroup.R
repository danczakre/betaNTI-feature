### Script to combine bNTI(feature) nulls and compare to observed values

Sample_Name = "bNTI_Feat_ASV_withConsp" # Input sample name
Factor = "Type"
Level = "DNA"
Sub.Factor = "Cumulative.Treatment"
Sub.Level = "Inundated"

# Switches for script behaviors
rm.conspec = F # Remove conspecifics
abund.weig = T # Weight values by relative abundances
type = F # This has limited functionality specific to my project, can be safely ignored
rm.tax = T # This configures whether the assemblage data has taxonomy

### Load in necessary libraries
require(Rfast) # For faster variant of finding column minimum
require(dplyr) # For inner joining (faster than merge)
require(abind) # For list->array function
require(picante) # For match.phylo.data
require(phytools) # For midpoint.root
require(ggtree); require(ggplot2); require(reshape2) # For plotting the results

# Combine matrices as array
acomb = function(...) abind(..., along = 3)


# ################## #
#### Load in data ####
# ################## #

# setwd("~/Documents/bNTI Feature Manuscript/FTICR Data/")
# data = read.csv("Processed_ECA_8ppm_Data.csv", row.names = 1) # Load in assemblage
# tree = read.tree("ECA_8ppm_MCD_UPGMA.tre") # Load in tree
# meta = read.csv("../FTICR_Metadata.csv")

setwd("~/Documents/bNTI Feature Manuscript/16S ASV Data/")
data = read.table("ECA_feature-table-with-taxonomy.tsv", skip = 1, comment.char = "", sep = "\t", row.names = 1, header = T)
tree = read.tree("ECA_tree.nwk")
meta = read.csv("ASV_Metadata.csv")

# #################### #
#### Pre-processing ####
# #################### #

# If there is a taxonomy column, it needs to be removed
if(rm.tax == T){
  data = data[,-which(colnames(data) %in% "taxonomy")]
}

# Removing low sequence data
data = data[,-which(colSums(data) < 10000)]
data = data[-which(rowSums(data) == 0),]

# If not abundance weighted, setting data to presence/absence
if(abund.weig == F){
  data[data > 0] = 1
}

# Selecting a subset of the data if specified
if(!type == F){
  data = data[,grep(type, colnames(data))]
  data = data[-which(rowSums(data) == 0),]
}

### Subsetting data based upon the group
# Ensure meta and data are in the same order
data = data[,which(colnames(data) %in% meta$Sample.ID)]
meta = meta[which(meta$Sample.ID %in% colnames(data)),]

if(!identical(meta$Sample.ID, colnames(data))){
  print("Your metadata/factor sheet doesn't match the provided data, attepmting to fix it")
  data = data[,meta$Sample.ID]
  
  if(!identical(meta$Sample.ID, colnames(data))){
    stop("Your sample names couldn't be easily fixed, so I'm stopping to prevent damage")
  }
  
  print("Names were fixed")
}

# Selecting 1st level subset
fac.col = which(colnames(meta) %in% Factor)
fac.samp = meta$Sample.ID[which(meta[,fac.col] %in% Level)]
data = data[,which(colnames(data) %in% fac.samp)]
meta = meta[which(meta$Sample.ID %in% fac.samp),]

# Selecting 2nd level subset
fac.col = which(colnames(meta) %in% Sub.Factor)
fac.samp = meta$Sample.ID[which(meta[,fac.col] %in% Sub.Level)]
data = data[,which(colnames(data) %in% fac.samp)]

if(min(rowSums(data)) == 0){
  data = data[-which(rowSums(data) == 0),]
}

rm(fac.col, fac.samp, meta)

### Matching data and rooting tree
# tree = midpoint.root(tree) # Rooting the tree for consistent results
phylo = match.phylo.data(tree, data) # Matching ICR dataset to the tree

data = t(phylo$data)
tree = phylo$phy

rm("phylo")

# Storing dimensions
samp.num = dim(data)[1]
mem.num = dim(data)[2]

# ####################### #
#### Merging null reps ####
# ####################### #

# Merging the separate bMNTD files
files = list.files(path = paste0("/Users/danc783/Feature-level by Samp Null Reps/", Sample_Name, "/", Factor, "-", Level, "_", Sub.Factor, "-", Sub.Level),
                   pattern = "Feature_Samp_Null", full.names = T) # Listing files

null.by.samp = NULL # Dummy object

for(curr.file in files){
  temp = as.data.frame(read.csv(curr.file, row.names = 1))
  null.by.samp = c(null.by.samp, list(temp))
} # Merging individual 

null.by.samp = do.call(acomb, null.by.samp)
rm("curr.file")

# ###################################### #
#### Determine observed bNTI(feature) ####
# ###################################### #

# Measuring distances
coph = cophenetic(tree)

# Removing conspecifics, if desired
if(rm.conspec){
  coph[coph == 0] = NA
}

# Creating empty object 
min.neigh = NULL # Creating an empty matrix to store pairwise results
comp.names = NULL

# Running through the pairwise comparisons
for (i in 1:(samp.num - 1)) {
  for (j in (i + 1):samp.num) {
    
    # Selecting members in samples and subsetting coph. correspondingly
    samp1 = colnames(data[i, data[i, ] > 0, drop = FALSE]) # Members in the first sample
    samp2 = colnames(data[j, data[j, ] > 0, drop = FALSE]) # Members in the second sample
    
    pair.dist = coph[samp1, samp2, drop = FALSE]
    
    # First sample minimums
    if(length(which(is.na(pair.dist[,1]))) > 0){
      min.dist1 = apply(pair.dist, 1, min, na.rm = T)
    } else {
      min.dist1 = rowMins(pair.dist, value = T)
    } # There is a bug with Rfast minimum calculations that causes it to report NA if it is the first value
    
    names(min.dist1) = row.names(pair.dist)
    
    # Second sample minimums
    if(length(which(is.na(pair.dist[1,]))) > 0){
      min.dist2 = apply(pair.dist, 2, min, na.rm = T)
    } else {
      min.dist2 = colMins(pair.dist, value = T)
    } # There is a bug with Rfast minimum calculations that causes it to report NA if it is the first value
    
    names(min.dist2) = colnames(pair.dist)
    
    # If desired, weighting distances by relative abundance
    if(abund.weig){
      min.dist1 = min.dist1*(data[i, samp1]/sum(data[i, samp1]))
      min.dist2 = min.dist2*(data[j, samp2]/sum(data[j, samp2]))
    }
    
    min.dist = c(min.dist1, min.dist2) # Combining those minimum distances
    min.dist = data.frame(Names = names(min.dist), Dist = min.dist) # Converting to a data frame for easier aggregation
    
    if(rm.conspec){
      # If there are duplicates when conspecifics have been removed, these are averaged as they will have different values
      min.dist = aggregate(Dist~Names, min.dist, FUN = mean) 
    } else {
      # If conspecifics haven't been removed, duplicates will by necessity be 0 - if a community member 
      # is in both samples, it has to be its own nearest neighbor
      if(length(which(duplicated(min.dist$Names) %in% TRUE)) > 0){
        if(mean(min.dist$Dist[duplicated(min.dist$Names)]) == 0){
          min.dist = min.dist[!duplicated(min.dist$Names),]
        } else {
          # If the duplicates aren't 0, something is wrong in the data - stopping the script
          stop("Something odd is happening with your duplicated values. Check that out.")
        }
      }
    } # Checking to ensure duplicates exist at all
    
    # Creating an object with all commiunity members to merge in those which were present in these two samples
    merge.dist = data.frame(Names = colnames(data))

    # Adding observed community member minimum distances to 
    merge.dist = left_join(x = merge.dist, y = min.dist, by = "Names")
    row.names(merge.dist) = merge.dist$Names
    
    # Removing names column
    merge.dist = merge.dist[,-which(colnames(merge.dist) %in% "Names"),drop = F]
    
    # Adding this pairwise comparison to the overall object
    min.neigh = cbind(min.neigh, as.matrix(merge.dist))
    comp.names = c(comp.names, paste0(row.names(data)[i], "-", row.names(data)[j]))
    
    # Clocking
    print(c(i, j, date()))
  }
} # Loop works through the pairwise comparisons

# Clean-up
rm(i, j, samp1, samp2, min.dist, min.dist1, min.dist2, merge.dist)

# Setting column names
colnames(min.neigh) = comp.names

# Finding average minimum neighbor distance by sample
min.by.samp = matrix(data = NA, nrow = nrow(min.neigh), ncol = nrow(data))
row.names(min.by.samp) = row.names(min.neigh)
colnames(min.by.samp) = row.names(data)

for(i in 1:samp.num){
  # Selecting current sample
  curr.samp = grep(paste0("^",row.names(data)[i], "-|", "-", row.names(data)[i], "$"), colnames(min.neigh))
  
  # Selecting members in current sample; min.neigh row.names and column names on data are identical
  curr.mem = which(row.names(min.by.samp) %in% names(data[i,which(data[i,] > 0)]))
  
  # Creating temp object
  temp = min.neigh[,curr.samp]
  
  # Setting all values not for the current sample to NA
  temp[-curr.mem,] = NA
  
  # Adding to final output matrix
  min.by.samp[,i] = rowMeans(temp, na.rm = F)
  
  # Clean-up
  rm(temp, curr.mem, curr.samp)
} 


# ############################# #
#### Calculate bNTI(feature) ####
# ############################# #

# Measuring the difference between observed and null distances
if(!identical(row.names(min.by.samp), row.names(null.by.samp))){
  stop("Your null and observed results do not match. Please double check them.")
}

feat.by.samp = matrix(c(NA), nrow = nrow(null.by.samp), ncol = ncol(null.by.samp))
row.names(feat.by.samp) = row.names(min.by.samp)
colnames(feat.by.samp) = colnames(min.by.samp)

for(i in 1:ncol(null.by.samp)){
  for(j in 1:nrow(null.by.samp)){
    m = null.by.samp[j,i,] # Just setting all the randomizations for a given comparison to a matrix
    feat.by.samp[j,i] = ((min.by.samp[j,i]-mean(m))/sd(m)) # The bNTI calculation
  }
}

write.csv(feat.by.samp, paste0(Sample_Name, "_", Factor, "-", Level, "_", Sub.Factor, "-", Sub.Level, "_bNTI_feature_by_samp_", length(files), "rep.csv"), quote = F)


# ################## #
#### Plot results ####
# ################## #

# Creating data frame to make plotting easier
feat.samp.frame = data.frame(Member = row.names(feat.by.samp), Value = feat.by.samp[,2], Direction = "Insignificant", stringsAsFactors = F)
feat.samp.frame$Direction[which(feat.samp.frame$Value <= -1)] = "Low"
feat.samp.frame$Direction[which(feat.samp.frame$Value >= 1)] = "High"
feat.samp.frame$Direction[which(feat.samp.frame$Value <= -2)] = "Sig. Low"
feat.samp.frame$Direction[which(feat.samp.frame$Value >= 2)] = "Sig. High"

feat.samp.frame$Direction = factor(feat.samp.frame$Direction, levels = c("Sig. High", "High", "Insignificant", "Low", "Sig. Low"))

# Plotting a standard tree with a point graph to see the scale of differences
p = ggtree(tree)
p = facet_plot(p, panel = "Sample Nearest Neighbor", data = feat.samp.frame, 
               geom = geom_point, mapping = aes(x = Value, y = y, color = Direction))
p + scale_color_manual(values = c("firebrick4", "firebrick1", "goldenrod1", "dodgerblue1", "dodgerblue4"), drop = F) + 
  theme_tree2()+theme(legend.position = "left")

# Plotting variable members
feat.order = order(apply(feat.by.samp, 1, mean, na.rm = F), decreasing = F) # Needed to find most variable members
feat.melt = melt(feat.by.samp[feat.order[1:5],]) # Melting most variable members for ggplot

ggplot(data = feat.melt, aes(x = Var2, y = value, group = as.character(Var1)))+
  geom_point(aes(color = as.character(Var1)))+
  geom_line(aes(color = as.character(Var1)))+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(colour = "black", angle = -40, vjust = 1, hjust = 0),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank())

