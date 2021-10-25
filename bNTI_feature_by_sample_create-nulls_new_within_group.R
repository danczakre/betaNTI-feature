### R-script to generate individual null bNTI(feature) reps
# Adapted from the "comdistnt" in the picante package by RE Danczak

Sample_Name = "bNTI_Feat_ASV_withConsp" # Input sample name
Factor = "Type"
Level = "DNA"

# Switches for script behaviors
rm.conspec = F # Remove conspecifics
abund.weig = T # Weight values by relative abundances
type = F # This has limited functionality specific to my project, can be safely ignored
rm.tax = T # This configures whether the assemblage data has taxonomy
noise = T # This adds "noise" to the nulls to ensure 0's don't exist in the null
range = 11:50 # Setting replicate numbers (I don't recommend more than 99 at the moment, time-consuming)


### Load in necessary libraries
require(Rfast) # For faster variant of finding column minimum
require(dplyr) # For left joining
library(picante) # For match.phylo.data
require(phytools) # For midpoint.root

print(date())


# ################## #
#### Load in data ####
# ################## #

# setwd("~/Documents/bNTI Feature Manuscript/FTICR Data/")
# data = read.csv("Processed_ECA_8ppm_Data.csv", row.names = 1) # Load in assemblage
# tree = read.tree("ECA_8ppm_MCD_UPGMA.tre") # Load in tree
# meta = read.csv("FTICR_Metadata.csv")

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
  
  # Removing low sequence data
  data = data[,-which(colSums(data) < 10000)]
  data = data[-which(rowSums(data) == 0),]
}

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

# Selecting sample subset
fac.col = which(colnames(meta) %in% Factor)
fac.samp = meta$Sample.ID[which(meta[,fac.col] %in% Level)]
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


# ##################################### #
#### Looping through null replicates ####
# ##################################### #

# Measuring distances
coph = cophenetic(tree)

# Removing conspecifics, if desired
if(rm.conspec){
  coph[coph == 0] = NA
}

# Creating noise object
if(noise){
  min.coph = min(coph[!(coph == 0)])
  noise.list = seq(from = min.coph*(10^-20), to = min.coph*(5*10^-20), by = min.coph*(10^-21))
}

# Creating empty null object
for(rep in range){
  print(paste0("Start of rep #", rep, " - ", date()))
  
  null = taxaShuffle(coph) # Randomizing cophenetic distances
  null.rep = NULL # Creating empty object
  comp.names = NULL # Creating an empty object to store comparison names
  
  for (i in 1:(samp.num - 1)) {
    for (j in (i + 1):samp.num) {  
      samp1 = colnames(data[i, data[i, ] > 0, drop = FALSE])
      samp2 = colnames(data[j, data[j, ] > 0, drop = FALSE])
      
      pair.null = null[samp1, samp2, drop = FALSE]
      
      # First sample minimums
      if(length(which(is.na(pair.null[,1]))) > 0){
        min.null1 = apply(pair.null, 1, min, na.rm = T)
      } else {
        min.null1 = rowMins(pair.null, value = T)
      } # There is a bug with Rfast minimum calculations that causes it to report NA if it is the first value
      
      names(min.null1) = row.names(pair.null)
      
      # Second sample minimums
      if(length(which(is.na(pair.null[1,]))) > 0){
        min.null2 = apply(pair.null, 2, min, na.rm = T)
      } else {
        min.null2 = colMins(pair.null, value = T)
      } # There is a bug with Rfast minimum calculations that causes it to report NA if it is the first value
      
      names(min.null2) = colnames(pair.null)
      
      # Abundance weighting, if set
      if(abund.weig){
        min.null1 = min.null1*(data[i, samp1]/sum(data[i, samp1]))
        min.null2 = min.null2*(data[j, samp2]/sum(data[j, samp2]))
      }
      
      # Combining minimum distances
      null.dist = c(min.null1, min.null2)
      null.dist = data.frame(Names = names(null.dist), Dist = null.dist)
      
      if(rm.conspec){
        null.dist = aggregate(Dist~Names, null.dist, FUN = mean) 
      } else {
        if(length(which(duplicated(null.dist$Names) %in% TRUE)) > 0){
          if(mean(null.dist$Dist[duplicated(null.dist$Names)]) == 0){
            null.dist = null.dist[!duplicated(null.dist$Names),]
          } else {
            stop("Something odd is happening with your duplicated values. Check that out.")
          }
        }
      } # Resolving repeats - if conspecifics are removed, averaging differences
      # If conspecifics weren't removed, then duplicated values should always be zero
      
      if(noise){
        null.dist$Dist[null.dist$Dist == 0] = sample(noise.list, size = length(null.dist$Dist[null.dist$Dist == 0]), replace = T)
      } # Divide-by-zero errors result in NaN values which leave the members untrackable
      
      # Creating an object with all commiunity members to merge in those which were present in these two samples
      merge.null = data.frame(Names = colnames(data))
      
      # Adding observed community member minimum distances to 
      merge.null = left_join(x = merge.null, y = null.dist, by = "Names")
      row.names(merge.null) = merge.null$Names
      
      # Removing names column
      merge.null = merge.null[,-which(colnames(merge.null) %in% "Names"),drop = F]
      
      # Combining pairwise comparisons
      null.rep = cbind(null.rep, as.matrix(merge.null))
      comp.names = c(comp.names, paste0(row.names(data)[i], "-", row.names(data)[j]))
      
    }
  }
  
  # Setting column names
  colnames(null.rep) = comp.names
  
  # Creating empty matrix to store by sample data
  null.by.samp = matrix(data = NA, nrow = nrow(null.rep), ncol = nrow(data))
  row.names(null.by.samp) = row.names(null.rep)
  colnames(null.by.samp) = row.names(data)
  
  for(i in 1:samp.num){
    # Selecting current sample
    curr.samp = grep(paste0("^",row.names(data)[i], "-|", "-", row.names(data)[i], "$"), colnames(null.rep))
    
    # Selecting members in current sample; null.rep row.names and column names on data are identical
    curr.mem = which(row.names(null.by.samp) %in% names(data[i,which(data[i,] > 0)]))
    
    # Creating temp object
    temp = null.rep[,curr.samp]
    
    # Setting all values not for the current sample to NA
    temp[-curr.mem,] = NA
    
    # Adding to final output matrix
    null.by.samp[,i] = rowMeans(temp, na.rm = F)
    
    # Clean-up
    rm(temp, curr.mem, curr.samp)
  }
  
  # Creating output directories
  if(!dir.exists(paste0("/Users/danc783/Feature-level by Samp Null Reps/", Sample_Name, "/"))){
    dir.create(paste0("/Users/danc783/Feature-level by Samp Null Reps/", Sample_Name, "/"))
  }
  
  if(!dir.exists(paste0("/Users/danc783/Feature-level by Samp Null Reps/", Sample_Name, "/", Factor, "-", Level))){
    dir.create(paste0("/Users/danc783/Feature-level by Samp Null Reps/", Sample_Name, "/", Factor, "-", Level))
  }
  
  write.csv(null.by.samp, paste0("/Users/danc783/Feature-level by Samp Null Reps/", Sample_Name, "/", Factor, "-", Level, "/", 
                                 Sample_Name, "-", Factor, "_", Level, "_Feature_Samp_Null_rep", rep, ".csv"), quote = F, row.names = T)
  
  print(paste0("End of rep #", rep, " - ", date()))
}
