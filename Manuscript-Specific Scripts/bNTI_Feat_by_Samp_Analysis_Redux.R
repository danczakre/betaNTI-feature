### Analyzing single sample bNTI feature

# Load in libraries
library(ggplot2); library(reshape2); library(gplots); library(ggpubr); library(ggtree)
library(plyr); library(dplyr); library(stringr)
library(vegan)
library(picante)
library(Hmisc)
library(igraph); library(qgraph); library(WGCNA)

# Focus sample
focus = "ECA_0Cyc_R2"

# Flags
plot.results = F
spear.cor = F
WGCNA = T


# ###################### #
#### Define Functions ####
# ###################### #

# ggplot theme
plot.theme = theme(axis.text = element_text(colour = "black", size = 12),
                   axis.title = element_text(colour = "black", size = 14),
                   panel.border = element_rect(size = 1, colour = "black"),
                   panel.grid = element_blank())

# Add direction to a frame
make.direction = function(feat.frame, val.col){
  
  # Empty object
  direction = rep("Insignificant", nrow(feat.frame))
  
  # Set direction based on bNTI value
  direction[which(feat.frame[,val.col] <= -1)] = "Low"
  direction[which(feat.frame[,val.col] >= 1)] = "High"
  direction[which(feat.frame[,val.col] <= -2)] = "Sig. Low"
  direction[which(feat.frame[,val.col] >= 2)] = "Sig. High"
  direction = factor(direction, levels = c("Sig. High", "High", "Insignificant", "Low", "Sig. Low"))
  
  return(direction)
}

# Rename function
rename.samples = function(bnti.feat, factors, type){
  samp.col = grep(type, colnames(factors))
  
  for(i in 1:ncol(bnti.feat)){
    w = which(meta[,samp.col] %in% colnames(bnti.feat)[i])
    colnames(bnti.feat)[i] = meta$SampleName[w]
  }
  
  bnti.feat = bnti.feat[,order(colnames(bnti.feat))]
  
  return(bnti.feat)
}

# Divergent plots
div.plot.by.factor = function(bnti.feat, factors, time.plot, adjust.p){
  
  # Missing flag
  if(missing(adjust.p)){
    adjust.p = T
  }
  
  if(missing(time.plot)){
    time.plot = F
  }
  
  # Empty objects
  treat.mwu = NULL
  cycle.ano = NULL
  mem.names = NULL
  
  # Loop through stats
  for(i in 1:nrow(bnti.feat)){
    
    # Temp dataset
    temp = data.frame(bnti_value = as.numeric(bnti.feat[i,]), Cycle = factors$Cycles, Treat = factors$Cumulative.Treatment)
    
    # Treatment comparison
    w = wilcox.test(bnti_value~Treat, data = temp)
    treat.mwu = rbind(treat.mwu, data.frame(MWU = w$statistic, p.value = w$p.value))
    
    # Cycle comparison
    ano = summary(aov(bnti_value~Cycle, data = temp)) # Across all cycles
    cycle.ano = rbind(cycle.ano, data.frame(ANOVA = ano[[1]]$`F value`[1], p.value = ano[[1]]$`Pr(>F)`[1]))
    
    # Member names
    mem.names = c(mem.names, row.names(bnti.feat)[i])
  }
  
  rm(ano, w, i, temp)
  
  # Assign names
  row.names(treat.mwu) = mem.names
  row.names(cycle.ano) = mem.names
  
  # p-value adjustmnet
  if(adjust.p){
    treat.mwu$p.value = p.adjust(treat.mwu$p.value, method = "fdr")
    cycle.ano$p.value = p.adjust(cycle.ano$p.value, method = "fdr")
  }
  
  ### Plotting most divergent across treatments
  # Creating object for plotting
  div.treat = as.data.frame(t(bnti.feat[which(treat.mwu$p.value < 0.05),]))
  div.treat$Treat = factors$Cumulative.Treatment
  div.treat = melt(div.treat, id.vars = "Treat")
  
  # Plotting boxplot
  print(
    ggplot(data = div.treat, aes(x = reorder(variable, value, FUN = median), y = value))+
      geom_boxplot(aes(fill = Treat)) + 
      geom_hline(yintercept = c(-2,2), lty = 2, color = "red") +
      xlab(NULL) + theme_bw() + plot.theme + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  )
  
  # Plotting most variable features by sample
  if(time.plot){
    if(length(which(colnames(factors) %in% "Type")) > 0){
      div.treat = as.data.frame(t(bnti.feat[which(treat.mwu$p.value < 0.05),]))
      most.var = apply(div.treat, 2, sd)
      most.var = names(most.var)[order(most.var, decreasing = T)][1:5]
      div.treat = div.treat[,most.var]
      div.treat$Samples = factors$SampleName
      div.treat$Type = factors$Type
      div.treat = melt(div.treat, id.vars = c("Samples", "Type"))
      
      print(
        ggplot(div.treat, aes(x = Samples, y = value))+
          geom_point()+
          geom_hline(yintercept = 0, color = "gray50", lwd = 0.5)+
          geom_hline(yintercept = c(-1,1), color = "blue", lwd = 0.5, lty = 2)+
          geom_hline(yintercept = c(-2,2), color = "red", lwd = 0.5, lty = 2)+
          facet_grid(variable~Type)+
          xlab(NULL) + theme_bw() + plot.theme +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
      )
    } else {
      div.treat = as.data.frame(t(bnti.feat))
      most.var = apply(div.treat, 2, sd)
      most.var = names(most.var)[order(most.var, decreasing = T)]
      div.treat = div.treat[,most.var]
      div.treat$Samples = factors$SampleName
      div.treat = melt(div.treat, id.vars = c("Samples"))
      
      print(
        ggplot(div.treat, aes(x = Samples, y = value))+
          geom_point()+
          geom_hline(yintercept = 0, color = "gray50", lwd = 0.5)+
          geom_hline(yintercept = c(-1,1), color = "blue", lwd = 0.5, lty = 2)+
          geom_hline(yintercept = c(-2,2), color = "red", lwd = 0.5, lty = 2)+
          facet_grid(variable~.)+
          xlab(NULL) + theme_bw() + plot.theme +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
      )
    }
  }
}

# Summarize data by groups
summarize.data = function(bnti.feat, group.data, group){
  
  bnti.feat$Members = row.names(bnti.feat) # Columns for left joining
  group.data$Members = row.names(group.data)
  
  w = which(colnames(group.data) %in% group) # Setting group name for easy restructuring
  colnames(group.data)[w] = "Group"
  
  bnti.feat = bnti.feat %>% left_join(group.data[,c("Members", "Group")], by = "Members") %>% 
    select(-c("Members")) %>% group_by(Group) %>% summarise_all(.funs = mean) # Summarizing
  
  bnti.feat = as.data.frame(bnti.feat) # Row name stuff
  row.names(bnti.feat) = bnti.feat$Group; bnti.feat = bnti.feat[,-1]
  
  return(bnti.feat)
}


# ############################# #
#### Load and clean the data ####
# ############################# #

setwd("~/Documents/bNTI Feature Manuscript/")

# Sample-resolved bNTI feature
dna.bnti = read.csv("16S ASV Data/bNTI_Feat_ASV_withConsp_bNTI_feature_by_pair_DNA5_999rep.csv", row.names = 1)
rna.bnti = read.csv("16S ASV Data/bNTI_Feat_ASV_withConsp_bNTI_feature_by_pair_cDNA5_999rep.csv", row.names = 1)
icr.bnti = read.csv("FTICR Data//bNTI_Feat_ICR_withConsp_bNTI_feature_by_pair_Stegen_EC_05_0cycles_01Oct19_Alder_Infuse_p08_1_01_49067_999rep.csv", row.names = 1)

# Load in trees/dendrogram
asv.tree = read.tree("16S ASV Data/ECA_tree.nwk")
icr.tree = read.tree("FTICR Data/ECA_8ppm_MCD_UPGMA.tre")

# Taxonomic information
data = read.table("16S ASV Data/ECA_feature-table-with-taxonomy.tsv", skip = 1, comment.char = "", sep = "\t", row.names = 1, header = T)
row.tax = row.names(data)
tax = data[,"taxonomy"]

tax = strsplit(tax, "; ")
tax = ldply(tax, rbind)
colnames(tax) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

for(i in 2:ncol(tax)){
  tax[is.na(tax[,i]), i] = tax[is.na(tax[,i]), (i-1)]
}

row.names(tax) = row.tax

rm(i,row.tax)

# Molecular information
mol = read.csv("FTICR Data/Processed_ECA_8ppm_Mol.csv", row.names = 1)

# Load in metadata
meta = read.csv("Metadata.csv")
meta$Sample_ID_OTU_Table_gDNA = gsub("g", "", meta$Sample_ID_OTU_Table_gDNA)
meta = meta[order(meta$SampleName),]

# Removing focus sample name
colnames(dna.bnti) = gsub("^DNA5.", "", colnames(dna.bnti))
colnames(rna.bnti) = gsub("^cDNA5.", "", colnames(rna.bnti))
colnames(icr.bnti) = gsub("^Stegen_EC_05_0cycles_01Oct19_Alder_Infuse_p08_1_01_49067.", "", colnames(icr.bnti))

# Ensuring data with meta-data is selected
dna.bnti = dna.bnti[,which(colnames(dna.bnti) %in% c(meta$Sample_ID_OTU_Table_gDNA, meta$Sample_ID_OTU_Table_cDNA))]
rna.bnti = rna.bnti[,which(colnames(rna.bnti) %in% c(meta$Sample_ID_OTU_Table_gDNA, meta$Sample_ID_OTU_Table_cDNA))]

# Removing members without bNTI values
dna.bnti = dna.bnti[-which(is.na(rowSums(dna.bnti))),]
rna.bnti = rna.bnti[-which(is.na(rowSums(rna.bnti))),]
icr.bnti = icr.bnti[-which(is.na(rowSums(icr.bnti))),]


# ###################### #
#### Renaming samples ####
# ###################### #

### Creating factor sheet for DNA bNTI
# Empty object
dna.factors = NULL

# Grabbing metadata
for(i in 1:ncol(dna.bnti)){
  if(length(grep("^DNA", colnames(dna.bnti)[i])) > 0){
    w = which(meta$Sample_ID_OTU_Table_gDNA %in% colnames(dna.bnti)[i])
    dna.factors = rbind(dna.factors, meta[w,])
  } else if(length(grep("cDNA", colnames(dna.bnti)[i])) > 0){
    w = which(meta$Sample_ID_OTU_Table_cDNA %in% colnames(dna.bnti)[i])
    dna.factors = rbind(dna.factors, meta[w,])
  }
}

# Matching sample IDs
dna.factors$Sample.ID = colnames(dna.bnti)

# Adding type column
dna.factors$Type = ifelse(grepl("c", dna.factors$Sample.ID), "RNA", "DNA")

# DNA taxonomy
dna.tax = tax[which(row.names(tax) %in% row.names(dna.bnti)),]
dna.tax = dna.tax[match(row.names(dna.bnti), row.names(dna.tax)),]
dna.tax$Members = row.names(dna.tax)

### Creating factor sheet for DNA bNTI
# Empty object
rna.factors = NULL

# Grabbing metadata
for(i in 1:ncol(rna.bnti)){
  if(length(grep("^DNA", colnames(rna.bnti)[i])) > 0){
    w = which(meta$Sample_ID_OTU_Table_gDNA %in% colnames(rna.bnti)[i])
    rna.factors = rbind(rna.factors, meta[w,])
  } else if(length(grep("cDNA", colnames(rna.bnti)[i])) > 0){
    w = which(meta$Sample_ID_OTU_Table_cDNA %in% colnames(rna.bnti)[i])
    rna.factors = rbind(rna.factors, meta[w,])
  }
}

# Matching sample IDs
rna.factors$Sample.ID = colnames(rna.bnti)

# Adding type column
rna.factors$Type = ifelse(grepl("c", rna.factors$Sample.ID), "RNA", "DNA")

# DNA taxonomy
rna.tax = tax[which(row.names(tax) %in% row.names(rna.bnti)),]
rna.tax = rna.tax[match(row.names(rna.bnti), row.names(rna.tax)),]
rna.tax$Members = row.names(rna.tax)


### Rename FTICR
# Renaming ICR bNTI
icr.bnti = rename.samples(icr.bnti, meta, "FTICR")

# ICR factors
icr.factors = meta[which(meta$SampleName %in% colnames(icr.bnti)),]

# Matching molecular information
mol = mol[which(row.names(mol) %in% row.names(icr.bnti)),]
mol = mol[match(row.names(icr.bnti), row.names(mol)),]
mol$Mass = row.names(mol)

rm(i,w)


# #################################### #
#### Generating some standard plots ####
# #################################### #

if(plot.results){
  # Comparison boxplots
  div.plot.by.factor(dna.bnti, dna.factors)
  div.plot.by.factor(rna.bnti, rna.factors)
  div.plot.by.factor(icr.bnti, icr.factors, adjust.p = F)
  
  # ASV tree
  dna.mean = data.frame(Members = row.names(dna.bnti), DNA_means = rowMeans(dna.bnti))
  dna.mean$Direction = make.direction(dna.mean, 2)
  rna.mean = data.frame(Members = row.names(rna.bnti), RNA_means = rowMeans(rna.bnti))
  rna.mean$Direction = make.direction(rna.mean, 2)
  
  p = ggtree(asv.tree)
  p = facet_plot(p, panel = "DNA bNTI Values", data = dna.mean, geom = geom_point,
                 mapping = aes(x = DNA_means, y = y, color = as.factor(Direction)))
  p = facet_plot(p, panel = "RNA bNTI Values", data = rna.mean, geom = geom_point,
                 mapping = aes(x = RNA_means, y = y, color = as.factor(Direction)))
  p + scale_color_manual(values = c("firebrick4", "firebrick1", "goldenrod1", "dodgerblue1", "dodgerblue4"), drop = F) + 
    theme_tree2()+theme(legend.position = "left")
  
  rm(dna.mean, rna.mean)
  
  # ICR dendrogram
  icr.mean = data.frame(Members = row.names(icr.bnti), ICR_means = rowMeans(icr.bnti))
  icr.mean$Direction = make.direction(icr.mean, 2)
  
  p = ggtree(icr.tree)
  p = facet_plot(p, panel = "FTICR-MS bNTI Values", data = icr.mean, geom = geom_point,
                 mapping = aes(x = ICR_means, y = y, color = as.factor(Direction)))
  p + scale_color_manual(values = c("firebrick4", "firebrick1", "goldenrod1", "dodgerblue1", "dodgerblue4"), drop = F) + 
    theme_tree2()+theme(legend.position = "left")
  
  rm(icr.mean, p)
}

# ############################################ #
#### Looking at taxonomic/functional groups ####
# ############################################ #

if(plot.results){
  ### Summarizing bNTI by taxonomy/groups
  # DNA stuff
  dna.bnti.by.tax = summarize.data(dna.bnti, dna.tax, "Family")
  
  # RNA stuff
  rna.bnti.by.tax = summarize.data(rna.bnti, rna.tax, "Family")
  
  # ICR stuff
  icr.bnti.by.ele = summarize.data(icr.bnti, mol, "El_comp")
  icr.bnti.by.class = summarize.data(icr.bnti, mol, "Class")
  
  ### Doing the plotting
  div.plot.by.factor(dna.bnti.by.tax, dna.factors, time.plot = T)
  div.plot.by.factor(rna.bnti.by.tax, rna.factors, time.plot = T)
  div.plot.by.factor(icr.bnti.by.ele, icr.factors, time.plot = T, adjust.p = F)
  
  rm(dna.bnti.by.tax, rna.bnti.by.tax, icr.bnti.by.ele, icr.bnti.by.class, dna.factors, rna.factors, icr.factors)
}
# ################################ #
#### Correlation-based analyses ####
# ################################ #

### Breaking data up (with awful object names to boot)
# DNA bNTI
dna.bnti.dna = dna.bnti[,grep("^DNA", colnames(dna.bnti))]; dna.bnti.dna = rename.samples(dna.bnti.dna, meta, "gDNA")
dna.bnti.rna = dna.bnti[,grep("cDNA", colnames(dna.bnti))]; dna.bnti.rna = rename.samples(dna.bnti.rna, meta, "cDNA")

# RNA bNTI
rna.bnti.dna = rna.bnti[,grep("^DNA", colnames(rna.bnti))]; rna.bnti.dna = rename.samples(rna.bnti.dna, meta, "gDNA")
rna.bnti.rna = rna.bnti[,grep("cDNA", colnames(rna.bnti))]; rna.bnti.rna = rename.samples(rna.bnti.rna, meta, "cDNA")

# Remove focus and non-ICR samples (for matching purposes)
to.remove = unique(c(meta$SampleName[is.na(meta$SampleID_FTICR)], focus))
dna.bnti.dna = dna.bnti.dna[,-which(colnames(dna.bnti.dna) %in% to.remove)]; dna.bnti.rna = dna.bnti.rna[,-which(colnames(dna.bnti.rna) %in% to.remove)]
rna.bnti.dna = rna.bnti.dna[,-which(colnames(rna.bnti.dna) %in% to.remove)]; rna.bnti.rna = rna.bnti.rna[,-which(colnames(rna.bnti.rna) %in% to.remove)]
meta = meta[-which(meta$SampleName %in% to.remove),]

# Matching minimum sample number
dna.bnti.dna = dna.bnti.dna[,which(colnames(dna.bnti.dna) %in% colnames(rna.bnti.rna))]
rna.bnti.dna = rna.bnti.dna[,which(colnames(rna.bnti.dna) %in% colnames(rna.bnti.rna))]
icr.bnti = icr.bnti[,which(colnames(icr.bnti) %in% colnames(rna.bnti.rna))]
meta = meta[which(meta$SampleName %in% colnames(rna.bnti.rna)),]

rm(to.remove, focus)

### DNA -> ICR correlations
if(spear.cor){
  # Correlations and select significance
  dna.icr.cor = rcorr(as.matrix(t(dna.bnti.dna)), as.matrix(t(icr.bnti)), type = "spearman") # Correlation matrix
  # dna.icr.cor$r = dna.icr.cor$r[which(row.names(dna.icr.cor$r) %in% row.names(dna.bnti.dna)),
  #                               which(colnames(dna.icr.cor$r) %in% row.names(icr.bnti))]
  # dna.icr.cor$P = dna.icr.cor$P[which(row.names(dna.icr.cor$P) %in% row.names(dna.bnti.dna)),
  #                               which(colnames(dna.icr.cor$P) %in% row.names(icr.bnti))]
  
  # Removing duplicate comparisons
  dna.icr.cor$r[upper.tri(dna.icr.cor$r, diag = T)] = NA
  dna.icr.cor$P[upper.tri(dna.icr.cor$P, diag = T)] = NA
  
  # Adjusting p-values
  cor.p = melt(dna.icr.cor$P)
  cor.p = cor.p[!is.na(cor.p$value),]
  # cor.p$value = p.adjust(cor.p$value, method = "fdr")
  
  # Removing insignificant correlations
  cor.r = melt(dna.icr.cor$r); colnames(cor.r)[1:2] = c("Members", "Mass") # Melted and named
  cor.r = cor.r[!is.na(cor.r$value),]
  cor.r = cor.r[which(cor.p$value < 0.0001),]
  cor.r$Mass = as.character(cor.r$Mass)
  
  # Adding in taxonomic/molecular data
  cor.r = cor.r %>% left_join(dna.tax[,c("Members", "Order")], by = "Members") %>% 
    left_join(mol[,c("Mass", "El_comp", "bs1_class")], by = "Mass")
  
  # Summarizing data
  avg.cor = cor.r %>% group_by(Order, El_comp) %>% summarise(mean_rho = mean(value), median_rho = median(value), rho_count = n())
  avg.cor %>% group_by(Order) %>% filter(max(rho_count) > 500) %>% 
    ggplot(aes(Order, median_rho)) + geom_bar(stat = "identity", aes(fill = median_rho)) + facet_grid(El_comp~.) + 
    scale_fill_gradient2(high = "firebrick2", mid = "goldenrod2", low = "dodgerblue2", midpoint = 0) + 
    theme_bw() + plot.theme + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  ### Testing network analysis
  d = with(cor.r, data.frame(to = Members, from = Mass, weight = value))
  vertices = with(d, data.frame(id = unique(c(to, from)), type = "Microbes"), stringsAsFactors = F)
  vertices$type[which(vertices$id %in% row.names(icr.bnti))] = "Metabolites"
  
  net = graph_from_data_frame(d = d, vertices = vertices, directed = F)
  plot(net)
  
  # write.csv(d, "DNA_Edge_Object.csv", quote = F, row.names = F)
  # write.csv(vertices, "DNA_Node_Object.csv", quote = F, row.names = F)
}

### RNA -> ICR comparison
if(spear.cor){
  # Correlations and select significance
  rna.icr.cor = rcorr(as.matrix(t(rna.bnti.rna)), as.matrix(t(icr.bnti)), type = "spearman") # Correlation matrix
  rna.icr.cor$r = rna.icr.cor$r[which(row.names(rna.icr.cor$r) %in% row.names(dna.bnti.dna)),
                                which(colnames(rna.icr.cor$r) %in% row.names(icr.bnti))]
  rna.icr.cor$P = rna.icr.cor$P[which(row.names(rna.icr.cor$P) %in% row.names(dna.bnti.dna)),
                                which(colnames(rna.icr.cor$P) %in% row.names(icr.bnti))]
  
  # Adjusting p-values
  cor.p = melt(rna.icr.cor$P)
  # cor.p$value = p.adjust(cor.p$value, method = "fdr")
  
  # Removing insignificant correlations
  cor.r = melt(rna.icr.cor$r); colnames(cor.r)[1:2] = c("Members", "Mass") # Melted and named
  cor.r = cor.r[which(cor.p$value < 0.001),]
  cor.r$Mass = as.character(cor.r$Mass)
  
  # Adding in taxonomic/molecular data
  cor.r = cor.r %>% left_join(rna.tax[,c("Members", "Order")], by = "Members") %>% 
    left_join(mol[,c("Mass", "El_comp", "bs1_class")], by = "Mass")
  
  # Summarizing data
  avg.cor = cor.r %>% group_by(Order, El_comp) %>% summarise(mean_rho = mean(value), median_rho = median(value), rho_count = n())
  avg.cor %>% group_by(Order) %>% 
    ggplot(aes(Order, median_rho)) + geom_bar(stat = "identity", aes(fill = median_rho)) + facet_grid(El_comp~.) + 
    scale_fill_gradient2(high = "firebrick2", mid = "goldenrod2", low = "dodgerblue2", midpoint = 0) + 
    theme_bw() + plot.theme + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  ### Testing network analysis
  d = with(cor.r, data.frame(to = Members, from = Mass, weight = value))
  vertices = with(cor.r, data.frame(id = unique(Members), type = "Microbes"))
  vertices$type[which(vertices$id %in% row.names(icr.bnti))] = "Metabolites"
  
  net = graph_from_data_frame(d = d, vertices = vertices, directed = F)
  plot(net)
  
  # write.csv(d, "RNA_Edge_Object.csv", quote = F, row.names = F)
  # write.csv(vertices, "RNA_Node_Object.csv", quote = F, row.names = F)
}

# ########### #
#### WGCNA ####
# ########### #

if(WGCNA){ # Switch controlling whether WCGNA will be run or not
  input = t(rbind(rna.bnti.rna, icr.bnti))
  thresh = 0.1
  sign = "unsigned"
  
  powers = c(1:25) # Generates a series of powers to test for future soft thresholding; soft thresholding allows for continuous data to be described
  sft = pickSoftThreshold(input, powerVector = powers, verbose = 5) # Tests the powers against the data to see which power fits the best
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2]) # Plots the test to see if they worked (where it flattens is where it 'worked')
  
  p = -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
  p = sft$fitIndices[,1][which(p %in% max(-sign(sft$fitIndices[,3])*sft$fitIndices[,2]))] # Automatically generating the most 'optimal' power
  
  
  block = blockwiseModules(input, power = p, TOMType = sign, minModuleSize = 5,
                           numericLabels = TRUE) # Assigns the data to modules by transforming the data with the calculated power
  
  moduleLabels = block$colors # Saving labels from the modules
  moduleColors = labels2colors(block$colors) # Saves colors of the modules
  MEs = block$MEs # Saving module eigenvalues which contain the information used to designate modules
  taxTree = block$dendrograms[[1]] # Saving the dendrograms which were generated from TOM data created during module assignment
  
  plotDendroAndColors(taxTree, moduleColors[block$blockGenes[[1]]], "Module colors",
                      dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05) # Plots modules and dendrograms
  
  plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", plotDendrograms = FALSE, xLabelsAngle = 90) # Plots relationships between eigenvectors to see what modules are most closely related
  
  plot(hclust(as.dist(1-cor(MEs)), method = "average")) # Just a plot which can help to see if any module merging can occur
  
  # TOM calculations
  TOM = TOMsimilarity(adjacency(input, power = p, type = sign)) # Calculates the "Topographcial Overlap Matrix"; a measure of how 'related' two nodes are, what their connections are, how adjacent, etc.
  dissTOM = 1 - TOM # Changes similarity matrix to dissimilarity matrix
  plotTOM = dissTOM^7 # Raises the dis. matrix to a power in order to bring out the effects of 'moderate' effects
  diag(plotTOM) = NA # Sets the diagonal of the data to NA's so that there isn't any weird scaling (they are normally one due to self-comparisons)
  
  # # The TOM plot isn't really necessary
  # TOMplot(plotTOM, taxTree, moduleColors) # Generates the plot which looks at relationships between each constituent
  
  # Meta-data for exporting
  w = apply(input, 2, median)
  w = cbind(w, moduleColors)
  
  dimnames(TOM) = list(colnames(input), colnames(input)) # Sets the dimension names to over this matrix to the taxa it was derived from
  cyt = exportNetworkToCytoscape(TOM, edgeFile = paste0("Cytoscape Files/RNA-ICR_Edges_", thresh, "thresh_", sign, ".txt"), 
                                 nodeFile = paste0("Cytoscape Files/RNA-ICR_Nodes_", thresh, "thresh_", sign, ".txt"), 
                                 threshold = thresh, nodeAttr = w) # Threshold is the minimum allowable value for in the TOM matrix
  
  # Creating a module list
  colors = data.frame(Members = colnames(input), colors = moduleColors, hexColor = col2hex(moduleColors), type = "Microbes")
  colors$type[which(colors$Members %in% row.names(icr.bnti))] = "Metabolites"
  color.mol = mol[,c("Mass", "El_comp", "bs1_class", "bs3_class", "MolForm", "NOSC", "AI_Mod", "NtoC_ratio")]; colnames(color.mol)[1] = "Members"
  color.tax = data.frame(Members = row.names(tax), Order = tax$Order)
  colors = colors %>% left_join(color.tax, by = "Members") %>% left_join(color.mol, by = "Members")
  write.csv(colors, paste0("Cytoscape Files/RNA-ICR_Modules_", thresh, "thresh_", sign, ".csv"), row.names = F, quote = F)
}
