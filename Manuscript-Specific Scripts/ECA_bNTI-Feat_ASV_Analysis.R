### Analyzing ASV bNTI Feature

# Load in libraries
library(ggplot2); library(reshape2); library(ggpubr); library(ggtree); library(ggstance)
library(plyr); library(dplyr); library(stringr); library(tidyr)
library(vegan)
library(picante)

# Switch
fill = "Alt"


# ###################### #
#### Define Functions ####
# ###################### #

# Theme elements used in each plot
plot.theme = theme(axis.text = element_text(colour = "black", size = 12),
                   axis.title = element_text(colour = "black", size = 14),
                   panel.border = element_rect(size = 1, colour = "black"),
                   panel.grid = element_blank())

# Plot bNTI by taxa/taxa by bNTI
taxa.feat.plot = function(feat, abund, taxa, taxa.col, suppress.plot){
  
  if(missing(suppress.plot)){
    suppress.plot = F
  }
  
  ### Average relative abundance by bNTI feature classification
  feat.samp.frame = data.frame(Member = row.names(feat), Value = rowMeans(feat, na.rm = T), Direction = "Insignificant", stringsAsFactors = F)
  feat.samp.frame$Direction[which(feat.samp.frame$Value <= -1)] = "Low"
  feat.samp.frame$Direction[which(feat.samp.frame$Value >= 1)] = "High"
  feat.samp.frame$Direction[which(feat.samp.frame$Value <= -2)] = "Sig. Low"
  feat.samp.frame$Direction[which(feat.samp.frame$Value >= 2)] = "Sig. High"
  feat.samp.frame$Direction = factor(feat.samp.frame$Direction, levels = c("Sig. High", "High", "Insignificant", "Low", "Sig. Low"))
  
  # Adding extra data
  feat.samp.frame$Abundance = rowSums(abund)
  feat.samp.frame = feat.samp.frame %>% group_by(Direction) %>% mutate(RelAbund = Abundance/sum(Abundance))
  feat.samp.frame = cbind(feat.samp.frame, Taxonomy = taxa[,taxa.col])
  
  if(!suppress.plot){
    # Plotting stacked bar plot
    print(
      ggplot(feat.samp.frame, aes(x = Direction, y = RelAbund))+
        geom_bar(stat = "identity", aes(fill = Taxonomy))+
        ggtitle(deparse(substitute(feat))) + theme_bw() + plot.theme
    )
    
    # Boxplot to compare bNTI feat by taxa
    print(
      ggplot(feat.samp.frame, aes(x = Taxonomy, y = Value))+
        geom_boxplot(aes(color = Taxonomy)) + xlab(NULL) + ylab("bNTI-feature") + 
        geom_hline(yintercept = c(-2,2), lty = 2, color = "red") + 
        geom_hline(yintercept = c(-1,1), lty = 2, color = "blue") + ggtitle(deparse(substitute(feat))) + 
        theme_bw() + plot.theme + theme(legend.position = "none",
                                        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    )
    
    # Average barchart by taxonomy
    print(
      feat.samp.frame %>% group_by(Taxonomy) %>% summarise(Average_bNTI = mean(Value)) %>%
        ggplot(aes(x = Taxonomy, y = Average_bNTI)) + geom_bar(stat = "identity", aes(fill = Taxonomy)) + xlab(NULL) + 
        ylab("bNTI-feature") + geom_hline(yintercept = c(-2,2), lty = 2, color = "red") + 
        geom_hline(yintercept = c(-1,1), lty = 2, color = "blue") + ggtitle(deparse(substitute(feat))) + 
        theme_bw() + plot.theme + theme(legend.position = "none",
                                        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    )
  }
  
  return(feat.samp.frame)
}

# Parsing data based on group
parse.by.group = function(feat, data, meta, tax, tax.col, factor, level){
  # Identify column with the factor and select the appropriate level
  fact.col = which(colnames(meta) %in% factor)
  fact.loc = which(meta[,fact.col] %in% level)
  feat.fact = feat[,fact.loc]
  data.fact = data[,fact.loc]
  
  # Remove missing formula
  tax.fact = tax[-which(rowSums(is.na(feat.fact)) == ncol(feat.fact)),]
  data.fact = data[-which(rowSums(is.na(feat.fact)) == ncol(feat.fact)),]
  feat.fact = feat.fact[-which(rowSums(is.na(feat.fact)) == ncol(feat.fact)),]
  
  # Create and return new matrix
  fact.frame = data.frame(Members = row.names(feat.fact), Average = rowMeans(feat.fact, na.rm = T), Direction = "Insignificant")
  fact.frame$Direction[which(fact.frame$Average <= -1)] = "Low"
  fact.frame$Direction[which(fact.frame$Average >= 1)] = "High"
  fact.frame$Direction[which(fact.frame$Average <= -2)] = "Sig. Low"
  fact.frame$Direction[which(fact.frame$Average >= 2)] = "Sig. High"
  fact.frame$Direction = factor(fact.frame$Direction, levels = c("Sig. High", "High", "Insignificant", "Low", "Sig. Low"))
  
  fact.frame$Abundance = rowSums(data.fact)
  fact.frame = fact.frame %>% group_by(Direction) %>% mutate(RelAbund = Abundance/sum(Abundance))
  fact.frame = cbind(fact.frame, Taxonomy = tax.fact[,tax.col])
  
  fact.frame$Type = level
  
  return(fact.frame)
}

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


# ############################### #
#### Loading and cleaning data ####
# ############################### #

# Load in data
setwd("~/Documents/bNTI Feature Manuscript/16S ASV Data/")
data = read.table("ECA_feature-table-with-taxonomy.tsv", skip = 1, comment.char = "", sep = "\t", row.names = 1, header = T)
meta = read.csv("ASV_Metadata.csv")

all.bnti = read.csv("bNTI_Feat_ASV_withConsp_bNTI_feature_by_samp_999rep.csv", row.names = 1)
dna.bnti = read.csv("bNTI_Feat_ASV_withConsp_bNTI_feature_by_samp_Type-DNA_999rep.csv", row.names = 1)
rna.bnti = read.csv("bNTI_Feat_ASV_withConsp_bNTI_feature_by_samp_Type-RNA_999rep.csv", row.names = 1)
tree = read.tree("ECA_tree.nwk")

# Removing taxonomy from data and creating a taxonomy object
tax = data$taxonomy
data = data[,-which(colnames(data) %in% "taxonomy")]

tax = strsplit(tax, "; ")
tax = ldply(tax, rbind)
colnames(tax) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
row.names(tax) = row.names(data)

if(fill == "Copy"){
  for(i in 2:ncol(tax)){
    tax[is.na(tax[,i]), i] = tax[is.na(tax[,i]), (i-1)]
  }
} else if(fill == "Alt"){
  tax$Phylum[is.na(tax$Phylum)] = "p__?"
  tax$Class[is.na(tax$Class)] = "c__?"
  tax$Order[is.na(tax$Order)] = "o__?"
  tax$Family[is.na(tax$Family)] = "f__?"
  tax$Genus[is.na(tax$Genus)] = "g__?"
  tax$Species[is.na(tax$Species)] = "s__?"
  
  for(i in 2:ncol(tax)){
    tax[,i] = paste0(tax[,(i-1)], "; ", tax[,i])
  }
}
  
rm(i)

# Filtering data
data = data[,-which(colSums(data) < 10000)] # Remove low sequence samples
data = data[-which(rowSums(data) == 0),]

# Running phylo match
phylo = match.phylo.data(tree, data) # Matching ICR dataset to the tree

data = phylo$data
tree = phylo$phy

rm("phylo")

# Matching datasets
if(!identical(row.names(all.bnti), tree$tip.label)){
  stop("Either the wrong feature or tree file were loaded.")
}

tax = tax[which(row.names(tax) %in% row.names(all.bnti)),]
tax = tax[order(match(row.names(tax), row.names(all.bnti))),]

# Setting taxa below a certain threshold as low abundance
tax.data = as.data.frame(apply(data, 2, function(x) (x/sum(x))*100))
tax.col = which(colnames(tax) %in% "Family")

add.tax = tax[,tax.col]
tax.data$taxonomy = add.tax
tax.data = aggregate(.~taxonomy, tax.data, FUN = sum)

add.tax = tax.data$taxonomy
tax.data = tax.data[,-1]
low.tax = add.tax[which(apply(tax.data, 1, max) < 0.1)]

tax[which(tax[,tax.col] %in% low.tax),tax.col] = "Low Abundance"

rm(tax.data, add.tax, low.tax)

# Matching metadata
meta = meta[which(meta$Sample.ID %in% colnames(data)),]
data = data[,which(colnames(data) %in% meta$Sample.ID)]
all.bnti = all.bnti[,which(colnames(all.bnti) %in% colnames(data))]

# Cleaning up some missing data
all.bnti = all.bnti[-which(rowSums(data) == 0),]
tax = tax[-which(rowSums(data) == 0),]
data = data[-which(rowSums(data) == 0),]

# Creating DNA/RNA data
dna.data = data[,grep("^DNA", colnames(data))]
dna.data = dna.data[-which(rowSums(dna.data) == 0),]
dna.tax = tax[row.names(dna.data),]

rna.data = data[,grep("cDNA", colnames(data))]
rna.data = rna.data[-which(rowSums(rna.data) == 0),]
rna.tax = tax[row.names(rna.data),]


# ############################## #
#### Basic bNTI feature plots ####
# ############################## #

### Stacked barplots, bNTI boxplots, and bNTI barplots
all.feat.frame = taxa.feat.plot(all.bnti, data, tax, tax.col, suppress.plot = T)
dna.feat.frame = taxa.feat.plot(dna.bnti, dna.data, dna.tax, tax.col, suppress.plot = T)
rna.feat.frame = taxa.feat.plot(rna.bnti, rna.data, rna.tax, tax.col, suppress.plot = T)

### Combining data
# Some data rearranging
all.feat.frame = with(all.feat.frame, data.frame(Member = Member, Taxonomy = Taxonomy, 
                                                 All = Value))
dna.feat.frame = with(dna.feat.frame, data.frame(Member = Member, DNA = Value))
rna.feat.frame = with(rna.feat.frame, data.frame(Member = Member, RNA = Value))

# Merging
combined.feat = all.feat.frame %>% left_join(dna.feat.frame, by = "Member") %>% left_join(rna.feat.frame, by = "Member")

# Melting
combined.melt = melt(combined.feat, id.vars = c("Member", "Taxonomy"))
combined.melt = combined.melt %>% filter(!is.na(value))

# Plotting combined boxplot
plot1 = combined.melt %>% 
  group_by(Taxonomy, variable) %>% 
  mutate(mean_value = mean(value, na.rm = T)) %>%
  group_by(Taxonomy) %>%
  filter(max(abs(mean_value)) > 2) %>%
  filter(median(abs(mean_value)) > 1) %>%
  filter(!is.na(value)) %>%
  mutate(Trunc_Tax = str_trunc(Taxonomy, 40, side = "left")) %>%
  ggplot(aes(x = Trunc_Tax, y = value))+
  geom_boxplot(aes(fill = variable)) +  xlab(NULL) + ylab("bNTI-feature") + 
  geom_hline(yintercept = c(-2,2), lty = 2, color = "red") + 
  geom_hline(yintercept = c(-1,1), lty = 2, color = "blue") + theme_bw() +
  scale_fill_manual(values = c("goldenrod1", "purple", "forestgreen"))+
  plot.theme + theme(legend.position = c(0.8, 0.8),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Selecting taxonomies and ploting average relative abundances
selected.tax = combined.melt %>% 
  group_by(Taxonomy, variable) %>% 
  mutate(mean_value = mean(value, na.rm = T)) %>%
  group_by(Taxonomy) %>%
  filter(max(abs(mean_value)) > 2) %>%
  filter(median(abs(mean_value)) > 1) %>%
  select(Taxonomy) # Identify taxa from plot 1

rel.data = apply(data, 2, function(x) (x/sum(x))*100)
selected.data = as.data.frame(rel.data[which(tax$Family %in% selected.tax$Taxonomy),])
selected.data$Taxonomy = tax$Family[which(tax$Family %in% selected.tax$Taxonomy)]

selected.data %>% group_by(Taxonomy) %>%
  dplyr::summarise_all(.funs = sum) %>%
  gather(Sample, Rel_Abund, -Taxonomy) %>%
  mutate(Trun_Tax = str_trunc(Taxonomy, 40, side = "left")) %>%
  ggplot(aes(x = Trun_Tax, y = Rel_Abund)) +
  geom_boxplot() + xlab("Taxonomy") + ylab("Relative Abundance (%)") +
  theme_bw() + plot.theme + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# Plotting combined barplot
plot2 = combined.melt %>% group_by(Taxonomy, variable) %>% 
  dplyr::summarise(mean_value = mean(value, na.rm = T)) %>% 
  group_by(Taxonomy) %>%
  filter(max(abs(mean_value)) > 2) %>% 
  filter(median(abs(mean_value)) > 1) %>%
  ggplot(aes(x = Taxonomy, y = mean_value))+
  geom_bar(stat = "identity", aes(fill = Taxonomy)) + facet_wrap(.~variable, ncol = 1)+
  xlab(NULL) + ylab("bNTI-feature") + 
  geom_hline(yintercept = c(-2,2), lty = 2, color = "red") + 
  geom_hline(yintercept = 0, color = "black") + theme_bw() + 
  plot.theme + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Alternative barplot
plot3 = combined.melt %>% group_by(Taxonomy, variable) %>% 
  dplyr::summarise(mean_value = mean(value, na.rm = T)) %>% 
  group_by(Taxonomy) %>%
  filter(max(abs(mean_value)) > 2) %>% 
  filter(median(abs(mean_value)) > 1) %>%
  mutate(Trunc_Tax = str_trunc(Taxonomy, 40, side = "left")) %>%
  ggplot(aes(x = Trunc_Tax, y = mean_value))+
  geom_bar(stat = "identity", position = "dodge", aes(fill = variable))+
  xlab(NULL) + ylab("bNTI-feature") + 
  geom_hline(yintercept = c(-2,2), lty = 2, color = "red") + 
  geom_hline(yintercept = 0, color = "black") + theme_bw() + 
  plot.theme + theme(legend.position = c(0.8, 0.2), 
                     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Combining plots
print(
  ggarrange(plot1, plot3, labels = c("A", "B"), ncol = 1 )
)

### Plotting proportions by taxa
# Assigning direction
combined.melt$Direction = make.direction(combined.melt, 4)

# Difference based on proportion
taxa.diff = combined.melt %>% group_by(Taxonomy, Direction, variable) %>% 
  summarise(Count = n()) %>% 
  group_by(Taxonomy, variable) %>%
  summarise(Direction = Direction, Rel_Prop = (Count/sum(Count)*100)) %>%
  filter(Direction == "Insignificant") %>%
  filter(variable %in% c("DNA", "RNA")) %>% group_by(Taxonomy) %>% 
  spread(variable, Rel_Prop) %>% replace(is.na(.), 0) %>% 
  summarise(diff = DNA - RNA) # Relative proportion differences of insignificant values

prop.diff = c(taxa.diff$Taxonomy[order(taxa.diff$diff)][1:8],
              taxa.diff$Taxonomy[order(taxa.diff$diff, decreasing = T)][1:8]) # Selecting the largest differences

# Difference based on average
taxa.diff = combined.melt %>% group_by(Taxonomy, variable) %>% 
  summarise(mean_value = mean(value, na.rm = T)) %>%
  spread(variable, mean_value) %>% summarise(diff = abs(DNA) - abs(RNA))

mean.diff = c(taxa.diff$Taxonomy[order(taxa.diff$diff)][1:8],
              taxa.diff$Taxonomy[order(taxa.diff$diff, decreasing = T)][1:8]) # Selecting the largest differences

# Doing the plotting things
plot2 = combined.melt %>% group_by(Taxonomy, Direction, variable) %>% 
  summarise(Count = n()) %>% 
  group_by(Taxonomy, variable) %>%
  summarise(Direction = Direction, Rel_Prop = (Count/sum(Count)*100)) %>% 
  filter(Taxonomy %in% prop.diff) %>%
  ggplot(aes(variable, Rel_Prop)) + facet_wrap(.~Taxonomy)+
  geom_bar(stat = "identity", aes(fill = Direction)) + 
  xlab(NULL) + ylab("Relative Proporiton (%)") +
  scale_fill_manual(values = c("firebrick4", "firebrick1", "goldenrod1", 
                               "dodgerblue1", "dodgerblue4"), drop = F) +
  theme_bw() + plot.theme + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

plot3 = ggarrange(plot1, plot2)
plot3

rm(combined.melt, taxa.diff, mean.diff, prop.diff)


# ############################ #
### Creating a combined tree ###
# ############################ #

# Adding direction to individual feat frames
all.feat.frame$Direction = make.direction(all.feat.frame, 3)
dna.feat.frame$Direction = make.direction(dna.feat.frame, 2)
rna.feat.frame$Direction = make.direction(rna.feat.frame, 2)

# Creating difference object
combined.diff = combined.feat
combined.diff[is.na(combined.diff)] = 0
combined.diff$DNA = scale(combined.diff$DNA)
combined.diff = with(combined.diff, data.frame(Member = Member, Difference = c(abs(DNA) - abs(RNA)), Favored = "DNA favored"))
combined.diff$Favored[combined.diff$Difference < 0] = "RNA favored"
combined.diff$Favored[combined.diff$Difference == 0] = "Neither"
combined.diff$Favored = factor(combined.diff$Favored, levels = c("DNA favored", "Neither", "RNA favored"))

if(identical(row.names(tax), combined.diff$Member)){
  combined.diff$Class = tax$Class
} else {
  stop()
}

# Plotting tree
p = ggtree(tree)
p = facet_plot(p, panel = "Overall bNTI-feature", data = all.feat.frame, geom = geom_point,
               mapping = aes(x = All, y = y, color = as.factor(Direction)))
p = facet_plot(p, panel = "DNA vs. RNA", data = combined.diff, geom = geom_barh,
               mapping = aes(x = Difference, fill = as.factor(Favored)), stat = 'identity')
# p = facet_plot(p, panel = "This is for text", data = combined.diff, geom = geom_label,
#                mapping = aes(x = Difference, y = y, label = Class))
p + scale_color_manual(values = c("firebrick4", "firebrick1", "goldenrod1", "dodgerblue1", "dodgerblue4"), drop = F) + 
  scale_fill_manual(values = c("purple", "goldenrod1", "forestgreen"), drop = F) + labs(color = "Significance", fill = "Type Pref.")+
  theme_tree2()+theme(legend.position = "left")

# ggsave("Large_Tree_Figure_w_Labels.pdf", device = "pdf",
#        width = 25, height = 2000, limitsize = F)
