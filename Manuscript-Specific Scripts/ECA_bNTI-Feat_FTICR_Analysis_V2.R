### Analyzing FTICR bNTI Feature

# Load in libraries
library(ggplot2); library(reshape2); library(ggpubr); library(ggstance); library(ggtree)
library(plyr); library(dplyr); library(tidyr); library(stringr)
library(vegan)
library(picante)

# ###################### #
#### Define Functions ####
# ###################### #

plot.theme = theme(axis.text = element_text(colour = "black", size = 12),
                   axis.title = element_text(colour = "black", size = 14),
                   panel.border = element_rect(size = 1, colour = "black"),
                   panel.grid = element_blank())

# Add direction to a frame
make.direction = function(feat.frame){
  
  # Value column
  val.col = grep("value", colnames(feat.frame), ignore.case = T)
  
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

# Creating standard plots
group.feat.plot = function(feat, mol, suppress.plot){
  
  if(missing(suppress.plot)){
    suppress.plot = F
  }
  
  # Density plot
  feat.melt = melt(as.matrix(feat))
  feat.melt = feat.melt[!is.na(feat.melt$value),]
  
  ggplot(feat.melt, aes(x = value))+
    geom_density(fill = "steelblue2")+xlab("bNTI")+
    geom_vline(xintercept = c(-2,2), lty = 2, color = "red")+
    geom_vline(xintercept = c(-1,1), lty = 2, color = "blue")+
    theme_bw() + plot.theme
  
  # Average bNTI-Feat dendrogram
  feat.samp.frame = data.frame(Member = row.names(feat), Value = rowMeans(feat, na.rm = T), stringsAsFactors = F)
  feat.samp.frame$Direction = make.direction(feat.samp.frame)
  
  feat.samp.frame = cbind(feat.samp.frame, mol)
  
  if(!suppress.plot){
    # Plotting a standard tree with a point graph to see the scale of differences
    # p = ggtree(tree)
    # p = facet_plot(p, panel = "Nearest Neighbor", data = feat.samp.frame, 
    #                geom = geom_point, mapping = aes(x = Value, y = y, color = Direction))
    # p = facet_plot(p, panel = "AI_Mod", data = feat.samp.frame, 
    #                geom = geom_point, mapping = aes(x = AI_Mod, y = y))
    # p = facet_plot(p, panel = "DBE", data = feat.samp.frame, 
    #                geom = geom_point, mapping = aes(x = DBE, y = y))
    # p = facet_plot(p, panel = "NOSC", data = feat.samp.frame, 
    #                geom = geom_point, mapping = aes(x = NOSC, y = y))
    # p = facet_plot(p, panel = "N:C", data = feat.samp.frame, 
    #                geom = geom_point, mapping = aes(x = NtoC_ratio, y = y))
    # 
    # print(
    #   p + scale_color_manual(values = c("firebrick4", "firebrick1", "goldenrod1", "dodgerblue1", "dodgerblue4"), drop = F) + 
    #     theme_tree2()+theme(legend.position = "left")
    # )
    
    # Boxplot to compare characteristics by direction
    feat.melt = feat.samp.frame[,c("Direction", "AI_Mod", "DBE", "NOSC", "OtoC_ratio", "HtoC_ratio", "NtoC_ratio")]
    feat.melt = melt(feat.melt, id.vars = "Direction")
    
    print(
      ggplot(data = feat.melt, aes(x = Direction, y = value))+
        geom_boxplot()+facet_wrap(.~variable, ncol = 1, scales = "free_y")+
        xlab(NULL) + theme_bw() + plot.theme
    )
    
    ### Elemental Composition plots
    # Analyzing elemental composition by direction
    feat.melt = feat.samp.frame %>% group_by(Direction) %>% count(El_comp) %>% mutate(RelAbund = (n/sum(n)*100))
    
    plot1 = ggplot(data = feat.melt, aes(x = Direction, y = RelAbund))+
      geom_bar(stat = "identity", aes(fill = El_comp))+
      xlab(NULL) + ylab("Relative Abundance (%)") + theme_bw() + plot.theme
    
    # Boxplot to compare bNTI feat by El_comp
    plot2 = ggplot(feat.samp.frame, aes(x = El_comp, y = Value))+
      geom_boxplot(aes(color = El_comp)) + xlab(NULL) + ylab("bNTI-feature") + 
      geom_hline(yintercept = c(-2,2), lty = 2, color = "red") + 
      geom_hline(yintercept = c(-1,1), lty = 2, color = "blue") + 
      theme_bw() + plot.theme + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    # Merging plots
    print(ggarrange(plot1, plot2, common.legend = T))
    
    ### BS1 plots
    # Analyzing boundary sets by direction
    feat.melt = feat.samp.frame %>% group_by(Direction) %>% count(bs1_class) %>% mutate(RelAbund = (n/sum(n)*100))
    
    plot1 = ggplot(data = feat.melt, aes(x = Direction, y = RelAbund))+
      geom_bar(stat = "identity", aes(fill = bs1_class))+
      xlab(NULL) + ylab("Relative Abundance (%)") + theme_bw() + plot.theme
    
    # Boxplot to compare bNTI feat by El_comp
    plot2 = ggplot(feat.samp.frame, aes(x = bs1_class, y = Value))+
      geom_boxplot(aes(color = bs1_class)) + xlab(NULL) + ylab("bNTI-feature") + 
      geom_hline(yintercept = c(-2,2), lty = 2, color = "red") + 
      geom_hline(yintercept = c(-1,1), lty = 2, color = "blue") + 
      theme_bw() + plot.theme + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    # Merging plots
    print(ggarrange(plot1, plot2, common.legend = T))

  }
  
  return(feat.samp.frame)
}

# Organize data by group (may be outdated)
parse.by.group = function(feat, mol, meta, factor, level){
  # Identify column with the factor and select the appropriate level
  fact.col = which(colnames(meta) %in% factor)
  fact.loc = which(meta[,fact.col] %in% level)
  feat.fact = feat[,fact.loc]
  
  # Remove missing formula
  mol.fact = mol[-which(rowSums(feat.fact, na.rm = T) == 0),]
  feat.fact = feat.fact[-which(rowSums(feat.fact, na.rm = T) == 0),]
  
  # Create and return new matrix
  fact.frame = data.frame(Members = row.names(feat.fact), Average = rowMeans(feat.fact, na.rm = T), Direction = "Insignificant")
  fact.frame$Direction[which(fact.frame$Average <= -1)] = "Low"
  fact.frame$Direction[which(fact.frame$Average >= 1)] = "High"
  fact.frame$Direction[which(fact.frame$Average <= -2)] = "Sig. Low"
  fact.frame$Direction[which(fact.frame$Average >= 2)] = "Sig. High"
  fact.frame$Direction = factor(fact.frame$Direction, levels = c("Sig. High", "High", "Insignificant", "Low", "Sig. Low"))
  
  fact.frame = cbind(fact.frame, mol.fact)
  fact.frame$Type = level
  
  return(fact.frame)
}


# ############################### #
#### Loading and cleaning data ####
# ############################### #

# Load in data
setwd("~/Documents/bNTI Feature Manuscript/FTICR Data/")
data = read.csv("Processed_ECA_8ppm_Data.csv", row.names = 1)
mol = read.csv("Processed_ECA_8ppm_Mol.csv", row.names = 1)
meta = read.csv("FTICR_Metadata.csv")

all.bnti = read.csv("bNTI_Feat_ICR_withConsp_bNTI_feature_by_samp_999rep.csv", row.names = 1)
dry.bnti = read.csv("bNTI_Feat_ICR_withConsp_bNTI_feature_by_samp_Cumulative.Treatment-Dry_999rep.csv", row.names = 1)
inun.bnti = read.csv("bNTI_Feat_ICR_withConsp_bNTI_feature_by_samp_Cumulative.Treatment-Inundated_999rep.csv", row.names = 1)
tree = read.tree("ECA_8ppm_MCD_UPGMA.tre")

# Matching datasets
if(!identical(row.names(all.bnti), tree$tip.label)){
  stop("Either the wrong feature or tree file were loaded.")
}

data = data[which(row.names(data) %in% row.names(all.bnti)),] # Matching data/mol to the feature data
mol = mol[which(row.names(mol) %in% row.names(all.bnti)),]

all.bnti = all.bnti[,meta$Sample.ID] # Matching column order
data = data[,meta$Sample.ID]

data = data[row.names(all.bnti),] # Matching row order
mol = mol[row.names(all.bnti),]

# Some data cleanup
mol$bs1_class[grep(";", mol$bs1_class)] = "Multi-Class"

# Creating objects for groups
dry.mol = mol[row.names(dry.bnti),]
inun.mol = mol[row.names(inun.bnti),]


# ############################### #
#### Generating standard plots ####
# ############################### #

### Stacked barplots, bNTI boxplots, and bNTI barplots
all.feat.frame = group.feat.plot(all.bnti, mol, suppress.plot = T)
dry.feat.frame = group.feat.plot(dry.bnti, dry.mol, suppress.plot = T)
inun.feat.frame = group.feat.plot(inun.bnti, inun.mol, suppress.plot = T)

### Combining data
# Some data rearranging
colnames(all.feat.frame)[2] = "All"
dry.feat.frame = with(dry.feat.frame, data.frame(Member = Member, Dry = Value))
inun.feat.frame = with(inun.feat.frame, data.frame(Member = Member, Inundated = Value))

# Merging
combined.feat = all.feat.frame %>% left_join(dry.feat.frame, by = "Member") %>% left_join(inun.feat.frame, by = "Member")

# More rearranging
combined.feat = melt(combined.feat, id.vars = colnames(all.feat.frame)[-2])
combined.feat = combined.feat[!is.na(combined.feat$value),]

# Combined elem. comp. boxplot
plot1 = ggplot(combined.feat, aes(x = El_comp, y = value))+
  geom_boxplot(aes(fill = variable)) +  xlab(NULL) + ylab("bNTI-feature") + 
  geom_hline(yintercept = c(-2,2), lty = 2, color = "red") + 
  geom_hline(yintercept = c(-1,1), lty = 2, color = "blue") + theme_bw() + 
  plot.theme + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Combined bs1 boxplot
plot2 = ggplot(combined.feat, aes(x = bs1_class, y = value))+
  geom_boxplot(aes(fill = variable)) +  xlab(NULL) + ylab("bNTI-feature") + 
  geom_hline(yintercept = c(-2,2), lty = 2, color = "red") + 
  geom_hline(yintercept = c(-1,1), lty = 2, color = "blue") + theme_bw() + 
  plot.theme + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Merging boxplots
ggarrange(plot1, plot2)

# Combined elem. comp. barplot
plot1 = combined.feat %>% group_by(El_comp, variable) %>% summarise(mean_value = mean(value, na.rm = T)) %>% 
  ggplot(aes(x = El_comp, y = mean_value))+
  geom_bar(stat = "identity", aes(fill = El_comp)) + facet_wrap(.~variable, ncol = 1)+
  xlab(NULL) + ylab("bNTI-feature") + geom_hline(yintercept = c(-2,2), lty = 2, color = "red") + 
  geom_hline(yintercept = c(-1,1), lty = 2, color = "blue") +
  geom_hline(yintercept = c(0), color = "black") + theme_bw() + 
  plot.theme + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Combined bs1 barplot
plot2 = combined.feat %>% group_by(bs1_class, variable) %>% summarise(mean_value = mean(value, na.rm = T)) %>% 
  ggplot(aes(x = bs1_class, y = mean_value))+
  geom_bar(stat = "identity", aes(fill = bs1_class)) + facet_wrap(.~variable, ncol = 1)+
  xlab(NULL) + ylab("bNTI-feature") + geom_hline(yintercept = c(-2,2), lty = 2, color = "red") + 
  geom_hline(yintercept = c(-1,1), lty = 2, color = "blue") +
  geom_hline(yintercept = c(0), color = "black") + theme_bw() + 
  plot.theme + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggarrange(plot1, plot2)

### Plotting proportions by elem. comp. and bs1
# Assigning direction
combined.feat$Direction = make.direction(combined.feat)

# Doing the plotting things for elem. comp.
combined.feat %>% group_by(El_comp, Direction, variable) %>% 
  summarise(Count = n()) %>% 
  group_by(El_comp, variable) %>%
  mutate(Rel_Prop = (Count/sum(Count)*100)) %>% 
  ggplot(aes(variable, Rel_Prop)) + facet_wrap(.~El_comp)+
  geom_bar(stat = "identity", aes(fill = Direction)) + xlab(NULL) + ylab("Relative Proporiton (%)") +
  scale_fill_manual(values = c("firebrick4", "firebrick1", "goldenrod1", "dodgerblue1", "dodgerblue4"), drop = F) +
  theme_bw() + plot.theme + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Doing the plotting things for bs1
combined.feat %>% group_by(bs1_class, Direction, variable) %>% summarise(Count = n()) %>% group_by(bs1_class, variable) %>%
  mutate(Rel_Prop = (Count/sum(Count)*100)) %>% ggplot(aes(variable, Rel_Prop)) + facet_wrap(.~bs1_class)+
  geom_bar(stat = "identity", aes(fill = Direction)) + xlab(NULL) + ylab("Relative Proporiton (%)") +
  scale_fill_manual(values = c("firebrick4", "firebrick1", "goldenrod1", "dodgerblue1", "dodgerblue4"), drop = F) +
  theme_bw() + plot.theme + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


#### Generating bNTI feature trees
# Creating summary dataframe for plotting on top of trees
all.feat.frame = data.frame(Peaks = row.names(all.bnti), value = rowMeans(all.bnti, na.rm = T))
dry.feat.frame = data.frame(Peaks = row.names(dry.bnti), value = rowMeans(dry.bnti, na.rm = T))
inun.feat.frame = data.frame(Peaks = row.names(inun.bnti), value = rowMeans(inun.bnti, na.rm = T))

# Setting direction
all.feat.frame$Direction = make.direction(all.feat.frame)
dry.feat.frame$Direction = make.direction(dry.feat.frame)
inun.feat.frame$Direction = make.direction(inun.feat.frame)

# Plotting tree with information
p = ggtree(tree)
p = facet_plot(p, panel = "Overall Dataset Nearest Neighbor", data =  all.feat.frame,
               geom = geom_point, mapping = aes(x = value, y = y, color = Direction))
p = facet_plot(p, panel = "Dry Nearest Neighbor", data =  dry.feat.frame,
               geom = geom_point, mapping = aes(x = value, y = y, color = Direction))
p = facet_plot(p, panel = "Inundated Nearest Neighbor", data =  inun.feat.frame,
               geom = geom_point, mapping = aes(x = value, y = y, color = Direction))

p + scale_color_manual(values = c("firebrick4", "firebrick1", "goldenrod1", 
                                  "dodgerblue1", "dodgerblue4"), drop = F) + 
  theme_tree()+theme(legend.position = "left")


# Difference based on proportion
combined.feat = all.feat.frame %>% left_join(inun.feat.frame, by = "Peaks") %>% left_join(dry.feat.frame, by = "Peaks")
combined.feat = combined.feat[,-grep("Direction", colnames(combined.feat))]
colnames(combined.feat) = c("Peaks", "All", "Dry", "Inun")
combined.feat[is.na(combined.feat)] = 0
combined.feat$Dry = scale(combined.feat$Dry)
combined.feat$Inun = scale(combined.feat$Inun)
combined.feat = with(combined.feat, data.frame(Peaks = Peaks, 
                                               Difference = c(abs(Dry) - abs(Inun)), Favored = "Dry favored"))
combined.feat$Favored[combined.feat$Difference < 0] = "Inundated favored"
combined.feat$Favored[combined.feat$Difference == 0] = "Neither"
combined.feat$Favored = factor(combined.feat$Favored, levels = c("Dry favored", "Neither", "Inundated favored"))

mol$Peaks = row.names(mol)
combined.feat = combined.feat %>% left_join(mol, by = "Peaks")

# Plotting tree
p = ggtree(tree)
p = facet_plot(p, panel = "Overall bNTI-feature", data = all.feat.frame, geom = geom_point,
               mapping = aes(x = value, y = y, color = as.factor(Direction)))
p = facet_plot(p, panel = "Dry vs. Inundated", data = combined.feat, geom = geom_barh,
               mapping = aes(x = Difference, fill = as.factor(Favored)), stat = 'identity')
# p = facet_plot(p, panel = "This is for text", data = combined.feat, geom = geom_label,
#                mapping = aes(x = Difference, y = y, label = MolForm))
p + scale_color_manual(values = c("firebrick4", "firebrick1", "goldenrod1", "dodgerblue1", "dodgerblue4"), drop = F) + 
  scale_fill_manual(values = c("purple", "goldenrod1", "forestgreen"), drop = F) + 
  labs(color = "Significance", fill = "Type Pref.")+
  theme_tree2()+theme(legend.position = "left")

# ggsave("Large_Tree_Figure_w_Labels.pdf", device = "pdf",
#        width = 25, height = 2000, limitsize = F)

# Identifying molecular properties for divergent/convergent N/P formulas
all.melt = all.feat.frame %>% left_join(mol, by = "Peaks") %>% filter(Direction %in% c("Sig. High", "Sig. Low")) %>%
  filter(P > 0 | N > 0)
all.melt = melt(all.melt[,c("Direction", "El_comp", "NOSC", "AI_Mod", "DBE")], id.vars = c("Direction", "El_comp"))
all.melt$variable = gsub("AI_Mod", "Modified AI", all.melt$variable)
all.melt$Direction = gsub("High", "Diver.", all.melt$Direction)
all.melt$Direction = gsub("Low", "Conver.", all.melt$Direction)
all.melt$Direction = factor(all.melt$Direction, levels = c("Sig. Diver.", "Sig. Conver."))

ggplot(all.melt, aes(x = Direction, y = value))+
  geom_boxplot(aes(fill = Direction))+
  facet_grid(variable~El_comp, scales = "free_y")+
  stat_compare_means(label = "p.signif", method = "wilcox.test")+
  xlab("Direction") + ylab("Molecular Property Value")+
  scale_fill_manual(values = c("firebrick4", "dodgerblue4"))+
  theme_bw() + plot.theme + theme(legend.position = "none",
                                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

stats = rbind(data.frame(variable = "DBE", El_comp = unique(all.melt$El_comp), W.stat = NA, p.value = NA),
              data.frame(variable = "Modified AI", El_comp = unique(all.melt$El_comp), W.stat = NA, p.value = NA),
              data.frame(variable = "NOSC", El_comp = unique(all.melt$El_comp), W.stat = NA, p.value = NA))

for(i in 1:nrow(stats)){
  temp = all.melt[which(all.melt$variable %in% stats$variable[i] & all.melt$El_comp %in% stats$El_comp[i]), ]
  temp.stat = wilcox.test(value~Direction, data = temp)
  stats$W.stat[i] = temp.stat$statistic
  stats$p.value[i] = temp.stat$p.value
}

stats$p.value = p.adjust(stats$p.value, method = "fdr")
stats$round.p = round(stats$p.value, digits = 5)
