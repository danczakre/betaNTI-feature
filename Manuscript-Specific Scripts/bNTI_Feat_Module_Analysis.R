### Processing bNTI network information
library(reshape2); library(ggplot2); library(ggpubr)
library(forcats); library(tidyr); library(dplyr);
library(vegan)
library(Hmisc)

# Load in data
setwd("~/Documents/bNTI Feature Manuscript/Cytoscape Files/Module Analyses/")
mod = read.csv("RNA-ICR_0.1thresh_Processed_Nodes.csv")

### Define functions
# ggplot theme
plot.theme = theme(axis.text = element_text(colour = "black", size = 12),
                   axis.title = element_text(colour = "black", size = 14),
                   panel.border = element_rect(size = 1, colour = "black"),
                   panel.grid = element_blank())

### Removing parts of modules which were doublets
mod = mod[-which(mod$Degree == 1 & mod$NeighborhoodConnectivity == 1 & mod$NumberOfUndirectedEdges == 1),]

### Characterize each module based upon metabolite/microbe membership
# Determining membership
type.mod = mod %>% select(moduleColors, type) %>% group_by(moduleColors, type) %>% summarise(Count = n()) %>%
  mutate(RelAbund = (Count/sum(Count))*100)
type.mod = type.mod[-which(type.mod$RelAbund == 100),]
type.mod.count = type.mod %>% group_by(moduleColors) %>% summarise(Count = sum(Count), RelAbund = 100)

# Plotting modules membership
type.mod %>% group_by(type) %>% mutate(moduleColors = fct_reorder(moduleColors, RelAbund)) %>% 
  ggplot(aes(moduleColors, RelAbund)) + geom_bar(stat = "identity", aes(fill = type)) +
  geom_text(data = type.mod.count, aes(x = moduleColors, y = 99, label = Count), size = 3) + 
  xlab(NULL) + ylab("Relative Abundance (%)") + theme_bw() + plot.theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

### Characterizing proportions by taxa type
# Taxonomy reshaping
tax.mod = mod[which(mod$type %in% "Microbes"),] %>% group_by(Order, moduleColors) %>% 
  summarise(Count = n()) %>% ungroup() %>% group_by(moduleColors) %>% 
  mutate(RelCount = (Count/sum(Count))*100)
tax.mod = tax.mod[which(tax.mod$moduleColors %in% unique(type.mod$moduleColors)),] %>% 
  left_join(type.mod[which(type.mod$type %in% "Microbes"),c("moduleColors", "RelAbund")], by = "moduleColors") %>%
  mutate(Adj_RelAbund = (RelCount*(RelAbund/100)))

# Creating a sorting object
sort.mod = tax.mod[,c("moduleColors", "RelAbund")]; sort.mod = sort.mod[!duplicated(sort.mod$moduleColors),]
sort.mod = sort.mod[order(sort.mod$RelAbund),] # Need a separate sorting object because things were being silly
tax.mod$moduleColors = factor(tax.mod$moduleColors, levels = sort.mod$moduleColors)

# Plotting stacked bar chart (not recommended)
ggplot(tax.mod, aes(moduleColors, Adj_RelAbund)) + geom_bar(stat = "identity", aes(fill = Order)) + 
  xlab(NULL) + ylab("Relative Abundance (%)") + theme_bw() + plot.theme +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Playing with diversity instead
tax.mod = tax.mod %>% spread(moduleColors, Count, fill = 0) %>% select(-c("RelCount", "RelAbund", "Adj_RelAbund"))
tax.mod = aggregate(.~Order, data = tax.mod, FUN = "sum")
mod.div = data.frame(modules = colnames(tax.mod)[-1], Shannon_H = vegan::diversity(t(tax.mod[,-1]), index = "shannon")) %>% 
  mutate(exp_H = exp(Shannon_H)) %>%
  mutate(Pielou_J = Shannon_H/log(specnumber(t(tax.mod[,-1]))), No_of_Orders = specnumber(t(tax.mod[,-1])),
         Richness = colSums(tax.mod[,-1])) # Diversity indices
mod.div = melt(mod.div, id.vars = "modules")
mod.div$modules = factor(mod.div$modules, levels = sort.mod$moduleColors) # Order based on rel. abund
mod.div = mod.div[!is.nan(mod.div$value),]

# Plotting the diversity indices
div.plot = ggplot(mod.div, aes(modules, value)) + geom_bar(stat = "identity", aes(fill = variable)) + 
  facet_grid(variable~., scales = "free_y") + 
  xlab(NULL) + ylab("Diversity values") + theme_bw() + plot.theme +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

### Characterizing elemental composition of modules
# Rehaping elemental composition
ele.mod = mod[which(mod$type %in% "Metabolites"),] %>% group_by(El_comp, moduleColors) %>% 
  summarise(Count = n()) %>% ungroup() %>% group_by(moduleColors) %>% 
  mutate(RelCount = (Count/sum(Count))*100)
ele.mod = ele.mod[which(ele.mod$moduleColors %in% unique(type.mod$moduleColors)),] %>% 
  left_join(type.mod[which(type.mod$type %in% "Metabolites"),c("moduleColors", "RelAbund")], by = "moduleColors") %>%
  mutate(Adj_RelAbund = (RelCount*(RelAbund/100)))

# Sorting modules
ele.mod$moduleColors =  factor(ele.mod$moduleColors, levels = sort.mod$moduleColors)

# Plotting stacked bar chart (not recommended)
ele.plot = ggplot(ele.mod, aes(moduleColors, RelCount)) + 
  geom_bar(stat = "identity", aes(fill = El_comp)) + facet_grid(El_comp~.) + 
  xlab(NULL) + ylab("Relative Abundance (%)") + theme_bw() + plot.theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Combine plots
ggarrange(div.plot, ele.plot, labels = c("A", "B"))