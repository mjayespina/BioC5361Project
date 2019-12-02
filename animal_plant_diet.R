# BioC5361 Project
# Title: Animal vs. Plant Diet
# Group : Group 2
# filename: animal_plant_diet.R
# Author: Mary Jane Espina, Claire Peichel , Sarah Lichtenberger


#### Data importing #####

# Make sure you all the packages installed for the following libraries

# Libraries
library(tidyverse) # install.packages("tidyverse")--for data wrangling
library(reshape2) # install.packages("reshape2")--to convert data from wide to long format
library(vegan) # install.packages("vegan") -- for making alpha diversity
library(ggplot2) # install.packages("ggplot2") -- to make plots
library(ggpubr) #install.packages ("ggpubr") -- to put together plots in one page
library(ape) # install.packages("ape") -- to make hierarchal clustering
library(dendextend) # install.pacakges ("dendoextend") --to make dendogram


# Import data --> Taxable from Green Genes Database
Taxa <- read.delim("https://knights-lab.github.io/MLRepo/datasets/david/gg/taxatable.txt", header = TRUE, row.names = 1)
dim(Taxa) #get the dimension of data
glimpse(Taxa) #glimpse of the data set


# Import original metadata from David 2014
Meta <- read.delim("https://raw.githubusercontent.com/knights-lab/MLRepo/master/datasets/david/mapping-orig.txt", header = TRUE)

# Subset only relevant information from the metadata
Meta <-select(Meta, X.SampleID, SubjectFood, Diet, Day)



#### Data Processing ####

# Summary Stats for the Data

Sum <- colSums(Taxa) # get the column sums of the Taxa data frame

view(Sum) # view the Sum

summary(Sum) # get summary statistics of the colSums

# Quality control of the total reads
# Only samples with more than 20,000 reads were included

TaxaQC <- Taxa[colSums(Taxa)>= 20045]

# Convert to relative abundance
TaxaNorm <- sweep(TaxaQC, 2, colSums(TaxaQC), FUN = "/")

# Transpose TaxaNorm data to easily merge two files
TaxaT <- t(TaxaNorm) 

row.names(Meta) <- Meta$X.SampleID

Merge <- merge(x=Meta,y=TaxaT,by=0) # merge two data together



#### Plotting alpha diversity ####

## Alpha Diversity

# Compute for shannon-alpha diversity
TaxaAlpha <-  diversity(TaxaNorm,index = "shannon",MARGIN = 2 )
View(TaxaAlpha)


# Merge Alpha and Metadata into one file

AlphaMerge <- merge(x=Meta,y=TaxaAlpha,by=0) %>% # merge the alpha diversity and metadata
  filter(!Diet == "NA")  %>% # filter out diet in the data frame
  select(-Row.names) %>% # remove Row.names since it has the same info as X.SampleID
  rename(Subject=SubjectFood) %>% # rename SubjectFood to Subject
  group_by(Day) %>% # Grooup data by day
  melt(id.vars=c("X.SampleID", 
                 "Subject","Day", "Diet"),
       variable.name=c("y"), value.name=c("Alpha")) # transform data into long format

# Plotting alpha diversity, grouped by Day using ggplot

AlphaPlot1 <- ggplot(AlphaMerge, aes(x=Day, y=Alpha, group=Day)) +
  geom_boxplot() + # make a boxplot
  facet_wrap (~ Diet) + # facet based on diet
  labs(x="Day", y="Alpha Diversity (Shannon)", # some plot aesthetics and formatting
       title = "Within-species alpha diversity by day for animal vs plant diet",
       subtitle = "There is decline in within-species alpha diversity on the 4th day of plant-based diet.") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 10, color = "darkslategrey", margin = margin(b = 25)),
        plot.caption = element_text(size = 8, margin = margin(t = 10), color = "grey70", hjust = 0))

AlphaPlot1

# Plotting alpha diversity, grouped by Subject using ggplot

AlphaPlot2 <- ggplot(AlphaMerge, aes(x=reorder(Subject, -Alpha), y=Alpha, fill=Subject)) +
  geom_boxplot() + # make a boxplot
  facet_wrap (~ Diet) + # facet based on diet
  labs(x="Subject", y="Alpha Diversity (Shannon)", # some plot aesthetics and formatting
       title = "Within-species alpha diversity by Subject for animal vs plant diet",
       subtitle = "Each subject has different within-species alpha diversity but pattern remains the same in both diet arms.") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction="horizontal",
        plot.title = element_text(size = 16, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 10, color = "darkslategrey", margin = margin(b = 25)),
        plot.caption = element_text(size = 8, margin = margin(t = 10), color = "grey70", hjust = 0))

AlphaPlot2

# Saving  the plot into pdf using ggsave

alpha <-ggarrange(AlphaPlot1, AlphaPlot2, nrow = 2, ncol = 1) %>% # put together two plots into a page
  ggsave( file="~/Desktop/alpha.pdf", width = 7, height = 10, dpi=300) # save plot into pdf



#### Dendogram ####

# Subset into Baseline Samples

Baseline <- Merge %>%
  filter(Day == -4) %>% # filter only samples on Day -4
  unite(Subject_Diet, SubjectFood:Diet, sep ="-Subject_") %>% # unite subject and diet column
  select (-Row.names, -X.SampleID, -Day) %>% # remove day and other data
  column_to_rownames(var = "Subject_Diet") %>% # make column into rownames
  t() %>% # transpose the matrix
  cor(method = "spearman") # calculate for spearman correlation


# Subset into Diet Samples
Diet <- Merge %>%
  filter(Day == 4) %>% # filter only samples on Day -4
  unite(Subject_Diet, SubjectFood:Diet, sep ="-Subject_") %>% # unite subject and diet column
  select (-Row.names, -X.SampleID, -Day) %>% # remove day and other data
  column_to_rownames(var = "Subject_Diet") %>% # make column into rownames
  t() %>% # transpose the matrix
  cor(method = "spearman") # calculate for spearman correlation

# Making heirarchal clustering for baseline and diet

base <- Baseline %>% # baseline data
  dist %>% # calculate a distance matrix, 
  hclust(method = "complete") %>% # hierarchical clustering 
  as.dendrogram # turn the object into a dendrogram.


diet<- Diet %>%# diet data
  dist %>% # calculate a distance matrix, 
  hclust(method = "complete") %>% # hierarchical clustering 
  as.dendrogram # turn the object into a dendrogram.

## Make a phylogentic tree for d baseline and diet samples

# Saving the dendogram  into one pdf file using base R
pdf(file="~/Desktop/phylo.pdf",paper="letter")

# Setting the graphical parameters
par(mfrow =c(2,1), oma=c(0.75,1,0.75,1)) # setting overall margins

par(mar=c(4,2.5,1,7)) # plot margins

# Baseline dendogram
base %>% set("leaves_pch", 19) %>% #set the leaf shape
  set("leaves_cex", 1.2) %>% # set leaf size
  set("leaves_col",  value = c("purple", "purple", "purple", "purple",
                               "salmon","purple", "salmon","purple", 
                               "purple", "purple", "purple", "purple",
                               "salmon", "purple")) %>% # manually assigning colors
  set("labels_cex", 0.75) %>% # font size of the label
  plot(main = "16s rRNA DNA-Seq Day -4 (Baseline)", # header of the plot
       xlab= "Spearman distance", horiz=TRUE) # label of x-axis

par(mar=c(4,2.5,2.5,7)) # plot margins

# Diet dendogram
diet %>% set("leaves_pch", 19) %>% # set the leaf shape
  set("leaves_cex", 1.2) %>% # set leaf size
  set("leaves_col",  value = c("salmon","salmon", "salmon", "salmon", "salmon", "salmon",
                               "salmon", "purple", "purple", "purple", "purple", "purple",
                               "purple", "purple")) %>% # manually assigning colors
  set("labels_cex", 0.75) %>% # font size of the label
  plot(main = "16s rRNA DNA-Seq Day 4 (Diet)", # header of the plot
       xlab= "Spearman distance", horiz=TRUE) # label of x-axis

dev.off() # turning device-off
