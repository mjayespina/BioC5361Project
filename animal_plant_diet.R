# BioC5361 Project
# Title: Animal vs. Plant Diet
# Group : Group 2
# filename: animal_plant_diet.R
# Author: Mary Jane Espina, Claire Peichel , Sarah Lichtenberger

#### Data importing #####

# Libraries
library(tidyverse) # install.packages("tidyverse")
library(reshape2) # install.packages("reshape2")
library(ape) # install.packages("ape")
library(vegan) # install.packages("vegan")
library(ggplot2) # install.packages("ggplot2")
library(gridExtra) # install.packages("gridExtra")

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

Sum <- colSums(Taxa) #get the column sums of the Taxa data frame

view(Sum) #view the Sum

summary(Sum) #get summary statistics of the colSums

# Quality control of the total reads
# Only samples with more than 20,000 reads were included

TaxaQC <- Taxa[colSums(Taxa)>= 20045]

# Convert to relative abundance
TaxaNorm <- sweep(TaxaQC, 2, colSums(TaxaQC), FUN = "/")

# Transpose TaxaNorm data to easily merge two files
TaxaT <- t(TaxaNorm) 

row.names(Meta) <- Meta$X.SampleID

Merge <- merge(x=Meta,y=TaxaT,by=0) #merge two data together


#### Plotting alpha diversity ####

# Alpha Diversity

# Compute for shannon alpha diversity
TaxaAlpha <-  diversity(TaxaNorm,index = "shannon",MARGIN = 2 )
View(TaxaAlpha)


# Merge Alpha and Metadata into one file

AlphaMerge <- merge(x=Meta,y=TaxaAlpha,by=0) %>% #merge the alpha diversity and metadata
  filter(!Diet == "NA")  %>% #filter out diet in the data frame
  select(-Row.names) %>% #remove Row.names since it has the same info as X.SampleID
  rename(Subject=SubjectFood) %>% #rename SubjectFood to Subject
  group_by(Day) %>% # Grooup data by day
  melt(id.vars=c("X.SampleID", 
                 "Subject","Day", "Diet"),
       variable.name=c("y"), value.name=c("Alpha")) #transform data into long format

# Plotting alpha diversity, grouped by Day

AlphaPlot1 <- ggplot(AlphaMerge, aes(x=Day, y=Alpha, group=Day)) +
  geom_boxplot() + 
  
  facet_wrap (~ Diet) + #facet based on diet
  labs(x="Day", y="Alpha Diversity (Shannon)",
       title = "Within-species alpha diversity by day for animal vs plant diet",
       subtitle = "There is decline in within-species alpha diversity on the 4th day of plant-based diet arm.") +
  theme_bw() +
  theme(text = element_text(family = "Georgia"), 
        plot.title = element_text(size = 18, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 12, color = "darkslategrey", margin = margin(b = 25)),
        plot.caption = element_text(size = 8, margin = margin(t = 10), color = "grey70", hjust = 0))

AlphaPlot1

# Plotting alpha diversity, grouped by Subject

AlphaPlot2 <- ggplot(AlphaMerge, aes(x=reorder(Subject, -Alpha), y=Alpha, fill=Subject)) +
  geom_boxplot() +
  facet_wrap (~ Diet) +
  labs(x="Subject", y="Alpha Diversity (Shannon)",
       title = "Within-species alpha diversity by Subject for animal vs plant diet",
       subtitle = "Each subject has different within-species alpha diversity but pattern remains the same in both diet arms.") +
  theme_bw() +
  theme(text = element_text(family = "Arial"), 
        legend.position = "bottom",
        legend.direction="horizontal",
        plot.title = element_text(size = 18, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 12, color = "darkslategrey", margin = margin(b = 25)),
        plot.caption = element_text(size = 8, margin = margin(t = 10), color = "grey70", hjust = 0))

AlphaPlot2

#### Dendogram ####

# Libraries
library(ggdendro)
library(ggraph)

# Subset into Baseline Samples

Baseline <- Merge %>%
  filter(Day == -4) %>% #filter only samples on Day -4
  unite(Subject_Diet, SubjectFood:Diet, sep ="-Subject_") %>% #Unite subject and diet column
  select (-Row.names, -X.SampleID, -Day) %>% #remove day and other data
  column_to_rownames(var = "Subject_Diet") %>% #make column into rownames
  t() %>% #transpose the matrix
  cor(method = "spearman") #calculate for spearman correlation

base <-hclust(dist(Baseline))

# Make a phylogentic tree for baseline samples
basephylo <-ggdendrogram(base, rotate = TRUE, size = 2.5, leaf_label=TRUE) +
  coord_flip() +
  labs(title = "16s Sequencing (Spearman distance)",
       subtitle = "Baseline Samples- Day -4")
basephylo

# Subset into Diet Samples
Diet <- Merge %>%
  filter(Day == 4) %>% #filter only samples on Day -4
  unite(Subject_Diet, SubjectFood:Diet, sep ="-Subject_") %>%
  select (-Row.names, -X.SampleID, -Day) %>%
  column_to_rownames(var = "Subject_Diet") %>%
  t() %>%
  cor(method = "spearman") 


# Make a phylogentic tree for diet samples
diet <-hclust(dist(Diet))
  
dietphylo <-ggdendrogram(diet, rotate = TRUE, size = 2.5, leaf_label=TRUE) +
  coord_flip() +
  labs(title = "16s Sequencing (Spearman distance)",
       subtitle = "Diet Samples - Day 4")

dietphylo

ggsave (basephylo, file="~/Desktop/base.pdf", width = 6, height = 5, dpi=300)
ggsave (dietphylo, file="~/Desktop/diet.pdf", width = 6, height = 5, dpi=300)








