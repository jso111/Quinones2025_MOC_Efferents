#title: Prestin immunolabeling stats
#author: John Oghalai
#date: 2/4/2025
#install.packages("tidyverse")
library(gsignal)
library(ggplot2)
library(gghighlight)
library(ggsignif)
library(gridExtra)
library(ggpubr)
library(plotrix)
library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)

#####
# Load data and clean it
#####
v1 <- tibble(read.csv('v1.csv'))
genotypes = c('WT','Alpha9KO','VGLUT3KO','double')
v1$genotype<-factor(v1$genotype,levels=genotypes)

# Run ANOVA
fit <- aov(prestin ~ genotype, data=v1)
summary(fit)


