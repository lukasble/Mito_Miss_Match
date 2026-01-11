install.packages("tidyverse")
install.packages("ggplot2")
install.packages("BioWorldR")
install.packages("corrplot")
install.packages("Hmisc")
install.packages("reshape2")
install.packages("factoextra")

library(factoextra)
library(tidyverse)
library(ggplot2)
library(BioWorldR)
library(ggplot2)
library(dplyr)
library(corrplot)
library(Hmisc)
library(reshape2)

#Set working directory
setwd("")
dataset <- read.csv("damaging_alleles_with_breed.csv")
syn_breed <- read.csv("syn_per_sample.tsv", sep="\t")
nonsyn_breed <- read.csv("nonsyn_per_sample.tsv", sep="\t")
mod_breed <- read.csv("modifier_per_sample.tsv", sep="\t")
low_breed <- read.csv("low_per_sample.tsv", sep="\t")
mid_breed <- read.csv("moderate_per_sample.tsv", sep="\t")
high_breed <- read.csv("high_per_sample.tsv", sep="\t")
count_breed <- read.csv("count_samp.tsv", sep="\t")

#Filling Dataset
dataset$GeneLoad <- mid_breed$LoadCount+high_breed$LoadCount #Summing count of HIGH and MID impact SNPs per dog sample ID
dataset["SynonymousCounts"] <- syn_breed$SynonymousCounts
dataset["NonsynonymousCounts"] <- nonsyn_breed$NonSynonymousCounts
dataset["ModifierCounts"] <- mod_breed$LoadCount
dataset["LowCounts"] <- low_breed$LoadCount
dataset["ModerateCounts"] <- mid_breed$LoadCount
dataset["HighCounts"] <- high_breed$LoadCount
dataset["TotalSNPs"] <- count_breed$snp_count

#Calculating and adding ratios to dataset
dataset["dNdS"] <- as.numeric(dataset$NonsynonymousCounts)/as.numeric(dataset$SynonymousCounts)
dataset["SNPProp"] <- as.numeric(dataset$TotalSNPs)/801164
dataset["ModifierSNPratio"] <- as.numeric(dataset$ModifierCounts)/as.numeric(dataset$TotalSNPs)
dataset["LowSNPratio"] <- as.numeric(dataset$LowCounts)/as.numeric(dataset$TotalSNPs)
dataset["ModerateSNPratio"] <- as.numeric(dataset$ModerateCounts)/as.numeric(dataset$TotalSNPs)
dataset["HighSNPratio"] <- as.numeric(dataset$HighCounts)/as.numeric(dataset$TotalSNPs)
dataset["Weight"] <- as.numeric(dataset$weight_value)

#Classifying dog breed by weight class
dataset["Breed_size"] <- ifelse(dataset$Weight < 5, "Toy",
                                ifelse(dataset$Weight < 10, "Small",
                                       ifelse(dataset$Weight < 20, "Medium",
                                              ifelse(dataset$Weight < 40, "Large",
                                                     "Giant"))))

dataset <- dataset %>%
  filter(!is.na(Breed_size)) %>%    #Filtering rows with breed size value NA
  filter(!grepl("CLUP", breed_code)) %>%   #Filtering rows with wolf data by sample ID code CLUP
  filter(!grepl("VILL", breed_code))       #Filtering rows with Asian Village Dog data by sample ID code VILL

#Ordering columns and filtering dataset
dataset<- dataset[,c(1,3,26,25,18,11,19,12:17,20:24)]

#Total dog sample number and per weight category
nrow(dataset)  #1602 total dog samples
nrow(dataset[dataset$Breed_size=="Toy",])
nrow(dataset[dataset$Breed_size=="Small",])
nrow(dataset[dataset$Breed_size=="Medium",])
nrow(dataset[dataset$Breed_size=="Large",])
nrow(dataset[dataset$Breed_size=="Giant",])

#Metric distributions
par(mfrow=c(2,2))
par(cex.lab = 1.8)
par(cex.axis = 1.6)
hist(dataset$Weight,
     main = "", ylab = "", xlab = "Weight (Kg)", 
     col="lightblue", ylim = c(0,250))
grid()
hist(dataset$SynonymousCounts, 
     main = "", ylab = "",  xlab = "Synonymous SNPs", 
     col="lightblue", ylim = c(0,1000))
grid()
hist(dataset$NonsynonymousCounts, 
     main = "", ylab = "",  xlab = "Nonsynonymous SNPs", 
     col="lightblue", ylim = c(0,1300))
grid()
hist(dataset$dNdS, 
     main = "", ylab = "",  xlab = "dN/dS SNP Ratio", 
     col="lightblue", ylim = c(0,850))
grid()

#Pearson correlation test dataset
cor.test(dataset$Weight, dataset$GeneLoad)
cor.test(dataset$Weight, dataset$SynonymousCounts)
cor.test(dataset$Weight, dataset$NonsynonymousCounts)
cor.test(dataset$Weight, dataset$dNdS)
cor.test(dataset$Weight, dataset$TotalSNPs)
cor.test(dataset$Weight, dataset$NonsynonymousCounts)
cor.test(dataset$Weight, dataset$dNdS)

ncol(dataset)
#Correlation analysis visualization
par(mfrow=c(1,1))
matrix <- cor(dataset[,c(4:7,14:18)], method = c("pearson"))
corrplot(matrix,
         method = "circle",
         type = "lower",
         diag = FALSE,
         tl.cex = 1.2,
         tl.offset = 1,   # vertical offset for top labels
         tl.col = "black",   # label color
         tl.srt = 30)        # rotate labels

#Boxplots per dog category of SNPs
par(mfrow=c(1,2))
par(cex.main = 2)
par(cex.lab = 1.2)
par(cex.axis = 1.6)
boxplot(TotalSNPs ~ Breed_size,
        data = dataset,
        col = "lightblue",
        main = "Total SNPs",
        xlab = "",
        ylab = "",
        outline = TRUE)
grid()

boxplot(GeneLoad ~ Breed_size,
        data = dataset,
        col = "lightblue",
        main = "High Effect Mutations",
        xlab = "",
        ylab = "",
        outline = TRUE)
grid()

#Scatterplots for trends and fitting linear regression line
par(mfrow=c(1,2))
par(cex.main = 2)
par(cex.lab = 1.2)
par(cex.axis = 1.6)
plot(dataset$Weight, dataset$TotalSNPs,
     xlab = "Weight (kg)",
     ylab = "",
     main = "Total SNPs vs Body Weight",
     pch = 19,
     col = "steelblue")

abline(lm(TotalSNPs ~ Weight, data = dataset),
       col = "red", lwd = 2)

plot(dataset$Weight, dataset$SNPProp,
     xlab = "Weight (kg)",
     ylab = "",
     main = "Gene Load proportion vs Body Weight",
     pch = 19,
     col = "steelblue")
abline(lm(SNPProp ~ Weight, data = dataset),
       col = "red", lwd = 2)


#PCA plot - dimensionality reduction to find predictors that drive most variation in data - excluding the Weight
pca <- prcomp(breed_means_avg[4:16], scale = TRUE)

#Calculate total variance explaned by each component in the pca and create screeplot
var_pc <-  round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

#display principal components and plot
pca$rotation

par(mfrow=c(1,1))
biplot(pca, scale = 0, cex = 0)
fviz_pca_biplot(
  pca,
  repel = TRUE,
  label = "var"
) +
  labs(
    x = paste0("PC1 (", var_pc[1], "%)"),
    y = paste0("PC2 (", var_pc[2], "%)")
  )

#Scree plot
scree_df <- data.frame(
  PC = factor(paste0("PC", seq_along(var_pc)),
              levels = paste0("PC", seq_along(var_pc))),
  Variance = var_pc
)
#Scree plot of variance explained by PCs
ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Scree Plot: Variance Explained by Principal Components",
    x = "Principal Component",
    y = "Proportion of Variance Explained"
  ) +
  ylim(0, 1) +
  theme_minimal()