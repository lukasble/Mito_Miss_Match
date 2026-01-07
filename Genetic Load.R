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

setwd("")  #Set your working directory
dataset <- read.csv("damaging_alleles_with_breed.csv")
syn_breed <- read.csv("syn_per_sample.tsv", sep="\t")
nonsyn_breed <- read.csv("nonsyn_per_sample.tsv", sep="\t")
mod_breed <- read.csv("modifier_per_sample.tsv", sep="\t")
low_breed <- read.csv("low_per_sample.tsv", sep="\t")
mid_breed <- read.csv("moderate_per_sample.tsv", sep="\t")
high_breed <- read.csv("high_per_sample.tsv", sep="\t")
#total_breed <- read.csv("snp_counts_per_sample.tsv", sep="\t")
count_breed <- read.csv("count_samp.tsv", sep="\t")

#Filling Dataset
dataset$DamagingAltAlleleCount <- mid_breed$LoadCount+high_breed$LoadCount #Summing count of HIGH and MID impact SNPs per dog sample ID
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

#Average of samples per dog breed and bining dog breed by weight
#Classifying dog breed by weight
dataset["Breed_size"] <- ifelse(dataset$weight_value < 5, "Toy",
                                ifelse(dataset$weight_value < 11, "Small",
                                       ifelse(dataset$weight_value < 25, "Medium",
                                              ifelse(dataset$weight_value < 39, "Large",
                                                     "Giant"))))

#Calculate average per breed for each metric (according to breed code) and merge them in a new breed-only set
breed_means_avg <- dataset %>%
  filter(!is.na(Breed_size)) %>%    #Filtering rows with breed size value NA
  filter(!grepl("CLUP", breed_code)) %>%   #Filtering rows with wolf data by sample ID code CLUP
  filter(!grepl("VILL", breed_code)) %>%   #Filtering rows with Asian Village Dog data by sample ID code VILL
  group_by(breed_code, Breed_size) %>%
  summarise(
    Weight = mean(weight_value, na.rm = TRUE),
    Dmg = mean(DamagingAltAlleleCount, na.rm = TRUE),
    Syn = mean(SynonymousCounts, na.rm = TRUE),
    Non = mean(NonsynonymousCounts, na.rm = TRUE),
    dNdS = mean(dNdS, na.rm = TRUE),
    Tot = mean(TotalSNPs, na.rm = TRUE),
    SNPProp = mean(SNPProp, na.rm = TRUE),
    Modi = mean(ModifierCounts, na.rm = TRUE),
    Low = mean(LowCounts, na.rm = TRUE),
    Mode = mean(ModerateCounts, na.rm = TRUE),
    High = mean(HighCounts, na.rm = TRUE),
    ModiRatio = mean(ModifierSNPratio, na.rm = TRUE),
    LowRatio = mean(LowSNPratio, na.rm = TRUE),
    ModeRatio = mean(ModerateSNPratio, na.rm = TRUE),
    HighRatio = mean(HighSNPratio, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

nrow(breed_means_avg) #349 remaining dog breeds
summary(breed_means_avg) #Summary statistics per dog breed
par(mfrow=c(2,2))
par(cex.lab = 1.2)
par(cex.axis = 1.6)
hist(breed_means_avg$Weight, main = "", xlab = "Weights")
hist(breed_means_avg$Syn, main = "", xlab = "Synonymous Counts")
hist(breed_means_avg$Non, main = "", xlab = "Nonsynonymous Counts")
hist(breed_means_avg$dNdS, main = "", xlab = "Total dN/dS Ratio")

#Pearson correlation test 
cor.test(breed_means_avg$Weight, breed_means_avg$Dmg)
cor.test(breed_means_avg$Weight, breed_means_avg$Syn)
cor.test(breed_means_avg$Weight, breed_means_avg$Non)
cor.test(breed_means_avg$Weight, breed_means_avg$dNdS)
cor.test(breed_means_avg$Weight, breed_means_avg$Tot)
cor.test(breed_means_avg$Weight, breed_means_avg$High)
cor.test(breed_means_avg$Weight, breed_means_avg$LowRatio)
cor.test(breed_means_avg$Weight, breed_means_avg$ModiRatio)
cor.test(breed_means_avg$Weight, breed_means_avg$ModeRatio)
cor.test(breed_means_avg$Weight, breed_means_avg$HighRatio)

#Correlation analysis visualization
par(mfrow=c(1,1))
matrix <- cor(breed_means_avg[,c(3:4,7:8,14:17)], method = c("pearson"))
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
boxplot(Tot ~ Breed_size,
        data = breed_means_avg,
        col = "lightblue",
        main = "Total Mutations",
        xlab = "",
        ylab = "",
        outline = TRUE)

boxplot(SNPProp ~ Breed_size,
        data = breed_means_avg,
        col = "lightblue",
        main = "High Effect Mutations",
        xlab = "",
        ylab = "",
        outline = TRUE)

boxplot(Mode ~ Breed_size,
        data = breed_means_avg,
        col = "lightblue",
        main = "Moderate Effect Mutations",
        xlab = "",
        ylab = "",
        outline = TRUE)

boxplot(High ~ Breed_size,
        data = breed_means_avg,
        col = "lightblue",
        main = "High Effect Mutations",
        xlab = "",
        ylab = "",
        outline = TRUE)

#Scatterplots for trends and fitting linear regression line
par(mfrow=c(1,2))
plot(breed_means_avg$Weight, breed_means_avg$Tot,
     xlab = "Body Weight (kg)",
     ylab = "Total SNPs",
     main = "Total SNPs vs Body weight",
     pch = 19,
     col = "steelblue")

abline(lm(Tot ~ Weight, data = breed_means_avg),
       col = "red", lwd = 2)

plot(breed_means_avg$Weight, breed_means_avg$SNPProp,
     xlab = "Body Weight (kg)",
     ylab = "ynonymous Counts",
     main = "Syn. SNPs vs Body Weight",
     pch = 19,
     col = "steelblue")
abline(lm(SNPProp ~ Weight, data = breed_means_avg),
       col = "red", lwd = 2)

plot(breed_means_avg$Weight, breed_means_avg$HighRatio,
     xlab = "Average Breed Weight",
     ylab = "Average High Effect SNPs",
     main = "High Impact Ratio vs Body Weight",
     pch = 19,
     col = "steelblue")
abline(lm(HighRatio ~ Weight, data = breed_means_avg),
       col = "red", lwd = 2)

plot(breed_means_avg$Weight, breed_means_avg$dNdS,
     xlab = "Avg Breed Weight",
     ylab = "Avg dNdS ratio",
     main = "dN/dS Ratio vs Body Weight",
     pch = 19,
     col = "steelblue")
abline(lm(dNdS ~ Weight, data = breed_means_avg),
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