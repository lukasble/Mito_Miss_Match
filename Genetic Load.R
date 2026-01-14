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
dataset$GeneticLoad <- mid_breed$LoadCount+high_breed$LoadCount #Summing count of HIGH and MID impact SNPs per dog sample ID
dataset["SynonymousCounts"] <- syn_breed$SynonymousCounts
dataset["NonsynonymousCounts"] <- nonsyn_breed$NonSynonymousCounts
dataset["ModifierCounts"] <- mod_breed$LoadCount
dataset["LowCounts"] <- low_breed$LoadCount
dataset["ModerateCounts"] <- mid_breed$LoadCount
dataset["HighCounts"] <- high_breed$LoadCount
dataset["TotalSNPs"] <- count_breed$snp_count

#Calculating and adding ratios to dataset
dataset["NS_Ratio"] <- as.numeric(dataset$NonsynonymousCounts)/as.numeric(dataset$SynonymousCounts)
dataset["GeneticLoadRatio"] <- as.numeric(dataset$GeneticLoad)/as.numeric(dataset$TotalSNPs)
dataset["SNP_Ratio"] <- as.numeric(dataset$TotalSNPs)/801164
dataset["ModifierRatio"] <- as.numeric(dataset$ModifierCounts)/as.numeric(dataset$TotalSNPs)
dataset["LowRatio"] <- as.numeric(dataset$LowCounts)/as.numeric(dataset$TotalSNPs)
dataset["ModerateRatio"] <- as.numeric(dataset$ModerateCounts)/as.numeric(dataset$TotalSNPs)
dataset["HighRatio"] <- as.numeric(dataset$HighCounts)/as.numeric(dataset$TotalSNPs)
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
dataset<- dataset[,c(1,3,27,26,18,11,21,20,12,13,19, 14:17,22:25)]

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
     main = "", ylab = "Number of individuals", xlab = "Weight (Kg)", 
     col="lightblue", ylim = c(0,250))
grid()
hist(dataset$SynonymousCounts, 
     main = "", ylab = "Number of individuals", xlab = "Synonymous SNP Count", 
     col="lightblue", ylim = c(0,1000))
grid()
hist(dataset$NonsynonymousCounts, 
     main = "", ylab = "Number of individuals",  xlab = "Nonsynonymous SNP Count", 
     col="lightblue", ylim = c(0,1300))
grid()
hist(dataset$NS_Ratio, 
     main = "", ylab = "Number of individuals",  xlab = "Nonsyn./Syn. Ratio", 
     col="lightblue", ylim = c(0,850))
grid()

#Pearson correlation test dataset
cor.test(dataset$Weight, dataset$TotalSNPs)
cor.test(dataset$Weight, dataset$GeneticLoad)
cor.test(dataset$Weight, dataset$GeneticLoadRatio)
cor.test(dataset$Weight, dataset$NS_Ratio)
cor.test(dataset$Weight, dataset$ModifierRatio)
cor.test(dataset$Weight, dataset$LowRatio)
cor.test(dataset$Weight, dataset$ModerateRatio)
cor.test(dataset$Weight, dataset$HighRatio)

ncol(dataset)
#Correlation analysis visualization
par(mfrow=c(1,1))
matrix <- cor(dataset[,c(4:6,8,11,16:19)], method = c("pearson"))
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
        main = "Total SNP Number across size classes",
        xlab = "",
        ylab = "Count",
        outline = TRUE)
grid()

boxplot(GeneticLoadRatio ~ Breed_size,
        data = dataset,
        col = "lightblue",
        main = "Genetic Load Proportion across size classes",
        xlab = "",
        ylab = "Proportion",
        outline = TRUE)
grid()


#Scatter plots for visualization of data dispersion and fitting of linear regression model
par(mfrow=c(1,2))
par(cex.main = 2)
par(cex.lab = 1.6)
par(cex.axis = 1.6)

#Total SNP count
model <- lm(TotalSNPs ~ Weight, data = dataset)
xlim <- range(dataset$Weight)
ylim <- range(dataset$TotalSNPs)

x_pos <- xlim[1] + 0.7 * diff(xlim)
y_pos <- ylim[1] + 0.5 * diff(ylim)

eqn <- bquote(
  italic(y) == .(round(coef(model)[1], 2)) +
    .(round(coef(model)[2], 2)) %.% italic(x)
)

stats <- bquote(
  R^2 == .(round(summary(model)$r.squared, 3)) * "," ~
    italic(p) == .(signif(summary(model)$coefficients[2, 4], 3))
)
plot(dataset$Weight, dataset$TotalSNPs,
     xlab = "Weight (kg)",
     ylab = "Count",
     main = "Total SNP Number vs Body Weight",
     pch = 19,
     col = "steelblue")

abline(model, col = "red", lwd = 2)
text(x = x_pos,
     y = y_pos,
     labels = as.expression(eqn),
     pos = 4,
     cex =1.3,
     col = "black")

text(x = x_pos,
     y = y_pos - diff(range(dataset$TotalSNPs)) * 0.035,
     labels = as.expression(stats),
     pos = 4,
     cex =1.3,
     col = "black")
grid()

#Genetic Load Ratio
model <- lm(GeneticLoadRatio ~ Weight, data = dataset)
xlim <- range(dataset$Weight)
ylim <- range(dataset$GeneticLoadRatio)

x_pos <- xlim[1] + 0.7 * diff(xlim)
y_pos <- ylim[1] + 0.45 * diff(ylim)

eqn <- bquote(
  italic(y) == .(round(coef(model)[1], 2)) +
    .(round(coef(model)[2], 6)) %.% italic(x)
)

stats <- bquote(
  R^2 == .(round(summary(model)$r.squared, 3)) * "," ~
    italic(p) == .(signif(summary(model)$coefficients[2, 4], 3))
)

plot(dataset$Weight, dataset$GeneticLoadRatio,
     xlab = "Weight (kg)",
     ylab = "Proportion",
     main = "Genetic Load Proportion vs Body Weight",
     pch = 19,
     col = "steelblue")
abline(model, col = "red", lwd = 2)
text(x = x_pos,
     y = y_pos,
     labels = as.expression(eqn),
     pos = 4,
     cex =1.3,
     col = "black")

text(x = x_pos,
     y = y_pos - diff(range(dataset$GeneticLoadRatio)) * 0.035,
     labels = as.expression(stats),
     pos = 4,
     cex =1.3,
     col = "black")
grid()

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