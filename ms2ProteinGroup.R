rm(list = ls());

# setting working directory - 
setwd("C:/Users/Fang-Ke/Desktop/Github");
# souring customerized functions
source('C:/Users/Fang-Ke/Desktop/Github/SILAC_Analysis_Functions.R');

# reading proteinGroup.txt file
AcetylK <- read.table.proteinGroup(30000);
filtered.AcetylK <- filter.proteinGroup.TMT(AcetylK);
quantified.AcetylK <- TMT.cleanup.10plex(filtered.AcetylK);

# calculating the median ratio of each channel vs mean
AcetylK.ratios <- get.TMT.ratios.10plex2(quantified.AcetylK);
# check.log2.ratios.distribution.10plex2(AcetylK.ratios, "Protein");

# calculating log2 and log10 intensities
AcetylK.intensities <- get.TMT.intensities.10plex(AcetylK.ratios);
# check.log10.intensity.distribution.10plex(AcetylK.intensities, "Protein");
# Ratios.vs.Intensity.10plex2(AcetylK.intensities, "Protein");
# Ratios.vs.normalized.Intensity.10plex2(AcetylK.intensities, "Protein");

Genename <- get.first.IDs(AcetylK.intensities[, "Gene.names"]);
# Genename <- get.first.IDs(AcetylK.intensities[, "Fasta.headers"]);
# Position.within.protein <- as.numeric(Position.within.protein);
# site.notation <- paste(Genename, AcetylK.intensities$Amino.acid, sep = ".");
# site.notation <- paste(site.notation, as.character(Position.within.protein), sep = "");
site.notation <- Genename;
AcetylK.intensities <- cbind(AcetylK.intensities, site.notation);

AcetylK.intensities <- AcetylK.intensities[!duplicated(site.notation),]
rownames(AcetylK.intensities) <- AcetylK.intensities$site.notation
intensities <- AcetylK.intensities[, 137:146]
colnames(intensities) <- c("AU565.1", "MDA-MB-231", "SKBR3/CAMA1", "T47D.1", "HCC1954.1",
                           "HCC1500.1", "AU565.2", "T47D.2", "HCC1954.2", "HCC1500.2")
library(ALL)
library(limma)
library(gplots)
intensities <- as.matrix(intensities)
# cluster.heatmap(intensities)

intensity.variance <- apply(intensities, 1, var)
intensity.mean <- apply(intensities, 1, mean)
# plot(intensity.mean, intensity.variance)
variance.cutoff <- quantile(intensity.variance, probs = seq(0, 1, 0.1))
# abline(h = variance.cutoff[10], col = "red")
varied.intensities <- intensities[unname(which(intensity.variance >= variance.cutoff[10])), ]


library(gtools)
# heatmap.2(intensities, col=colorpanel(75, low = "#00FF00", mid = "#FFFF00", high = "#EA0000"), key=TRUE, symkey=FALSE, 
#           density.info="none", trace="none", cexRow=0.5, 
#           hclustfun = function(x) hclust(x, method = "ward.D2"))
heatmap.2(varied.intensities, col=colorpanel(75, low = "#00FF00", mid = "#FFFF00", high = "#EA0000"), key=TRUE, symkey=FALSE, 
          density.info="none", trace="none", cexRow=0.5, 
          hclustfun = function(x) hclust(x, method = "ward.D2"))


library(Biobase)
# minimalSet <- ExpressionSet(assayData = intensities)

pData <- read.csv("C:/Users/Fang-Ke/Desktop/Github/phenotypicinfo.csv")
rownames(pData) <- as.character(pData[, 1])
pData <- pData[, 2:9]
summary(pData)
all(rownames(pData) == colnames(exprs))

metadata <- data.frame(labelDescription=c("Gene.Cluster", "ER +/-", "PR +/-",  "HER2 +/-",
                                          "Source", "Tumor.type", "Age", "Ethnicity"), 
                       row.names = c("Gene.Cluster","ER","PR","HER2","Source","Tumor.type","Age","Ethnicity"));
phenoData <- new("AnnotatedDataFrame", 
                 data = pData, varMetadata = metadata)

exampleSet <- ExpressionSet(assayData = intensities,
                            phenoData = phenoData);

f <- factor(as.character(exampleSet$HER2))
design <- model.matrix(~ f - 1)
colnames(design) <- c("Her2", "NEG")
fit <- eBayes(lmFit(exampleSet,design))
contrast.matrix <- makeContrasts("Her2-NEG", levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
statistics <- topTable(fit, number = nrow(intensities))
AcetylK.intensities <- merge(AcetylK.intensities, statistics, by = "row.names")
library(dplyr)
AcetylK.intensities <- mutate(AcetylK.intensities, minus.log10.adj.P.Val = -log10(adj.P.Val))
write.csv(AcetylK.intensities, "protein_entiretable_ms2.csv")

significant.list <- volcano.plot(AcetylK.intensities, cutoff = 2.3, filename = "protein",
                                 colselect = "minus.log10.adj.P.Val", colplotx = "logFC",
                                 colploty = "minus.log10.adj.P.Val", x_lab = "Log 2 Fold Change Her2+/Her2-",
                                 y_lab = "- log10 adjust P Value", 
                                 text_label = "Row.names")

table.s1 <- AcetylK.intensities[, c("Protein.names", "Gene.names", "Fasta.headers", 
                                    "Number.of.proteins", "Unique.peptides", "Unique.sequence.coverage....",
                                    "Mol..weight..kDa.", "Sequence.length", "normalized.Reporter.intensity.0", "normalized.Reporter.intensity.1",
                                    "normalized.Reporter.intensity.2", "normalized.Reporter.intensity.3", "normalized.Reporter.intensity.4", "normalized.Reporter.intensity.5",
                                    "normalized.Reporter.intensity.6", "normalized.Reporter.intensity.7", "normalized.Reporter.intensity.8", "normalized.Reporter.intensity.9",
                                    "Reporter.intensity.count.0", "Reporter.intensity.count.1", "Reporter.intensity.count.2", "Reporter.intensity.count.3", 
                                    "Reporter.intensity.count.4", "Reporter.intensity.count.5", "Reporter.intensity.count.6", "Reporter.intensity.count.7", 
                                    "Reporter.intensity.count.8", "Reporter.intensity.count.9", "logFC", "P.Value", "adj.P.Val")]
write.csv(table.s1, "Table.S1.Proteins identified and quantified by MS2 method.csv", row.names = FALSE)



nrow(significant.list)/nrow(AcetylK.intensities)

# f <- factor(as.character(exampleSet$HER2))
# design <- model.matrix(~f)
# fit <- eBayes(lmFit(exampleSet,design))
# topTable(fit, coef=2)
# selected  <- p.adjust(fit$p.value[, 2]) <0.05
# exampleSet.selected <- exampleSet[selected, ]
# significant.list <- exprs(exampleSet.selected);
# significant.genename <- get.first.IDs(rownames(significant.list), mark = ".");
# write.csv(significant.genename, "genelist_protein.csv");
# heatmap(exprs(exampleSet.selected))

#total table
selected.genes <- exprs(exampleSet)
selected.genes.p.values <- topTable(fit, coef=2, number = nrow(selected.genes))
selected.genes <- merge(selected.genes, selected.genes.p.values, by = "row.names")

#significnat table
significant.genes <- exprs(exampleSet.selected)
significant.genes.p.values <- selected.genes.p.values[rownames(significant.genes), ]
significant.genes <- merge(significant.genes, significant.genes.p.values, by = "row.names")

significant.genes.ms2 <- significant.genes
quantified.AcetylK75.ms2 <- AcetylK.intensities;


log10.AcetylK.intensities <- as.data.frame(AcetylK.intensities[, 127:136]);
colnames(log10.AcetylK.intensities) <- c("126", "127N", "127C", "128N", "128C",
                                         "129N", "129C", "130N", "130C", "131");
pairs.plot(log10.AcetylK.intensities, "Protein");
plot.correlation(log10.AcetylK.intensities, "126", "129C", "Protein");
plot.correlation(log10.AcetylK.intensities, "128N", "130N", "Protein");
plot.correlation(log10.AcetylK.intensities, "128C", "130C", "Protein");
plot.correlation(log10.AcetylK.intensities, "129N", "131", "Protein");



setwd("R:/neubertlabspace/Rawdata_Neubertlab/FangkeHuang/search results/CPTAC2016GrantApplication/Pro&Phospho_CPTAC_MS3/combined/txt");

AcetylK <- read.table.proteinGroup(30000);
filtered.AcetylK <- filter.proteinGroup.TMT(AcetylK);
quantified.AcetylK <- TMT.cleanup.10plex(filtered.AcetylK);
# AcetylK75 <- subset(filtered.AcetylK, Localization.prob >= 0.75);
# quantified.AcetylK75 <- TMT.cleanup.10plex(AcetylK75);

AcetylK.ratios <- get.TMT.ratios.10plex2(quantified.AcetylK);
check.log2.ratios.distribution.10plex2(AcetylK.ratios, "Protein");
AcetylK.intensities <- get.TMT.intensities.10plex(AcetylK.ratios);
check.log10.intensity.distribution.10plex(AcetylK.intensities, "Protein");
Ratios.vs.Intensity.10plex2(AcetylK.intensities, "Protein");
Ratios.vs.normalized.Intensity.10plex2(AcetylK.intensities, "Protein");

# Position.within.protein <- get.first.IDs(AcetylK.intensities[, "Positions.within.proteins"]);
# Genename <- get.first.IDs(AcetylK.intensities[, "Majority.protein.IDs"]);
Genename <- get.first.IDs(AcetylK.intensities[, "Gene.names"]);
# Position.within.protein <- as.numeric(Position.within.protein);
# site.notation <- paste(Genename, AcetylK.intensities$Amino.acid, sep = ".");
# site.notation <- paste(site.notation, as.character(Position.within.protein), sep = "");
site.notation <- Genename;
AcetylK.intensities <- cbind(AcetylK.intensities, site.notation);

AcetylK.intensities <- AcetylK.intensities[!duplicated(site.notation),]
rownames(AcetylK.intensities) <- AcetylK.intensities$site.notation
intensities <- AcetylK.intensities[, 137:146]
colnames(intensities) <- c("AU565.1", "MDA-MB-231", "SKBR3/CAMA1", "T47D.1", "HCC1954.1",
                           "HCC1500.1", "AU565.2", "T47D.2", "HCC1954.2", "HCC1500.2")
library(ALL)
library(limma)
library("gplots")
intensities <- as.matrix(intensities)
cluster.heatmap(intensities)

library(Biobase)
# minimalSet <- ExpressionSet(assayData = intensities)

pData <- read.csv("R:/neubertlabspace/Rawdata_Neubertlab/FangkeHuang/search results/CPTAC_final/combined/pYIP-txt/phenotypicinfo.csv")
rownames(pData) <- as.character(pData[, 1])
pData <- pData[, 2:9]
summary(pData)
all(rownames(pData) == colnames(exprs))

metadata <- data.frame(labelDescription=c("Gene.Cluster", "ER +/-", "PR +/-",  "HER2 +/-",
                                          "Source", "Tumor.type", "Age", "Ethnicity"), 
                       row.names = c("Gene.Cluster","ER","PR","HER2","Source","Tumor.type","Age","Ethnicity"));
phenoData <- new("AnnotatedDataFrame", 
                 data = pData, varMetadata = metadata)

exampleSet <- ExpressionSet(assayData = intensities,
                            phenoData = phenoData);

f <- factor(as.character(exampleSet$HER2))
design <- model.matrix(~ f - 1)
colnames(design) <- c("Her2", "NEG")
fit <- eBayes(lmFit(exampleSet,design))
contrast.matrix <- makeContrasts("Her2-NEG", levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
statistics <- topTable(fit, number = nrow(intensities))
AcetylK.intensities <- merge(AcetylK.intensities, statistics, by = "row.names")
library(dplyr)
AcetylK.intensities <- mutate(AcetylK.intensities, minus.log10.adj.P.Val = -log10(adj.P.Val))
write.csv(AcetylK.intensities, "protein_entiretable_ms2.csv")
significant.list <- volcano.plot(AcetylK.intensities, cutoff = 2.3, filename = "protein",
                                 colselect = "minus.log10.adj.P.Val", colplotx = "logFC",
                                 colploty = "minus.log10.adj.P.Val", x_lab = "Log 2 Fold Change Her2+/Her2-",
                                 y_lab = "- log10 adjust P Value", 
                                 text_label = "Row.names")

nrow(significant.list)/nrow(AcetylK.intensities)

# f <- factor(as.character(exampleSet$HER2))
# design <- model.matrix(~f)
# fit <- eBayes(lmFit(exampleSet,design))
# topTable(fit, coef=2)
# selected  <- p.adjust(fit$p.value[, 2]) <0.05
# exampleSet.selected <- exampleSet[selected, ]
# significant.list <- exprs(exampleSet.selected);
# significant.genename <- get.first.IDs(rownames(significant.list), mark = ".");
# write.csv(significant.genename, "genelist_protein.csv");
# heatmap(exprs(exampleSet.selected))

#total table
selected.genes <- exprs(exampleSet)
selected.genes.p.values <- topTable(fit, coef=2, number = nrow(selected.genes))
selected.genes <- merge(selected.genes, selected.genes.p.values, by = "row.names")

#significnat table
significant.genes <- exprs(exampleSet.selected)
significant.genes.p.values <- selected.genes.p.values[rownames(significant.genes), ]
significant.genes <- merge(significant.genes, significant.genes.p.values, by = "row.names")

significant.genes.ms2 <- significant.genes
quantified.AcetylK75.ms3 <- AcetylK.intensities;


log10.AcetylK.intensities <- as.data.frame(AcetylK.intensities[, 127:136]);
colnames(log10.AcetylK.intensities) <- c("126", "127N", "127C", "128N", "128C",
                                         "129N", "129C", "130N", "130C", "131");
pairs.plot(log10.AcetylK.intensities, "Protein");
plot.correlation(log10.AcetylK.intensities, "126", "129C", "Protein");
plot.correlation(log10.AcetylK.intensities, "128N", "130N", "Protein");
plot.correlation(log10.AcetylK.intensities, "128C", "130C", "Protein");
plot.correlation(log10.AcetylK.intensities, "129N", "131", "Protein");




#overlap between the two
AcetylK.overlap <- intersect(quantified.AcetylK75.ms2$site.notation,
                             quantified.AcetylK75.ms3$site.notation);
significant.genes.overlap <- intersect(significant.genes.ms2$Row.names, 
                                       significant.genes.ms3$Row.names)
#merge two table together
AcetylK.merge <- merge(quantified.AcetylK75.ms2, 
                       quantified.AcetylK75.ms3, 
                       by = "site.notation");
write.csv(AcetylK.merge, "protein.merge.csv");

#calculate ratios: channel 1 over channel 2; #plot the correlation of the ratios
ms2.ratio.1.2 <- AcetylK.merge$normalized.Reporter.intensity.1.x / AcetylK.merge$normalized.Reporter.intensity.2.x;
ms3.ratio.1.2 <- AcetylK.merge$normalized.Reporter.intensity.1.y / AcetylK.merge$normalized.Reporter.intensity.2.y;
AcetylK.merge <- cbind(AcetylK.merge, cbind(ms2.ratio.1.2, ms3.ratio.1.2));
plot.correlation(AcetylK.merge, "ms2.ratio.1.2", "ms3.ratio.1.2", "Protein.ratio.1.2", 
                 lab.note = "", diagonal = FALSE, log.2 = FALSE, smoothplot = FALSE);
# plot.ratio.ratio.intensity(AcetylK.merge, 
#                            "normalized.Reporter.intensity.1.x", 
#                            "normalized.Reporter.intensity.2.x",
#                            "normalized.Reporter.intensity.1.y",
#                            "normalized.Reporter.intensity.1.y",
#                            "ms2.ratio.1.2", "ms3.ratio.1.2",
#                            "ms2 ratio over ms3 ratio (log2)",
#                            dynamic.x.range = c(-4, 4));

ms2.ratio.06.59 <- (AcetylK.merge$normalized.Reporter.intensity.0.x + AcetylK.merge$normalized.Reporter.intensity.6.x) / (AcetylK.merge$normalized.Reporter.intensity.5.x + AcetylK.merge$normalized.Reporter.intensity.9.x);
ms3.ratio.06.59 <- (AcetylK.merge$normalized.Reporter.intensity.0.y + AcetylK.merge$normalized.Reporter.intensity.6.y) / (AcetylK.merge$normalized.Reporter.intensity.5.y + AcetylK.merge$normalized.Reporter.intensity.9.y);
AcetylK.merge <- cbind(AcetylK.merge, cbind(ms2.ratio.06.59, ms3.ratio.06.59));


plot.correlation(AcetylK.merge, "ms2.ratio.06.59", "ms3.ratio.06.59", "Protein.ratio.06.59", 
                 lab.note = "", diagonal = FALSE, log.2 = TRUE, smoothplot = TRUE);
# plot.ratio.ratio.intensity(AcetylK.merge, 
#                            "normalized.Reporter.intensity.0.x", 
#                            "normalized.Reporter.intensity.5.x",
#                            "normalized.Reporter.intensity.0.y",
#                            "normalized.Reporter.intensity.5.y",
#                            "ms2.ratio.06.59", "ms3.ratio.06.59",
#                            "ms2 ratio over ms3 ratio (log2)",
#                            dynamic.x.range = c(-4, 4));

ms2.ratio.37.48 <- (AcetylK.merge$normalized.Reporter.intensity.3.x + AcetylK.merge$normalized.Reporter.intensity.7.x) / (AcetylK.merge$normalized.Reporter.intensity.4.x + AcetylK.merge$normalized.Reporter.intensity.8.x);
ms3.ratio.37.48 <- (AcetylK.merge$normalized.Reporter.intensity.3.y + AcetylK.merge$normalized.Reporter.intensity.7.y) / (AcetylK.merge$normalized.Reporter.intensity.4.y + AcetylK.merge$normalized.Reporter.intensity.8.y);
AcetylK.merge <- cbind(AcetylK.merge, cbind(ms2.ratio.37.48, ms3.ratio.37.48));
plot.correlation(AcetylK.merge, "ms2.ratio.37.48", "ms3.ratio.37.48", "Protein.ratio.37.48", 
                 lab.note = "", diagonal = FALSE, log.2 = TRUE, smoothplot = TRUE);
# plot.ratio.ratio.intensity(AcetylK.merge, 
#                            "normalized.Reporter.intensity.3.x", 
#                            "normalized.Reporter.intensity.4.x",
#                            "normalized.Reporter.intensity.3.y",
#                            "normalized.Reporter.intensity.4.y",
#                            "ms2.ratio.37.48", "ms3.ratio.37.48",
#                            "ms2 ratio over ms3 ratio (log2)",
#                            dynamic.x.range = c(-4, 4));


AcetylK <- read.table.proteinGroup(30000);
filtered.AcetylK <- filter.proteinGroup.TMT(AcetylK);
quantified.AcetylK <- TMT.cleanup.10plex(filtered.AcetylK);
# AcetylK75 <- subset(filtered.AcetylK, Localization.prob >= 0.75);
# quantified.AcetylK75 <- TMT.cleanup.10plex(AcetylK75);

AcetylK.ratios <- get.TMT.ratios.10plex2(quantified.AcetylK);
check.log2.ratios.distribution.10plex2(AcetylK.ratios, "Protein");
AcetylK.intensities <- get.TMT.intensities.10plex(AcetylK.ratios);
check.log10.intensity.distribution.10plex(AcetylK.intensities, "Protein");
Ratios.vs.Intensity.10plex2(AcetylK.intensities, "Protein");
Ratios.vs.normalized.Intensity.10plex2(AcetylK.intensities, "Protein");

# Position.within.protein <- get.first.IDs(AcetylK.intensities[, "Positions.within.proteins"]);
# Genename <- get.first.IDs(AcetylK.intensities[, "Majority.protein.IDs"]);
Genename <- get.first.IDs(AcetylK.intensities[, "Gene.names"]);
# Position.within.protein <- as.numeric(Position.within.protein);
# site.notation <- paste(Genename, AcetylK.intensities$Amino.acid, sep = ".");
# site.notation <- paste(site.notation, as.character(Position.within.protein), sep = "");
site.notation <- Genename;
AcetylK.intensities <- cbind(AcetylK.intensities, site.notation);

AcetylK.intensities <- AcetylK.intensities[!duplicated(site.notation),]
rownames(AcetylK.intensities) <- AcetylK.intensities$site.notation
intensities <- AcetylK.intensities[, 137:146]
colnames(intensities) <- c("AU565.1", "MDA-MB-231", "SKBR3/CAMA1", "T47D.1", "HCC1954.1",
                           "HCC1500.1", "AU565.2", "T47D.2", "HCC1954.2", "HCC1500.2")
library(ALL)
library(limma)
library("gplots")
intensities <- as.matrix(intensities)
cluster.heatmap(intensities)

library(Biobase)
# minimalSet <- ExpressionSet(assayData = intensities)

pData <- read.csv("R:/neubertlabspace/Rawdata_Neubertlab/FangkeHuang/search results/CPTAC_final/combined/pYIP-txt/phenotypicinfo.csv")
rownames(pData) <- as.character(pData[, 1])
pData <- pData[, 2:9]
summary(pData)
all(rownames(pData) == colnames(exprs))

metadata <- data.frame(labelDescription=c("Gene.Cluster", "ER +/-", "PR +/-",  "HER2 +/-",
                                          "Source", "Tumor.type", "Age", "Ethnicity"), 
                       row.names = c("Gene.Cluster","ER","PR","HER2","Source","Tumor.type","Age","Ethnicity"));
phenoData <- new("AnnotatedDataFrame", 
                 data = pData, varMetadata = metadata)

exampleSet <- ExpressionSet(assayData = intensities,
                            phenoData = phenoData);

f <- factor(as.character(exampleSet$HER2))
design <- model.matrix(~ f - 1)
colnames(design) <- c("Her2", "NEG")
fit <- eBayes(lmFit(exampleSet,design))
contrast.matrix <- makeContrasts("Her2-NEG", levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
statistics <- topTable(fit, number = nrow(intensities))
AcetylK.intensities <- merge(AcetylK.intensities, statistics, by = "row.names")
library(dplyr)
AcetylK.intensities <- mutate(AcetylK.intensities, minus.log10.adj.P.Val = -log10(adj.P.Val))
write.csv(AcetylK.intensities, "protein_entiretable_ms3.csv")
