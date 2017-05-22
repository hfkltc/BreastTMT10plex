rm(list = ls());

setwd("C:/Users/Fang-Ke/Desktop/Github/ms2pSTY");
source('C:/Users/Fang-Ke/Desktop/Github/SILAC_Analysis_Functions.R');

AcetylK <- read.table.PhosphoSites(30000);
filtered.AcetylK <- filter.proteinGroup.TMT(AcetylK);
quantified.AcetylK <- TMT.cleanup.10plex(filtered.AcetylK);
AcetylK75 <- subset(filtered.AcetylK, Localization.prob >= 0.75);
quantified.AcetylK75 <- TMT.cleanup.10plex(AcetylK75);

# table(filtered.AcetylK$Amino.acid)
# table(filtered.AcetylK$Amino.acid)/nrow(filtered.AcetylK)

# table(quantified.AcetylK$Amino.acid)
# table(quantified.AcetylK$Amino.acid)/nrow(quantified.AcetylK)

# table(quantified.AcetylK75$Amino.acid)
# table(quantified.AcetylK75$Amino.acid)/nrow(quantified.AcetylK75)

# table(AcetylK75$Amino.acid)
# table(AcetylK75$Amino.acid)/nrow(AcetylK75)

AcetylK.ratios <- get.TMT.ratios.10plex2(quantified.AcetylK);
# check.log2.ratios.distribution.10plex2(AcetylK.ratios, "PhosphoSTY");
AcetylK.intensities <- get.TMT.intensities.10plex(AcetylK.ratios);
# check.log10.intensity.distribution.10plex(AcetylK.intensities, "PhosphoSTY");
# Ratios.vs.Intensity.10plex2(AcetylK.intensities, "PhosphoSTY");
# Ratios.vs.normalized.Intensity.10plex2(AcetylK.intensities, "PhosphoSTY");

Position.within.protein <- get.first.IDs(AcetylK.intensities[, "Positions.within.proteins"]);
Genename <- get.first.IDs(AcetylK.intensities[, "Gene.names"]);
# Genename <- AcetylK.intensities$Protein
Position.within.protein <- as.numeric(Position.within.protein);
site.notation <- paste(Genename, AcetylK.intensities$Amino.acid, sep = ".");
site.notation <- paste(site.notation, as.character(Position.within.protein), sep = "");
AcetylK.intensities <- cbind(AcetylK.intensities, site.notation);

# remove duplications
# AcetylK.intensities <- AcetylK.intensities[!duplicated(site.notation),]
# temp <- AcetylK.intensities[(duplicated(site.notation) | duplicated(site.notation, fromLast = TRUE)), ]
# rownames(AcetylK.intensities) <- AcetylK.intensities$site.notation
rownames(AcetylK.intensities) <- 1:nrow(AcetylK.intensities)
intensities <- AcetylK.intensities[, 152:161]
colnames(intensities) <- c("AU565.1", "MDA-MB-231", "SKBR3/CAMA1", "T47D.1", "HCC1954.1",
                           "HCC1500.1", "AU565.2", "T47D.2", "HCC1954.2", "HCC1500.2")
library(ALL)
intensities <- as.matrix(intensities)
# write.csv(intensities, "pYIP_intensities.csv")
# cluster.heatmap(intensities)

library(gplots)
library(gtools)
heatmap.2(intensities, col=colorpanel(75, low = "#00FF00", mid = "#FFFF00", high = "#EA0000"), key=TRUE, symkey=FALSE, 
          density.info="none", trace="none", cexRow=0.5, 
          hclustfun = function(x) hclust(x, method = "ward.D2"))

library(Biobase)
library(ALL)
library(limma)
library("gplots")
# minimalSet <- ExpressionSet(assayData = intensities)

pData <- read.csv("phenotypicinfo.csv")
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
write.csv(AcetylK.intensities, "pYIP_entiretable.csv")
significant.list <- volcano.plot(AcetylK.intensities, cutoff = 2.30103, filename = "pYIP",
                                 colselect = "minus.log10.adj.P.Val", colplotx = "logFC",
                                 colploty = "minus.log10.adj.P.Val", x_lab = "Log 2 Ratio Her2+/Her2-",
                                 y_lab = "- log10 adjust P Value", text_label = FALSE,
                                 is.title = FALSE)

# table.s2 <- AcetylK.intensities[, c("Proteins", "Leading.proteins", "Protein.names", "Gene.names", "Fasta.headers", "Localization.prob", "PEP",
#                                     "Number.of.Phospho..STY.", "Amino.acid", "Sequence.window", "Phospho..STY..Probabilities", "Position.in.peptide", "Charge", 
#                                     "normalized.Reporter.intensity.0", "normalized.Reporter.intensity.1",
#                                     "normalized.Reporter.intensity.2", "normalized.Reporter.intensity.3", "normalized.Reporter.intensity.4", "normalized.Reporter.intensity.5",
#                                     "normalized.Reporter.intensity.6", "normalized.Reporter.intensity.7", "normalized.Reporter.intensity.8", "normalized.Reporter.intensity.9",
#                                     "Reporter.intensity.count.0", "Reporter.intensity.count.1", "Reporter.intensity.count.2", "Reporter.intensity.count.3", 
#                                     "Reporter.intensity.count.4", "Reporter.intensity.count.5", "Reporter.intensity.count.6", "Reporter.intensity.count.7", 
#                                     "Reporter.intensity.count.8", "Reporter.intensity.count.9", "logFC", "P.Value", "adj.P.Val")]
# write.csv(table.s2, "8. Table.S2.Phosphopeptides identified and quantified after TiO2 enrichment by MS2 method.csv", row.names = FALSE)
# 
# 
# #total table
# selected.genes <- exprs(exampleSet)
# selected.genes.p.values <- topTable(fit, coef=2, number = nrow(selected.genes))
# selected.genes <- merge(selected.genes, selected.genes.p.values, by = "row.names")
# 
# #significnat table
# significant.genes <- exprs(exampleSet.selected)
# significant.genes.p.values <- selected.genes.p.values[rownames(significant.genes), ]
# significant.genes <- merge(significant.genes, significant.genes.p.values, by = "row.names")
# 
# #volcano plot
# statistics <- topTable(fit, number = nrow(intensities))
# AcetylK.intensities <- merge(AcetylK.intensities, statistics, by = "row.names")
# significant.list <- volcano.plot(AcetylK.intensities, cutoff = 2.30103, filename = "pYIP",
#                                  colselect = "minus.log10.adj.P.Val", colplotx = "logFC",
#                                  colploty = "minus.log10.adj.P.Val", x_lab = "Log 2 Fold Change Her2+/Her2-",
#                                  y_lab = "- log10 adjust P Value", 
#                                  text_label = "Row.names")
# 
# quantified.AcetylK75.ms2 <- AcetylK.intensities;
# 
# 
# log10.AcetylK.intensities <- as.data.frame(AcetylK.intensities[, 142:151]);
# colnames(log10.AcetylK.intensities) <- c("126", "127N", "127C", "128N", "128C",
#                                          "129N", "129C", "130N", "130C", "131");
# pairs.plot(log10.AcetylK.intensities, "PhosphoSTY");
# plot.correlation(log10.AcetylK.intensities, "126", "129C", "PhosphoSTY");
# plot.correlation(log10.AcetylK.intensities, "128N", "130N", "PhosphoSTY");
# plot.correlation(log10.AcetylK.intensities, "128C", "130C", "PhosphoSTY");
# plot.correlation(log10.AcetylK.intensities, "129N", "131", "PhosphoSTY");
# 
# 
# 
# setwd("R:/neubertlabspace/Rawdata_Neubertlab/FangkeHuang/search results/CPTAC2016GrantApplication/Pro&Phospho_CPTAC_MS3/combined/txt");
# 
# AcetylK <- read.table.PhosphoSites(30000);
# filtered.AcetylK <- filter.proteinGroup.TMT(AcetylK);
# quantified.AcetylK <- TMT.cleanup.10plex(filtered.AcetylK);
# AcetylK75 <- subset(filtered.AcetylK, Localization.prob >= 0.75);
# quantified.AcetylK75 <- TMT.cleanup.10plex(AcetylK75);
# 
# AcetylK.ratios <- get.TMT.ratios.10plex2(quantified.AcetylK75);
# check.log2.ratios.distribution.10plex2(AcetylK.ratios, "PhosphoSTY");
# AcetylK.intensities <- get.TMT.intensities.10plex(AcetylK.ratios);
# check.log10.intensity.distribution.10plex(AcetylK.intensities, "PhosphoSTY");
# Ratios.vs.Intensity.10plex2(AcetylK.intensities, "PhosphoSTY");
# Ratios.vs.normalized.Intensity.10plex2(AcetylK.intensities, "PhosphoSTY");
# 
# Position.within.protein <- get.first.IDs(AcetylK.intensities[, "Positions.within.proteins"]);
# Genename <- get.first.IDs(AcetylK.intensities[, "Gene.names"]);
# Position.within.protein <- as.numeric(Position.within.protein);
# site.notation <- paste(Genename, AcetylK.intensities$Amino.acid, sep = ".");
# site.notation <- paste(site.notation, as.character(Position.within.protein), sep = "");
# AcetylK.intensities <- cbind(AcetylK.intensities, site.notation);
# 
# AcetylK.intensities <- AcetylK.intensities[!duplicated(site.notation),]
# rownames(AcetylK.intensities) <- AcetylK.intensities$site.notation
# intensities <- AcetylK.intensities[, 122:131]
# colnames(intensities) <- c("AU565.1", "MDA-MB-231", "SKBR3/CAMA1", "T47D.1", "HCC1954.1",
#                            "HCC1500.1", "AU565.2", "T47D.2", "HCC1954.2", "HCC1500.2")
# library(ALL)
# intensities <- as.matrix(intensities)
# cluster.heatmap(intensities)
# 
# 
# quantified.AcetylK75.ms3 <- AcetylK.intensities;
# 
# log10.AcetylK.intensities <- as.data.frame(AcetylK.intensities[, 142:151]);
# colnames(log10.AcetylK.intensities) <- c("a126", "n127", "c127", "n128", "c128",
#                                          "n129", "c129", "n130", "c130", "a131");
# pairs.plot(log10.AcetylK.intensities, "PhosphoSTY");
# plot.correlation(log10.AcetylK.intensities, "a126", "c129", "PhosphoSTY");
# plot.correlation(log10.AcetylK.intensities, "n128", "n130", "PhosphoSTY");
# plot.correlation(log10.AcetylK.intensities, "c128", "c130", "PhosphoSTY");
# plot.correlation(log10.AcetylK.intensities, "n129", "a131", "PhosphoSTY");
# 
# #overlap between the two
# AcetylK.overlap <- intersect(quantified.AcetylK75.ms2$site.notation,
#                              quantified.AcetylK75.ms3$site.notation);
# #merge two table together
# AcetylK.merge <- merge(quantified.AcetylK75.ms2, 
#                        quantified.AcetylK75.ms3, 
#                        by = "site.notation");
# write.csv(AcetylK.merge, "AcetylK.merge.csv");
# 
# #calculate ratios: channel 1 over channel 2; #plot the correlation of the ratios
# ms2.ratio.1.2 <- AcetylK.merge$normalized.Reporter.intensity.1.x / AcetylK.merge$normalized.Reporter.intensity.2.x;
# ms3.ratio.1.2 <- AcetylK.merge$normalized.Reporter.intensity.1.y / AcetylK.merge$normalized.Reporter.intensity.2.y;
# AcetylK.merge <- cbind(AcetylK.merge, cbind(ms2.ratio.1.2, ms3.ratio.1.2));
# plot.correlation(AcetylK.merge, "ms2.ratio.1.2", "ms3.ratio.1.2", "AcetylK.ratio.1.2", 
#                  lab.note = "", diagonal = FALSE, log.2 = FALSE);
# # plot.ratio.ratio.intensity(AcetylK.merge, 
# #                            "normalized.Reporter.intensity.1.x", 
# #                            "normalized.Reporter.intensity.2.x",
# #                            "normalized.Reporter.intensity.1.y",
# #                            "normalized.Reporter.intensity.1.y",
# #                            "ms2.ratio.1.2", "ms3.ratio.1.2",
# #                            "ms2 ratio over ms3 ratio (log2)",
# #                            dynamic.x.range = c(-4, 4));
# 
# ms2.ratio.06.59 <- (AcetylK.merge$normalized.Reporter.intensity.0.x + AcetylK.merge$normalized.Reporter.intensity.6.x) / (AcetylK.merge$normalized.Reporter.intensity.5.x + AcetylK.merge$normalized.Reporter.intensity.9.x);
# ms3.ratio.06.59 <- (AcetylK.merge$normalized.Reporter.intensity.0.y + AcetylK.merge$normalized.Reporter.intensity.6.y) / (AcetylK.merge$normalized.Reporter.intensity.5.y + AcetylK.merge$normalized.Reporter.intensity.9.y);
# AcetylK.merge <- cbind(AcetylK.merge, cbind(ms2.ratio.06.59, ms3.ratio.06.59));
# plot.correlation(AcetylK.merge, "ms2.ratio.06.59", "ms3.ratio.06.59", "AcetylK.ratio.06.59", 
#                  lab.note = "", diagonal = FALSE, log.2 = FALSE);
# # plot.ratio.ratio.intensity(AcetylK.merge, 
# #                            "normalized.Reporter.intensity.0.x", 
# #                            "normalized.Reporter.intensity.5.x",
# #                            "normalized.Reporter.intensity.0.y",
# #                            "normalized.Reporter.intensity.5.y",
# #                            "ms2.ratio.06.59", "ms3.ratio.06.59",
# #                            "ms2 ratio over ms3 ratio (log2)",
# #                            dynamic.x.range = c(-4, 4));
# 
# ms2.ratio.37.48 <- (AcetylK.merge$normalized.Reporter.intensity.3.x + AcetylK.merge$normalized.Reporter.intensity.7.x) / (AcetylK.merge$normalized.Reporter.intensity.4.x + AcetylK.merge$normalized.Reporter.intensity.8.x);
# ms3.ratio.37.48 <- (AcetylK.merge$normalized.Reporter.intensity.3.y + AcetylK.merge$normalized.Reporter.intensity.7.y) / (AcetylK.merge$normalized.Reporter.intensity.4.y + AcetylK.merge$normalized.Reporter.intensity.8.y);
# AcetylK.merge <- cbind(AcetylK.merge, cbind(ms2.ratio.37.48, ms3.ratio.37.48));
# plot.correlation(AcetylK.merge, "ms2.ratio.37.48", "ms3.ratio.37.48", "AcetylK.ratio.37.48", 
#                  lab.note = "", diagonal = FALSE, log.2 = FALSE);
# # plot.ratio.ratio.intensity(AcetylK.merge, 
# #                            "normalized.Reporter.intensity.3.x", 
# #                            "normalized.Reporter.intensity.4.x",
# #                            "normalized.Reporter.intensity.3.y",
# #                            "normalized.Reporter.intensity.4.y",
# #                            "ms2.ratio.37.48", "ms3.ratio.37.48",
# #                            "ms2 ratio over ms3 ratio (log2)",
# #                            dynamic.x.range = c(-4, 4));
# 
# write.csv(AcetylK.merge, "merged.phosphosites.csv")
# 
