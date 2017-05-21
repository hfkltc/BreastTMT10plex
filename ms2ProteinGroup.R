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