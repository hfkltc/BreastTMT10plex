#SILAC Analysis Functions

#read Phospho (STY)Sites.txt file
read.table.PhosphoSites <- function (n_row) {
  mydata <- read.table("Phospho (STY)Sites.txt", header = TRUE, sep="\t", 
                       quote="\"", nrows = n_row, fill = TRUE, comment.char="")
  mydata
}

#read peptides.txt file
read.table.peptides <- function (n_row) {
  mydata <- read.table("peptides.txt", header = TRUE, sep="\t", 
                       quote="\"", nrows = n_row, fill = TRUE, comment.char="")
  mydata
}

#read evidance.txt file
read.table.evidence <- function (n_row) {
  mydata <- read.table("evidence.txt", header = TRUE, sep="\t", 
                       quote="\"", nrows = n_row, fill = TRUE, comment.char="")
  mydata
}

#read evidance.txt file
read.table.proteinGroup <- function (n_row) {
  mydata <- read.table("proteinGroups.txt", header = TRUE, sep="\t", 
                       quote="\"", nrows = n_row, fill = TRUE, comment.char="")
  mydata
}

read.table.AcetylK.Group <- function (n_row) {
  mydata <- read.table("Acetyl (K)Sites.txt", header = TRUE, sep="\t", 
                       quote="\"", nrows = n_row, fill = TRUE, comment.char="")
  mydata
}

#read msms.txt file
read.table.msms <- function (n_row) {
  mydata <- read.table("msms.txt", header = TRUE, sep="\t", 
                       quote="\"", nrows = n_row, fill = TRUE, comment.char="")
  mydata
}

#filter data
#remove contaminants, remoe reverse, filter score

ratioCountForwardReverseScore <- function(df, scoreForRev = 2) {
  ratioCount <- df[, c("Ratio.H.M.count.Forward", "Ratio.H.M.count.Reverse")]
  countScore <- ratioCount$Ratio.H.M.count.Forward + ratioCount$Ratio.H.M.count.Reverse
  df <- cbind(df, countScore)
  df <- df[which(df$countScore > scoreForRev),]
  df
}

ratioCountForwardReverseScore2 <- function(df) {
  df <- df[which(df$Ratio.H.M.count.Forward > 1),]
  df <- df[which(df$Ratio.H.M.count.Reverse > 1),]
  df
}

ratioVariability <- function(df, variability = 50) {
  df <- df[which(df$Ratio.H.M.variability.....Forward <= variability),]
  df <- df[which(df$Ratio.H.M.variability.....Reverse <= variability),]
  df
}

ratioVariability2 <- function(df, variability = 50) {
  df <- df[which(df$Ratio.H.M.variability.....Forward1 <= variability),]
  df <- df[which(df$Ratio.H.M.variability.....Reverse1 <= variability),]
  df <- df[which(df$Ratio.H.M.variability.....Forward2 <= variability),]
  df <- df[which(df$Ratio.H.M.variability.....Reverse2 <= variability),]
  df
}

filter.PhosphoSites <- function (mydata, score = 80) {
  mydata <- mydata[which(mydata$Reverse != "+"),]
  mydata <- mydata[which(mydata$Potential.contaminant != "+"),]
  mydata <- mydata[which(mydata$Score >= score),]
  mydata
}

filter.TMT.PhosphoSites <- function (mydata, score = 0) {
  mydata <- mydata[which(mydata$Reverse != "+"),]
  mydata <- mydata[which(mydata$Potential.contaminant != "+"),]
  mydata <- mydata[which(mydata$Score >= score),]
  mydata
}

filter.peptides <- function (mydata, score = 0) {
  mydata <- mydata[which(mydata$Reverse != "+"),]
  mydata <- mydata[which(mydata$Contaminant != "+"),]
  mydata <- mydata[which(mydata$Score >= score),]
  mydata
}

filter.TMT.peptides <- function (mydata, score = 0) {
  mydata <- mydata[which(mydata$Reverse != "+"),]
  mydata <- mydata[which(mydata$Potential.contaminant != "+"),]
  mydata <- mydata[which(mydata$Score >= score),]
  mydata
}

filter.evidence <- function (mydata, score = 0) {
  mydata <- mydata[which(mydata$Reverse != "+"),]
  mydata <- mydata[which(mydata$Contaminant != "+"),]
  mydata <- mydata[which(mydata$Score >= score),]
  mydata
}

filter.TMT.evidence <- function (mydata) {
  mydata <- mydata[which(mydata$Reverse != "+"),]
  mydata <- mydata[which(mydata$Potential.contaminant != "+"),]
  mydata
}

filter.proteinGroup <- function (mydata) {
  mydata <- mydata[which(mydata$Reverse != "+"),]
  mydata <- mydata[which(mydata$Potential.contaminant != "+"),]
  mydata <- mydata[which(mydata$Only.identified.by.site != "+"),]
  mydata
}

filter.proteinGroup.TMT <- function (mydata) {
  if ("+" %in% mydata$Only.identified.by.site) {
    mydata <- mydata[which(mydata$Only.identified.by.site != "+"),]};
  if ("+" %in% mydata$Reverse) {
    mydata <- mydata[which(mydata$Reverse != "+"),]};
  if ("+" %in% mydata$Potential.contaminant) {
    mydata <- mydata[which(mydata$Potential.contaminant != "+"),]}
  mydata
}

filter.phosphosites.TMT <- function (mydata) {
  if ("+" %in% mydata$Only.identified.by.site) {
    mydata <- mydata[which(mydata$Only.identified.by.site != "+"),]};
  if ("+" %in% mydata$Reverse) {
    mydata <- mydata[which(mydata$Reverse != "+"),]};
  if ("+" %in% mydata$Potential.contaminant) {
    mydata <- mydata[which(mydata$Potential.contaminant != "+"),]}
  mydata
}

filter.proteinGroup2 <- function (mydata) {
  mydata <- mydata[which(mydata$Reverse != "+"),]
  mydata <- mydata[which(mydata$Contaminant != "+"),]
  mydata <- mydata[which(mydata$Only.identified.by.site != "+"),]
  mydata
}

#transform ratios and intensities
transform.ratio.intensity <- function (mydata) {
  
  ratio_col <- grep("Ratio.+", colnames(mydata), value = TRUE)
  ratio_col <- ratio_col[!grepl("*localized*", ratio_col)]
  
  forward_ratio_col <- grep("*Forward*", ratio_col, value = TRUE)
  reverse_ratio_col <- grep("*Reverse*", ratio_col, value = TRUE)
  
  mydata[forward_ratio_col] <- log2(mydata[forward_ratio_col])
  mydata[reverse_ratio_col] <- log2(mydata[reverse_ratio_col])
  
  inten_col <- grep("Intensity.+", colnames(mydata), value = TRUE)
  mydata[inten_col] <- log10(mydata[inten_col])
  
  mydata
}

transform.tmt.intensity <- function (mydata) {
  
  tmt_col <- grep("Reporter+", colnames(mydata), value = TRUE);
  
  mydata[tmt_col] <- log2(mydata[tmt_col]);
  mydata
}

# generate a numeric IDs from a IDs from Maxquant file
# for example, input IDs are Phospho(STY)Site id column
get.first.IDs <- function(IDs, mark = ";") {
  new.IDs <- NULL;
  n <- length(IDs);
  for (i in 1:n) {
    if (grepl(mark, IDs[i])) {
      new.IDs <- c(new.IDs, strsplit(toString(IDs[i]), mark, fixed = TRUE)[[1]][1]);
    }
    else {
      new.IDs <- c(new.IDs, toString(IDs[i]));
    }
  }
  new.IDs
}

get.unique.numeric.IDs <- function (IDs) {
  IDs <- unique(IDs)
  n <- length(IDs)
  new.IDs <- numeric()
  for (i in 1:n) {
    if (grepl(";", IDs[i])) {
      new.IDs <- c(new.IDs, as.numeric(unlist(strsplit(toString(IDs[i]), ";"))))
    }
    else {
      new.IDs <- c(new.IDs, as.numeric(IDs[i]))
    }
  }
  new.IDs <- unique(new.IDs)
  new.IDs <- sort(new.IDs)
}

get.unique.gene.names <- function (gene.names) {
  gene.names <- unique(gene.names)
  n <- length(gene.names)
  gene.names.david <- vector()
  for (i in 1:n) {
    if (grepl(";", gene.names[i])) {
      each.entry <- unlist(strsplit(toString(gene.names[i]), ";"))
      each.entry <- unique(each.entry)
      gene.names.david <- c(gene.names.david, each.entry)
    }
    else {
      each.entry <- toString(gene.names[i])
      gene.names.david <- c(gene.names.david, each.entry)
    }
  }
  gene.names.david <- unique(gene.names.david)
  cat("Number of unique genes is ", length(gene.names.david), "\n")
  gene.names.david <- sort(gene.names.david)
}

get.unique.gene.names2 <- function (gene.names) {
#   gene.names <- unique(gene.names)
  n <- length(gene.names)
  gene.names.david <- vector()
  for (i in 1:n) {
    if (grepl(";", gene.names[i])) {
      each.entry <- unlist(strsplit(toString(gene.names[i]), ";"))[1]
#       each.entry <- unique(each.entry)
      gene.names.david <- c(gene.names.david, each.entry)
    }
    else {
      each.entry <- toString(gene.names[i])
      gene.names.david <- c(gene.names.david, each.entry)
    }
  }
#   gene.names.david <- unique(gene.names.david)
  cat("Number of unique genes is ", length(gene.names.david), "\n")
#   gene.names.david <- sort(gene.names.david)
  gene.names.david
}

get.fraction.no <- function(a_str, n) {
  a_str <- toString(a_str)
  fraction <- substr(a_str, nchar(a_str) - n + 1, nchar(a_str))
#   if (grepl("_", fraction) || grepl("F", fraction)) {
#     fraction <- substr(fraction, 2, 2)
#   }
#   else {fraction}
}

phospho.fraction.table <- function (mydata) {
  phos.frac.table <- table(mydata$Raw.file, mydata$Phospho..STY.) 
  raw_files <- rownames(phos.frac.table)
  fractions <- lapply(raw_files, get.fraction.no, 2)
  total.evidence <- table(mydata$Raw.file)
  #   rownames(phos.frac.table) <- fractions
  phosphopep <- total.evidence - phos.frac.table[,1]
  phos.frac.table <- cbind(phosphopep, phos.frac.table)
  colnames(phos.frac.table)[1] <- "Phosphopeptide"
  
  phos.frac.table <- cbind(total.evidence, phos.frac.table)
  colnames(phos.frac.table)[1] <- "Total.Peptides"
  
  phos.frac.table <- cbind(as.numeric(fractions), phos.frac.table)
  colnames(phos.frac.table)[1] <- "Fraction.No"
  
  phos.frac.table <- as.table(phos.frac.table)
  phos.frac.table <- phos.frac.table[order(phos.frac.table[, 1]), ]
  rownames(phos.frac.table) <- as.vector(phos.frac.table[, 1])
  
  phos.frac.table
}

# plot nonphosphopeptides versus phosphopeptides in different fractions
plot.phos.nonphos <- function (to.plot.table, filename = "") {
  
  filename_tiff <- paste(filename, ".tiff", sep = "")
  tiff(file = filename_tiff)
  barplot(to.plot.table, col = c("deepskyblue3", "goldenrod1"), 
          beside=TRUE, cex.names = 0.7, space = c(0, 0.5), 
          main = "Phospho- and Nonphospho- peptides in each fraction",
          xlab = "Fractions", ylab = "Frequency")
  legend("topright", legend = rownames(to.plot.table), 
         col = c("deepskyblue3", "goldenrod1"),
         fill = c("deepskyblue3", "goldenrod1"))
  dev.off()
  
  filename_png <- paste(filename, ".png", sep = "")
  png(file = filename_png)
  barplot(to.plot.table, col = c("deepskyblue3", "goldenrod1"), 
          beside=TRUE, cex.names = 0.7, space = c(0, 0.5), 
          main = "Phospho- and Nonphospho- peptides in each fraction",
          xlab = "Fractions", ylab = "Frequency")
  legend("topright", legend = rownames(to.plot.table), 
         col = c("deepskyblue3", "goldenrod1"), 
         fill = c("deepskyblue3", "goldenrod1"))
  dev.off()
  
  filename_eps <- paste(filename, ".eps", sep = "")
  setEPS(horizontal = FALSE, paper = "special", width = 8, height = 8)
  postscript(file = filename_eps)
  barplot(to.plot.table, col = c("deepskyblue3", "goldenrod1"), 
          beside=TRUE, cex.names = 0.8, space = c(0, 0.5), 
          main = "Phospho- and Nonphospho- peptides in each fraction",
          xlab = "Fractions", ylab = "Frequency")
  legend("topright", legend = rownames(to.plot.table), 
         col = c("deepskyblue3", "goldenrod1"), 
         fill = c("deepskyblue3", "goldenrod1"))
  dev.off()
}

plot.phos.percent <- function (to.plot.table, filename = "") {
    
  filename_tiff <- paste(filename, ".tiff", sep = "")
  tiff(file = filename_tiff)
  barplot(to.plot.table, col = c("deepskyblue3", "goldenrod1"), 
          beside=TRUE, cex.names = 0.7, space = c(0, 0.5), 
          main = "Total and phospho- peptides in each fraction",
          xlab = "Fractions", ylab = "Frequency")
  legend("topright", legend = rownames(to.plot.table), 
         col = c("deepskyblue3", "goldenrod1"),
         fill = c("deepskyblue3", "goldenrod1"))
  dev.off()
  
  filename_png <- paste(filename, ".png", sep = "")
  png(file = filename_png)
  barplot(to.plot.table, col = c("deepskyblue3", "goldenrod1"), 
          beside=TRUE, cex.names = 0.7, space = c(0, 0.5), 
          main = "Total and phospho- peptides in each fraction",
          xlab = "Fractions", ylab = "Frequency")
  legend("topright", legend = rownames(to.plot.table), 
         col = c("deepskyblue3", "goldenrod1"), 
         fill = c("deepskyblue3", "goldenrod1"))
  dev.off()
}

plot.phos.loca.prob <- function (PhosphoSites) {
  tiff(file = "phos_loca_prob.tiff");
  hist(PhosphoSites$Localization.prob,
       xlab = "Probability",
       main = "PhosphoSite Localization Probabilities Distribution");
  dev.off();
  
  png(file = "phos_loca_prob.png")
  hist(PhosphoSites$Localization.prob, 
       xlab = "Probability",
       main = "PhosphoSite Localization Probabilities Distribution")
  dev.off()
  
  setEPS(horizontal = FALSE, paper = "special", width = 8, height = 8)
  postscript(file = "phos_loca_prob.eps")
  hist(PhosphoSites$Localization.prob, 
       xlab = "Probability",
       main = "PhosphoSite Localization Probabilities Distribution")
  dev.off()
}


plot.evi.peptide <- function (to.plot.table, filename = "") {
  
  filename_tiff <- paste(filename, ".tiff", sep = "")
  tiff(file = filename_tiff)
  barplot(to.plot.table, col = c("deepskyblue3", "goldenrod1"), 
          beside=TRUE, cex.names = 0.7, space = c(0, 0.5), 
          main = "Number of evidence and peptides",
          xlab = "Factors", ylab = "Frequency")
  legend("topright", legend = rownames(to.plot.table), 
         col = c("deepskyblue3", "goldenrod1"),
         fill = c("deepskyblue3", "goldenrod1"))
  dev.off()
  
  filename_png <- paste(filename, ".png", sep = "")
  png(file = filename_png)
  barplot(to.plot.table, col = c("deepskyblue3", "goldenrod1"), 
          beside=TRUE, cex.names = 0.7, space = c(0, 0.5), 
          main = "Number of evidence and peptides",
          xlab = "Factors", ylab = "Frequency")
  legend("topright", legend = rownames(to.plot.table), 
         col = c("deepskyblue3", "goldenrod1"), 
         fill = c("deepskyblue3", "goldenrod1"))
  dev.off()
}

expand_data_frame <- function(df, col_fix, col_expand) {
  df <- df[, c(col_fix, col_expand)]
  expand_df <- data.frame()
  df_colnames <- colnames(df)
  for (i in 1:nrow(df)) {
    if (grepl(";", df[, col_expand][i])) {
      proID <- unlist(strsplit(toString(df[, col_expand][i]), ";"))
      frac_proID <- as.data.frame(cbind(toString(df[, col_fix][i]), proID))
      colnames(frac_proID) <- df_colnames
      expand_df <- rbind(expand_df, frac_proID)
    }
    else {
      expand_df <- rbind(expand_df, df[i, ])
    }
  }
  expand_df
}

compress_pro_frac <- function(frac_pro) {
  n <- 1
  ProGrID <- frac_pro[1, 1]
  pro_frac <- data.frame()
  m <- 0
  fraction <- ""
  for (i in 1:(nrow(frac_pro) + 1)) {
    #cat("i is", i, "n is", n, "ProGrID is", ProGrID, "frac_pro[i, 1] is", frac_pro[i, 1], "\n")
    if (frac_pro[i, 1] == ProGrID && i != nrow(frac_pro) + 1) {
      m <- m + 1
      fraction <- paste(fraction, frac_pro[i, 2], sep = ",")
    }
    else if (frac_pro[i, 1] != ProGrID && i != nrow(frac_pro) + 1) {
      pro_frac[n, 1] <- ProGrID
      pro_frac[n, 2] <- substr(fraction, 2, nchar(fraction))
      pro_frac[n, 3] <- m
      ProGrID <- frac_pro[i, 1]
      m <- 1
      fraction <- ""
      fraction <- paste(fraction, frac_pro[i, 2], sep = ",")
      n <- n + 1
    }
    else if (i == nrow(frac_pro) + 1) {
      #cat("I'm here")
      pro_frac[n, 1] <- ProGrID
      pro_frac[n, 2] <- substr(fraction, 2, nchar(fraction))
      pro_frac[n, 3] <- m      
    }
  }
  colnames(pro_frac) <- c("ProteinGroupID", "FractionCounts", "Counts")
  pro_frac[,c(2,3)]
}


get_evi_pep_pho_split <- function(split_evi, split_f, filename) {
  evi_n <- vector()
  pep_n <- vector()
  phopep_n <- vector()
    
  for (i in 1:length(split_f)) {
    var_name <- toString(split_f[i])
    cat(var_name, "\n")
    assign(var_name, as.data.frame(split_evi[i]))
    evi_n[i] <- nrow(get(var_name))
    pep_n[i] <- length(unique(get(var_name)[,6]))
    phopep_n[i] <- length(unique(get(var_name)[which(get(var_name)[,1] >0), ][,6]))
  }
  
  evi_per_pep <- evi_n/pep_n
  pho_percent <- phopep_n/pep_n
  evi_pep_pho <- cbind(evi_n, pep_n, phopep_n, evi_per_pep, pho_percent)
  rownames(evi_pep_pho) <- as.vector(split_f)
  colnames(evi_pep_pho) <- c("Evidence", "Peptides", "Phosphopeptides", 
                             "Evidence.Per.Peptides", "Enrichment.Efficiency")
  
  write.csv(format(evi_pep_pho, digits = 2, scientific = FALSE), 
            filename, row.names = TRUE)
  
  evi_pep_pho
}

plot.density.comparison <- function(var1, by.var2, df) {
  by.var2.f <- factor(df[,c(by.var2)]);
  filename <- paste(var1, " Disctribution by ", sep = "");
  filename <- paste(filename, by.var2, sep = "");
  png(file = paste(filename, ".png", sep = ""));
  sm.density.compare(df[, c(var1)], df[, c(by.var2)], xlab=var1, col = c("deepskyblue3", "goldenrod1"));
  title(main=filename);
  legend("topright", legend = levels(by.var2.f), 
         fill=c("deepskyblue3", "goldenrod1"), 
         col = c("deepskyblue3", "goldenrod1"));
  dev.off();
}

countCharOccurrences <- function(str, char) {
   str2 <- gsub(char, "", str);
   return (nchar(str) - nchar(str2) + 1)
}

saveGeneEntries <- function(df, geneName) {
  geneEntries <- df[which(df$Gene.names == geneName), ]
  #geneEntries <- df[which(grep(df$Gene.names, geneName)),]
  filename <- paste(geneName, ".csv", sep = "")
  geneEntries <- geneEntries[, c("Ratio.H.M.normalized.Forward", 
                                 "Ratio.H.M.normalized.Reverse", 
                                 "Leading.proteins", "Protein.names", 
                                 "Gene.names", "Score.diff", "Score", 
                                 "id", "Number.of.Phospho..STY.", "Amino.acid", 
                                 "Modified.sequence", 
                                 "Phospho..STY..Probabilities", "Positions")];
  write.csv(geneEntries, filename, row.names = FALSE)
}

plotRatioCorrelation <- function(df) {
  ratios <- df[complete.cases(df[, c("Ratio.H.M.normalized.Forward", 
                                     "Ratio.H.M.normalized.Reverse")]), ];
  ratios[, c("Ratio.H.M.normalized.Reverse")] <- 1/ratios[, c("Ratio.H.M.normalized.Reverse")]
  tf.df <- transform.ratio.intensity(ratios);
  ratios.df <- tf.df[, c("Ratio.H.M.normalized.Forward", 
                                 "Ratio.H.M.normalized.Reverse")];
  completeRatios.df <- ratios.df[complete.cases(ratios.df), ];
  cat("Number of total ratios is ", nrow(completeRatios.df), "\n");
  cat("The Pearson correlation is ", cor(completeRatios.df)[1,2], "\n");
  colnames(completeRatios.df) <- c("F", "R");
  par(mfcol = c(1, 1), mar = c(5, 5, 4, 4));
  png(filename = "ratio_correlation.png", width = 480, height = 500);
  plot(completeRatios.df$F, completeRatios.df$R, 
       xlab = "Forward Raio (log2)", ylab = "Reverse Raio (log2)", pch = 16, 
       main = "Forward vs Reverse", cex.main = 2, cex.lab = 1.5, cex.axis = 1.5,
       xlim=c(-6, 6), ylim=c(-6, 6))
  dev.off()
}

plotRatioCorrelation2 <- function(df, filename, cutoff) {
  ratios <- df[complete.cases(df[, c("Ratio.H.M.normalized.Forward", 
                                             "Ratio.H.M.normalized.Reverse")]), ];
  ratios[, c("Ratio.H.M.normalized.Reverse")] <- 1/ratios[, c("Ratio.H.M.normalized.Reverse")]
  tf.df <- transform.ratio.intensity(ratios);
  if (filename == "Phospho (STY)Sites") {
    ratios <- ratios[, c("Ratio.H.M.normalized.Forward", 
                     "Ratio.H.M.normalized.Reverse", 
                     "Ratio.H.M.variability.....Forward",
                     "Ratio.H.M.variability.....Reverse", 
                     "Ratio.H.M.count.Forward", "Ratio.H.M.count.Reverse", 
                     "Leading.proteins", "Protein.names", "Gene.names", 
                     "Score", "id", "Number.of.Phospho..STY.", 
                     "Amino.acid", "Modified.sequence", 
                     "Phospho..STY..Probabilities", "Positions")];
  }
  else if (filename == "proteinGroups") {
    ratios <- ratios[, c("Ratio.H.M.normalized.Forward", 
                         "Ratio.H.M.normalized.Reverse", 
                         "Ratio.H.M.variability.....Forward",
                         "Ratio.H.M.variability.....Reverse", 
                         "Ratio.H.M.count.Forward", "Ratio.H.M.count.Reverse", 
                         "Protein.names", "Gene.names", "Unique.peptides", 
                         "Sequence.coverage....", "id")];
  }
  
  else if (filename == "peptides") {
    ratios <- ratios[, c("Ratio.H.M.normalized.Forward", 
                         "Ratio.H.M.normalized.Reverse", 
                         "Ratio.H.M.variability.....Forward",
                         "Ratio.H.M.variability.....Reverse", 
                         "Ratio.H.M.count.Forward", "Ratio.H.M.count.Reverse", 
                         "Protein.names", "Gene.names", "Score", "Intensity",
                         "Unique..Groups.", "id")];
  }
  
  ratios.df <- tf.df[, c("Ratio.H.M.normalized.Forward", 
                         "Ratio.H.M.normalized.Reverse")];
  colnames(ratios.df) <- c("F", "R");
  ratioDistance <- sqrt(ratios.df$F^2 + ratios.df$R^2);
  countLargeChange <- length(ratioDistance[ratioDistance >= 1]);
  cat("Number of total large changes is ", countLargeChange, "\n");
  
    
  ratios <- cbind(ratios.df, ratios)
  write.csv(ratios, "ratiosAndTransformedRatios.csv", row.names = FALSE)
  
  par(mfcol = c(1, 1), mar = c(5, 5, 4, 4));
  png(filename = "ratio_correlation.png", width = 480, height = 500);
  plot(ratios$F, ratios$R, 
       xlab = "Forward Raio (log2)", ylab = "Reverse Raio (log2)", pch = 16, 
       main = "Forward vs Reverse", cex.main = 2, cex.lab = 1.5, cex.axis = 1.5,
       xlim=c(-6, 6), ylim=c(-6, 6), col = "deepskyblue4");
  up.PhosphoSites <- ratios[which(ratios$F >= cutoff & ratios$R >= cutoff), ];
  down.PhosphoSites <- ratios[which(ratios$F <= -cutoff & ratios$R <= -cutoff), ];
  consistant.PhosphoSites <- rbind(up.PhosphoSites, down.PhosphoSites);
  cat("Number of consistant changes is ", nrow(consistant.PhosphoSites), "\n");
  write.csv(consistant.PhosphoSites, "SILAC_consistant_changes.csv",
            row.names = FALSE)
  points(consistant.PhosphoSites$F, consistant.PhosphoSites$R, pch = 16, col = "firebrick1");
  
  upForDownRev <- ratios[which(ratios$F >= cutoff & ratios$R <= -cutoff), ]
  upRevDownFor <- ratios[which(ratios$F <= -cutoff & ratios$R >= cutoff), ]
  inconsistant.PhosphoSites <- rbind(upForDownRev, upRevDownFor)
  write.csv(inconsistant.PhosphoSites, "SILAC_inconsistant_changes.csv", 
            row.names = FALSE)
#   points(inconsistant.PhosphoSites$F, inconsistant.PhosphoSites$R, pch = 16, col = "goldenrod1");
#   
  dev.off()
}

plotRatioCorrelation2.pdf <- function(df, filename, cutoff) {
  ratios <- df[complete.cases(df[, c("Ratio.H.M.normalized.Forward", 
                                     "Ratio.H.M.normalized.Reverse")]), ];
  ratios[, c("Ratio.H.M.normalized.Reverse")] <- 1/ratios[, c("Ratio.H.M.normalized.Reverse")]
  tf.df <- transform.ratio.intensity(ratios);
  if (filename == "Phospho (STY)Sites") {
    ratios <- ratios[, c("Ratio.H.M.normalized.Forward", 
                         "Ratio.H.M.normalized.Reverse", 
                         "Ratio.H.M.variability.....Forward",
                         "Ratio.H.M.variability.....Reverse", 
                         "Ratio.H.M.count.Forward", "Ratio.H.M.count.Reverse", 
                         "Leading.proteins", "Protein.names", "Gene.names", 
                         "Score", "id", "Number.of.Phospho..STY.", 
                         "Amino.acid", "Modified.sequence", 
                         "Phospho..STY..Probabilities", "Positions")];
  }
  else if (filename == "proteinGroups") {
    ratios <- ratios[, c("Ratio.H.M.normalized.Forward", 
                         "Ratio.H.M.normalized.Reverse", 
                         "Ratio.H.M.variability.....Forward",
                         "Ratio.H.M.variability.....Reverse", 
                         "Ratio.H.M.count.Forward", "Ratio.H.M.count.Reverse", 
                         "Protein.names", "Gene.names", "Unique.peptides", 
                         "Sequence.coverage....", "id")];
  }
  
  else if (filename == "peptides") {
    ratios <- ratios[, c("Ratio.H.M.normalized.Forward", 
                         "Ratio.H.M.normalized.Reverse", 
                         "Ratio.H.M.variability.....Forward",
                         "Ratio.H.M.variability.....Reverse", 
                         "Ratio.H.M.count.Forward", "Ratio.H.M.count.Reverse", 
                         "Protein.names", "Gene.names", "Score", "Intensity",
                         "Unique..Groups.", "id")];
  }
  
  ratios.df <- tf.df[, c("Ratio.H.M.normalized.Forward", 
                         "Ratio.H.M.normalized.Reverse")];
  colnames(ratios.df) <- c("F", "R");
  ratioDistance <- sqrt(ratios.df$F^2 + ratios.df$R^2);
  countLargeChange <- length(ratioDistance[ratioDistance >= 1]);
  cat("Number of total large changes is ", countLargeChange, "\n");
  
  
  ratios <- cbind(ratios.df, ratios)
  write.csv(ratios, "ratiosAndTransformedRatios.csv", row.names = FALSE)
  
#   par(mfcol = c(1, 1), mar = c(5, 5, 4, 4));
#   pdf(filename = "ratio_correlation.pdf", width = 480, height = 500);
  pdf("ratio_correlation.pdf", width = 480, height = 500);
  par(mfcol = c(1, 1), mar = c(5, 5, 4, 4));
  plot(ratios$F, ratios$R, 
       xlab = "Forward Raio (log2)", ylab = "Reverse Raio (log2)", pch = 16, 
       main = "Forward vs Reverse", cex.main = 2, cex.lab = 1.5, cex.axis = 1.5,
       xlim=c(-6, 6), ylim=c(-6, 6), col = "deepskyblue4");
  up.PhosphoSites <- ratios[which(ratios$F >= cutoff & ratios$R >= cutoff), ];
  down.PhosphoSites <- ratios[which(ratios$F <= -cutoff & ratios$R <= -cutoff), ];
  consistant.PhosphoSites <- rbind(up.PhosphoSites, down.PhosphoSites);
  cat("Number of consistant changes is ", nrow(consistant.PhosphoSites), "\n");
  write.csv(consistant.PhosphoSites, "SILAC_consistant_changes.csv",
            row.names = FALSE)
  points(consistant.PhosphoSites$F, consistant.PhosphoSites$R, pch = 16, col = "firebrick1");
#   print(plot1)
#   print(plot2)
  dev.off()
  
  upForDownRev <- ratios[which(ratios$F >= cutoff & ratios$R <= -cutoff), ]
  upRevDownFor <- ratios[which(ratios$F <= -cutoff & ratios$R >= cutoff), ]
  inconsistant.PhosphoSites <- rbind(upForDownRev, upRevDownFor)
  write.csv(inconsistant.PhosphoSites, "SILAC_inconsistant_changes.csv", 
            row.names = FALSE)
  #   points(inconsistant.PhosphoSites$F, inconsistant.PhosphoSites$R, pch = 16, col = "goldenrod1");
  #   
  dev.off()
}

plotRatioCorrelation3 <- function(df, var.forward, var.reverse) {
  ratios <- df[complete.cases(df[, c("Ratio.H.M.normalized.Forward", 
                                     "Ratio.H.M.normalized.Reverse")]), ];
  ratios[, c("Ratio.H.M.normalized.Reverse")] <- 1/ratios[, c("Ratio.H.M.normalized.Reverse")]
  ratios[, c("Ratio.H.M.normalized.Forward")] <- log2(ratios[, c("Ratio.H.M.normalized.Forward")]);
  ratios[, c("Ratio.H.M.normalized.Reverse")] <- log2(ratios[, c("Ratio.H.M.normalized.Reverse")]);
  cat("Number of total ratios is ", nrow(ratios), "\n");
  
  count1 <- subset(ratios, Ratio.H.M.count.Forward == 1 | Ratio.H.M.count.Reverse == 1);
  cat("Number of either ratio count = 1 is ", nrow(count1), "\n");
  ratios <- subset(ratios, Ratio.H.M.count.Forward > 1 & Ratio.H.M.count.Reverse > 1);
  cat("Number of ratio left is ", nrow(ratios), "\n");
  par(mfcol = c(1, 1), mar = c(5, 5, 4, 4));
  png(filename = "ratio_correlation.png", width = 480, height = 500);
  plot(count1$Ratio.H.M.normalized.Forward, count1$Ratio.H.M.normalized.Reverse, 
       xlab = "Forward Raio (log2)", ylab = "Reverse Raio (log2)", pch = 16, 
       main = "Forward vs Reverse", cex.main = 2, cex.lab = 1.5, cex.axis = 1.5,
       xlim=c(-6, 6), ylim=c(-6, 6), col = "deepskyblue4");
  
  small_var <- subset(ratios, Ratio.H.M.variability.....Forward <= var.forward & Ratio.H.M.variability.....Reverse <= var.reverse);
  cat("Number of both small variability is ", nrow(small_var), "\n");
#   write.csv(consistant.PhosphoSites, "SILAC_consistant_changes.csv",
#             row.names = FALSE)  
  points(small_var$Ratio.H.M.normalized.Forward, small_var$Ratio.H.M.normalized.Reverse, pch = 16, col = "firebrick1");
  
  big_var <- subset(ratios, Ratio.H.M.variability.....Forward > var.forward | Ratio.H.M.variability.....Reverse > var.reverse);
  cat("Number of either big variability is ", nrow(big_var), "\n");
  points(big_var$Ratio.H.M.normalized.Forward, big_var$Ratio.H.M.normalized.Reverse, pch = 16, col = "goldenrod1");
  
  dev.off()
  
  png(filename = "ratio_correlation_count1.png", width = 480, height = 500);
  plot(count1$Ratio.H.M.normalized.Forward, count1$Ratio.H.M.normalized.Reverse, 
       xlab = "Forward Raio (log2)", ylab = "Reverse Raio (log2)", pch = 16, 
       main = "Forward vs Reverse", cex.main = 2, cex.lab = 1.5, cex.axis = 1.5,
       xlim=c(-6, 6), ylim=c(-6, 6), col = "deepskyblue4");
  dev.off()

  png(filename = "ratio_correlation_small_variability.png", width = 480, height = 500);
  plot(small_var$Ratio.H.M.normalized.Forward, small_var$Ratio.H.M.normalized.Reverse, 
       xlab = "Forward Raio (log2)", ylab = "Reverse Raio (log2)", pch = 16, 
       main = "Forward vs Reverse", cex.main = 2, cex.lab = 1.5, cex.axis = 1.5,
       xlim=c(-6, 6), ylim=c(-6, 6), col = "firebrick1");
  dev.off()

  png(filename = "ratio_correlation_big_variability.png", width = 480, height = 500);
  plot(big_var$Ratio.H.M.normalized.Forward, big_var$Ratio.H.M.normalized.Reverse, 
       xlab = "Forward Raio (log2)", ylab = "Reverse Raio (log2)", pch = 16, 
       main = "Forward vs Reverse", cex.main = 2, cex.lab = 1.5, cex.axis = 1.5,
       xlim=c(-6, 6), ylim=c(-6, 6), col = "goldenrod1");
  dev.off()
}

plotRatioCorrelation4 <- function(ratios, sig_updown) {
  ratios[, c("Ratio.H.M.normalized.Forward")] <- log2(ratios[, c("Ratio.H.M.normalized.Forward")]);
  ratios[, c("Ratio.H.M.normalized.Reverse")] <- log2(ratios[, c("Ratio.H.M.normalized.Reverse")]);
  cat("Number of total ratios is ", nrow(ratios), "\n");
  
  sig_updown[, c("Ratio.H.M.normalized.Forward")] <- log2(sig_updown[, c("Ratio.H.M.normalized.Forward")]);
  sig_updown[, c("Ratio.H.M.normalized.Reverse")] <- log2(sig_updown[, c("Ratio.H.M.normalized.Reverse")]);
  cat("Number of updown ratios is ", nrow(sig_updown), "\n");
  
  par(mfcol = c(1, 1), mar = c(5, 5, 4, 4));
  png(filename = "ratio_correlation.png", width = 480, height = 500);
  plot(ratios$Ratio.H.M.normalized.Forward, ratios$Ratio.H.M.normalized.Reverse, 
       xlab = "Forward Raio (log2)", ylab = "Reverse Raio (log2)", pch = 16, 
       main = "Forward vs Reverse", cex.main = 2, cex.lab = 1.5, cex.axis = 1.5,
       xlim=c(-6, 6), ylim=c(-6, 6), col = "deepskyblue4");
  
  points(sig_updown$Ratio.H.M.normalized.Forward, sig_updown$Ratio.H.M.normalized.Reverse, pch = 16, col = "firebrick1");
  
  dev.off()
}

vennDiagram <- function(filename, n_row, idOrQuant) {
  library(plyr); 
  library(venneuler); 
  library(VennDiagram);
  library(rJava); 
  library(grid); 
  library(knitr);
  
  if (filename == "Phospho (STY)Sites" && idOrQuant == "id") {
    df <- read.table.PhosphoSites(n_row);
    df <- filter.PhosphoSites(df, 0);
    df <- df[, c("Localization.prob.Forward", "Localization.prob.Reverse")];
    cat("This is ", filename, " identification\n\n");
  }
  
  else if (filename == "Phospho (STY)Sites" && idOrQuant == "Quant") {
    df <- read.table.PhosphoSites(n_row);
    df <- filter.PhosphoSites(df, 0);
    df <- df[, c("Ratio.H.M.normalized.Forward", "Ratio.H.M.normalized.Reverse")];
    cat("This is ", filename, " Quantitation\n\n");
  }
  
  else if (filename == "proteinGroups" && idOrQuant == "Quant") {
    df <- read.table.proteinGroup(n_row);
    df <- filter.proteinGroup(df);
    df <- df[, c("Ratio.H.M.normalized.Forward", "Ratio.H.M.normalized.Reverse")];
    cat("This is ", filename, " Quantitation\n\n");
  }
  
  else if (filename == "proteinGroups" && idOrQuant == "id") {
    cat("Can't count proteinGroups identification right now!");
  }
  
  df <- as.data.frame(!is.na(df)); 
  colnames(df) <- c("F", "R");
  vennCounts <- ddply(df, .(df$F, df$R), nrow);
  colnames(vennCounts) <- c("F", "R", "Freq");
  cat("This is Venn Diagram counts\n");
  print(vennCounts);
  cat("\n");
  
  F_count <- sum(vennCounts[which(vennCounts$F == TRUE), 3]);
  R_count <- sum(vennCounts[which(vennCounts$R == TRUE), 3]);
  totalCount <- sum(vennCounts[, 3]);
  countTable <- rbind(F_count, R_count, totalCount);
  colnames(countTable) <- c("Freq");
  cat("This is from each file\n");
  print(countTable);
}

get.motif.x.input <- function(sequence.probability) {
  #getting the parentheses and its content
  sequence.probability <- toString(sequence.probability);
  probability <- gsub("[\\(\\)]", "", regmatches(sequence.probability, gregexpr("\\(.*?\\)", sequence.probability))[[1]]);
  probability <- as.numeric(probability);
  probability <- probability > 0.75;
  
  index <- gregexpr("\\(.*?\\)", sequence.probability);
  insert.index <- unlist(index)[probability];
  
  start <- c(1, insert.index);
  end <- c(insert.index - 1, nchar(sequence.probability));
  sel <- cbind(start, end);
  
  strings <- apply(sel, 1, function(x) substr(sequence.probability, x[1], x[2]));
  strings <- paste(strings, collapse='*')
  
  strings <- gsub("\\(.*?\\)", '', strings);
  strings
}

TMT.label.efficiency <- function(reporter, n_row) {
  evidence.le <- read.table.evidence(n_row);
  evidence.le <- filter.TMT.evidence(evidence.le);
  
  le.reporter <- paste("LE", reporter, sep = "");
  evidence <- evidence.le[which(evidence.le$Experiment == le.reporter),];
  
  l.reporter <- paste("TMT6plex.Lys", reporter, sep="");
  l.reporter <- paste(l.reporter, ".Modification", sep="");
  n.reporter <- paste("TMT6plex.Nter", reporter, sep="");
  n.reporter <- paste(n.reporter, ".Modification", sep="");
  evidence <- evidence[, c("Sequence", "Length", "Acetyl..Protein.N.term.",
                                   l.reporter, n.reporter,
                                   "Missed.cleavages", "Proteins", "Leading.Proteins",
                                   "Leading.Razor.Protein", "Gene.Names", 
                                   "Protein.Names", "Type", "Raw.file", "Experiment",
                                   "Mass", "Protein.group.IDs", "Peptide.ID",
                                   "Mod..peptide.ID")];
  
  evidence.sequence <- sort(unique(evidence$Sequence));
  
  peptides <- NULL;
  for (i in 1:length(evidence.sequence)) {
    cat ("The calculation reaches ", i/length(evidence.sequence), "percent\n");
    each.sequence <- toString(evidence.sequence[i]);
    evidences.each.sequence <- evidence[which(evidence$Sequence == each.sequence),];
    if (nrow(evidences.each.sequence) == 1) {
      peptides <- rbind(peptides, evidences.each.sequence);
    }
    if (nrow(evidences.each.sequence) > 1) {
      evidences.each.sequence[1, c("Acetyl..Protein.N.term.")] <- sum(evidences.each.sequence[, c("Acetyl..Protein.N.term.")]);
      evidences.each.sequence[1, l.reporter] <- sum(evidences.each.sequence[, l.reporter]);
      evidences.each.sequence[1, n.reporter] <- sum(evidences.each.sequence[, n.reporter]);
      peptides <- rbind(peptides, evidences.each.sequence[1,]);
    } 
  }
  
  total.peptides <- nrow(peptides);
  cat("The total number of peptide with", reporter, "is", total.peptides, "\n");
  
  peptides.Nacetyl <- peptides[which(peptides$Acetyl..Protein.N.term. > 0),];
  cat("The number of peptide with N-acetyl is", nrow(peptides.Nacetyl), "\n");
  cat("It is ", nrow(peptides.Nacetyl)/total.peptides, "of total peptides\n\n")
  
  peptides.Lreporter <- peptides[which(peptides[, l.reporter] > 0),];
  cat("The number of peptide with Lys modification is", nrow(peptides.Lreporter), "\n");
  cat("It is ", nrow(peptides.Lreporter)/total.peptides, "of total peptides\n\n")
  
  peptides.Nreporter <- peptides[which(peptides[, n.reporter] > 0),];
  cat("The number of peptide with Nter modification is", nrow(peptides.Nreporter), "\n");
  cat("It is ", nrow(peptides.Nreporter)/total.peptides, "of total peptides\n\n")
  
  peptides.LNreporter <- peptides[which(peptides[, l.reporter] > 0 | 
                                        peptides[, n.reporter] > 0), ];
  cat("The number of peptide with either Nter or Lys modification is", nrow(peptides.LNreporter), "\n");
  cat("It is ", nrow(peptides.LNreporter)/total.peptides, "of total peptides\n\n")
  
}

TMT.label.efficiency.monkey <- function(reporter, n_row, nc) {
  evidence.le <- read.table.evidence(n_row);
  evidence.le <- filter.TMT.evidence(evidence.le);
    
  le.reporter <- paste("L", reporter, sep = "");
  le.reporter <- paste(le.reporter, nc, sep = "");
  # evidence <- evidence.le[which(evidence.le$Experiment == le.reporter),];
  evidence <- evidence.le[which(evidence.le$Experiment == "3ul_2"),];
  
#   cat ("The number of ");
  
  l.reporter <- paste("TMT.Lys.", nc, sep="");
  l.reporter <- paste(l.reporter, reporter, sep="");
  # l.reporter <- "TMT6plex.Lys126.Modification"
#   l.reporter <- paste(l.reporter, ".Modification", sep="");
  n.reporter <- paste("TMT.Nter.", nc, sep="");
  n.reporter <- paste(n.reporter, reporter, sep="");
  # n.reporter <- "TMT6plex.Nter126.Modification"
#   n.reporter <- paste(n.reporter, ".Modification", sep="");
  evidence <- evidence[, c("Sequence", "Length", "Acetyl..Protein.N.term.",
                           l.reporter, n.reporter,
                           "Missed.cleavages", "Proteins", "Leading.Proteins",
                           "Leading.Razor.Protein", "Type", "Raw.file", "Experiment",
                           "Mass", "Protein.group.IDs", "Peptide.ID",
                           "Mod..peptide.ID")];
  cat("The total number of evidence is", nrow(evidence), "\n");
  
  evidence.LNreporter <- evidence[which(evidence[, l.reporter] > 0 | 
                                        evidence[, n.reporter] > 0), ];
  cat("The number of evidence with either Nter or Lys modification is", nrow(evidence.LNreporter), "\n");
  cat("It is ", nrow(evidence.LNreporter)/nrow(evidence), "of total evidence\n\n")
  
  evidence.Lreporter <- evidence[which(evidence[, l.reporter] > 0), ];
  evidence.Nreporter <- evidence[which(evidence[, n.reporter] > 0), ];
  
#   evidence.sequence <- sort(unique(evidence$Sequence));  
#   
#   peptides <- NULL;
#   for (i in 1:length(evidence.sequence)) {
#     cat ("The calculation reaches ", i/length(evidence.sequence), "percent\n");
#     each.sequence <- toString(evidence.sequence[i]);
#     evidences.each.sequence <- evidence[which(evidence$Sequence == each.sequence),];
#     if (nrow(evidences.each.sequence) == 1) {
#       peptides <- rbind(peptides, evidences.each.sequence);
#     }
#     if (nrow(evidences.each.sequence) > 1) {
#       evidences.each.sequence[1, c("Acetyl..Protein.N.term.")] <- sum(evidences.each.sequence[, c("Acetyl..Protein.N.term.")]);
#       evidences.each.sequence[1, l.reporter] <- sum(evidences.each.sequence[, l.reporter]);
#       evidences.each.sequence[1, n.reporter] <- sum(evidences.each.sequence[, n.reporter]);
#       peptides <- rbind(peptides, evidences.each.sequence[1,]);
#     } 
#   }
#   
#   total.peptides <- nrow(peptides);
#   cat("The total number of peptide with", reporter, "is", total.peptides, "\n");
#   
#   peptides.Nacetyl <- peptides[which(peptides$Acetyl..Protein.N.term. > 0),];
#   cat("The number of peptide with N-acetyl is", nrow(peptides.Nacetyl), "\n");
#   cat("It is ", nrow(peptides.Nacetyl)/total.peptides, "of total peptides\n\n")
#   
#   peptides.Lreporter <- peptides[which(peptides[, l.reporter] > 0),];
#   cat("The number of peptide with Lys modification is", nrow(peptides.Lreporter), "\n");
#   cat("It is ", nrow(peptides.Lreporter)/total.peptides, "of total peptides\n\n")
#   
#   peptides.Nreporter <- peptides[which(peptides[, n.reporter] > 0),];
#   cat("The number of peptide with Nter modification is", nrow(peptides.Nreporter), "\n");
#   cat("It is ", nrow(peptides.Nreporter)/total.peptides, "of total peptides\n\n")
#   
#   peptides.LNreporter <- peptides[which(peptides[, l.reporter] > 0 | 
#                                           peptides[, n.reporter] > 0), ];
#   cat("The number of peptide with either Nter or Lys modification is", nrow(peptides.LNreporter), "\n");
#   cat("It is ", nrow(peptides.LNreporter)/total.peptides, "of total peptides\n\n")
  
}

TMT.peptides.cleanup <- function(nrow) {
  peptides <- read.table.peptides(nrow);
  #removing reverse and potential contaminants
  peptides <- filter.TMT.peptides(peptides);
  #calculate how many positive reporter ion intensities
  total.peptides <- nrow(peptides);
  reporter0.count <- nrow(peptides[which(peptides$Reporter.intensity.0 > 0),]);
  reporter1.count <- nrow(peptides[which(peptides$Reporter.intensity.1 > 0),]);
  reporter2.count <- nrow(peptides[which(peptides$Reporter.intensity.2 > 0),]);
  reporter3.count <- nrow(peptides[which(peptides$Reporter.intensity.3 > 0),]);
  reporter4.count <- nrow(peptides[which(peptides$Reporter.intensity.4 > 0),]);
  reporter5.count <- nrow(peptides[which(peptides$Reporter.intensity.5 > 0),]);
  
  cat("Total peptide count is", total.peptides, "\n");
  cat("Reporter 126 count is", reporter0.count, "\n");
  cat("Reporter 127 count is", reporter1.count, "\n");
  cat("Reporter 128 count is", reporter2.count, "\n");
  cat("Reporter 129 count is", reporter3.count, "\n");
  cat("Reporter 130 count is", reporter4.count, "\n");
  cat("Reporter 131 count is", reporter5.count, "\n");
  
  peptides <- peptides[which(peptides$Reporter.intensity.0 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.1 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.2 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.3 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.4 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.5 > 0),];
  
  cat("Total non-zero peptide count is", nrow(peptides), "\n");
  
  peptide.count <- c(total.peptides, reporter0.count, reporter1.count, 
                     reporter2.count, reporter3.count, reporter4.count,
                     reporter5.count, nrow(peptides));
  
  output <- list(peptides, peptide.count);
  output;
}

TMT.peptides.cleanup.10plex <- function(nrow) {
  peptides <- read.table.peptides(nrow);
  #removing reverse and potential contaminants
  peptides <- filter.TMT.peptides(peptides);
  #calculate how many positive reporter ion intensities
  total.peptides <- nrow(peptides);
  reporter0.count <- nrow(peptides[which(peptides$Reporter.intensity.0 > 0),]);
  reporter1.count <- nrow(peptides[which(peptides$Reporter.intensity.1 > 0),]);
  reporter2.count <- nrow(peptides[which(peptides$Reporter.intensity.2 > 0),]);
  reporter3.count <- nrow(peptides[which(peptides$Reporter.intensity.3 > 0),]);
  reporter4.count <- nrow(peptides[which(peptides$Reporter.intensity.4 > 0),]);
  reporter5.count <- nrow(peptides[which(peptides$Reporter.intensity.5 > 0),]);
  reporter6.count <- nrow(peptides[which(peptides$Reporter.intensity.6 > 0),]);
  reporter7.count <- nrow(peptides[which(peptides$Reporter.intensity.7 > 0),]);
  reporter8.count <- nrow(peptides[which(peptides$Reporter.intensity.8 > 0),]);
  reporter9.count <- nrow(peptides[which(peptides$Reporter.intensity.9 > 0),]);
  
  cat("Total peptide count is", total.peptides, "\n");
  cat("Reporter 126 count is", reporter0.count, "\n");
  cat("Reporter 127N count is", reporter1.count, "\n");
  cat("Reporter 127C count is", reporter2.count, "\n");
  cat("Reporter 128N count is", reporter3.count, "\n");
  cat("Reporter 128C count is", reporter4.count, "\n");
  cat("Reporter 129N count is", reporter5.count, "\n");
  cat("Reporter 129C count is", reporter2.count, "\n");
  cat("Reporter 130N count is", reporter3.count, "\n");
  cat("Reporter 130C count is", reporter4.count, "\n");
  cat("Reporter 131 count is", reporter5.count, "\n");
  
  peptides <- peptides[which(peptides$Reporter.intensity.0 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.1 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.2 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.3 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.4 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.5 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.6 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.7 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.8 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.9 > 0),];
  
  cat("Total non-zero peptide count is", nrow(peptides), "\n");
  
  peptide.count <- c(total.peptides, reporter0.count, reporter1.count, 
                     reporter2.count, reporter3.count, reporter4.count,
                     reporter5.count, reporter6.count, reporter7.count,
                     reporter8.count, reporter9.count, nrow(peptides));
  
  output <- list(peptides, peptide.count);
  output;
}

TMT.phosphosites.cleanup.10plex <- function(nrow) {
  peptides <- read.table.PhosphoSites(nrow);
  #removing reverse and potential contaminants
  peptides <- filter.TMT.PhosphoSites(peptides);
  #calculate how many positive reporter ion intensities
  total.peptides <- nrow(peptides);
  reporter0.count <- nrow(peptides[which(peptides$Reporter.intensity.0 > 0),]);
  reporter1.count <- nrow(peptides[which(peptides$Reporter.intensity.1 > 0),]);
  reporter2.count <- nrow(peptides[which(peptides$Reporter.intensity.2 > 0),]);
  reporter3.count <- nrow(peptides[which(peptides$Reporter.intensity.3 > 0),]);
  reporter4.count <- nrow(peptides[which(peptides$Reporter.intensity.4 > 0),]);
  reporter5.count <- nrow(peptides[which(peptides$Reporter.intensity.5 > 0),]);
  reporter6.count <- nrow(peptides[which(peptides$Reporter.intensity.6 > 0),]);
  reporter7.count <- nrow(peptides[which(peptides$Reporter.intensity.7 > 0),]);
  reporter8.count <- nrow(peptides[which(peptides$Reporter.intensity.8 > 0),]);
  reporter9.count <- nrow(peptides[which(peptides$Reporter.intensity.9 > 0),]);
  
  cat("Total phosphosites count is", total.peptides, "\n");
  cat("Reporter 126 count is", reporter0.count, "\n");
  cat("Reporter 127N count is", reporter1.count, "\n");
  cat("Reporter 127C count is", reporter2.count, "\n");
  cat("Reporter 128N count is", reporter3.count, "\n");
  cat("Reporter 128C count is", reporter4.count, "\n");
  cat("Reporter 129N count is", reporter5.count, "\n");
  cat("Reporter 129C count is", reporter2.count, "\n");
  cat("Reporter 130N count is", reporter3.count, "\n");
  cat("Reporter 130C count is", reporter4.count, "\n");
  cat("Reporter 131 count is", reporter5.count, "\n");
  
  peptides <- peptides[which(peptides$Reporter.intensity.0 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.1 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.2 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.3 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.4 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.5 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.6 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.7 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.8 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.9 > 0),];
  
  cat("Total non-zero phosphosites count is", nrow(peptides), "\n");
  
  peptide.count <- c(total.peptides, reporter0.count, reporter1.count, 
                     reporter2.count, reporter3.count, reporter4.count,
                     reporter5.count, reporter6.count, reporter7.count,
                     reporter8.count, reporter9.count, nrow(peptides));
  
  output <- list(peptides, peptide.count);
  output;
}

TMT.proteins.cleanup <- function (nrow) {
  proteins <- read.table.proteinGroup(nrow);
  proteins <- filter.proteinGroup.TMT(proteins);
  
  total.proteins <- nrow(proteins);
  reporter0.count <- nrow(proteins[which(proteins$Reporter.intensity.0 > 0),]);
  reporter1.count <- nrow(proteins[which(proteins$Reporter.intensity.1 > 0),]);
  reporter2.count <- nrow(proteins[which(proteins$Reporter.intensity.2 > 0),]);
  reporter3.count <- nrow(proteins[which(proteins$Reporter.intensity.3 > 0),]);
  reporter4.count <- nrow(proteins[which(proteins$Reporter.intensity.4 > 0),]);
  reporter5.count <- nrow(proteins[which(proteins$Reporter.intensity.5 > 0),]);
  
  cat("Total peptide count is", total.proteins, "\n");
  cat("Reporter 126 count is", reporter0.count, "\n");
  cat("Reporter 127 count is", reporter1.count, "\n");
  cat("Reporter 128 count is", reporter2.count, "\n");
  cat("Reporter 129 count is", reporter3.count, "\n");
  cat("Reporter 130 count is", reporter4.count, "\n");
  cat("Reporter 131 count is", reporter5.count, "\n");
  
  proteins <- proteins[which(proteins$Reporter.intensity.0 > 0),];
  proteins <- proteins[which(proteins$Reporter.intensity.1 > 0),];
  proteins <- proteins[which(proteins$Reporter.intensity.2 > 0),];
  proteins <- proteins[which(proteins$Reporter.intensity.3 > 0),];
  proteins <- proteins[which(proteins$Reporter.intensity.4 > 0),];
  proteins <- proteins[which(proteins$Reporter.intensity.5 > 0),];
  
  cat("Total non-zero protein count is", nrow(proteins), "\n");
  
  proteins
}

TMT.proteins.cleanup.10plex <- function (nrow) {
  proteins <- read.table.proteinGroup(nrow);
  proteins <- filter.proteinGroup.TMT(proteins);
  
  total.proteins <- nrow(proteins);
  reporter0.count <- nrow(proteins[which(proteins$Reporter.intensity.0 > 0),]);
  reporter1.count <- nrow(proteins[which(proteins$Reporter.intensity.1 > 0),]);
  reporter2.count <- nrow(proteins[which(proteins$Reporter.intensity.2 > 0),]);
  reporter3.count <- nrow(proteins[which(proteins$Reporter.intensity.3 > 0),]);
  reporter4.count <- nrow(proteins[which(proteins$Reporter.intensity.4 > 0),]);
  reporter5.count <- nrow(proteins[which(proteins$Reporter.intensity.5 > 0),]);
  reporter6.count <- nrow(proteins[which(proteins$Reporter.intensity.6 > 0),]);
  reporter7.count <- nrow(proteins[which(proteins$Reporter.intensity.7 > 0),]);
  reporter8.count <- nrow(proteins[which(proteins$Reporter.intensity.8 > 0),]);
  reporter9.count <- nrow(proteins[which(proteins$Reporter.intensity.9 > 0),]);
  
  cat("Total peptide count is", total.proteins, "\n");
  cat("Reporter 126 count is", reporter0.count, "\n");
  cat("Reporter 127n count is", reporter1.count, "\n");
  cat("Reporter 127c count is", reporter2.count, "\n");
  cat("Reporter 128n count is", reporter3.count, "\n");
  cat("Reporter 128c count is", reporter4.count, "\n");
  cat("Reporter 129n count is", reporter5.count, "\n");
  cat("Reporter 129c count is", reporter6.count, "\n");
  cat("Reporter 130n count is", reporter7.count, "\n");
  cat("Reporter 130c count is", reporter8.count, "\n");
  cat("Reporter 131 count is", reporter9.count, "\n");
  
  proteins <- proteins[which(proteins$Reporter.intensity.0 > 0),];
  proteins <- proteins[which(proteins$Reporter.intensity.1 > 0),];
  proteins <- proteins[which(proteins$Reporter.intensity.2 > 0),];
  proteins <- proteins[which(proteins$Reporter.intensity.3 > 0),];
  proteins <- proteins[which(proteins$Reporter.intensity.4 > 0),];
  proteins <- proteins[which(proteins$Reporter.intensity.5 > 0),];
  proteins <- proteins[which(proteins$Reporter.intensity.6 > 0),];
  proteins <- proteins[which(proteins$Reporter.intensity.7 > 0),];
  proteins <- proteins[which(proteins$Reporter.intensity.8 > 0),];
  proteins <- proteins[which(proteins$Reporter.intensity.9 > 0),];
  
  cat("Total non-zero protein count is", nrow(proteins), "\n");
  
  proteins
}

TMT.phosphosites.cleanup.10plex <- function (nrow) {
  phosphosites <- read.table.PhosphoSites(nrow);
  phosphosites <- filter.phosphosites.TMT(phosphosites);
  
  total.phosphosites <- nrow(phosphosites);
  reporter0.count <- nrow(phosphosites[which(phosphosites$Reporter.intensity.0 > 0),]);
  reporter1.count <- nrow(phosphosites[which(phosphosites$Reporter.intensity.1 > 0),]);
  reporter2.count <- nrow(phosphosites[which(phosphosites$Reporter.intensity.2 > 0),]);
  reporter3.count <- nrow(phosphosites[which(phosphosites$Reporter.intensity.3 > 0),]);
  reporter4.count <- nrow(phosphosites[which(phosphosites$Reporter.intensity.4 > 0),]);
  reporter5.count <- nrow(phosphosites[which(phosphosites$Reporter.intensity.5 > 0),]);
  reporter6.count <- nrow(phosphosites[which(phosphosites$Reporter.intensity.6 > 0),]);
  reporter7.count <- nrow(phosphosites[which(phosphosites$Reporter.intensity.7 > 0),]);
  reporter8.count <- nrow(phosphosites[which(phosphosites$Reporter.intensity.8 > 0),]);
  reporter9.count <- nrow(phosphosites[which(phosphosites$Reporter.intensity.9 > 0),]);
  
  cat("Total phosphosites count is", total.phosphosites, "\n");
  cat("Reporter 126 count is", reporter0.count, "\n");
  cat("Reporter 127n count is", reporter1.count, "\n");
  cat("Reporter 127c count is", reporter2.count, "\n");
  cat("Reporter 128n count is", reporter3.count, "\n");
  cat("Reporter 128c count is", reporter4.count, "\n");
  cat("Reporter 129n count is", reporter5.count, "\n");
  cat("Reporter 129c count is", reporter6.count, "\n");
  cat("Reporter 130n count is", reporter7.count, "\n");
  cat("Reporter 130c count is", reporter8.count, "\n");
  cat("Reporter 131 count is", reporter9.count, "\n");
  
  phosphosites <- phosphosites[which(phosphosites$Reporter.intensity.0 > 0),];
  phosphosites <- phosphosites[which(phosphosites$Reporter.intensity.1 > 0),];
  phosphosites <- phosphosites[which(phosphosites$Reporter.intensity.2 > 0),];
  phosphosites <- phosphosites[which(phosphosites$Reporter.intensity.3 > 0),];
  phosphosites <- phosphosites[which(phosphosites$Reporter.intensity.4 > 0),];
  phosphosites <- phosphosites[which(phosphosites$Reporter.intensity.5 > 0),];
  phosphosites <- phosphosites[which(phosphosites$Reporter.intensity.6 > 0),];
  phosphosites <- phosphosites[which(phosphosites$Reporter.intensity.7 > 0),];
  phosphosites <- phosphosites[which(phosphosites$Reporter.intensity.8 > 0),];
  phosphosites <- phosphosites[which(phosphosites$Reporter.intensity.9 > 0),];
  
  cat("Total non-zero protein count is", nrow(phosphosites), "\n");
  
  phosphosites
}

get.TMT.ratios <- function(peptides.cleanup) {
  ratios.126.128 <- peptides.cleanup[, c("Reporter.intensity.0")]/peptides.cleanup[, c("Reporter.intensity.2")];
  ratios.127.128 <- peptides.cleanup[, c("Reporter.intensity.1")]/peptides.cleanup[, c("Reporter.intensity.2")];
  ratios.128.128 <- peptides.cleanup[, c("Reporter.intensity.2")]/peptides.cleanup[, c("Reporter.intensity.2")];
  ratios.129.128 <- peptides.cleanup[, c("Reporter.intensity.3")]/peptides.cleanup[, c("Reporter.intensity.2")];
  ratios.130.128 <- peptides.cleanup[, c("Reporter.intensity.4")]/peptides.cleanup[, c("Reporter.intensity.2")];
  ratios.131.128 <- peptides.cleanup[, c("Reporter.intensity.5")]/peptides.cleanup[, c("Reporter.intensity.2")];
  
  ratios.129.129 <- peptides.cleanup[, c("Reporter.intensity.3")]/peptides.cleanup[, c("Reporter.intensity.3")];
  ratios.130.129 <- peptides.cleanup[, c("Reporter.intensity.4")]/peptides.cleanup[, c("Reporter.intensity.3")];
  ratios.131.129 <- peptides.cleanup[, c("Reporter.intensity.5")]/peptides.cleanup[, c("Reporter.intensity.3")];
  
  log2.ratios.126.128 <- log2(ratios.126.128);
  log2.ratios.127.128 <- log2(ratios.127.128);
  log2.ratios.128.128 <- log2(ratios.128.128);
  log2.ratios.129.128 <- log2(ratios.129.128);
  log2.ratios.130.128 <- log2(ratios.130.128);
  log2.ratios.131.128 <- log2(ratios.131.128);
  
  log2.ratios.129.129 <- log2(ratios.129.129);
  log2.ratios.130.129 <- log2(ratios.130.129);
  log2.ratios.131.129 <- log2(ratios.131.129);
  
  TMT.ratios <- data.frame(ratios.126.128, ratios.127.128, ratios.128.128, 
                           ratios.129.128, ratios.130.128, ratios.131.128, 
                           ratios.129.129, ratios.130.129, ratios.131.129, 
                           log2.ratios.126.128, log2.ratios.127.128, log2.ratios.128.128, 
                           log2.ratios.129.128, log2.ratios.130.128, log2.ratios.131.128, 
                           log2.ratios.129.129, log2.ratios.130.129, log2.ratios.131.129 
                           );
  
  peptides.cleanup <- cbind(peptides.cleanup, TMT.ratios);
  
  peptides.cleanup;
}

get.TMT.ratios2 <- function(peptides.cleanup) {
  ratios.126.126 <- peptides.cleanup[, c("Reporter.intensity.0")]/peptides.cleanup[, c("Reporter.intensity.0")];
  ratios.127.126 <- peptides.cleanup[, c("Reporter.intensity.1")]/peptides.cleanup[, c("Reporter.intensity.0")];
  ratios.128.126 <- peptides.cleanup[, c("Reporter.intensity.2")]/peptides.cleanup[, c("Reporter.intensity.0")];
  ratios.129.126 <- peptides.cleanup[, c("Reporter.intensity.3")]/peptides.cleanup[, c("Reporter.intensity.0")];
  ratios.130.126 <- peptides.cleanup[, c("Reporter.intensity.4")]/peptides.cleanup[, c("Reporter.intensity.0")];
  ratios.131.126 <- peptides.cleanup[, c("Reporter.intensity.5")]/peptides.cleanup[, c("Reporter.intensity.0")];
  
  log2.ratios.126.126 <- log2(ratios.126.126);
  log2.ratios.127.126 <- log2(ratios.127.126);
  log2.ratios.128.126 <- log2(ratios.128.126);
  log2.ratios.129.126 <- log2(ratios.129.126);
  log2.ratios.130.126 <- log2(ratios.130.126);
  log2.ratios.131.126 <- log2(ratios.131.126);
  
  TMT.ratios <- data.frame(ratios.126.126, ratios.127.126, ratios.128.126, 
                           ratios.129.126, ratios.130.126, ratios.131.126, 
                           log2.ratios.126.126, log2.ratios.127.126, log2.ratios.128.126, 
                           log2.ratios.129.126, log2.ratios.130.126, log2.ratios.131.126
                           );
  
  peptides.cleanup <- cbind(peptides.cleanup, TMT.ratios);
  
  peptides.cleanup;
}

#with total intensity normalization
get.TMT.ratios3 <- function(peptides.cleanup) {
  sum.126 <- sum(peptides.cleanup[, c("Reporter.intensity.0")]);
  sum.127 <- sum(peptides.cleanup[, c("Reporter.intensity.1")]);
  sum.128 <- sum(peptides.cleanup[, c("Reporter.intensity.2")]);
  sum.129 <- sum(peptides.cleanup[, c("Reporter.intensity.3")]);
  sum.130 <- sum(peptides.cleanup[, c("Reporter.intensity.4")]);
  sum.131 <- sum(peptides.cleanup[, c("Reporter.intensity.5")]);
  
  ratios.126.126 <- peptides.cleanup[, c("Reporter.intensity.0")]/peptides.cleanup[, c("Reporter.intensity.0")];
  ratios.127.126 <- peptides.cleanup[, c("Reporter.intensity.1")]/peptides.cleanup[, c("Reporter.intensity.0")];
  ratios.128.126 <- peptides.cleanup[, c("Reporter.intensity.2")]/peptides.cleanup[, c("Reporter.intensity.0")];
  ratios.129.126 <- peptides.cleanup[, c("Reporter.intensity.3")]/peptides.cleanup[, c("Reporter.intensity.0")];
  ratios.130.126 <- peptides.cleanup[, c("Reporter.intensity.4")]/peptides.cleanup[, c("Reporter.intensity.0")];
  ratios.131.126 <- peptides.cleanup[, c("Reporter.intensity.5")]/peptides.cleanup[, c("Reporter.intensity.0")];
  
  normalized.ratios.126.126 <- ratios.126.126 * sum.126 / sum.126;
  normalized.ratios.127.126 <- ratios.127.126 * sum.126 / sum.127;
  normalized.ratios.128.126 <- ratios.128.126 * sum.126 / sum.128;
  normalized.ratios.129.126 <- ratios.129.126 * sum.126 / sum.129;
  normalized.ratios.130.126 <- ratios.130.126 * sum.126 / sum.130;
  normalized.ratios.131.126 <- ratios.131.126 * sum.126 / sum.131;
  
  log2.ratios.126.126 <- log2(ratios.126.126);
  log2.ratios.127.126 <- log2(ratios.127.126);
  log2.ratios.128.126 <- log2(ratios.128.126);
  log2.ratios.129.126 <- log2(ratios.129.126);
  log2.ratios.130.126 <- log2(ratios.130.126);
  log2.ratios.131.126 <- log2(ratios.131.126);
  
  log2.normalized.ratios.126.126 <- log2(normalized.ratios.126.126);
  log2.normalized.ratios.127.126 <- log2(normalized.ratios.127.126);
  log2.normalized.ratios.128.126 <- log2(normalized.ratios.128.126);
  log2.normalized.ratios.129.126 <- log2(normalized.ratios.129.126);
  log2.normalized.ratios.130.126 <- log2(normalized.ratios.130.126);
  log2.normalized.ratios.131.126 <- log2(normalized.ratios.131.126);
  
  TMT.ratios <- data.frame(ratios.126.126, ratios.127.126, ratios.128.126, 
                           ratios.129.126, ratios.130.126, ratios.131.126, 
                           normalized.ratios.126.126, normalized.ratios.127.126,
                           normalized.ratios.128.126, normalized.ratios.129.126,
                           normalized.ratios.130.126, normalized.ratios.131.126,
                           log2.ratios.126.126, log2.ratios.127.126, log2.ratios.128.126, 
                           log2.ratios.129.126, log2.ratios.130.126, log2.ratios.131.126,
                           log2.normalized.ratios.126.126, log2.normalized.ratios.127.126,
                           log2.normalized.ratios.128.126, log2.normalized.ratios.129.126,
                           log2.normalized.ratios.130.126, log2.normalized.ratios.131.126
                           );
  
  peptides.cleanup <- cbind(peptides.cleanup, TMT.ratios);
  
  peptides.cleanup;
}

#total intensity normalization 
normalize.TMT.Intensity.6plex <- function(proteins) {
  sum.126 <- sum(proteins[, c("Reporter.intensity.0")]);
  sum.127 <- sum(proteins[, c("Reporter.intensity.1")]);
  sum.128 <- sum(proteins[, c("Reporter.intensity.2")]);
  sum.129 <- sum(proteins[, c("Reporter.intensity.3")]);
  sum.130 <- sum(proteins[, c("Reporter.intensity.4")]);
  sum.131 <- sum(proteins[, c("Reporter.intensity.5")]);
  
  normalized.126 <- proteins[, c("Reporter.intensity.0")] * sum.126 / sum.126;
  normalized.127 <- proteins[, c("Reporter.intensity.1")] * sum.126 / sum.127;
  normalized.128 <- proteins[, c("Reporter.intensity.2")] * sum.126 / sum.128;
  normalized.129 <- proteins[, c("Reporter.intensity.3")] * sum.126 / sum.129;
  normalized.130 <- proteins[, c("Reporter.intensity.4")] * sum.126 / sum.130;
  normalized.131 <- proteins[, c("Reporter.intensity.5")] * sum.126 / sum.131;
  
  log2.intensity.126 <- log2(normalized.126);
  log2.intensity.127 <- log2(normalized.127);
  log2.intensity.128 <- log2(normalized.128);
  log2.intensity.129 <- log2(normalized.129);
  log2.intensity.130 <- log2(normalized.130);
  log2.intensity.131 <- log2(normalized.131);
  log2.intensity <- log2(proteins$Intensity);
  
  normalized.intensity <- data.frame(normalized.126, 
                                     normalized.127,
                                     normalized.128,
                                     normalized.129,
                                     normalized.130,
                                     normalized.131,
                                     log2.intensity.126,
                                     log2.intensity.127,
                                     log2.intensity.128,
                                     log2.intensity.129,
                                     log2.intensity.130,
                                     log2.intensity.131)
  
  proteins <- cbind(proteins, normalized.intensity);
  
  proteins;
  
}


normalize.TMT.Intensity.10plex <- function(proteins) {
  sum.126 <- sum(proteins[, c("Reporter.intensity.0")]);
  sum.127n <- sum(proteins[, c("Reporter.intensity.1")]);
  sum.127c <- sum(proteins[, c("Reporter.intensity.2")]);
  sum.128n <- sum(proteins[, c("Reporter.intensity.3")]);
  sum.128c <- sum(proteins[, c("Reporter.intensity.4")]);
  sum.129n <- sum(proteins[, c("Reporter.intensity.5")]);
  sum.129c <- sum(proteins[, c("Reporter.intensity.6")]);
  sum.130n <- sum(proteins[, c("Reporter.intensity.7")]);
  sum.130c <- sum(proteins[, c("Reporter.intensity.8")]);
  sum.131 <- sum(proteins[, c("Reporter.intensity.9")]);
  
  normalized.126 <- proteins[, c("Reporter.intensity.0")] * sum.126 / sum.126;
  normalized.127n <- proteins[, c("Reporter.intensity.1")] * sum.126 / sum.127n;
  normalized.127c <- proteins[, c("Reporter.intensity.2")] * sum.126 / sum.127c;
  normalized.128n <- proteins[, c("Reporter.intensity.3")] * sum.126 / sum.128n;
  normalized.128c <- proteins[, c("Reporter.intensity.4")] * sum.126 / sum.128c;
  normalized.129n <- proteins[, c("Reporter.intensity.5")] * sum.126 / sum.129n;
  normalized.129c <- proteins[, c("Reporter.intensity.6")] * sum.126 / sum.129c;
  normalized.130n <- proteins[, c("Reporter.intensity.7")] * sum.126 / sum.130n;
  normalized.130c <- proteins[, c("Reporter.intensity.8")] * sum.126 / sum.130c;
  normalized.131 <- proteins[, c("Reporter.intensity.9")] * sum.126 / sum.131;
  
  log10.intensity.126 <- log10(normalized.126);
  log10.intensity.127n <- log10(normalized.127n);
  log10.intensity.127c <- log10(normalized.127c);
  log10.intensity.128n <- log10(normalized.128n);
  log10.intensity.128c <- log10(normalized.128c);
  log10.intensity.129n <- log10(normalized.129n);
  log10.intensity.129c <- log10(normalized.129c);
  log10.intensity.130n <- log10(normalized.130n);
  log10.intensity.130c <- log10(normalized.130c);
  log10.intensity.131 <- log10(normalized.131);
  log10.intensity <- log10(proteins$Intensity);
  
  normalized.intensity <- data.frame(normalized.126, 
                                     normalized.127n, normalized.127c,
                                     normalized.128n, normalized.128c,
                                     normalized.129n, normalized.129c,
                                     normalized.130n, normalized.130c,
                                     normalized.131,
                                     log10.intensity.126,
                                     log10.intensity.127n, log10.intensity.127c,
                                     log10.intensity.128n, log10.intensity.128c,
                                     log10.intensity.129n, log10.intensity.129c,
                                     log10.intensity.130n, log10.intensity.130c,
                                     log10.intensity.131)
  
  proteins <- cbind(proteins, normalized.intensity);
  
  proteins;
                           
}

get.TMT.ratios.10plex <- function(peptides.cleanup) {
  sum.126 <- sum(peptides.cleanup[, c("Reporter.intensity.0")]);
  sum.127n <- sum(peptides.cleanup[, c("Reporter.intensity.1")]);
  sum.127c <- sum(peptides.cleanup[, c("Reporter.intensity.2")]);
  sum.128n <- sum(peptides.cleanup[, c("Reporter.intensity.3")]);
  sum.128c <- sum(peptides.cleanup[, c("Reporter.intensity.4")]);
  sum.129n <- sum(peptides.cleanup[, c("Reporter.intensity.5")]);
  sum.129c <- sum(peptides.cleanup[, c("Reporter.intensity.6")]);
  sum.130n <- sum(peptides.cleanup[, c("Reporter.intensity.7")]);
  sum.130c <- sum(peptides.cleanup[, c("Reporter.intensity.8")]);
  sum.131 <- sum(peptides.cleanup[, c("Reporter.intensity.9")]);
  reporters.intensity.0.9 <- peptides.cleanup[, c("Reporter.intensity.0", "Reporter.intensity.1",
                                                  "Reporter.intensity.2", "Reporter.intensity.3",
                                                  "Reporter.intensity.4", "Reporter.intensity.5",
                                                  "Reporter.intensity.6", "Reporter.intensity.7",
                                                  "Reporter.intensity.8", "Reporter.intensity.9")];
  average.intensity <- apply(reporters.intensity.0.9, 1, mean);
  
  ratios.126.126 <- peptides.cleanup[, c("Reporter.intensity.0")]/peptides.cleanup[, c("Reporter.intensity.0")];
  ratios.127n.126 <- peptides.cleanup[, c("Reporter.intensity.1")]/peptides.cleanup[, c("Reporter.intensity.0")];
  ratios.127c.126 <- peptides.cleanup[, c("Reporter.intensity.2")]/peptides.cleanup[, c("Reporter.intensity.0")];
  ratios.128n.126 <- peptides.cleanup[, c("Reporter.intensity.3")]/peptides.cleanup[, c("Reporter.intensity.0")];
  ratios.128c.126 <- peptides.cleanup[, c("Reporter.intensity.4")]/peptides.cleanup[, c("Reporter.intensity.0")];
  ratios.129n.126 <- peptides.cleanup[, c("Reporter.intensity.5")]/peptides.cleanup[, c("Reporter.intensity.0")];
  ratios.129c.126 <- peptides.cleanup[, c("Reporter.intensity.6")]/peptides.cleanup[, c("Reporter.intensity.0")];
  ratios.130n.126 <- peptides.cleanup[, c("Reporter.intensity.7")]/peptides.cleanup[, c("Reporter.intensity.0")];
  ratios.130c.126 <- peptides.cleanup[, c("Reporter.intensity.8")]/peptides.cleanup[, c("Reporter.intensity.0")];
  ratios.131.126 <- peptides.cleanup[, c("Reporter.intensity.9")]/peptides.cleanup[, c("Reporter.intensity.0")];
  
  cat("The median ratio of 127n/126 is", median(ratios.127n.126), "\n");
  cat("The median ratio of 127c/126 is", median(ratios.127c.126), "\n");
  cat("The median ratio of 128n/126 is", median(ratios.128n.126), "\n");
  cat("The median ratio of 128c/126 is", median(ratios.128c.126), "\n");
  cat("The median ratio of 129n/126 is", median(ratios.129n.126), "\n");
  cat("The median ratio of 129c/126 is", median(ratios.129c.126), "\n");
  cat("The median ratio of 130n/126 is", median(ratios.130n.126), "\n");
  cat("The median ratio of 130c/126 is", median(ratios.130c.126), "\n");
  cat("The median ratio of 131/126 is", median(ratios.131.126), "\n");
  
  normalized.ratios.126.126 <- ratios.126.126 * sum.126 / sum.126;
  normalized.ratios.127n.126 <- ratios.127n.126 * sum.126 / sum.127n;
  normalized.ratios.127c.126 <- ratios.127c.126 * sum.126 / sum.127c;
  normalized.ratios.128n.126 <- ratios.128n.126 * sum.126 / sum.128n;
  normalized.ratios.128c.126 <- ratios.128c.126 * sum.126 / sum.128c;
  normalized.ratios.129n.126 <- ratios.129n.126 * sum.126 / sum.129n;
  normalized.ratios.129c.126 <- ratios.129c.126 * sum.126 / sum.129c;
  normalized.ratios.130n.126 <- ratios.130n.126 * sum.126 / sum.130n;
  normalized.ratios.130c.126 <- ratios.130c.126 * sum.126 / sum.130c;
  normalized.ratios.131.126 <- ratios.131.126 * sum.126 / sum.131;
  
  log2.ratios.126.126 <- log2(ratios.126.126);
  log2.ratios.127n.126 <- log2(ratios.127n.126);
  log2.ratios.127c.126 <- log2(ratios.127c.126);
  log2.ratios.128n.126 <- log2(ratios.128n.126);
  log2.ratios.128c.126 <- log2(ratios.128c.126);
  log2.ratios.129n.126 <- log2(ratios.129n.126);
  log2.ratios.129c.126 <- log2(ratios.129c.126);
  log2.ratios.130n.126 <- log2(ratios.130n.126);
  log2.ratios.130c.126 <- log2(ratios.130c.126);
  log2.ratios.131.126 <- log2(ratios.131.126);
  
  log2.normalized.ratios.126.126 <- log2(normalized.ratios.126.126);
  log2.normalized.ratios.127n.126 <- log2(normalized.ratios.127n.126);
  log2.normalized.ratios.127c.126 <- log2(normalized.ratios.127c.126);
  log2.normalized.ratios.128n.126 <- log2(normalized.ratios.128n.126);
  log2.normalized.ratios.128c.126 <- log2(normalized.ratios.128c.126);
  log2.normalized.ratios.129n.126 <- log2(normalized.ratios.129n.126);
  log2.normalized.ratios.129c.126 <- log2(normalized.ratios.129c.126);
  log2.normalized.ratios.130n.126 <- log2(normalized.ratios.130n.126);
  log2.normalized.ratios.130c.126 <- log2(normalized.ratios.130c.126);
  log2.normalized.ratios.131.126 <- log2(normalized.ratios.131.126);
  
  TMT.ratios <- data.frame(ratios.126.126, ratios.127n.126, ratios.127c.126, 
                           ratios.128n.126, ratios.128c.126, ratios.129n.126,
                           ratios.129c.126, ratios.130n.126, ratios.130c.126,
                           ratios.131.126, 
                           normalized.ratios.126.126, normalized.ratios.127c.126,
                           normalized.ratios.127n.126, normalized.ratios.128n.126,
                           normalized.ratios.128c.126, normalized.ratios.129n.126,
                           normalized.ratios.129c.126, normalized.ratios.130n.126,
                           normalized.ratios.130c.126, normalized.ratios.131.126,
                           log2.ratios.126.126, log2.ratios.127n.126, log2.ratios.127c.126, 
                           log2.ratios.128n.126, log2.ratios.128c.126, log2.ratios.129n.126,
                           log2.ratios.129c.126, log2.ratios.130n.126, log2.ratios.130c.126, 
                           log2.ratios.131.126,
                           log2.normalized.ratios.126.126, log2.normalized.ratios.127c.126,
                           log2.normalized.ratios.127n.126, log2.normalized.ratios.128n.126,
                           log2.normalized.ratios.128c.126, log2.normalized.ratios.129n.126, 
                           log2.normalized.ratios.129c.126, log2.normalized.ratios.130n.126,
                           log2.normalized.ratios.130c.126, log2.normalized.ratios.131.126
  );
  
  peptides.cleanup <- cbind(peptides.cleanup, TMT.ratios);
  
  peptides.cleanup;
}

#ratios were calculated using the average intensity for each line
get.TMT.ratios.10plex2 <- function(peptides.cleanup) {
  
  sum.126 <- sum(peptides.cleanup[, c("Reporter.intensity.0")]);
  sum.127n <- sum(peptides.cleanup[, c("Reporter.intensity.1")]);
  sum.127c <- sum(peptides.cleanup[, c("Reporter.intensity.2")]);
  sum.128n <- sum(peptides.cleanup[, c("Reporter.intensity.3")]);
  sum.128c <- sum(peptides.cleanup[, c("Reporter.intensity.4")]);
  sum.129n <- sum(peptides.cleanup[, c("Reporter.intensity.5")]);
  sum.129c <- sum(peptides.cleanup[, c("Reporter.intensity.6")]);
  sum.130n <- sum(peptides.cleanup[, c("Reporter.intensity.7")]);
  sum.130c <- sum(peptides.cleanup[, c("Reporter.intensity.8")]);
  sum.131 <- sum(peptides.cleanup[, c("Reporter.intensity.9")]);
  reporters.intensity.0.9 <- peptides.cleanup[, c("Reporter.intensity.0", "Reporter.intensity.1",
                                                  "Reporter.intensity.2", "Reporter.intensity.3",
                                                  "Reporter.intensity.4", "Reporter.intensity.5",
                                                  "Reporter.intensity.6", "Reporter.intensity.7",
                                                  "Reporter.intensity.8", "Reporter.intensity.9")];
  average.intensity <- apply(reporters.intensity.0.9, 1, mean);
  sum.average <- sum(average.intensity);
  
  ratios.126 <- peptides.cleanup[, c("Reporter.intensity.0")]/average.intensity;
  ratios.127n <- peptides.cleanup[, c("Reporter.intensity.1")]/average.intensity;
  ratios.127c <- peptides.cleanup[, c("Reporter.intensity.2")]/average.intensity;
  ratios.128n <- peptides.cleanup[, c("Reporter.intensity.3")]/average.intensity;
  ratios.128c <- peptides.cleanup[, c("Reporter.intensity.4")]/average.intensity;
  ratios.129n <- peptides.cleanup[, c("Reporter.intensity.5")]/average.intensity;
  ratios.129c <- peptides.cleanup[, c("Reporter.intensity.6")]/average.intensity;
  ratios.130n <- peptides.cleanup[, c("Reporter.intensity.7")]/average.intensity;
  ratios.130c <- peptides.cleanup[, c("Reporter.intensity.8")]/average.intensity;
  ratios.131 <- peptides.cleanup[, c("Reporter.intensity.9")]/average.intensity;
  
  cat("The median ratio of 126/mean is", median(ratios.126), "\n");
  cat("The median ratio of 127n/mean is", median(ratios.127n), "\n");
  cat("The median ratio of 127c/mean is", median(ratios.127c), "\n");
  cat("The median ratio of 128n/mean is", median(ratios.128n), "\n");
  cat("The median ratio of 128c/mean is", median(ratios.128c), "\n");
  cat("The median ratio of 129n/mean is", median(ratios.129n), "\n");
  cat("The median ratio of 129c/mean is", median(ratios.129c), "\n");
  cat("The median ratio of 130n/mean is", median(ratios.130n), "\n");
  cat("The median ratio of 130c/mean is", median(ratios.130c), "\n");
  cat("The median ratio of 131/mean is", median(ratios.131), "\n");
  
  normalized.ratios.126 <- ratios.126 * sum.average / sum.126;
  normalized.ratios.127n <- ratios.127n * sum.average / sum.127n;
  normalized.ratios.127c <- ratios.127c * sum.average / sum.127c;
  normalized.ratios.128n <- ratios.128n * sum.average / sum.128n;
  normalized.ratios.128c <- ratios.128c * sum.average / sum.128c;
  normalized.ratios.129n <- ratios.129n * sum.average / sum.129n;
  normalized.ratios.129c <- ratios.129c * sum.average / sum.129c;
  normalized.ratios.130n <- ratios.130n * sum.average / sum.130n;
  normalized.ratios.130c <- ratios.130c * sum.average / sum.130c;
  normalized.ratios.131 <- ratios.131 * sum.average / sum.131;
  
  normalized.Reporter.intensity.0 <- peptides.cleanup[, c("Reporter.intensity.0")] * sum.average / sum.126;
  normalized.Reporter.intensity.1 <- peptides.cleanup[, c("Reporter.intensity.1")] * sum.average / sum.127n;
  normalized.Reporter.intensity.2 <- peptides.cleanup[, c("Reporter.intensity.2")] * sum.average / sum.127c;
  normalized.Reporter.intensity.3 <- peptides.cleanup[, c("Reporter.intensity.3")] * sum.average / sum.128n;
  normalized.Reporter.intensity.4 <- peptides.cleanup[, c("Reporter.intensity.4")] * sum.average / sum.128c;
  normalized.Reporter.intensity.5 <- peptides.cleanup[, c("Reporter.intensity.5")] * sum.average / sum.129n;
  normalized.Reporter.intensity.6 <- peptides.cleanup[, c("Reporter.intensity.6")] * sum.average / sum.129c;
  normalized.Reporter.intensity.7 <- peptides.cleanup[, c("Reporter.intensity.7")] * sum.average / sum.130n;
  normalized.Reporter.intensity.8 <- peptides.cleanup[, c("Reporter.intensity.8")] * sum.average / sum.130c;
  normalized.Reporter.intensity.9 <- peptides.cleanup[, c("Reporter.intensity.9")] * sum.average / sum.131;
  
  cat("The normalized median ratio of 126/mean is", median(normalized.ratios.126), "\n");
  cat("The normalized median ratio of 127n/mean is", median(normalized.ratios.127n), "\n");
  cat("The normalized median ratio of 127c/mean is", median(normalized.ratios.127c), "\n");
  cat("The normalized median ratio of 128n/mean is", median(normalized.ratios.128n), "\n");
  cat("The normalized median ratio of 128c/mean is", median(normalized.ratios.128c), "\n");
  cat("The normalized median ratio of 129n/mean is", median(normalized.ratios.129n), "\n");
  cat("The normalized median ratio of 129c/mean is", median(normalized.ratios.129c), "\n");
  cat("The normalized median ratio of 130n/mean is", median(normalized.ratios.130n), "\n");
  cat("The normalized median ratio of 130c/mean is", median(normalized.ratios.130c), "\n");
  cat("The normalized median ratio of 131/mean is", median(normalized.ratios.131), "\n");
  
  log2.ratios.126 <- log2(ratios.126);
  log2.ratios.127n <- log2(ratios.127n);
  log2.ratios.127c <- log2(ratios.127c);
  log2.ratios.128n <- log2(ratios.128n);
  log2.ratios.128c <- log2(ratios.128c);
  log2.ratios.129n <- log2(ratios.129n);
  log2.ratios.129c <- log2(ratios.129c);
  log2.ratios.130n <- log2(ratios.130n);
  log2.ratios.130c <- log2(ratios.130c);
  log2.ratios.131 <- log2(ratios.131);
  
  log2.normalized.ratios.126 <- log2(normalized.ratios.126);
  log2.normalized.ratios.127n <- log2(normalized.ratios.127n);
  log2.normalized.ratios.127c <- log2(normalized.ratios.127c);
  log2.normalized.ratios.128n <- log2(normalized.ratios.128n);
  log2.normalized.ratios.128c <- log2(normalized.ratios.128c);
  log2.normalized.ratios.129n <- log2(normalized.ratios.129n);
  log2.normalized.ratios.129c <- log2(normalized.ratios.129c);
  log2.normalized.ratios.130n <- log2(normalized.ratios.130n);
  log2.normalized.ratios.130c <- log2(normalized.ratios.130c);
  log2.normalized.ratios.131 <- log2(normalized.ratios.131);
  
  TMT.ratios <- data.frame(ratios.126, ratios.127n, ratios.127c, 
                           ratios.128n, ratios.128c, ratios.129n,
                           ratios.129c, ratios.130n, ratios.130c,
                           ratios.131, 
                           normalized.ratios.126, normalized.ratios.127n,
                           normalized.ratios.127c, normalized.ratios.128n,
                           normalized.ratios.128c, normalized.ratios.129n,
                           normalized.ratios.129c, normalized.ratios.130n,
                           normalized.ratios.130c, normalized.ratios.131,
                           log2.ratios.126, log2.ratios.127n, log2.ratios.127c, 
                           log2.ratios.128n, log2.ratios.128c, log2.ratios.129n,
                           log2.ratios.129c, log2.ratios.130n, log2.ratios.130c, 
                           log2.ratios.131,
                           log2.normalized.ratios.126, log2.normalized.ratios.127c,
                           log2.normalized.ratios.127n, log2.normalized.ratios.128n,
                           log2.normalized.ratios.128c, log2.normalized.ratios.129n, 
                           log2.normalized.ratios.129c, log2.normalized.ratios.130n,
                           log2.normalized.ratios.130c, log2.normalized.ratios.131,
                           normalized.Reporter.intensity.0, normalized.Reporter.intensity.1,
                           normalized.Reporter.intensity.2, normalized.Reporter.intensity.3,
                           normalized.Reporter.intensity.4, normalized.Reporter.intensity.5,
                           normalized.Reporter.intensity.6, normalized.Reporter.intensity.7,
                           normalized.Reporter.intensity.8, normalized.Reporter.intensity.9);
  
  peptides.cleanup <- cbind(peptides.cleanup, TMT.ratios);
  
  peptides.cleanup;
}


check.log2.ratios.distribution <- function(peptide.ratios, peporpro) {
  file.name <- paste(peporpro, "log2_ratios_distribution.png", sep = "_");
  
  png(filename = file.name, width = 1000, height = 1000);
  par(mfrow=c(3,3), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  hist(peptide.ratios$log2.ratios.126.128, main = NULL, 
       xlab = "log2 126/128", col = "dodgerblue4");
  hist(peptide.ratios$log2.ratios.127.128, main = NULL, 
       xlab = "log2 127/128", col = "dodgerblue4");
  hist(peptide.ratios$log2.ratios.128.128, main = NULL, 
       xlab = "log2 128/128", col = "dodgerblue4");
  hist(peptide.ratios$log2.ratios.129.128, main = NULL, 
       xlab = "log2 129/128", col = "dodgerblue4");
  hist(peptide.ratios$log2.ratios.130.128, main = NULL, 
       xlab = "log2 130/128", col = "dodgerblue4");
  hist(peptide.ratios$log2.ratios.131.128, main = NULL, 
       xlab = "log2 131/128", col = "dodgerblue4");
  hist(peptide.ratios$log2.ratios.129.129, main = NULL, 
       xlab = "log2 129/129", col = "dodgerblue4");
  hist(peptide.ratios$log2.ratios.130.129, main = NULL, 
       xlab = "log2 130/129", col = "dodgerblue4");
  hist(peptide.ratios$log2.ratios.131.129, main = NULL, 
       xlab = "log2 131/129", col = "dodgerblue4");
  
  dev.off()
}

check.log2.ratios.distribution2 <- function(peptide.ratios, peporpro) {
  file.name <- paste(peporpro, "log2_ratios_distribution2.png", sep = "_");
  
  png(filename = file.name, width = 900, height = 600);
  par(mfrow=c(2,3), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  hist(peptide.ratios$log2.ratios.126.126, main = NULL, 
       xlab = "log2 126/126", col = "dodgerblue4");
  hist(peptide.ratios$log2.ratios.127.126, main = NULL, 
       xlab = "log2 127/126", col = "dodgerblue4");
  hist(peptide.ratios$log2.ratios.128.126, main = NULL, 
       xlab = "log2 128/126", col = "dodgerblue4");
  hist(peptide.ratios$log2.ratios.129.126, main = NULL, 
       xlab = "log2 129/126", col = "dodgerblue4");
  hist(peptide.ratios$log2.ratios.130.126, main = NULL, 
       xlab = "log2 130/126", col = "dodgerblue4");
  hist(peptide.ratios$log2.ratios.131.126, main = NULL, 
       xlab = "log2 131/126", col = "dodgerblue4");
    
  dev.off()
}

check.log2.ratios.distribution3 <- function(peptide.ratios, peporpro) {
  file.name <- paste(peporpro, "log2_ratios_distribution3.png", sep = "_");
  
  png(filename = file.name, width = 900, height = 600);
  par(mfrow=c(2,3), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  hist(peptide.ratios$log2.normalized.ratios.126.126, main = NULL, 
       xlab = "Normalized log2 126/126", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.127.126, main = NULL, 
       xlab = "Normalized log2 127/126", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.128.126, main = NULL, 
       xlab = "Normalized log2 128/126", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.129.126, main = NULL, 
       xlab = "Normalized log2 129/126", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.130.126, main = NULL, 
       xlab = "Normalized log2 130/126", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.131.126, main = NULL, 
       xlab = "Normalized log2 131/126", col = "dodgerblue4");
  
  dev.off()
}

check.log2.ratios.distribution.10plex <- function(peptide.ratios, peporpro) {
  file.name <- paste(peporpro, "log2_ratios_distribution3.png", sep = "_");
  
  png(filename = file.name, width = 900, height = 900);
  par(mfrow=c(3,3), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  hist(peptide.ratios$log2.normalized.ratios.127n.126, main = NULL, 
       xlab = "Normalized log2 127n/126", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.127c.126, main = NULL, 
       xlab = "Normalized log2 127c/126", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.128n.126, main = NULL, 
       xlab = "Normalized log2 128n/126", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.128c.126, main = NULL, 
       xlab = "Normalized log2 128c/126", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.129n.126, main = NULL, 
       xlab = "Normalized log2 129n/126", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.129c.126, main = NULL, 
       xlab = "Normalized log2 129c/126", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.130n.126, main = NULL, 
       xlab = "Normalized log2 130n/126", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.130c.126, main = NULL, 
       xlab = "Normalized log2 130c/126", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.131.126, main = NULL, 
       xlab = "Normalized log2 131/126", col = "dodgerblue4");
  
  dev.off()
}

#based on ratio of each reporter and average intensity
check.log2.ratios.distribution.10plex2 <- function(peptide.ratios, peporpro) {
  file.name <- paste(peporpro, "log2_ratios_distribution_average.png", sep = "_");
  
  png(filename = file.name, width = 1500, height = 600);
  par(mfrow=c(2,5), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  hist(peptide.ratios$log2.normalized.ratios.126, main = NULL, 
       xlab = "Normalized log2 126/mean", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.127n, main = NULL, 
       xlab = "Normalized log2 127n/mean", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.127c, main = NULL, 
       xlab = "Normalized log2 127c/mean", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.128n, main = NULL, 
       xlab = "Normalized log2 128n/mean", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.128c, main = NULL, 
       xlab = "Normalized log2 128c/mean", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.129n, main = NULL, 
       xlab = "Normalized log2 129n/mean", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.129c, main = NULL, 
       xlab = "Normalized log2 129c/mean", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.130n, main = NULL, 
       xlab = "Normalized log2 130n/mean", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.130c, main = NULL, 
       xlab = "Normalized log2 130c/mean", col = "dodgerblue4");
  hist(peptide.ratios$log2.normalized.ratios.131, main = NULL, 
       xlab = "Normalized log2 131/mean", col = "dodgerblue4");
  
  dev.off()
}

check.ratios.intensity.counts <- function(peptide.ratios, method.name, peporpro) {
  file.name <- paste(peporpro, "Ratios_IntensityCount.png", sep = "_");
  
  png(filename = file.name, width = 500, height = 500);
  par(mfrow=c(1,1), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  method.name <- paste("Reporter.intensity.count.2.", method.name, sep="");
  boxplot(peptide.ratios$log2.ratios.126.128 ~ peptide.ratios[, c(method.name)], 
           xlab = "128 Reporter Intensity Count",
           ylab = "log2 126/128 ratios",
           main = NULL, col = NULL);
  dev.off()
  
  boxplot.data <- boxplot(peptide.ratios$log2.ratios.126.128 ~ peptide.ratios[, c(method.name)],
                          plot = FALSE);
  summary.data <- rbind(round(boxplot.data$n), round(boxplot.data$stats[3,], 2));
  colnames(summary.data) <- boxplot.data$names;
  rownames(summary.data) <- c("# of observations", "Median")
  summary.data
}

get.TMT.intensities <- function(peptide.ratios) {
  log10.intensity.126 <- log10(peptide.ratios$Reporter.intensity.0);
  log10.intensity.127 <- log10(peptide.ratios$Reporter.intensity.1);
  log10.intensity.128 <- log10(peptide.ratios$Reporter.intensity.2);
  log10.intensity.129 <- log10(peptide.ratios$Reporter.intensity.3);
  log10.intensity.130 <- log10(peptide.ratios$Reporter.intensity.4);
  log10.intensity.131 <- log10(peptide.ratios$Reporter.intensity.5);
  log10.intensity <- log10(peptide.ratios$Intensity);
  
  log10.intensity.126.128 <- log10.intensity.126 + log10.intensity.128;
  log10.intensity.127.128 <- log10.intensity.127 + log10.intensity.128;
  log10.intensity.130.129 <- log10.intensity.130 + log10.intensity.129;
  log10.intensity.131.129 <- log10.intensity.131 + log10.intensity.129;
  
  TMT.log10.intensities <- data.frame(log10.intensity.126,
                                      log10.intensity.127,
                                      log10.intensity.128,
                                      log10.intensity.129,
                                      log10.intensity.130,
                                      log10.intensity.131,
                                      log10.intensity,
                                      log10.intensity.126.128,
                                      log10.intensity.127.128,
                                      log10.intensity.130.129,
                                      log10.intensity.131.129
                                      );
  
  peptide.intensities <- cbind(peptide.ratios, TMT.log10.intensities);
  peptide.intensities
}

get.TMT.intensities.10plex <- function(peptide.ratios) {
  log10.intensity.126 <- log10(peptide.ratios$Reporter.intensity.0);
  log10.intensity.127n <- log10(peptide.ratios$Reporter.intensity.1);
  log10.intensity.127c <- log10(peptide.ratios$Reporter.intensity.2);
  log10.intensity.128n <- log10(peptide.ratios$Reporter.intensity.3);
  log10.intensity.128c <- log10(peptide.ratios$Reporter.intensity.4);
  log10.intensity.129n <- log10(peptide.ratios$Reporter.intensity.5);
  log10.intensity.129c <- log10(peptide.ratios$Reporter.intensity.6);
  log10.intensity.130n <- log10(peptide.ratios$Reporter.intensity.7);
  log10.intensity.130c <- log10(peptide.ratios$Reporter.intensity.8);
  log10.intensity.131 <- log10(peptide.ratios$Reporter.intensity.9);
  log10.intensity <- log10(peptide.ratios$Intensity);
  
  log10.normalized.intensity.126 <- log10(peptide.ratios$normalized.Reporter.intensity.0);
  log10.normalized.intensity.127n <- log10(peptide.ratios$normalized.Reporter.intensity.1);
  log10.normalized.intensity.127c <- log10(peptide.ratios$normalized.Reporter.intensity.2);
  log10.normalized.intensity.128n <- log10(peptide.ratios$normalized.Reporter.intensity.3);
  log10.normalized.intensity.128c <- log10(peptide.ratios$normalized.Reporter.intensity.4);
  log10.normalized.intensity.129n <- log10(peptide.ratios$normalized.Reporter.intensity.5);
  log10.normalized.intensity.129c <- log10(peptide.ratios$normalized.Reporter.intensity.6);
  log10.normalized.intensity.130n <- log10(peptide.ratios$normalized.Reporter.intensity.7);
  log10.normalized.intensity.130c <- log10(peptide.ratios$normalized.Reporter.intensity.8);
  log10.normalized.intensity.131 <- log10(peptide.ratios$normalized.Reporter.intensity.9);
  
  log2.normalized.intensity.126 <- log2(peptide.ratios$normalized.Reporter.intensity.0);
  log2.normalized.intensity.127n <- log2(peptide.ratios$normalized.Reporter.intensity.1);
  log2.normalized.intensity.127c <- log2(peptide.ratios$normalized.Reporter.intensity.2);
  log2.normalized.intensity.128n <- log2(peptide.ratios$normalized.Reporter.intensity.3);
  log2.normalized.intensity.128c <- log2(peptide.ratios$normalized.Reporter.intensity.4);
  log2.normalized.intensity.129n <- log2(peptide.ratios$normalized.Reporter.intensity.5);
  log2.normalized.intensity.129c <- log2(peptide.ratios$normalized.Reporter.intensity.6);
  log2.normalized.intensity.130n <- log2(peptide.ratios$normalized.Reporter.intensity.7);
  log2.normalized.intensity.130c <- log2(peptide.ratios$normalized.Reporter.intensity.8);
  log2.normalized.intensity.131 <- log2(peptide.ratios$normalized.Reporter.intensity.9);
  
#   log10.intensity.126.128 <- log10.intensity.126 + log10.intensity.128;
#   log10.intensity.127.128 <- log10.intensity.127 + log10.intensity.128;
#   log10.intensity.130.129 <- log10.intensity.130 + log10.intensity.129;
#   log10.intensity.131.129 <- log10.intensity.131 + log10.intensity.129;
  
  TMT.log10.intensities <- data.frame(log10.intensity.126,
                                      log10.intensity.127n,
                                      log10.intensity.127c,
                                      log10.intensity.128n,
                                      log10.intensity.128c,
                                      log10.intensity.129n,
                                      log10.intensity.129c,
                                      log10.intensity.130n,
                                      log10.intensity.130c,
                                      log10.intensity.131,
                                      log10.normalized.intensity.126,
                                      log10.normalized.intensity.127n,
                                      log10.normalized.intensity.127c,
                                      log10.normalized.intensity.128n,
                                      log10.normalized.intensity.128c,
                                      log10.normalized.intensity.129n,
                                      log10.normalized.intensity.129c,
                                      log10.normalized.intensity.130n,
                                      log10.normalized.intensity.130c,
                                      log10.normalized.intensity.131,
                                      log2.normalized.intensity.126,
                                      log2.normalized.intensity.127n,
                                      log2.normalized.intensity.127c,
                                      log2.normalized.intensity.128n,
                                      log2.normalized.intensity.128c,
                                      log2.normalized.intensity.129n,
                                      log2.normalized.intensity.129c,
                                      log2.normalized.intensity.130n,
                                      log2.normalized.intensity.130c,
                                      log2.normalized.intensity.131,
                                      log10.intensity
  );
  
  peptide.intensities <- cbind(peptide.ratios, TMT.log10.intensities);
  peptide.intensities
}

check.log10.intensity.distribution <- function(peptide.intensities, peporpro) {
  file.name <- paste(peporpro, "Ratios_IntensityCount.png", sep = "_");
  
  png(filename = file.name, width = 1000, height = 500);
  par(mfrow=c(2,3), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  hist(peptide.intensities$log10.intensity.126, main = NULL, 
       xlab = "log10 126 Intensity", col = "dodgerblue4");
  hist(peptide.intensities$log10.intensity.127, main = NULL, 
       xlab = "log10 127 Intensity", col = "dodgerblue4");
  hist(peptide.intensities$log10.intensity.128, main = NULL, 
       xlab = "log10 128 Intensity", col = "dodgerblue4");
  hist(peptide.intensities$log10.intensity.129, main = NULL, 
       xlab = "log10 129 Intensity", col = "dodgerblue4");
  hist(peptide.intensities$log10.intensity.130, main = NULL, 
       xlab = "log10 130 Intensity", col = "dodgerblue4");
  hist(peptide.intensities$log10.intensity.131, main = NULL, 
       xlab = "log10 131 Intensity", col = "dodgerblue4");  
  dev.off()
}

check.log10.intensity.distribution.10plex <- function(peptide.intensities, peporpro) {
  file.name <- paste(peporpro, "Ratios_IntensityCount.png", sep = "_");
  
  png(filename = file.name, width = 1500, height = 600);
  par(mfrow=c(2,5), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  hist(peptide.intensities$log10.intensity.126, main = NULL, 
       xlab = "log10 126 Intensity", col = "dodgerblue4");
  hist(peptide.intensities$log10.intensity.127n, main = NULL, 
       xlab = "log10 127n Intensity", col = "dodgerblue4");
  hist(peptide.intensities$log10.intensity.127c, main = NULL, 
       xlab = "log10 127c Intensity", col = "dodgerblue4");
  hist(peptide.intensities$log10.intensity.128n, main = NULL, 
       xlab = "log10 128n Intensity", col = "dodgerblue4");
  hist(peptide.intensities$log10.intensity.128c, main = NULL, 
       xlab = "log10 128c Intensity", col = "dodgerblue4");
  hist(peptide.intensities$log10.intensity.129n, main = NULL, 
       xlab = "log10 129n Intensity", col = "dodgerblue4");
  hist(peptide.intensities$log10.intensity.129c, main = NULL, 
       xlab = "log10 129c Intensity", col = "dodgerblue4");
  hist(peptide.intensities$log10.intensity.130n, main = NULL, 
       xlab = "log10 130n Intensity", col = "dodgerblue4");
  hist(peptide.intensities$log10.intensity.130c, main = NULL, 
       xlab = "log10 130c Intensity", col = "dodgerblue4");
  hist(peptide.intensities$log10.intensity.131, main = NULL, 
       xlab = "log10 131 Intensity", col = "dodgerblue4");  
  dev.off()
}

Ratios.vs.Intensity <- function(peptide.intensities, peporpro) {
  file.name <- paste(peporpro, "Ratios_vs_Intensity.png", sep = "_");
  
  png(filename = file.name, width = 1000, height = 1000);
  par(mfrow=c(3,3), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  plot(peptide.intensities$log10.intensity.128, peptide.intensities$log2.ratios.126.128, 
       xlab = "Log10 128 Intensity", ylab = "Log2 126/128 Ratio", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log10.intensity.128, peptide.intensities$log2.ratios.127.128, 
       xlab = "Log10 128 Intensity", ylab = "Log2 127/128 Ratio", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log10.intensity.128, peptide.intensities$log2.ratios.128.128, 
       xlab = "Log10 128 Intensity", ylab = "Log2 128/128 Ratio", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log10.intensity.128, peptide.intensities$log2.ratios.129.128, 
       xlab = "Log10 128 Intensity", ylab = "Log2 129/128 Ratio", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log10.intensity.128, peptide.intensities$log2.ratios.130.128, 
       xlab = "Log10 128 Intensity", ylab = "Log2 130/128 Ratio", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log10.intensity.128, peptide.intensities$log2.ratios.131.128, 
       xlab = "Log10 128 Intensity", ylab = "Log2 131/128 Ratio", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log10.intensity.129, peptide.intensities$log2.ratios.129.129, 
       xlab = "Log10 129 Intensity", ylab = "Log2 129/129 Ratio", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log10.intensity.129, peptide.intensities$log2.ratios.130.129, 
       xlab = "Log10 129 Intensity", ylab = "Log2 130/129 Ratio", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log10.intensity.129, peptide.intensities$log2.ratios.131.129, 
       xlab = "Log10 129 Intensity", ylab = "Log2 131/129 Ratio", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  dev.off()
}

Ratios.vs.Intensity2 <- function(peptide.intensities, peporpro) {
  file.name <- paste(peporpro, "Ratios_vs_Intensity2.png", sep = "_");
  
  png(filename = file.name, width = 900, height = 600);
  par(mfrow=c(2,3), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  plot(peptide.intensities$log2.ratios.126.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 126/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log2.ratios.127.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 127/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log2.ratios.128.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 128/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log2.ratios.129.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 129/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log2.ratios.130.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 130/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log2.ratios.131.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 131/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  dev.off()
}

Ratios.vs.Intensity.10plex <- function(peptide.intensities, peporpro) {
  file.name <- paste(peporpro, "Ratios_vs_Intensity.png", sep = "_");
  
  png(filename = file.name, width = 900, height = 900);
  par(mfrow=c(3,3), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  plot(peptide.intensities$log2.ratios.127n.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 127n/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.ratios.127c.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 127c/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.ratios.128n.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 128n/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.ratios.128c.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 128c/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.ratios.129n.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 129n/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.ratios.129c.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 129c/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.ratios.130n.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 130n/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.ratios.130c.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 130c/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.ratios.131.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 131/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  dev.off()
}

Ratios.vs.Intensity.10plex2 <- function(peptide.intensities, peporpro) {
  file.name <- paste(peporpro, "AverageRatios_vs_Intensity.png", sep = "_");
  
  png(filename = file.name, width = 1500, height = 600);
  par(mfrow=c(2,5), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  plot(peptide.intensities$log2.ratios.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 126/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.ratios.127n, peptide.intensities$log10.intensity, 
       xlab = "Log2 127n/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.ratios.127c, peptide.intensities$log10.intensity, 
       xlab = "Log2 127c/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.ratios.128n, peptide.intensities$log10.intensity, 
       xlab = "Log2 128n/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.ratios.128c, peptide.intensities$log10.intensity, 
       xlab = "Log2 128c/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.ratios.129n, peptide.intensities$log10.intensity, 
       xlab = "Log2 129n/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.ratios.129c, peptide.intensities$log10.intensity, 
       xlab = "Log2 129c/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.ratios.130n, peptide.intensities$log10.intensity, 
       xlab = "Log2 130n/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.ratios.130c, peptide.intensities$log10.intensity, 
       xlab = "Log2 130c/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.ratios.131, peptide.intensities$log10.intensity, 
       xlab = "Log2 131/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  dev.off()
}

#plot the normalized ratios
Ratios.vs.Intensity3 <- function(peptide.intensities, peporpro) {
  file.name <- paste(peporpro, "Ratios_vs_Intensity3.png", sep = "_");
  
  png(filename = file.name, width = 900, height = 600);
  par(mfrow=c(2,3), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  plot(peptide.intensities$log2.normalized.ratios.126.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 126/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log2.normalized.ratios.127.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 127/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log2.normalized.ratios.128.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 128/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log2.normalized.ratios.129.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 129/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log2.normalized.ratios.130.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 130/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log2.normalized.ratios.131.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 131/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  dev.off()
}

#plot the normalized ratios
Ratios.vs.normalized.Intensity.10plex <- function(peptide.intensities, peporpro) {
  file.name <- paste(peporpro, "Ratios_vs_Normalized_Intensity.png", sep = "_");
  
  png(filename = file.name, width = 900, height = 900);
  par(mfrow=c(3,3), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  plot(peptide.intensities$log2.normalized.ratios.127n.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 127n/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.normalized.ratios.127c.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 127c/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.normalized.ratios.128n.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 128n/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.normalized.ratios.128c.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 128c/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.normalized.ratios.129n.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 129n/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.normalized.ratios.129c.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 129c/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.normalized.ratios.130n.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 130n/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.normalized.ratios.130c.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 130c/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.normalized.ratios.131.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 131/126 Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  dev.off()
}

Ratios.vs.normalized.Intensity.10plex2 <- function(peptide.intensities, peporpro) {
  file.name <- paste(peporpro, "AverageRatios_vs_Normalized_Intensity.png", sep = "_");
  
  png(filename = file.name, width = 1500, height = 600);
  par(mfrow=c(2,5), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  plot(peptide.intensities$log2.normalized.ratios.126, peptide.intensities$log10.intensity, 
       xlab = "Log2 126/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.normalized.ratios.127n, peptide.intensities$log10.intensity, 
       xlab = "Log2 127n/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.normalized.ratios.127c, peptide.intensities$log10.intensity, 
       xlab = "Log2 127c/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.normalized.ratios.128n, peptide.intensities$log10.intensity, 
       xlab = "Log2 128n/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.normalized.ratios.128c, peptide.intensities$log10.intensity, 
       xlab = "Log2 128c/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.normalized.ratios.129n, peptide.intensities$log10.intensity, 
       xlab = "Log2 129n/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.normalized.ratios.129c, peptide.intensities$log10.intensity, 
       xlab = "Log2 129c/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.normalized.ratios.130n, peptide.intensities$log10.intensity, 
       xlab = "Log2 130n/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.normalized.ratios.130c, peptide.intensities$log10.intensity, 
       xlab = "Log2 130c/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  plot(peptide.intensities$log2.normalized.ratios.131, peptide.intensities$log10.intensity, 
       xlab = "Log2 131/mean Ratio", ylab = "Log10 Intensity", pch = 16, 
       col = "dodgerblue4", cex = 0.6, xlim=c(-4, 4));
  dev.off()
}


#R-I plot
Ratios.vs.Intensity4 <- function(peptide.intensities, peporpro, cutoff = 0, reporter = "128") {
  base.reporter <- paste("log10.intensity.", reporter, sep = "");
  peptide.intensities <- peptide.intensities[which(peptide.intensities[, c(base.reporter)] > log10(cutoff)), ];
  file.name <- paste(peporpro, "Ratios_vs_Intensity4", sep = "_");
  file.name <- paste(file.name, toString(cutoff), sep = "");
  file.name <- paste(file.name, ".png", sep = "");
  
  png(filename = file.name, width = 600, height = 600);
  par(mfrow=c(2,2), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  plot(peptide.intensities$log10.intensity.126.128, 
       peptide.intensities$log2.ratios.126.128,
       xlab = "Log10 Intensity (126*128)", ylab = "Log2 (126/128)", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log10.intensity.127.128, 
       peptide.intensities$log2.ratios.127.128,
       xlab = "Log10 Intensity (127*128)", ylab = "Log2 (127/128)", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log10.intensity.130.129, 
       peptide.intensities$log2.ratios.130.129,
       xlab = "Log10 Intensity (130*129)", ylab = "Log2 (130/129)", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  plot(peptide.intensities$log10.intensity.131.129, 
       peptide.intensities$log2.ratios.131.129,
       xlab = "Log10 Intensity (131*129)", ylab = "Log2 (131/129)", pch = 16, 
       col = "dodgerblue4", cex = 0.6);
  
  dev.off()
}


means.sds.plot <- function(peptide.intensities, cutoff, reporter, peporpro) {
  base.reporter <- paste("log10.intensity.", reporter, sep = "");
  peptide.intensities <- peptide.intensities[which(peptide.intensities[, c(base.reporter)] > log10(cutoff)), ];
  cat("# of peptides above ", cutoff, "is", nrow(peptide.intensities), "\n\n");
  
  ratio.means <- c(mean(peptide.intensities$ratios.126.128), 
                   mean(peptide.intensities$ratios.127.128), 
#                    mean(peptide.intensities$ratios.128.128),
#                    mean(peptide.intensities$ratios.129.128),
#                    mean(peptide.intensities$ratios.130.128),
#                    mean(peptide.intensities$ratios.131.128),
#                    mean(peptide.intensities$ratios.129.129),
                   mean(peptide.intensities$ratios.130.129),
                   mean(peptide.intensities$ratios.131.129));
  
  ratio.sds <- c(sd(peptide.intensities$ratios.126.128),
                 sd(peptide.intensities$ratios.127.128),
#                  sd(peptide.intensities$ratios.128.128),
#                  sd(peptide.intensities$ratios.129.128),
#                  sd(peptide.intensities$ratios.130.128),
#                  sd(peptide.intensities$ratios.131.128),
#                  sd(peptide.intensities$ratios.129.129),
                 sd(peptide.intensities$ratios.130.129),
                 sd(peptide.intensities$ratios.131.129));
  
  ratio.medians <- c(median(peptide.intensities$ratios.126.128),
                     median(peptide.intensities$ratios.127.128),
#                      median(peptide.intensities$ratios.128.128),
#                      median(peptide.intensities$ratios.129.128),
#                      median(peptide.intensities$ratios.130.128),
#                      median(peptide.intensities$ratios.131.128),
#                      median(peptide.intensities$ratios.129.129),
                     median(peptide.intensities$ratios.130.129),
                     median(peptide.intensities$ratios.131.129));
  
  ratio.list <- c(peptide.intensities$ratios.126.128,
                  peptide.intensities$ratios.127.128,
                  peptide.intensities$ratios.128.128,
                  peptide.intensities$ratios.129.128,
                  peptide.intensities$ratios.130.128,
                  peptide.intensities$ratios.131.128);
  
  ratio.factor <- c(rep("126", length(peptide.intensities$ratios.126.128)),
                    rep("127", length(peptide.intensities$ratios.127.128)),
                    rep("128", length(peptide.intensities$ratios.128.128)),
                    rep("129", length(peptide.intensities$ratios.129.128)),
                    rep("130", length(peptide.intensities$ratios.130.128)),
                    rep("131", length(peptide.intensities$ratios.131.128)));
  
  ratio.by.group <- data.frame(ratio.list, ratio.factor);
  
  file.name <- paste(peporpro, "_means_sds_", sep = "");
  file.name <- paste(file.name, toString(cutoff), sep = "");
  file.name <- paste(file.name, ".png", sep = "");
  png(filename = file.name, width = 1000, height = 500);
  par(mfrow=c(1,2), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  bargraph.CI(x.factor = ratio.factor, response = ratio.list, 
              data = ratio.by.group, col = "dodgerblue4", ylab = "Fold",
              ci.fun = function(x) {c(mean(x) - sd(x), 
                                      mean(x) + sd(x))});
  
  ratio.list <- c(peptide.intensities$ratios.129.129,
                  peptide.intensities$ratios.130.129,
                  peptide.intensities$ratios.131.129);
  
  ratio.factor <- c(rep("129", length(peptide.intensities$ratios.129.129)),
                    rep("130", length(peptide.intensities$ratios.130.129)),
                    rep("131", length(peptide.intensities$ratios.131.129)));
  
  ratio.by.group <- data.frame(ratio.list, ratio.factor);
  
  bargraph.CI(x.factor = ratio.factor, response = ratio.list, 
              data = ratio.by.group, col = "dodgerblue4", ylab = "Fold",
              ci.fun = function(x) {c(mean(x) - sd(x), 
                                      mean(x) + sd(x))});
  dev.off()
  
  means.sds <- c(ratio.means, ratio.sds, ratio.medians);
  means.sds <- round(means.sds, 2);
#   colnames(means.sds) <- c("126/128", "127/128", "128/128",
#                            "129/128", "130/128", "131/128",
#                            "129/129", "130/129", "131/129");
#   rownames(means.sds) <- c("Mean", "S.D.", "Median");
  means.sds
}

TMT.peptide.analysis <- function (method1, method2) {
  work.directory <- paste("C:/Users/huangf01/Documents/software/MaxQuant_1.5.2.8_SearchFolder/TMT evaluation 2/combined/txt fraction d7 ", 
                          method1, sep = "");
  setwd(work.directory);
  
  peptides.info <- TMT.peptides.cleanup(3000);
  pc.method <- paste("pc.", method1, sep ="");
  assign(pc.method, as.numeric(unlist(peptides.info[2])));
  peptides <- as.data.frame(peptides.info[1]);
  
  peptide.ratios <- get.TMT.ratios(peptides);
  check.log2.ratios.distribution(peptide.ratios, "peptide");
  check.ratios.intensity.counts(peptide.ratios, method2, "peptide");
  peptide.intensities <- get.TMT.intensities(peptide.ratios);
  check.log10.intensity.distribution(peptide.intensities, "peptide");
  Ratios.vs.Intensity(peptide.intensities, "peptide");
  msm.method <- paste("msm.", method1, sep = "");
  assign(msm.method, means.sds.plot(peptide.intensities, cutoff = 10, reporter = "128", "peptide"));
  
  output <- list(get(pc.method), get(msm.method));
  output;
}

count.zero.values <- function(proteins) {
  for (i in 1:ncol(proteins)) {
    cat(colnames(proteins)[i], 'has', sum(proteins[, i] == 0), 'values\n');
  }
}

replace.zero.values <- function(proteins, colrange, replicate = 3) {
  second.min.col <- NULL;
  for (n in colrange) {
     intensity <- proteins[, n];
     intensity <- unique(intensity);
     second.min <- sort(intensity)[2];
     second.min.col <- c(second.min.col, second.min);
  }
  cat("The miminum intensity is ", min(second.min.col), "\n");
  cat ("Press [enter] to continue");
  line <- readline();
  
  one.measurement <- NULL;
  two.measurement <- NULL;
  three.measurement <- NULL;
  two.measurement.spread <- NULL;
  
  for (j in 1:(length(colrange)/replicate)) {
    cat(j, "\n");
#     cat ("Press [enter] to continue");
#     line <- readline();
    for (i in 1:nrow(proteins)) {
      cat(i, "\n");
      col.range <- (min(colrange)+(j-1)*replicate): (min(colrange)+(j-1)*replicate+replicate-1);
      if (sum(proteins[i, col.range] == 0) == 0) {
        three.measurement <- c(three.measurement, unlist(proteins[i, col.range]));
      }
      else if (sum(proteins[i, col.range] == 0) == 1) {
        two.values <- unlist(proteins[i, col.range][proteins[i, col.range] > 0]);
        two.measurement <- c(two.measurement, two.values);
        two.measurement.spread <- c(two.measurement.spread,
                                    (max(two.values) - min(two.values))*2/(max(two.values) + min(two.values)));
        match.col <- match(0, proteins[i, col.range]);
        match.col <- min(col.range) + match.col - 1;
        proteins[i, match.col] <- min(proteins[i, col.range][proteins[i, col.range] > 0]);
      }
      else if (sum(proteins[i, col.range] == 0) == 2) {
        one.measurement <- c(one.measurement, unlist(max(proteins[i, col.range])));
        proteins[i, col.range] <- max(proteins[i, col.range]);
       }
      else if (sum(proteins[i, col.range] == 0) == 3) {
        proteins[i, col.range] <- min(second.min.col);
      }
    }
  }
  cat("Number of one measurement is ", length(one.measurement), "\n");
  cat("Number of two measurement is ", length(two.measurement)/2, "\n");
  cat("Number of three measurement is ", length(three.measurement)/3, "\n");

  file.name <- "IntensityDistribution_byMeasurement.png";

  png(filename = file.name, width = 300, height = 900);
  par(mfrow=c(3,1), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
    col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
    mgp = c(3, 1, 0));
  hist(log10(three.measurement), main = NULL, 
     xlab = "Three measurement intensity", col = "dodgerblue4");
  hist(log10(two.measurement), main = NULL, xlim = c(5,11),
     xlab = "Two measurement intensity", col = "dodgerblue4");
  hist(log10(one.measurement), main = NULL, xlim = c(5,11),
     xlab = "One measurement intensity", col = "dodgerblue4");
  dev.off();

  file.name <- "SpreadDistribution.png";

  png(filename = file.name, width = 300, height = 300);
  par(mfrow=c(1,1), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
    col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
    mgp = c(3, 1, 0));
  hist(two.measurement.spread, main = NULL, 
     xlab = "Relative two measurement spread", col = "dodgerblue4");
  dev.off();
  
  summary.mean <- c(mean(log10(one.measurement)), 
                    mean(log10(two.measurement)), 
                    mean(log10(three.measurement)));
  summary.median <- c(median(log10(one.measurement)), 
                      median(log10(two.measurement)), 
                      median(log10(three.measurement)));
  summary.sd <- c(sd(log10(one.measurement)), 
                      sd(log10(two.measurement)), 
                      sd(log10(three.measurement)));
  summary.statistics <- cbind(summary.mean, summary.median, summary.sd);
  rownames(summary.statistics) <- c("One", "Two", "Three");
  colnames(summary.statistics) <- c("Mean", "Median", "S.D.");

  write.csv(summary.statistics, "statistics.csv", row.names = TRUE);

#   cat("Summary statistics\n", summary.statistics, "\n");
#   cat("Summary statistics for two measurement", summary(two.measurement), "\n");
#   cat("Summary statistics for three measurement", summary(three.measurement), "\n");

  proteins
}

boxplot.panel <- function(df, parameter, file.name) {
  png(filename = file.name, width = 2500, height = 1500);
  par(mfrow=c(3,5), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  boxplot(df[, c("SelectedNeuronCount")]~df[, parameter],data=df,
          ylab = "SelectedNeuronCount",
          main = "SelectedNeuronCount", col = NULL);
  boxplot(df[, c("NeuriteTotalCountPerNeuronCh2")]~df[, parameter],data=df,
          ylab = "NeuriteTotalCountPerNeuronCh2",
          main = NULL, col = NULL);
  boxplot(df[, c("NeuriteTotalLengthPerNeuronCh2")]~df[, parameter],data=df,
          ylab = "NeuriteTotalLengthPerNeuronCh2",
          main = NULL, col = NULL);
  boxplot(df[, c("NeuriteTotalLengthPerNeuriteCh2")]~df[, parameter],data=df,
          ylab = "NeuriteTotalLengthPerNeuriteCh2",
          main = NULL, col = NULL);
  boxplot(df[, c("BranchPointTotalCountPerNeuronCh2")]~df[, parameter],data=df,
          ylab = "BranchPointTotalCountPerNeuronCh2",
          main = NULL, col = NULL);
  boxplot(df[, c("BranchPointTotalCountPerNeuriteCh2")]~df[, parameter],data=df,
          ylab = "BranchPointTotalCountPerNeuriteCh2",
          main = NULL, col = NULL);
  boxplot(df[, c("BranchPointCountPerNeuriteLengthCh2")]~df[, parameter],data=df,
          ylab = "BranchPointCountPerNeuriteLengthCh2",
          main = NULL, col = NULL);
  boxplot(df[, c("MEAN_NeuriteTotalCountCh2")]~df[, parameter],data=df,
          ylab = "MEAN_NeuriteTotalCountCh2",
          main = NULL, col = NULL);
  boxplot(df[, c("MEAN_NeuriteTotalLengthCh2")]~df[, parameter],data=df,
          ylab = "MEAN_NeuriteTotalLengthCh2",
          main = NULL, col = NULL);
  boxplot(df[, c("MEAN_NeuriteAvgLengthCh2")]~df[, parameter],data=df,
          ylab = "MEAN_NeuriteAvgLengthCh2",
          main = NULL, col = NULL);
  boxplot(df[, c("MEAN_NeuriteMaxLengthWithBranchesCh2")]~df[, parameter],data=df,
          ylab = "MEAN_NeuriteMaxLengthWithBranchesCh2",
          main = NULL, col = NULL);
  boxplot(df[, c("MEAN_NeuriteMaxLengthWithoutBranchesCh2")]~df[, parameter],data=df,
          ylab = "MEAN_NeuriteMaxLengthWithoutBranchesCh2",
          main = NULL, col = NULL);
  boxplot(df[, c("MEAN_NeuriteTotalAreaCh2")]~df[, parameter],data=df,
          ylab = "MEAN_NeuriteTotalAreaCh2",
          main = NULL, col = NULL);
  boxplot(df[, c("MEAN_NeuriteWidthCh2")]~df[, parameter],data=df,
          ylab = "MEAN_NeuriteWidthCh2",
          main = NULL, col = NULL);
  boxplot(df[, c("MEAN_NeuriteAvgIntenCh2")]~df[, parameter],data=df,
          ylab = "MEAN_NeuriteAvgIntenCh2",
          main = NULL, col = NULL);
  dev.off()
}

boxplot.multi.factor <- function(df, file.name) {
  library(ggplot2);
  library(grid);
  library(gridExtra);
  p1 <- ggplot(data = df, aes(x = siRNA, y = SelectedNeuronCount)) + 
    geom_boxplot(aes(fill = Treatment), width = 0.8) + theme_bw() + 
    ggtitle("SelectedNeuronCount");
  p2 <- ggplot(data = df, aes(x = siRNA, y = NeuriteTotalCountPerNeuronCh2)) + 
    geom_boxplot(aes(fill = Treatment), width = 0.8) + theme_bw() + 
    ggtitle("NeuriteTotalCountPerNeuronCh2");
  p3 <- ggplot(data = df, aes(x = siRNA, y = NeuriteTotalLengthPerNeuronCh2)) + 
    geom_boxplot(aes(fill = Treatment), width = 0.8) + theme_bw() + 
    ggtitle("NeuriteTotalLengthPerNeuronCh2");
  p4 <- ggplot(data = df, aes(x = siRNA, y = NeuriteTotalLengthPerNeuriteCh2)) + 
    geom_boxplot(aes(fill = Treatment), width = 0.8) + theme_bw() + 
    ggtitle("NeuriteTotalLengthPerNeuriteCh2");
  p5 <- ggplot(data = df, aes(x = siRNA, y = BranchPointTotalCountPerNeuronCh2)) + 
    geom_boxplot(aes(fill = Treatment), width = 0.8) + theme_bw() + 
    ggtitle("BranchPointTotalCountPerNeuronCh2");
  p6 <- ggplot(data = df, aes(x = siRNA, y = BranchPointTotalCountPerNeuriteCh2)) + 
    geom_boxplot(aes(fill = Treatment), width = 0.8) + theme_bw() + 
    ggtitle("BranchPointTotalCountPerNeuriteCh2");
  p7 <- ggplot(data = df, aes(x = siRNA, y = BranchPointCountPerNeuriteLengthCh2)) + 
    geom_boxplot(aes(fill = Treatment), width = 0.8) + theme_bw() + 
    ggtitle("BranchPointCountPerNeuriteLengthCh2");
  p8 <- ggplot(data = df, aes(x = siRNA, y = MEAN_NeuriteTotalCountCh2)) + 
    geom_boxplot(aes(fill = Treatment), width = 0.8) + theme_bw() + 
    ggtitle("MEAN_NeuriteTotalCountCh2");
  p9 <- ggplot(data = df, aes(x = siRNA, y = MEAN_NeuriteTotalLengthCh2)) + 
    geom_boxplot(aes(fill = Treatment), width = 0.8) + theme_bw() + 
    ggtitle("MEAN_NeuriteTotalLengthCh2");
  p10 <- ggplot(data = df, aes(x = siRNA, y = MEAN_NeuriteAvgLengthCh2)) + 
    geom_boxplot(aes(fill = Treatment), width = 0.8) + theme_bw() + 
    ggtitle("MEAN_NeuriteAvgLengthCh2");
  p11 <- ggplot(data = df, aes(x = siRNA, y = MEAN_NeuriteMaxLengthWithBranchesCh2)) + 
    geom_boxplot(aes(fill = Treatment), width = 0.8) + theme_bw() + 
    ggtitle("MEAN_NeuriteMaxLengthWithBranchesCh2");
  p12 <- ggplot(data = df, aes(x = siRNA, y = MEAN_NeuriteMaxLengthWithoutBranchesCh2)) + 
    geom_boxplot(aes(fill = Treatment), width = 0.8) + theme_bw() + 
    ggtitle("MEAN_NeuriteMaxLengthWithoutBranchesCh2");
  p13 <- ggplot(data = df, aes(x = siRNA, y = MEAN_NeuriteTotalAreaCh2)) + 
    geom_boxplot(aes(fill = Treatment), width = 0.8) + theme_bw() + 
    ggtitle("MEAN_NeuriteTotalAreaCh2");
  p14 <- ggplot(data = df, aes(x = siRNA, y = MEAN_NeuriteWidthCh2)) + 
    geom_boxplot(aes(fill = Treatment), width = 0.8) + theme_bw() + 
    ggtitle("MEAN_NeuriteWidthCh2");
  p15 <- ggplot(data = df, aes(x = siRNA, y = MEAN_NeuriteAvgIntenCh2)) + 
    geom_boxplot(aes(fill = Treatment), width = 0.8) + theme_bw() + 
    ggtitle("MEAN_NeuriteAvgIntenCh2");
  
  png(filename = file.name, width = 1500, height = 1200);
  grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15, ncol = 3);
  dev.off();
  
  
#   grobula <- arrangeGrob(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,
#                          ncol=5)
#   
#   grobula <- arrangeGrob(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,
#                          ncol=5, nrow=3, 
#                          top = textGrob(file.name, gp = gpar(fontsize=18, fontface="bold.italic", fontsize=18))
#                          );
#   print(grobula)
#   ggsave <- ggplot2::ggsave; body(ggsave) <- body(ggplot2::ggsave)[-2]
#   ggsave(grobula, file = file.name, width = 40, height = 24, dpi = 300);
#   ggsave(filename = file.name, plot = grobula, width = 40, height = 24, dpi = 300);
  
}

create.protein.ArrayTrack.files <- function(proteins, Assay, Sample) {
  #need to load xlsx library
  library(xlsx);
  
  #creat array file
  
  chosen.columns <- c("id", "Fasta.headers", "Mol..weight..kDa.",
                      "Sequence.coverage....");
  
  array.file <- proteins[, c(chosen.columns, "Reporter.intensity.0",
                             "Reporter.intensity.1", "Reporter.intensity.2",
                             "Reporter.intensity.3", "Reporter.intensity.4",
                             "Reporter.intensity.5", "Reporter.intensity.6",
                             "Reporter.intensity.7", "Reporter.intensity.8",
                             "Reporter.intensity.9")];
  
  #simplify the protein IDs: this may change depending on the database you use
#   first.Protein.IDs <- get.first.IDs(proteins[, c("Protein.IDs")]);
#   proteins <- cbind(first.Protein.IDs, proteins);
#   array.file <- cbind(first.Protein.IDs, array.file);
  
  #create sample files: one file per sample!
  sample.126 <- proteins[, c(chosen.columns,
                             "Reporter.intensity.0")];  
  sample.127n <- proteins[, c(chosen.columns,
                              "Reporter.intensity.1")];  
  sample.127c <- proteins[, c(chosen.columns,
                              "Reporter.intensity.2")];  
  sample.128n <- proteins[, c(chosen.columns,
                              "Reporter.intensity.3")];  
  sample.128c <- proteins[, c(chosen.columns,
                              "Reporter.intensity.4")];  
  sample.129n <- proteins[, c(chosen.columns,
                              "Reporter.intensity.5")];  
  sample.129c <- proteins[, c(chosen.columns,
                              "Reporter.intensity.6")];  
  sample.130n <- proteins[, c(chosen.columns,
                              "Reporter.intensity.7")];  
  sample.130c <- proteins[, c(chosen.columns,
                              "Reporter.intensity.8")];  
  sample.131 <- proteins[, c(chosen.columns,
                             "Reporter.intensity.9")];  
  colnames(sample.126) <- c(chosen.columns,
                            "Intensity");
  colnames(sample.127n) <- c(chosen.columns,
                             "Intensity");
  colnames(sample.127c) <- c(chosen.columns,
                             "Intensity");
  colnames(sample.128n) <- c(chosen.columns,
                             "Intensity");
  colnames(sample.128c) <- c(chosen.columns,
                             "Intensity");
  colnames(sample.129n) <- c(chosen.columns,
                             "Intensity");
  colnames(sample.129c) <- c(chosen.columns,
                             "Intensity");
  colnames(sample.130n) <- c(chosen.columns,
                             "Intensity");
  colnames(sample.130c) <- c(chosen.columns,
                             "Intensity");
  colnames(sample.131) <- c(chosen.columns,
                            "Intensity");
  
  Intensity <- c(sum(sample.126[,c("Intensity")]), 
                 sum(sample.127n[,c("Intensity")]), sum(sample.127c[,c("Intensity")]),
                 sum(sample.128n[,c("Intensity")]), sum(sample.128c[,c("Intensity")]),
                 sum(sample.129n[,c("Intensity")]), sum(sample.129c[,c("Intensity")]),
                 sum(sample.130n[,c("Intensity")]), sum(sample.130c[,c("Intensity")]),
                 sum(sample.131[,c("Intensity")]));
  
  Filename <- c("protein_126.txt", 
                "protein_127n.txt", "protein_127c.txt",
                "protein_128n.txt", "protein_128c.txt",
                "protein_129n.txt", "protein_129c.txt",
                "protein_130n.txt", "protein_130c.txt",
                "protein_131.txt");
  
  hybridization.file <- cbind(Assay, Sample, Intensity, Filename);
  
  cat("Writing Array file\n");
  write.table(array.file, "protein_array.txt", sep="\t", row.names = FALSE, quote = FALSE);
  
  cat("Writing Sample 126 file\n");
  write.table(sample.126, "protein_126.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 127n file\n");
  write.table(sample.127n, "protein_127n.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 127c file\n");
  write.table(sample.127c, "protein_127c.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 128n file\n");
  write.table(sample.128n, "protein_128n.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 128c file\n");
  write.table(sample.128c, "protein_128c.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 129n file\n");
  write.table(sample.129n, "protein_129n.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 129c file\n");
  write.table(sample.129c, "protein_129c.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 130n file\n");
  write.table(sample.130n, "protein_130n.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 130c file\n");
  write.table(sample.130c, "protein_130c.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 131 file\n");
  write.table(sample.131, "protein_131.txt", sep="\t", row.names = FALSE, quote = FALSE);
  
  cat("Writing hybridization file\n");
  write.table(hybridization.file, "protein_hybridization.txt", sep="\t", row.names = FALSE, quote = FALSE);
}

creat.peptide.ArrayTrack.files <- function(peptides, Assay, Sample) {
  #need to load txt library
  library(xlsx);
  
  chosen.columns <- c("id", "Leading.razor.protein", "PEP", "Score", "Charges");
  
  array.file <- peptides[, c(chosen.columns, "Reporter.intensity.0",
                             "Reporter.intensity.1", "Reporter.intensity.2",
                             "Reporter.intensity.3", "Reporter.intensity.4",
                             "Reporter.intensity.5", "Reporter.intensity.6",
                             "Reporter.intensity.7", "Reporter.intensity.8",
                             "Reporter.intensity.9")];
  
#   #simplify the protein IDs: this may change depending on the database you use
#   first.Protein.IDs <- get.first.IDs(peptides[, c("Protein.IDs")]);
#   peptides <- cbind(first.Protein.IDs, peptides);
#   array.file <- cbind(first.Protein.IDs, array.file);
  
  #create sample files: one file per sample!
  sample.126 <- peptides[, c(chosen.columns, "Reporter.intensity.0")];
  sample.127n <- peptides[, c(chosen.columns, "Reporter.intensity.1")];
  sample.127c <- peptides[, c(chosen.columns, "Reporter.intensity.2")];
  sample.128n <- peptides[, c(chosen.columns, "Reporter.intensity.3")];
  sample.128c <- peptides[, c(chosen.columns, "Reporter.intensity.4")];
  sample.129n <- peptides[, c(chosen.columns, "Reporter.intensity.5")];
  sample.129c <- peptides[, c(chosen.columns, "Reporter.intensity.6")];
  sample.130n <- peptides[, c(chosen.columns, "Reporter.intensity.7")];
  sample.130c <- peptides[, c(chosen.columns, "Reporter.intensity.8")];
  sample.131 <- peptides[, c(chosen.columns, "Reporter.intensity.9")];

  colnames(sample.126) <- c(chosen.columns, "Intensity");
  colnames(sample.127n) <- c(chosen.columns, "Intensity");
  colnames(sample.127c) <- c(chosen.columns, "Intensity");
  colnames(sample.128n) <- c(chosen.columns, "Intensity");
  colnames(sample.128c) <- c(chosen.columns, "Intensity");
  colnames(sample.129n) <- c(chosen.columns, "Intensity");
  colnames(sample.129c) <- c(chosen.columns, "Intensity");
  colnames(sample.130n) <- c(chosen.columns, "Intensity");
  colnames(sample.130c) <- c(chosen.columns, "Intensity");
  colnames(sample.131) <- c(chosen.columns, "Intensity");
  
  
  Intensity <- c(sum(sample.126[,c("Intensity")]), 
                 sum(sample.127n[,c("Intensity")]), sum(sample.127c[,c("Intensity")]),
                 sum(sample.128n[,c("Intensity")]), sum(sample.128c[,c("Intensity")]),
                 sum(sample.129n[,c("Intensity")]), sum(sample.129c[,c("Intensity")]),
                 sum(sample.130n[,c("Intensity")]), sum(sample.130c[,c("Intensity")]),
                 sum(sample.131[,c("Intensity")]));
  
  Filename <- c("peptide_126.txt", 
                "peptide_127n.txt", "peptide_127c.txt",
                "peptide_128n.txt", "peptide_128c.txt",
                "peptide_129n.txt", "peptide_129c.txt",
                "peptide_130n.txt", "peptide_130c.txt",
                "peptide_131.txt");
  
  hybridization.file <- cbind(Assay, Sample, Intensity, Filename);
  
  cat("Writing Array file\n");
  write.table(array.file, "peptide_array.txt", sep="\t", row.names = FALSE, quote = FALSE);
  
  cat("Writing Sample 126 file\n");
  write.table(sample.126, "peptide_126.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 127n file\n");
  write.table(sample.127n, "peptide_127n.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 127c file\n");
  write.table(sample.127c, "peptide_127c.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 128n file\n");
  write.table(sample.128n, "peptide_128n.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 128c file\n");
  write.table(sample.128c, "peptide_128c.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 129n file\n");
  write.table(sample.129n, "peptide_129n.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 129c file\n");
  write.table(sample.129c, "peptide_129c.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 130n file\n");
  write.table(sample.130n, "peptide_130n.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 130c file\n");
  write.table(sample.130c, "peptide_130c.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 131 file\n");
  write.table(sample.131, "peptide_131.txt", sep="\t", row.names = FALSE, quote = FALSE);

  cat("Writing hybridization file\n");
  write.table(hybridization.file, "peptide_hybridization.txt", sep="\t", row.names = FALSE, quote = FALSE);
}

create.phosphosites.ArrayTrack.files <- function(phosphosites, Assay, Sample) {
  #need to load xlsx library
  library(xlsx);
  
  chosen.columns <- c("id", "Protein", "Fasta.headers", "Localization.prob",
                      "PEP", "Score", "Charge", "EntrezGeneID", "SWISSPROT_ENTRY_NAME",
                      "SWISSPROT_ACC_NUMBER", "DESCRIPTION", "GENENAME", "FUNCTIONS",
                      "SIMILARITY", "IMAGEID");
  
  array.file <- phosphosites[, c(chosen.columns, "Reporter.intensity.0",
                                 "Reporter.intensity.1", "Reporter.intensity.2",
                                 "Reporter.intensity.3", "Reporter.intensity.4",
                                 "Reporter.intensity.5", "Reporter.intensity.6",
                                 "Reporter.intensity.7", "Reporter.intensity.8",
                                 "Reporter.intensity.9")];
  
  #   #simplify the protein IDs: this may change depending on the database you use
  #   first.Protein.IDs <- get.first.IDs(phosphosites[, c("Protein.IDs")]);
  #   phosphosites <- cbind(first.Protein.IDs, phosphosites);
  #   array.file <- cbind(first.Protein.IDs, array.file);
  
  #create sample files: one file per sample!
  sample.126 <- phosphosites[, c(chosen.columns, "Reporter.intensity.0")];
  sample.127n <- phosphosites[, c(chosen.columns, "Reporter.intensity.1")];
  sample.127c <- phosphosites[, c(chosen.columns, "Reporter.intensity.2")];
  sample.128n <- phosphosites[, c(chosen.columns, "Reporter.intensity.3")];
  sample.128c <- phosphosites[, c(chosen.columns, "Reporter.intensity.4")];
  sample.129n <- phosphosites[, c(chosen.columns, "Reporter.intensity.5")];
  sample.129c <- phosphosites[, c(chosen.columns, "Reporter.intensity.6")];
  sample.130n <- phosphosites[, c(chosen.columns, "Reporter.intensity.7")];
  sample.130c <- phosphosites[, c(chosen.columns, "Reporter.intensity.8")];
  sample.131 <- phosphosites[, c(chosen.columns, "Reporter.intensity.9")];
  
  colnames(sample.126) <- c(chosen.columns, "Intensity");
  colnames(sample.127n) <- c(chosen.columns, "Intensity");
  colnames(sample.127c) <- c(chosen.columns, "Intensity");
  colnames(sample.128n) <- c(chosen.columns, "Intensity");
  colnames(sample.128c) <- c(chosen.columns, "Intensity");
  colnames(sample.129n) <- c(chosen.columns, "Intensity");
  colnames(sample.129c) <- c(chosen.columns, "Intensity");
  colnames(sample.130n) <- c(chosen.columns, "Intensity");
  colnames(sample.130c) <- c(chosen.columns, "Intensity");
  colnames(sample.131) <- c(chosen.columns, "Intensity");
  
  
  Intensity <- c(sum(sample.126[,c("Intensity")]), 
                 sum(sample.127n[,c("Intensity")]), sum(sample.127c[,c("Intensity")]),
                 sum(sample.128n[,c("Intensity")]), sum(sample.128c[,c("Intensity")]),
                 sum(sample.129n[,c("Intensity")]), sum(sample.129c[,c("Intensity")]),
                 sum(sample.130n[,c("Intensity")]), sum(sample.130c[,c("Intensity")]),
                 sum(sample.131[,c("Intensity")]));
  
  Filename <- c("phospho_126.txt", 
                "phospho_127n.txt", "phospho_127c.txt",
                "phospho_128n.txt", "phospho_128c.txt",
                "phospho_129n.txt", "phospho_129c.txt",
                "phospho_130n.txt", "phospho_130c.txt",
                "phospho_131.txt");
  
  hybridization.file <- cbind(Assay, Sample, Intensity, Filename);
  
  cat("Writing Array file\n");
  write.table(array.file, "phospho_array.txt", sep="\t", row.names = FALSE, quote = FALSE);
    
  cat("Writing Sample 126 file\n");
  write.table(sample.126, "phospho_126.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 127n file\n");
  write.table(sample.127n, "phospho_127n.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 127c file\n");
  write.table(sample.127c, "phospho_127c.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 128n file\n");
  write.table(sample.128n, "phospho_128n.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 128c file\n");
  write.table(sample.128c, "phospho_128c.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 129n file\n");
  write.table(sample.129n, "phospho_129n.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 129c file\n");
  write.table(sample.129c, "phospho_129c.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 130n file\n");
  write.table(sample.130n, "phospho_130n.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 130c file\n");
  write.table(sample.130c, "phospho_130c.txt", sep="\t", row.names = FALSE, quote = FALSE);
  cat("Writing Sample 131 file\n");
  write.table(sample.131, "phospho_131.txt", sep="\t", row.names = FALSE, quote = FALSE);
  
  
  cat("Writing hybridization file\n");
  write.table(hybridization.file, "phospho_hybridization.txt", sep="\t", row.names = FALSE, quote = FALSE);
}

#take the fasta headers and remove the protein sequences
grep.fasta.headers <- function(fasta) {
  headers <- NULL;
  for (i in 1:length(fasta)) {
    cat(i/length(fasta), "\n");
    if (substr(fasta[i], 0, 1) == ">") {
      headers <- c(headers, fasta[i])
    }
  }
  headers
}

grep.uniprot <- function(headers) {
  uniprot <- NULL;
  protein.names <- NULL;
  for (i in 1:length(headers)) {
    cat(i/length(headers), "\n");
    uniprot <- c(uniprot, 
                 unlist(strsplit(headers[i], "|", fixed = TRUE))[2]);
    name.start <- unlist(gregexpr(" ", headers[i]))[1] + 1;
    name.end <- unlist(gregexpr("OS=", headers[i]))[1] - 2;
    protein.name <- substr(headers[i], name.start, name.end);
    protein.names <- c(protein.names, protein.name);
  }
  uniprot <- cbind(uniprot, protein.names);
  uniprot
}

grep.Chlorocebus.sabeus.fasta <- function(headers) {
  uniprot <- NULL;
  protein.names <- NULL;
  for (i in 1:length(headers)) {
    cat(i/length(headers), "\n");
    each.uniprot <- unlist(strsplit(headers[i], " ", fixed = TRUE))[1];
    each.uniprot <- substr(each.uniprot, 2, nchar(each.uniprot));
    uniprot <- c(uniprot, each.uniprot);
    name.start <- unlist(gregexpr(": ", headers[i]))[1] + 2;
    name.end <- unlist(gregexpr(" [", headers[i], fixed = TRUE))[1] - 1;
    protein.name <- substr(headers[i], name.start, name.end);
    protein.names <- c(protein.names, protein.name);
  }
  uniprot <- cbind(uniprot, protein.names);
  uniprot
}

remove.characters.after.aString <- function(genelist, aString, col.name = "protein.names", adjust = 2) {
  short.protein.names <- NULL;
  for (i in 1:nrow(genelist)) {
    each.protein <- toString(unlist(genelist[i, col.name]));
    if (grepl(aString, each.protein)) {
      end.index <- unlist(regexpr(aString, each.protein)[1]);
      each.protein <- substr(each.protein, 1, end.index - adjust);
      short.protein.names <- c(short.protein.names, each.protein)
    }
    else {
      short.protein.names <- c(short.protein.names, each.protein)
    }
  }
  short.protein.names
}

remove.partial.endof.character <- function(character.vector, aString) {
  protein.names <- NULL;
  for (i in 1:length(character.vector)) {
    each.protein <- character.vector[i];
    if (grepl(aString, each.protein, fixed = TRUE)) {
      end.index <- unlist(regexpr(aString, each.protein)[1]);
      each.protein <- substr(each.protein, 1, end.index - 1);
      protein.names <- c(protein.names, each.protein)
    }
    else {
      protein.names <- c(protein.names, each.protein)
    }
  }
  protein.names;
}

remove.lowquality.infrontof.character <- function(character.vector, aString) {
  protein.names <- NULL;
  for (i in 1:length(character.vector)) {
    each.protein <- character.vector[i];
    if (grepl(aString, each.protein, fixed = TRUE)) {
      start.index <- nchar(aString);
      end.index <- nchar(each.protein);
      each.protein <- substr(each.protein, start.index + 1, end.index);
      protein.names <- c(protein.names, each.protein)
    }
    else {
      protein.names <- c(protein.names, each.protein)
    }
  }
  protein.names;
}

capitalize.first.letter.string <- function(aString) {
  if (is.character(substr(aString, 1, 1))) {
    aString <- paste(toupper(substr(aString, 1, 1)), substr(aString, 2, nchar(aString)), sep = "");
  }
  aString
}

match.swissprot.acc <- function(genelist, human.library.arraytrack) {
  protein.index <- NULL;
  for (i in 1: nrow(genelist)) {
    cat(i, "\n");
    protein.name <- toString(unlist(genelist[i, c("short.protein.names")]));
    protein.name <- paste(protein.name, " ", sep = "");
    count.match <- 0
    for (n in 1:nrow(human.library.arraytrack)) {
#       cat(n, "\n");
      protein.description <- toString(unlist(human.library.arraytrack[n, c("DESCRIPTION")]));
      protein.description <- paste(protein.description, " ", sep = "");
      if (grepl(protein.name, protein.description, fixed = TRUE) && substr(protein.name, 1, 1) == substr(protein.description, 1, 1)) {
        m <- c(i, n);
        protein.index <- rbind(protein.index, m);
        count.match <- count.match + 1;
      }
      if (n == nrow(human.library.arraytrack) && count.match == 0) {
        m <- c(i, NA);
        protein.index <- rbind(protein.index, m);
      }
    }
  }
  protein.index;
}

get.dotplot.data.frame <- function(annotated.proteins, gene.name, data.input.colnames, unique.id,
                                   gene.group = c("control", "control", "control", "diabetic", "diabetic",
                                                  "diabetic", "diabetic", "control", "control", "diabetic")) {
  dotplot.data.frame <- NULL;
  cat("Total number of gene is", length(gene.name), "\n");
  for (i in 1:length(gene.name)) {
    cat(i, "\n");
    gene.protein <- annotated.proteins[which(annotated.proteins[, unique.id] == gene.name[i]), ];
    gene.intensity <- t(gene.protein[, data.input.colnames]);
#     gene.group <- c("control", "control", "control", "diabetic", "diabetic",
#                     "diabetic", "diabetic", "control", "control", "diabetic");
    genename <- rep(gene.name[i], 10);
    df <- data.frame(gene.intensity, gene.group, genename);
    colnames(df) <- c("gene.intensity", "gene.group", "gene.name");
    dotplot.data.frame <- rbind(dotplot.data.frame, df);
  }
  dotplot.data.frame
}

draw.protein.dotplot <- function(dotplot.data.frame, genename, unique.id) {
  library(devtools);
  library(easyGgplot2);
  for (i in 1:length(genename)) {
    cat(i, "\n");
#   dotplot.data.frame <- get.dotplot.data.frame(annotated.proteins, genename);
    file.name <- paste(genename[i], "pro.png", sep = "_");
    df <- dotplot.data.frame[which(dotplot.data.frame[, unique.id] == genename[i]), ];
#   par(mfrow=c(1,1), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
#       col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
#       mgp = c(3, 1, 0));

    #geometric mean
    log2.intensity <- log2(df$gene.intensity);
    log2.mean.by.group <- tapply(log2.intensity, as.factor(df$gene.group), mean);
    diabetic.control.ratio <- 2^log2.mean.by.group["diabetic"]/2^log2.mean.by.group["control"];
    diabetic.control.ratio <- format(round(diabetic.control.ratio, 2), nsmall = 2);
    diabetic.control.ratio <- paste("(", diabetic.control.ratio, sep = "");
    diabetic.control.ratio <- paste(diabetic.control.ratio, ")", sep = "");

    title.txt <- paste(genename[i], diabetic.control.ratio, sep = " ");
    png(filename = file.name, width = 250, height = 300);
    p <- ggplot2.dotplot(data=df, xName='gene.group',yName='gene.intensity', 
                  groupName='gene.group', groupColors=c('#999999', '#56B4E9'), 
                  showLegend=FALSE, backgroundColor="white", ytitle="Normalized intensity", 
                  mainTitle=title.txt, xShowTitle=FALSE,
                  xTickLabelFont=c(14,"bold", "#003399"), ytitleFont=c(14,"bold", "#003399"),
                  yTickLabelFont=c(12,"bold", "#003399"), mainTitleFont=c(16,"bold", "#003399"),
                  dotsize=1.5,addMean=TRUE, meanPointShape=23, meanPointSize=4,
                  meanPointColor="black", meanPointFill="firebrick2");
    print(p);
    dev.off();
  }
}

draw.phosphosites.dotplot <- function(dotplot.data.frame, gene.names, unique.id) {
  library(devtools);
  library(easyGgplot2);
  
  for (i in 1:length(gene.names)) {
    cat(i, "\n");
    genename <- gene.names[i];
    df <- dotplot.data.frame[which(dotplot.data.frame$gene.name == genename), ];
    #geometric mean
    log2.intensity <- log2(df$gene.intensity);
    log2.mean.by.group <- tapply(log2.intensity, as.factor(df$gene.group), mean);
    diabetic.control.ratio <- 2^log2.mean.by.group["diabetic"]/2^log2.mean.by.group["control"];
    diabetic.control.ratio <- format(round(diabetic.control.ratio, 2), nsmall = 2);
    diabetic.control.ratio <- paste("(", diabetic.control.ratio, sep = "");
    diabetic.control.ratio <- paste(diabetic.control.ratio, ")", sep = "");
    
    #arithmetic mean
    #     mean.by.group <- tapply(df$gene.intensity, as.factor(df$gene.group), mean);
    #     diabetic.control.ratio <- mean.by.group["diabetic"]/mean.by.group["control"];
    #     diabetic.control.ratio <- format(round(diabetic.control.ratio, 2), nsmall = 2);
    #     diabetic.control.ratio <- paste("(", diabetic.control.ratio, sep = "");
    #     diabetic.control.ratio <- paste(diabetic.control.ratio, ")", sep = "");
    
    file.name <- paste(genename, "phospho.png", sep = "_");    
    png(filename = file.name, width = 250, height = 300);
    title.txt <- paste(genename, diabetic.control.ratio, sep = " ");
    p <- ggplot2.dotplot(data=df, xName='gene.group',yName='gene.intensity', 
                         groupName='gene.group', groupColors=c('#999999', '#56B4E9'), 
                         showLegend=FALSE, backgroundColor="white", ytitle="Normalized intensity", 
                         mainTitle=title.txt, xShowTitle=FALSE,
                         xTickLabelFont=c(14,"bold", "#003399"), ytitleFont=c(14,"bold", "#003399"),
                         yTickLabelFont=c(12,"bold", "#003399"), mainTitleFont=c(16,"bold", "#003399"),
                         dotsize=1.5,
                         addMean=TRUE, meanPointShape=23, meanPointSize=4,
                         meanPointColor="black", meanPointFill="firebrick2");
    print(p);
    dev.off()  
  }
}

ggplot2.dotplot.protein <- function(gene.names) {
  #1. read proteins data; this step removes all zero values. need to double check on this step
  # phosphosites <- TMT.phosphosites.cleanup.10plex(10000);
  cat("Reading annotated protein file\n");
  annotated.proteins <- read.table("protein_sam_nofilter_annotated_normalized_11062015_final.txt", header = TRUE, sep="\t", 
                                   quote="\"", fill = TRUE, comment.char="");
  
  ######################tell the program where to look for the intensity data#############
  data.input.colnames <- c("Normalized.0", "Normalized.1", "Normalized.2", 
                           "Normalized.3", "Normalized.4", "Normalized.5", 
                           "Normalized.6", "Normalized.7", "Normalized.8", 
                           "Normalized.9")
  
  ######################generate the data.frame for the plot##############################
  cat("Generating dataframe for the plot\n");
  dotplot.data.frame <- get.dotplot.data.frame(annotated.proteins, gene.names, data.input.colnames, unique.id = "GENENAME");
  
  ######################draw the individual plot##########################################
  cat("Plotting the dotplot\n");
  draw.protein.dotplot(dotplot.data.frame, gene.names, unique.id = "gene.name");  
}

ggplot2.dotplot.phosphosites <- function(gene.names) {
  #####################getting the raw intensity and phosphosite id#######################
  phosphosites.id <- TMT.phosphosites.cleanup.10plex(10000)[, c("id", "Protein.group.IDs")];
  colnames(phosphosites.id) <- c("phospho.id", "Protein.group.IDs");
  annotated.phosphosites <- read.table("phospho_sam_nofilter_annotated_11062015_final.txt", header = TRUE, sep="\t", 
                                       quote="\"", fill = TRUE, comment.char="");
  annotated.phosphosites <- merge(annotated.phosphosites, phosphosites.id, by = "phospho.id");
  genename.phosphoid <- paste(annotated.phosphosites$GENENAME, annotated.phosphosites$phospho.id, sep = "_");
  annotated.phosphosites <- cbind(annotated.phosphosites, genename.phosphoid);
  
  raw.intensity <- annotated.phosphosites[, c("Reporter.intensity.0", "Reporter.intensity.1",
                                              "Reporter.intensity.2", "Reporter.intensity.3",
                                              "Reporter.intensity.4", "Reporter.intensity.5",
                                              "Reporter.intensity.6", "Reporter.intensity.7",
                                              "Reporter.intensity.8", "Reporter.intensity.9")];
  #######################calculating the normalized intensity###############################
  cat("Calculating normalized phosphosite intensity\n");
  raw.median.intensity <- apply(raw.intensity, 2, median);
  raw.median.intensity.all <- median(unlist(raw.intensity));
  normalized.intensity <- sweep(raw.intensity*raw.median.intensity.all, MARGIN = 2, 
                                STATS = raw.median.intensity, FUN = "/");
  #   normalized.intensity <- cbind(normalized.intensity, diabetic.control.ratio);
  normalized.intensity <- cbind(normalized.intensity, annotated.phosphosites$GENENAME);
  normalized.intensity <- cbind(normalized.intensity, genename.phosphoid);
  
  #######################prepare the data table for a certain gene name##################
  cat("Preparing data frame for the plot\n");
  data.input.colnames <- c("Reporter.intensity.0", "Reporter.intensity.1", "Reporter.intensity.2", 
                           "Reporter.intensity.3", "Reporter.intensity.4", "Reporter.intensity.5",
                           "Reporter.intensity.6", "Reporter.intensity.7", "Reporter.intensity.8",
                           "Reporter.intensity.9")
  colnames(normalized.intensity)[11] <- "GENENAME";
  
  for (j in 1:length(gene.names)) {
    
    genename <- gene.names[j];
    cat("###############################################################################", genename, "\n\n");
    cat("Reading annotated phosphosites file\n");
    
    gene.phospho <- normalized.intensity[which(normalized.intensity$GENENAME == genename), ];
    phospho.names <- as.vector(unlist(gene.phospho$genename.phosphoid));
    dotplot.data.frame <- get.dotplot.data.frame(gene.phospho, phospho.names, data.input.colnames, unique.id = "genename.phosphoid");
    
    cat(gene.names, "\n");
    cat("Plotting the dotplot\n");
    draw.phosphosites.dotplot(dotplot.data.frame, phospho.names, unique.id = "gene.name");
    
    ######################plot the protein information#####################################
    ggplot2.dotplot.protein(genename);
  }
}

calculate.di.sam <- function(group.1, group.2, sorting = FALSE) {
  group1.col.mean <- apply(group.1, 1, mean);
  group2.col.mean <- apply(group.2, 1, mean);
  
  group1.minus.mean <- sweep(group.1, MARGIN = 1, group1.col.mean, "-");
  group1.minus.mean.square <- group1.minus.mean^2;
  group1.minus.mean.square.sum <- apply(group1.minus.mean.square, 1, sum);
  
  group2.minus.mean <- sweep(group.2, MARGIN = 1, group2.col.mean, "-");
  group2.minus.mean.square <- group2.minus.mean^2;
  group2.minus.mean.square.sum <- apply(group2.minus.mean.square, 1, sum);
  
  n1 <- ncol(group.1);
  n2 <- ncol(group.2);
  s.0 <- 0;
  s.factor <- (1/n1 + 1/n2)/(n1+n2-2);
  s.i <- sqrt(s.factor * (group1.minus.mean.square.sum + group2.minus.mean.square.sum));
  d.i <- (group1.col.mean - group2.col.mean)/(s.i + s.0);
  
  if (sorting) {
    d.i <- sort(d.i, decreasing = TRUE);
  }
    
  d.i
}

sam.analysis.v0 <- function(dataf = "proteins", n = 500000, intensity.column.names,
                            groups.id = 1:6, true.group = 1:3, n.group = 2,
                            filename.prefix = "group1vsgroup2") {
#     dataf = "proteins";
#     n = 500000;
#     groups.id = 1:6
#     n.group <- 2;
#     true.group <- 1:3;
#     filename.prefix = "group1vsgroup2"
  
  if (dataf == "proteins") {
    proteins <- read.table.proteinGroup(n);
    proteins <- filter.proteinGroup.TMT(proteins);
  }
  
  protein.id <- proteins[, "id"];
  raw.intensity <- proteins[, intensity.column.names[groups.id]];
  # raw.intensity[raw.intensity == 0] <- NA;
  
#   ######knn imputation 1############################################
#   if(exists(".Random.seed")) rm(.Random.seed)
#   raw.intensity.imputed <- kNN(as.matrix(raw.intensity));
#   raw.intensity.imputed[raw.intensity.imputed == 0] <- NA;
#   row_sub.imputed <- apply(raw.intensity.imputed, 1, function(row) all(!is.na(row)));
#   raw.intensity.imputed <- raw.intensity[row_sub.imputed, ];
#   protein.id.imputed <- protein.id[row_sub.imputed];
  
  #####removing rows containing zero values#######################
  row_sub <- apply(raw.intensity, 1, function(row) all(row != 0));
  # row_sub <- apply(raw.intensity, 1, function(row) all(!is.na(row)));
  raw.intensity <- raw.intensity[row_sub, ];
  protein.id <- protein.id[row_sub];
  
  #####data normalization and log conversion#####################
  raw.median.intensity <- apply(raw.intensity, 2, median);
  raw.median.intensity.all <- median(unlist(raw.intensity));
  normalized.intensity <- sweep(raw.intensity*raw.median.intensity.all, MARGIN = 2, 
                                STATS = raw.median.intensity, FUN = "/");
  log2.normalized.intensity <- log2(normalized.intensity);
  
  #####generate all possible permutation pattern
  n.sample <- ncol(log2.normalized.intensity);
  combn.groupings <- combn(1:n.sample, n.sample/n.group);
#   combn.groupings <- combn.groupings[, 1: (ncol(combn.groupings)/2)];
  
  d.is <- NULL;
  
  for (i in 1:ncol(combn.groupings)) {
    cat("i = ", i, "\n");
    group.1 <- log2.normalized.intensity[, combn.groupings[, i]];
    group.2 <- log2.normalized.intensity[, setdiff(1:n.sample, combn.groupings[, i])];  
    d.i <- calculate.di.sam(group.1, group.2, sorting = TRUE);
    d.is <- cbind(d.is, d.i);
  }
  
  d.i.e <- apply(d.is, 1, mean);
  d.is.delta <- sweep(d.is, MARGIN = 1, STATS = d.i.e, FUN = "-");
 
  group.1 <- log2.normalized.intensity[, true.group];
  group.2 <- log2.normalized.intensity[, setdiff(1:n.sample, true.group)];
  mean.group.1 <- apply(group.1, 1, mean);
  mean.group.2 <- apply(group.2, 1, mean);
  group1.over.group2 <- 2^mean.group.1/2^mean.group.2;
  abs.ratio <- sapply(group1.over.group2, absolute.fold.change);
  welth.t.test.p.value <- apply(cbind(group.1, group.2), 1, report.welch.t.test.p.value);
  
  #please note this function is using the sorting is false
  d.i.o <- calculate.di.sam(group.1, group.2, sorting = FALSE);
  
  #########################Plotting d.i.o vs d.i.e############################################
  file.name <- paste(filename.prefix, "dio_vs_die.png", sep = "_");
  png(filename = file.name, width = 500, height = 500);
  par(mfrow=c(1,1), cex.lab = 1.5, cex.axis = 1.5, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  plot(d.i.e, sort(d.i.o, decreasing = TRUE),
       xlab = "Expected d(i) values", ylab = "Observed d(i) values", pch = 1, 
       col = "dodgerblue4", cex = 0.6);
  abline(a = 0, b = 1);
  
  dev.off()
  
  #####################calculating delta and adding more information#################
  d.i.o.id <- cbind(d.i.o, protein.id);
  colnames(d.i.o.id) <- c("Observed.di", "Protein.id");
  d.i.o.id <- cbind(d.i.o.id, group1.over.group2);
  d.i.o.id <- cbind(d.i.o.id, abs.ratio);
  d.i.o.id <- cbind(d.i.o.id, welth.t.test.p.value);
  d.i.o.id <- cbind(d.i.o.id, group.1);
  d.i.o.id <- cbind(d.i.o.id, group.2);
  d.i.o.id <- d.i.o.id[order(-d.i.o.id[, "Observed.di"]), ];
  diff.dio.vs.die <- d.i.o.id[, "Observed.di"] - d.i.e;
  d.i.o.id <- cbind(d.i.o.id, d.i.e);
  d.i.o.id <- cbind(d.i.o.id, diff.dio.vs.die);
  abs.diff.dio.vs.die <- abs(diff.dio.vs.die)
  d.i.o.id <- cbind(d.i.o.id, abs.diff.dio.vs.die);
  
  protein.id.annotation <- proteins[, c("id", "Protein.IDs", "Protein.names", "Gene.names", "Fasta.headers")];
  colnames(protein.id.annotation) <- c("Protein.id", "Protein.IDs", "Protein.names", "Gene.names", "Fasta.headers");
  d.i.o.id <- merge(d.i.o.id, protein.id.annotation, by = "Protein.id");
  d.i.o.id <- d.i.o.id[order(-d.i.o.id[, "abs.diff.dio.vs.die"]), ];
  
  ################estimate FDR#####################################################
  fdr.sam <- sapply(d.i.o.id[, "abs.diff.dio.vs.die"], estimate.fdr.sam.test, 
                    stat.table = d.is.delta);
  d.i.o.id <- cbind(d.i.o.id, fdr.sam);

  file.name <- paste(filename.prefix, ".txt", sep = "");
  write.table(d.i.o.id, file.name, sep="\t", row.names = FALSE, quote = FALSE);
#   file.name <- paste(filename.prefix, "di_permutations.txt", sep = "_");
#   write.table(d.is, file.name, sep="\t", row.names = FALSE, quote = FALSE);
}

absolute.fold.change <- function (a.ratio) {
  if (a.ratio < 1) {
    abs.ratio <- 1/a.ratio;
  }
  else if (a.ratio >= 1) {
    abs.ratio <- a.ratio;
  }
  else {
    cat("Error! Input ratio cannot be negative or the input ratio is not numeric\n")
  }
  
  abs.ratio;
}

report.welch.t.test.p.value <- function(input, n = length(input)/2) {
  n.length <- length(input);
  x <- input[1:n];
  y <- input[(n+1):n.length];
  t.test(x, y)$p.value
}

estimate.fdr.sam.test <- function(stat.table, posi.delta.cutoff, 
                                  nega.delta.cutoff = -posi.delta.cutoff) {
  if (posi.delta.cutoff > 0 & nega.delta.cutoff <0) {
    significant.count <- apply(stat.table, 2, count.above.and.below, 
                               posi.delta.cutoff, nega.delta.cutoff);
    median(significant.count)/nrow(stat.table);
  }
  else {
    stop("Please check positive and negative delta cutoff\n")
  }
}

count.above.and.below <- function(x, above.s, below.t = - above.s) {
  sum(x >= above.s | x <= below.t);
}


#http://rstudio-pubs-static.s3.amazonaws.com/13301_6641d73cfac741a59c0a851feb99e98b.html
#https://cran.r-project.org/web/packages/VennDiagram/VennDiagram.pdf
my.draw.quad.venn <- function(df, choose.col = c("mch.mut1.sig", "mch.mut2.sig", "mch.mut3.sig", "mch.wt.sig"),
                              mark.sign = "+", 
                              set.names = c("mch.mut1", "mch.mut2", "mch.mut3", "mch.wt"), 
                              fill.colors = c("skyblue", "pink1", "mediumorchid", "orange"),
                              file.name = "test.png") {
  # require(lima)
  require(plyr)
  # require(lima)
  require(venneuler)
  require(knitr)
  require(VennDiagram)
  require(limma)
#   biocLite("statmod")
#   biocLite("limma")
  
  #   choose.col = c("mch.mut1.sig", "mch.mut2.sig", "mch.mut3.sig", "mch.wt.sig")
  #   mark.sign = "+"
  #   set.names = c("mch.mut1", "mch.mut2", "mch.mut3", "mch.wt")
  #   fill.colors = c("orange", "red", "green", "blue")
  #   file.name = "test.png"
  
  group1 <- test[, choose.col[1]];
  group1.positive <- group1 == mark.sign;
  
  group2 <- test[, choose.col[2]];
  group2.positive <- group2 == mark.sign;
  
  group3 <- test[, choose.col[3]];
  group3.positive <- group3 == mark.sign;
  
  group4 <- test[, choose.col[4]];
  group4.positive <- group4 == mark.sign;
  
  g4 <- cbind(group1.positive, group2.positive, group3.positive, group4.positive);
  
  a <- vennCounts(g4);
  a <- a[, "Counts"];
  test.area1 <- sum(a[9:16])
  test.area2 <- sum(a[c(5:8, 13:16)])
  test.area3 <- sum(a[c(3, 4, 7, 8, 11, 12, 15, 16)])
  test.area4 <- sum(a[c(2, 4, 6, 8, 10, 12, 14, 16)])
  test.n12 <- sum(a[c(13:16)]);
  test.n13 <- sum(a[c(11, 12, 15, 16)]);
  test.n14 <- sum(a[c(10, 12, 14, 16)]);
  test.n23 <- sum(a[c(7, 8, 15, 16)]);
  test.n24 <- sum(a[c(6, 8, 14, 16)]);
  test.n34 <- sum(a[c(4, 8, 12, 16)]);
  test.n123 <- sum(a[c(15, 16)]);
  test.n124 <- sum(a[c(14, 16)]);
  test.n134 <- sum(a[c(12, 16)]);
  test.n234 <- sum(a[c(8, 16)]);
  test.n1234 <- sum(a[16]);
  
  
  png(filename = file.name, width = 500, height = 500);
  venn.plot <- draw.quad.venn(test.area1, test.area2, test.area3, test.area4, 
                              test.n12, test.n13, test.n14, test.n23,
                              test.n24, test.n34, test.n123, test.n124,
                              test.n134, test.n234, test.n1234, 
                              category = set.names,
                              fill = fill.colors);
  grid.draw(venn.plot);
  dev.off();
}

ggplot2.dotplot.protein.general <- function(df, value.colnames, gene.names, unique.id, group.name.repeat) {
  #1. read proteins data; this step removes all zero values. need to double check on this step
  # phosphosites <- TMT.phosphosites.cleanup.10plex(10000);
  # df <- beforeskyline;
  
  
  cat("Reading annotated protein file\n");
  annotated.proteins <- df;
  gene.group <- c(rep(group.name.repeat[1], group.name.repeat[2]), 
                  rep(group.name.repeat[3], group.name.repeat[4]));
  
  ######################tell the program where to look for the intensity data#############
  data.input.colnames <- value.colnames;
  
  ######################generate the data.frame for the plot##############################
  cat("Generating dataframe for the plot\n");
  dotplot.data.frame <- get.dotplot.data.frame.general(annotated.proteins, gene.names, data.input.colnames, unique.id, gene.group);
  
  ######################draw the individual plot##########################################
  cat("Plotting the dotplot\n");
  draw.protein.dotplot.general(dotplot.data.frame, gene.names, unique.id, group.name.repeat);  
}

get.dotplot.data.frame.general <- function(annotated.proteins, gene.name, data.input.colnames, unique.id, gene.group) {
  dotplot.data.frame <- NULL;
  cat("Total number of gene is", length(gene.name), "\n");
  for (i in 1:length(gene.name)) {
    cat(i, "\n");
    gene.protein <- annotated.proteins[which(annotated.proteins[, unique.id] == gene.name[i]), ];
    gene.intensity <- t(gene.protein[, data.input.colnames]);
    genename <- rep(gene.name[i], length(gene.group));
    df <- data.frame(gene.intensity, gene.group, genename);
    colnames(df) <- c("gene.intensity", "gene.group", unique.id);
    dotplot.data.frame <- rbind(dotplot.data.frame, df);
  }
  dotplot.data.frame
}

draw.protein.dotplot.general <- function(dotplot.data.frame, genename, unique.id, group.name.repeat) {
  library(devtools);
  library(easyGgplot2);
  for (i in 1:length(genename)) {
    cat(i, "\n");
    #   dotplot.data.frame <- get.dotplot.data.frame(annotated.proteins, genename);
    file.name <- paste(genename[i], "pro.png", sep = "_");
    df <- dotplot.data.frame[which(dotplot.data.frame[, unique.id] == genename[i]), ];
    #   par(mfrow=c(1,1), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
    #       col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
    #       mgp = c(3, 1, 0));
    
    #geometric mean
    log2.intensity <- log2(df$gene.intensity);
    log2.mean.by.group <- tapply(log2.intensity, as.factor(df$gene.group), mean);
    diabetic.control.ratio <- 2^log2.mean.by.group[group.name.repeat[1]]/2^log2.mean.by.group[group.name.repeat[3]];
    diabetic.control.ratio <- format(round(diabetic.control.ratio, 2), nsmall = 2);
    diabetic.control.ratio <- paste("(", diabetic.control.ratio, sep = "");
    diabetic.control.ratio <- paste(diabetic.control.ratio, ")", sep = "");
    
    title.txt <- paste(genename[i], diabetic.control.ratio, sep = " ");
    png(filename = file.name, width = 250, height = 300);
    p <- ggplot2.dotplot(data=df, xName='gene.group',yName='gene.intensity', 
                         groupName='gene.group', groupColors=c('#999999', '#56B4E9'), 
                         showLegend=FALSE, backgroundColor="white", ytitle="Normalized intensity", 
                         mainTitle=title.txt, xShowTitle=FALSE,
                         xTickLabelFont=c(14,"bold", "#003399"), ytitleFont=c(14,"bold", "#003399"),
                         yTickLabelFont=c(12,"bold", "#003399"), mainTitleFont=c(16,"bold", "#003399"),
                         dotsize=1.5,addMean=TRUE, meanPointShape=23, meanPointSize=4,
                         meanPointColor="black", meanPointFill="firebrick2");
    print(p);
    dev.off();
  }
}

Qpro.from.precursorQ <- function(preQ) {
  #paste two column together with a sep sign
  preQ <- within(preQ, PepSeqSample <- paste(Peptide.Modified.Sequence, Replicate.Name, sep = "_"));
  #remove an entire column
  preQ$Precursor <- NULL;
  #Sum of rows based on column value
  #https://cran.r-project.org/web/packages/plyr/plyr.pdf
  # pepQ <- ddply(preQ, .(Peptide.Modified.Sequence, Replicate.Name), transform, Total.Area.MS1 = sum(Total.Area.MS1));
  # pepQ <- ddply(preQ, .(Peptide.Modified.Sequence, Replicate.Name), Total.Area.MS1 = sum(Total.Area.MS1));
  pepQ_sum <- ddply(preQ, .(PepSeqSample), summarize, Total.Area = sum(Total.Area.MS1));
  preQ$Total.Area.MS1 <- NULL;
  preQ <- unique(preQ);
  pepQ <- merge(preQ, pepQ_sum, by = "PepSeqSample");
  #calculate protein intensity
  pepQ <- within(pepQ, ProSample <- paste(Protein.Accession, Replicate.Name, sep = "_"));
  pepQ$PepSeqSample <- NULL;
  pepQ$Peptide.Modified.Sequence <- NULL;
  proQ_sum <- ddply(pepQ, .(ProSample), summarize, Total.Protein.Area = sum(Total.Area))
  pepQ$Total.Area <- NULL;
  pepQ <- unique(pepQ);
  proQ <- merge(pepQ, proQ_sum, by = "ProSample");
  #calculate student.t value and ratios
  # GroupNames <- rep(c(rep("Group1", 3), rep("Group2", 3)), 11);
  # proQ <- cbind(proQ, GroupNames);
  # proQ$ProSample <- NULL;
  # proQ$Replicate.Name <- NULL;
  # proQ_summary <- ddply(proQ, .(ProSample), summarize, Total.Protein.Area = sum(Total.Area))
  proIDs <- proQ[, c("Protein.Name", "Protein.Accession", "Protein.Gene")];
  proIDs <- unique(proIDs);
  proValues <- proQ[, c("Protein.Accession", "Replicate.Name", "Total.Protein.Area")];
  # proValues_melt <- melt(proValues, id = "Protein.Accession");
  #http://i2.wp.com/www.r-statistics.com/wp-content/uploads/2012/01/reshaping-data-using-melt-and-cast.png
  proValues <- cast(proValues, Protein.Accession~Replicate.Name);
  #calculate p-values and fold change
  
  group.1 <- log2(proValues[, 2:4]);
  group.2 <- log2(proValues[, 5:7]);
  mean.group.1 <- apply(group.1, 1, mean);
  mean.group.2 <- apply(group.2, 1, mean);
  group1.over.group2 <- 2^mean.group.1/2^mean.group.2;
  abs.ratio <- sapply(group1.over.group2, absolute.fold.change);
  welth.t.test.p.value <- apply(cbind(group.1, group.2), 1, report.welch.t.test.p.value);
  proValues <- cbind(proValues, group1.over.group2, abs.ratio, welth.t.test.p.value);
  proValues <- cbind(proIDs, proValues);
  write.csv(proValues, "ProteinQuantitation.csv");
  proValues
}

label.free.analysis <- function(proteins, s1, s2, filename, 
                                col.names = c("Protein.IDs", "Protein.names", "Gene.names",
                                               "Fasta.headers")
                                ) {
  raw.intensity <- proteins[, intensity.column.names[c(s1, s2)]];
  
  #remove all rows with zero values
  row_sub <- apply(raw.intensity, 1, function(row) all(row != 0));
  raw.intensity <- raw.intensity[row_sub, ];
  proteins <- proteins[row_sub, ];
  
  raw.median.intensity <- apply(raw.intensity, 2, median);
  raw.median.intensity.all <- median(unlist(raw.intensity));
  normalized.intensity <- sweep(raw.intensity*raw.median.intensity.all, MARGIN = 2, 
                                STATS = raw.median.intensity, FUN = "/");
  
  group.1 <- normalized.intensity[, 1:length(s1)];
  group.2 <- normalized.intensity[, (length(s1) + 1):(length(s1) + length(s2))];
  
  log2.normalized.intensity <- log2(normalized.intensity);
  log2.group.1 <- log2(group.1);
  log2.group.2 <- log2(group.2);
  log2.group.1.mean <- apply(log2.group.1, 1, mean);
  log2.group.2.mean <- apply(log2.group.2, 1, mean);
  
  group1.over.group2 <- 2^log2.group.1.mean/2^log2.group.2.mean;
  abs.ratio <- sapply(group1.over.group2, absolute.fold.change);
  welth.t.test.p.value <- apply(cbind(group.1, group.2), 1, report.welch.t.test.p.value);
  log.welth.p <- -log10(welth.t.test.p.value);
  log2.ratio <- log2(group1.over.group2);
  logp.multiply.log2.ratio <- log.welth.p * abs(log2.ratio);
  p.value.BH <- p.adjust(welth.t.test.p.value, method = "BH");
  
  ###########permutation-based FDR##########
  n.sample <- ncol(log2.normalized.intensity);
  combn.groupings <- combn(1:n.sample, length(s1));
  combn.groupings <- combn.groupings[, 2: (ncol(combn.groupings)/2)];
  welch.t.p.values <- NULL;
  for (i in 1:ncol(combn.groupings)) {
    cat("i = ", i, "\n");
    group.1 <- log2.normalized.intensity[, combn.groupings[, i]];
    group.2 <- log2.normalized.intensity[, setdiff(1:n.sample, combn.groupings[, i])];  
    welch.t.p.value <- apply(cbind(group.1, group.2), 1, report.welch.t.test.p.value);
    # welch.t.p.value <- sort(welch.t.p.value, decreasing = TRUE);
    welch.t.p.values <- cbind(welch.t.p.values, welch.t.p.value);
  }
  
  permutation.count.table <- NULL;
  for (i in 1:length(welth.t.test.p.value)) {
    cat("i = ", i, "\n");
    pvalue <- welth.t.test.p.value[i];
    each.pvalue.count <- apply(welch.t.p.values, 2, function(x) {sum(x <= pvalue)});
    permutation.count.table <- rbind(permutation.count.table, each.pvalue.count);
  }
  permutation.fdr <- apply(permutation.count.table, 1, mean);
  number.of.genes <- length(permutation.fdr);
  permutation.fdr <- permutation.fdr/number.of.genes;
  #############################################
  
  proIDs <- proteins[, col.names];
  
  normalized.intensity <- cbind(normalized.intensity, group1.over.group2, abs.ratio, 
                                welth.t.test.p.value, log.welth.p, p.value.BH, log2.ratio, 
                                logp.multiply.log2.ratio, permutation.fdr);
  proValues <- cbind(proIDs, normalized.intensity);
  write.csv(proValues, paste(filename, "ProteinsWelchttest.csv", sep="_"), row.names = FALSE);
  proValues
}

volcano.plot <- function(label.free.result, cutoff = 2, cutoff2 = FALSE, filename, 
                         colselect = "logp.multiply.log2.ratio",
                         colplotx = "log2.ratio", colploty = "log.welth.p", 
                         x_lab = "log2 ratio", y_lab = "log10 welth t p value", 
                         text_label = "Gene.names", is.title = TRUE) {
  proValues <- label.free.result;
  # significant.pro <- subset(proValues, logp.multiply.log2.ratio >= cutoff);
  significant.pro <- proValues[proValues[, colselect] >= cutoff, ]
  if (cutoff2 != FALSE) {
    significant.pro <- significant.pro[significant.pro[, "logFC"] >= cutoff2 | significant.pro[, "logFC"] <= -cutoff2, ]
  }
  
  filename <- paste(filename, toString(cutoff), sep = "_cf");
  png(filename = paste(filename, "volcano_plot.png", sep = "_"), 
      width = 500, height = 500);
  par(mfrow=c(1,1), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  plot(proValues[, colplotx], proValues[, colploty],
       xlab = x_lab, ylab = y_lab, 
       pch = 16, 
       col = "dodgerblue4", cex = 1);
  points(significant.pro[, colplotx], significant.pro[, colploty], pch = 16, 
         col = "red", cex = 1);
  abline(h = cutoff, col = "red", lty = 2)
  if (cutoff2 != FALSE) {
    abline(v = cutoff2, col = "red", lty = 2)
    abline(v = -cutoff2, col = "red", lty = 2)
  }
  if (text_label != FALSE) {
    text(significant.pro[, colplotx], significant.pro[, colploty], 
         labels = as.character(unlist(significant.pro[, text_label])), 
         cex = 0.8, col = "red", pos = 4);
  }
  # with(subset(proValues, log.welth.p > 2/abs(log2.ratio)), 
  #      points(log2.ratio, log.welth.p, pch = 16, 
  #             col = "red", cex = 0.6));
  if(is.title == TRUE) {
    title(main = filename, col = "dodgerblue4");
  }
  dev.off()
  significant.pro <- significant.pro[order(-significant.pro[, colselect]), ];
  write.csv(significant.pro, paste(filename, "significant.list.csv", sep="_"), row.names = FALSE);
  significant.pro
}

TMT.cleanup.10plex <- function(peptides) {
#   peptides <- read.table.peptides(nrow);
#   #removing reverse and potential contaminants
#   peptides <- filter.TMT.peptides(peptides);
  #calculate how many positive reporter ion intensities
  total.peptides <- nrow(peptides);
  reporter0.count <- nrow(peptides[which(peptides$Reporter.intensity.0 > 0),]);
  reporter1.count <- nrow(peptides[which(peptides$Reporter.intensity.1 > 0),]);
  reporter2.count <- nrow(peptides[which(peptides$Reporter.intensity.2 > 0),]);
  reporter3.count <- nrow(peptides[which(peptides$Reporter.intensity.3 > 0),]);
  reporter4.count <- nrow(peptides[which(peptides$Reporter.intensity.4 > 0),]);
  reporter5.count <- nrow(peptides[which(peptides$Reporter.intensity.5 > 0),]);
  reporter6.count <- nrow(peptides[which(peptides$Reporter.intensity.6 > 0),]);
  reporter7.count <- nrow(peptides[which(peptides$Reporter.intensity.7 > 0),]);
  reporter8.count <- nrow(peptides[which(peptides$Reporter.intensity.8 > 0),]);
  reporter9.count <- nrow(peptides[which(peptides$Reporter.intensity.9 > 0),]);
  
  cat("Total peptide count is", total.peptides, "\n");
  cat("Reporter 126 count is", reporter0.count, "\n");
  cat("Reporter 127N count is", reporter1.count, "\n");
  cat("Reporter 127C count is", reporter2.count, "\n");
  cat("Reporter 128N count is", reporter3.count, "\n");
  cat("Reporter 128C count is", reporter4.count, "\n");
  cat("Reporter 129N count is", reporter5.count, "\n");
  cat("Reporter 129C count is", reporter2.count, "\n");
  cat("Reporter 130N count is", reporter3.count, "\n");
  cat("Reporter 130C count is", reporter4.count, "\n");
  cat("Reporter 131 count is", reporter5.count, "\n");
  
  peptides <- peptides[which(peptides$Reporter.intensity.0 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.1 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.2 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.3 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.4 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.5 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.6 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.7 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.8 > 0),];
  peptides <- peptides[which(peptides$Reporter.intensity.9 > 0),];
  
  cat("Total non-zero count is", nrow(peptides), "\n");
  
#   peptide.count <- c(total.peptides, reporter0.count, reporter1.count, 
#                      reporter2.count, reporter3.count, reporter4.count,
#                      reporter5.count, reporter6.count, reporter7.count,
#                      reporter8.count, reporter9.count, nrow(peptides));
  
#   output <- list(peptides, peptide.count);
#   output;
  peptides;
}

cor.prob <- function(X, dfr = nrow(X) - 2) {
  R <- cor(X)
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr / (1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  
  cor.mat <- t(R)
  cor.mat[upper.tri(cor.mat)] <- NA
  cor.mat
}

#plot correlation of two variables with linear regression line in red;
#also label the R square and p values
#http://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
plot.correlation <- function(df, col1, col2, addnote, lab.note = " log10 intensity", 
                             diagonal = FALSE, log.2 = FALSE, smoothplot = FALSE) {
  filename <- paste(col1, col2, sep = ".");
  filename <- paste(filename, addnote, sep = ".");
  png(filename = paste(filename, ".png", sep = ""), width = 300, height = 300);
  par(mfrow=c(1,1), cex.lab = 1.2, cex.axis = 1.2, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  x <- df[, col1];
  y <- df[, col2];
  if (log.2) {
    x <- log2(x);
    y <- log2(y);
  }
  
  if (smoothplot) {
    smoothScatter(x, y,
                  xlab = paste(col1, lab.note, sep = ""), 
                  ylab = paste(col2, lab.note, sep = ""), 
                  pch = 16, col = "dodgerblue4", 
                  cex = 0.6)
  }
  else {
    plot(x, y, 
         xlab = paste(col1, lab.note, sep = ""), 
         ylab = paste(col2, lab.note, sep = ""), 
         pch = 16, col = "dodgerblue4", 
         cex = 0.6#, xlim = c(2, 7), ylim = c(2, 7)
    );
    
  }
  
  regl <- lm(y~x);
  cat("The linear regression model intercep is", regl[[1]][[1]], "\n")
  cat("The linear regression model slope is", regl[[1]][[2]], "\n")
  abline(regl, col="red");
  
  if (diagonal) {
    x.y.min <- min(min(x), min(y));
    x.y.max <- max(max(x), max(y));
    lines(x = c(x.y.min, x.y.max),
          y = c(x.y.min, x.y.max));
  }
  
  modsum <- summary(regl);
  r2 <- modsum$adj.r.squared;
  my.p = modsum$coefficients[2,4];
  rp = vector('expression',1)
  rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                     list(MYVALUE = format(r2,dig=3)))[2];
#   rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
#                      list(MYOTHERVALUE = format(my.p, digits = 5)))[2];
  legend('topleft', legend = rp, bty = 'n')
    
  dev.off();
}

pairs.plot <- function(df, addnote) {
  filename <- paste("PairedCorrelation", addnote, sep = ".");
  png(filename = paste(filename, ".png", sep = ""), width = 1500, height = 1500);
  par(mfrow=c(1,1), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  pairs(~a126+n127+c127+n128+c128+n129+c129+n130+c130+a131, 
        data=df,
        main=addnote);
  dev.off();
}

plot.ratio.intensity <- function(df, col1, col2, addnote, dynamic.x.range = FALSE) {
  filename <- paste(col1, col2, sep = ".");
  filename <- paste(filename, addnote, sep = ".");
  filename <- paste(filename, "RatioIntensity", sep = ".");
  sum.channel1 <- sum(df[, col1]);
  sum.channel2 <- sum(df[, col2]);
  ratios <- df[, col1]/df[, col2];
  normalized.ratios <- ratios; #* sum.channel2 / sum.channel1;
  log2.ratio <- log2(normalized.ratios);
  log10.intensity <- (log10(df[, col1]*df[, col2]))/2;
  # cat(log10.intensity);
  x.name <- paste(col1, col2, sep=".over.");
  png(filename = paste(filename, ".png", sep = ""), width = 300, height = 300);
  par(mfrow=c(1,1), cex.lab = 1.2, cex.axis = 1.2, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  if (dynamic.x.range) {
    plot(log2.ratio, log10.intensity, 
         xlab = paste(x.name, "Log2 Ratio", sep = " "), 
         ylab = "Log10 Intensity", pch = 16, 
         col = "dodgerblue4", cex = 0.6);
  }
  
  else {
    plot(log2.ratio, log10.intensity, 
         xlab = paste(x.name, "Log2 Ratio", sep = " "), 
         ylab = "Log10 Intensity", pch = 16, 
         col = "dodgerblue4", cex = 0.6, xlim = c(-4, 4));
  }
  
  dev.off();
}

plot.ratio.ratio.intensity <- function(df, inten.col1, inten.col2, inten.col3, inten.col4,
                                       ratio.col1, ratio.col2, addnote, dynamic.x.range = FALSE) {
  filename <- paste(inten.col1, inten.col2, sep = ".");
  # filename <- paste(filename, addnote, sep = ".");
  filename <- paste(filename, "R&I", sep = ".");
#   sum.channel1 <- sum(df[, col1]);
#   sum.channel2 <- sum(df[, col2]);
#   ratios <- df[, col1]/df[, col2];
#   normalized.ratios <- ratios; #* sum.channel2 / sum.channel1;
  log2.ratio <- log2(df[, ratio.col1]/df[, ratio.col2]);
  cat("The median ratio is", median(df[, ratio.col1]/df[, ratio.col2]), "/n/n");
  log10.intensity <- (log10(df[, inten.col1]*df[, inten.col2]*df[, inten.col3]*df[, inten.col4]))/4;
  # cat(log10.intensity);
  # x.name <- paste(col1, col2, sep=".over.");
  png(filename = paste(filename, ".png", sep = ""), width = 300, height = 300);
  par(mfrow=c(1,1), cex.lab = 1.2, cex.axis = 1.2, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  if (!dynamic.x.range) {
    plot(log2.ratio, log10.intensity, 
         xlab = addnote, 
         ylab = "Log10 Intensity", pch = 16, 
         col = "dodgerblue4", cex = 0.6);
  }
  
  else {
    plot(log2.ratio, log10.intensity, 
         xlab = addnote, 
         ylab = "Log10 Intensity", pch = 16, 
         col = "dodgerblue4", cex = 0.6, xlim = dynamic.x.range);
  }
  
  dev.off();
}


create.phosphosites.ArrayTrack.xlsxfiles <- function(phosphosites, Assay, Sample) {
  #need to load xlsx library
  library(xlsx);
  
  chosen.columns <- c("id", "Protein", "Fasta.headers", "Localization.prob",
                      "PEP", "Score", "Charge");
  
  array.file <- phosphosites[, c(chosen.columns, "Reporter.intensity.0",
                                 "Reporter.intensity.1", "Reporter.intensity.2",
                                 "Reporter.intensity.3", "Reporter.intensity.4",
                                 "Reporter.intensity.5", "Reporter.intensity.6",
                                 "Reporter.intensity.7", "Reporter.intensity.8",
                                 "Reporter.intensity.9")];
  
  #   #simplify the protein IDs: this may change depending on the database you use
  #   first.Protein.IDs <- get.first.IDs(phosphosites[, c("Protein.IDs")]);
  #   phosphosites <- cbind(first.Protein.IDs, phosphosites);
  #   array.file <- cbind(first.Protein.IDs, array.file);
  
  #create sample files: one file per sample!
  sample.126 <- phosphosites[, c(chosen.columns, "Reporter.intensity.0")];
  sample.127n <- phosphosites[, c(chosen.columns, "Reporter.intensity.1")];
  sample.127c <- phosphosites[, c(chosen.columns, "Reporter.intensity.2")];
  sample.128n <- phosphosites[, c(chosen.columns, "Reporter.intensity.3")];
  sample.128c <- phosphosites[, c(chosen.columns, "Reporter.intensity.4")];
  sample.129n <- phosphosites[, c(chosen.columns, "Reporter.intensity.5")];
  sample.129c <- phosphosites[, c(chosen.columns, "Reporter.intensity.6")];
  sample.130n <- phosphosites[, c(chosen.columns, "Reporter.intensity.7")];
  sample.130c <- phosphosites[, c(chosen.columns, "Reporter.intensity.8")];
  sample.131 <- phosphosites[, c(chosen.columns, "Reporter.intensity.9")];
  
  colnames(sample.126) <- c(chosen.columns, "Intensity");
  colnames(sample.127n) <- c(chosen.columns, "Intensity");
  colnames(sample.127c) <- c(chosen.columns, "Intensity");
  colnames(sample.128n) <- c(chosen.columns, "Intensity");
  colnames(sample.128c) <- c(chosen.columns, "Intensity");
  colnames(sample.129n) <- c(chosen.columns, "Intensity");
  colnames(sample.129c) <- c(chosen.columns, "Intensity");
  colnames(sample.130n) <- c(chosen.columns, "Intensity");
  colnames(sample.130c) <- c(chosen.columns, "Intensity");
  colnames(sample.131) <- c(chosen.columns, "Intensity");
  
  
  Intensity <- c(sum(sample.126[,c("Intensity")]), 
                 sum(sample.127n[,c("Intensity")]), sum(sample.127c[,c("Intensity")]),
                 sum(sample.128n[,c("Intensity")]), sum(sample.128c[,c("Intensity")]),
                 sum(sample.129n[,c("Intensity")]), sum(sample.129c[,c("Intensity")]),
                 sum(sample.130n[,c("Intensity")]), sum(sample.130c[,c("Intensity")]),
                 sum(sample.131[,c("Intensity")]));
  
  Filename <- c("phospho_126.xlsx", 
                "phospho_127n.xlsx", "phospho_127c.xlsx",
                "phospho_128n.xlsx", "phospho_128c.xlsx",
                "phospho_129n.xlsx", "phospho_129c.xlsx",
                "phospho_130n.xlsx", "phospho_130c.xlsx",
                "phospho_131.xlsx");
  
  hybridization.file <- cbind(Assay, Sample, Intensity, Filename);
  
  cat("Writing Array file\n");
  write.xlsx(array.file, "phospho_array.xlsx", row.names = FALSE);
  cat("Writing Sample 126 file\n");
  write.xlsx(sample.126, "phospho_126.xlsx", row.names = FALSE);
  cat("Writing Sample 127n file\n");
  write.xlsx(sample.127n, "phospho_127n.xlsx", row.names = FALSE);
  cat("Writing Sample 127c file\n");
  write.xlsx(sample.127c, "phospho_127c.xlsx", row.names = FALSE);
  cat("Writing Sample 128n file\n");
  write.xlsx(sample.128n, "phospho_128n.xlsx", row.names = FALSE);
  cat("Writing Sample 128c file\n");
  write.xlsx(sample.128c, "phospho_128c.xlsx", row.names = FALSE);
  cat("Writing Sample 129n file\n");
  write.xlsx(sample.129n, "phospho_129n.xlsx", row.names = FALSE);
  cat("Writing Sample 129c file\n");
  write.xlsx(sample.129c, "phospho_129c.xlsx", row.names = FALSE);
  cat("Writing Sample 130n file\n");
  write.xlsx(sample.130n, "phospho_130n.xlsx", row.names = FALSE);
  cat("Writing Sample 130c file\n");
  write.xlsx(sample.130c, "phospho_130c.xlsx", row.names = FALSE);
  cat("Writing Sample 131 file\n");
  write.xlsx(sample.131, "phospho_131.xlsx", row.names = FALSE);
  
  
  cat("Writing hybridization file\n");
  write.xlsx(hybridization.file, "phospho_hybridization.xlsx", row.names = FALSE);
}

cluster.heatmap <- function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, 
          distfun = dist, hclustfun = hclust, reorderfun = function(d, 
                                                                    w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv, 
                                                                                                                               "Rowv"), scale = c("row", "column", "none"), na.rm = TRUE, 
          margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 + 
            1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
          labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, 
          verbose = getOption("verbose"), ...) 
{
  scale <- if (symm && missing(scale)) 
    "none"
  else match.arg(scale)
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
    stop("'x' must be a numeric matrix")
  nr <- di[1L]
  nc <- di[2L]
  if (nr <= 1 || nc <= 1) 
    stop("'x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2L) 
    stop("'margins' must be a numeric vector of length 2")
  doRdend <- !identical(Rowv, NA)
  doCdend <- !identical(Colv, NA)
  if (!doRdend && identical(Colv, "Rowv")) 
    doCdend <- FALSE
  if (is.null(Rowv)) 
    Rowv <- rowMeans(x, na.rm = na.rm)
  if (is.null(Colv)) 
    Colv <- colMeans(x, na.rm = na.rm)
  if (doRdend) {
    if (inherits(Rowv, "dendrogram")) 
      ddr <- Rowv
    else {
      hcr <- hclustfun(distfun(x), method = "ward.D2")
      ddr <- as.dendrogram(hcr)
      if (!is.logical(Rowv) || Rowv) 
        ddr <- reorderfun(ddr, Rowv)
    }
    if (nr != length(rowInd <- order.dendrogram(ddr))) 
      stop("row dendrogram ordering gave index of wrong length")
  }
  else rowInd <- 1L:nr
  if (doCdend) {
    if (inherits(Colv, "dendrogram")) 
      ddc <- Colv
    else if (identical(Colv, "Rowv")) {
      if (nr != nc) 
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      ddc <- ddr
    }
    else {
      hcc <- hclustfun(distfun(if (symm) 
        x
        else t(x)), method = "ward.D2")
      ddc <- as.dendrogram(hcc)
      if (!is.logical(Colv) || Colv) 
        ddc <- reorderfun(ddc, Colv)
    }
    if (nc != length(colInd <- order.dendrogram(ddc))) 
      stop("column dendrogram ordering gave index of wrong length")
  }
  else colInd <- 1L:nc
  x <- x[rowInd, colInd]
  labRow <- if (is.null(labRow)) 
    if (is.null(rownames(x))) 
      (1L:nr)[rowInd]
  else rownames(x)
  else labRow[rowInd]
  labCol <- if (is.null(labCol)) 
    if (is.null(colnames(x))) 
      (1L:nc)[colInd]
  else colnames(x)
  else labCol[colInd]
  if (scale == "row") {
    x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
    sx <- apply(x, 1L, sd, na.rm = na.rm)
    x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
  }
  else if (scale == "column") {
    x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
    sx <- apply(x, 2L, sd, na.rm = na.rm)
    x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
  }
  lmat <- rbind(c(NA, 3), 2:1)
  lwid <- c(if (doRdend) 1 else 0.05, 4)
  lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 
            4)
  if (!missing(ColSideColors)) {
    if (!is.character(ColSideColors) || length(ColSideColors) != 
        nc) 
      stop("'ColSideColors' must be a character vector of length ncol(x)")
    lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
    lhei <- c(lhei[1L], 0.2, lhei[2L])
  }
  if (!missing(RowSideColors)) {
    if (!is.character(RowSideColors) || length(RowSideColors) != 
        nr) 
      stop("'RowSideColors' must be a character vector of length nrow(x)")
    lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
                                   1), lmat[, 2] + 1)
    lwid <- c(lwid[1L], 0.2, lwid[2L])
  }
  lmat[is.na(lmat)] <- 0
  if (verbose) {
    cat("layout: widths = ", lwid, ", heights = ", lhei, 
        "; lmat=\n")
    print(lmat)
  }
  dev.hold()
  on.exit(dev.flush())
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1L], 0, 0, 0.5))
    image(rbind(if (revC) 
      nr:1L
      else 1L:nr), col = RowSideColors[rowInd], axes = FALSE)
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2L]))
    image(cbind(1L:nc), col = ColSideColors[colInd], axes = FALSE)
  }
  par(mar = c(margins[1L], 0, 0, margins[2L]))
  if (!symm || scale != "none") 
    x <- t(x)
  if (revC) {
    iy <- nr:1
    if (doRdend) 
      ddr <- rev(ddr)
    x <- x[, iy]
  }
  else iy <- 1L:nr
  image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
          c(0, nr), axes = FALSE, xlab = "", ylab = "", ...)
  axis(1, 1L:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
       cex.axis = cexCol)
  if (!is.null(xlab)) 
    mtext(xlab, side = 1, line = margins[1L] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
       cex.axis = cexRow)
  if (!is.null(ylab)) 
    mtext(ylab, side = 4, line = margins[2L] - 1.25)
  if (!missing(add.expr)) 
    eval(substitute(add.expr))
  par(mar = c(margins[1L], 0, 0, 0))
  if (doRdend) 
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  else frame()
  par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2L]))
  if (doCdend) 
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  else if (!is.null(main)) 
    frame()
  if (!is.null(main)) {
    par(xpd = NA)
    title(main, cex.main = 1.5 * op[["cex.main"]])
  }
  invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
                                                              doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
}

intensity.distribution.10plex <- function(peptide.ratios, peporpro) {
  file.name <- paste(peporpro, "log10_intensity_distribution.png", sep = "_");
  
  png(filename = file.name, width = 1500, height = 600);
  par(mfrow=c(2,5), cex.lab = 1.8, cex.axis = 1.8, col = "dodgerblue4", 
      col.axis = "dodgerblue4", col.lab = "dodgerblue4", mar = c(5, 5, 1, 1),
      mgp = c(3, 1, 0));
  hist(peptide.ratios$log10.intensity.126, main = NULL, 
       xlab = "log10.intensity.126", col = "dodgerblue4");
  hist(peptide.ratios$log10.intensity.127n, main = NULL, 
       xlab = "log10.intensity.127n", col = "dodgerblue4");
  hist(peptide.ratios$log10.intensity.127c, main = NULL, 
       xlab = "log10.intensity.127c", col = "dodgerblue4");
  hist(peptide.ratios$log10.intensity.128n, main = NULL, 
       xlab = "log10.intensity.128n", col = "dodgerblue4");
  hist(peptide.ratios$log10.intensity.128c, main = NULL, 
       xlab = "log10.intensity.128c", col = "dodgerblue4");
  hist(peptide.ratios$log10.intensity.129n, main = NULL, 
       xlab = "log10.intensity.129n", col = "dodgerblue4");
  hist(peptide.ratios$log10.intensity.129c, main = NULL, 
       xlab = "log10.intensity.129c", col = "dodgerblue4");
  hist(peptide.ratios$log10.intensity.130n, main = NULL, 
       xlab = "log10.intensity.130n", col = "dodgerblue4");
  hist(peptide.ratios$log10.intensity.130c, main = NULL, 
       xlab = "log10.intensity.130c", col = "dodgerblue4");
  hist(peptide.ratios$log10.intensity.131, main = NULL, 
       xlab = "log10.intensity.131", col = "dodgerblue4");
  
  dev.off()
}

SILAC.ratio.intensity.plot <- function(df, ratio.col, intensity.col,fm) {
  df <- df[!is.na(df[, ratio.col]), ];
  par(mfcol = c(1, 1), mar = c(5, 5, 4, 4));
  png(filename = fm, width = 480, height = 500);
  plot(log2(df[, ratio.col]), log10(df[, intensity.col]), 
       xlab = "Log2 Raio", ylab = "Log10 Intensity", pch = 1, 
       main = "Ratio Intensity plot", cex.main = 2, cex.lab = 1.5, cex.axis = 1.5,
       # xlim=c(-6, 6), ylim=c(-6, 6), 
       col = "deepskyblue4");
  abline(v = median(log2(df[, ratio.col])), col = "red")
  dev.off()
}