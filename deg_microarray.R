# Differential gene expression analysis in R

# CEL files preparation 
untar('GSE28735_RAW.tar')
cels <- list.files("GSE28735/", pattern = "[gz]")
sapply(paste("GSE28735", cels, sep = "/"), gunzip)
# PT = pancreatic tumor samples
# NT = Non-tumor samples

# Load library
library(GEOquery)
library(limma)
library(affy)
library(hugene10sttranscriptcluster.db)
library(annotate)

setwd("~/GSE28735")
# Read CEL files with affy package
GSE28735 <- ReadAffy()

# Convert AffyBatch object to ExpressionSet Object
# rma():
#1. background correction
#2. quantile normalization
#3. summarization by median polish algorithm  
eset <- rma(GSE28735)

# Write normalized expression values into file for saving purpose
write.exprs(eset, file = "Normalized.txt")

# Identify differentially expressed genes using linear model
# Generate matrix model, specifying design
x <- c(rep(1,45), rep(2,45))
design <- model.matrix(~ 0 + factor(x))
colnames(design) <- c("NT","PT")

# Add gene annotation for probeset ID
ID <- featureNames(eset)
Symbol <- getSYMBOL(ID, "hugene10sttranscriptcluster.db")
fData(eset) <- data.frame(ID = ID, Symbol = Symbol)

# Fit linear model to the design 
fit <- lmFit(eset, design)

# Construct contrast matrix to be used in Empirical Bayes method (eBayes) 
cont.matrix <- makeContrasts(NT-PT, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# Extract a table of the top-ranked genes from linear model fit 
# fdr = False discovery rate 
# B value = Log odds (The higher the Log Odds for each gene,
# the higher the probability that the gene is differentially expressed
# and not a false positive.)
results <- topTable(fit2, adjust.method = "fdr", 
                    sort.by = "B", number = Inf)
results_down <- results[1:10,]
results_down
results_up <- topTable(fit2, adjust.method = "fdr", 
                    sort.by = "B", number = Inf, lfc = 2.5)
results_up[1:20,]

# Construct volcano plot 
# Summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method = "fdr", p.value = 0.05)
# List contrast names
colnames(fit2)
# Choose contrast of interest
ct <- 1       
volcanoplot(fit2, coef = ct, main = colnames(fit2)[ct], pch = 20,
            highlight = length(which(dT[,ct]!= 0)), 
            names = rep('.', nrow(fit2)))
