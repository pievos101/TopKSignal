# TCGA Cancer -- Workflow
# Required packages 
library(TopKSignal)
#library(ggplot2)
library(limma) # installed from Bioconductor

### READ IN DATA ############################################
#############################################################
# KIDNEY Gene expression values and clinical data can be downloaded from here:
# http://acgt.cs.tau.ac.il/multi_omic_benchmark/download.html.

# Source the function to read TCGA data matrices
# source("readTCGA.R")
# Note, the readTCGA function can be found at the very end of this file.
# It is needed to read in the Cancer Gene Expression Matrix as well as the clinical data (survival status)

# Location of the file 
loc="~/TCGA_data/NAR Data/"

# Set parameter for dimension reduction
n.patients <- 50 

cat("Reading the data ...\n")
# READ IN!
tcga       <- readTCGA("kidney", loc)

## Read in known genes
Known_Genes   <- read.table("multiomics_features.tsv", sep="\t")

### END OF READ IN DATA #####################################
#############################################################

cat("Remove genes with ZERO expression ...\n")

DATA_full_DEATH <- tcga$FULL_TCGA_DEATH # for death outcomes

# Remove genes with ZERO expression
ids             <- which(apply(DATA_full_DEATH,1,function(x){any(x==0)}))
DATA_full_DEATH <- DATA_full_DEATH[-ids,]
# -------------------------------

DATA_full_SURV <- tcga$FULL_TCGA_SURV   # for survived outcomes

# Remove genes with ZERO expression
ids            <- which(apply(DATA_full_SURV,1,function(x){any(x==0)}))
DATA_full_SURV <- DATA_full_SURV[-ids,]
# -------------------------------

# Take the genes which have no zero expression in one of the pops
common_GENES <- intersect(rownames(DATA_full_DEATH),rownames(DATA_full_SURV))

DATA_full_DEATH <- DATA_full_DEATH[common_GENES,]
DATA_full_SURV  <- DATA_full_SURV[common_GENES,]


cat("Dimensions are:\n")
print(dim(DATA_full_DEATH))
print(dim(DATA_full_SURV))


cat("Filter by known cancer genes ...\n")
known_genes   <- Known_Genes[,1]

current_genes <- rownames(DATA_full_DEATH)
current_genes <- strsplit(current_genes,split="|", fixed=TRUE)
current_genes <- sapply(current_genes,function(x){x[1]})
ids           <- match(current_genes, known_genes)
del.genes     <- which(is.na(ids))


DATA_full_DEATH <- DATA_full_DEATH[-del.genes,]
DATA_full_SURV  <- DATA_full_SURV[-del.genes,]


cat("Dimensions are:\n")
print(dim(DATA_full_DEATH))
print(dim(DATA_full_SURV))
### End of filter by known_genes


cat("Subsample patients ...\n")

# Subsample patients ---------------------------- # 
ids              <- sample(1:dim(DATA_full_DEATH)[2], n.patients, replace=FALSE)
DATA_full_DEATH  <- DATA_full_DEATH[,ids]


ids              <- sample(1:dim(DATA_full_SURV)[2], n.patients, replace=FALSE)
DATA_full_SURV   <- DATA_full_SURV[,ids]
# ------------------------------------------ #

save(DATA_full_DEATH, file="DATA_full_DEATH.RData")
save(DATA_full_SURV,  file="DATA_full_SURV.RData")

cat("Dimensions are:\n")
print(dim(DATA_full_DEATH))
print(dim(DATA_full_SURV))

cat("Voom transformation by limma ... \n")

COUNTS   <- cbind(DATA_full_DEATH, DATA_full_SURV)
V        <- voom(COUNTS)
COUNTS_norm <- V$E
VAR      <- apply(COUNTS_norm,1,var)


cat("Select genes by variances ... \n")

# 5% high, 5% low, 5% random

n.genes <- floor((length(VAR)/100)*5) # 10%

# top and bottom VARIANCE genes
topgenes1 <- names(head(sort(VAR, decreasing=TRUE), n.genes))
topgenes2 <- names(tail(sort(VAR, decreasing=TRUE), n.genes))

# random genes
randgenes <- sample(names(VAR),n.genes)

topgenes  <- unique(c(topgenes1, topgenes2, randgenes))

# Reduce data - VOOM
DATA_top_DEATH <- COUNTS_norm[topgenes, DEATH_patients] 
DATA_top_SURV  <- COUNTS_norm[topgenes, SURV_patients]

cat("Final dimensions are:\n")
print(dim(DATA_top_DEATH))
print(dim(DATA_top_SURV))

save(DATA_top_DEATH, file="DATA_top_DEATH.RData")
save(DATA_top_SURV, file="DATA_top_SURV.RData")


# ---------------------- CONSENSUS SIGNALS ------------------------------ #
###########################################################################
# ----------------------------------------------------------------------- #


## ---------------------- DEATH OUTCOMES -------------------------------- #
###########################################################################

cat("Signal estimation (non-survived outcomes)\n")

rankMatrix_top <- apply(DATA_top_DEATH, 2, function(x) rank(-x, ties.method = "random"))

rownames(rankMatrix_top) <- paste("r",1:dim(rankMatrix_top)[1], sep="")
colnames(rankMatrix_top) <- paste("c",1:dim(rankMatrix_top)[2], sep="")

# Peform estimation via TopKSignal 
res_DEATH  <- estimateTheta(R.input = rankMatrix_top, b = 0.1, num.boot = 500, solver = "gurobi",type ="restrictedQuadratic", bootstrap.type =  "poisson.bootstrap", nCore=6)
save(res_DEATH, file="res_DEATH.RData")

TKSrank          <- rank(-res_DEATH$estimation$signal.estimate)

names(TKSrank)   <- rownames(DATA_top_DEATH)

TKSrank_DEATH    <- TKSrank 

signal_DEATH                 <- res_DEATH$estimation$signal.estimate
signal_SE_DEATH              <- res_DEATH$estimation$SE
names(signal_DEATH)          <- rownames(DATA_top_DEATH)
names(signal_SE_DEATH)       <- rownames(DATA_top_DEATH)


signal_sorted_DEATH    <- sort(signal_DEATH, decreasing=TRUE)
ids                    <- match(names(signal_sorted_DEATH),names(signal_SE_DEATH))
signal_sorted_SE_DEATH <- signal_SE_DEATH[ids]

# Plots - Elbow Plot #####################
##########################################
## ggplot
#library(ggplot2)
#pdf("DEATH_ELBOW_quadratic.pdf")
#plot.elbow <- ggplot(data = data.frame(),aes(x = 1:length(signal_sorted_DEATH), y = signal_sorted_DEATH)) + geom_point() + geom_line() +
#        geom_errorbar(aes(ymin=signal_sorted_DEATH-signal_sorted_SE_DEATH, ymax=signal_sorted_DEATH + signal_sorted_SE_DEATH), width=.2,
#                 position=position_dodge(0.05),colour="black") + 
#        xlab("Genes (ordered)") +
#        ylab("Signal estimate (theta)") +
#        ylim(c(-2,4)) +
#	theme_light()
#plot.elbow #
#dev.off()

## ---------------------- SURVIVED OUTCOMES -------------------------------- #
###########################################################################

cat("Signal estimation (survived outcomes)\n")

rankMatrix_top <- apply(DATA_top_SURV, 2, function(x) rank(-x, ties.method = "random"))

rownames(rankMatrix_top) <- paste("r",1:dim(rankMatrix_top)[1], sep="")
colnames(rankMatrix_top) <- paste("c",1:dim(rankMatrix_top)[2], sep="")

# Peform estimation via TopKSignal 
res_SURV  <- estimateTheta(R.input = rankMatrix_top, b = 0.1, num.boot = 500, solver = "gurobi",type ="restrictedQuadratic", bootstrap.type =  "poisson.bootstrap", nCore=6)
save(res_SURV, file="res_SURV.RData")

TKSrank          <- rank(-res_SURV$estimation$signal.estimate)

names(TKSrank)   <- rownames(DATA_top_SURV)

TKSrank_SURV     <- TKSrank 

signal_SURV                 <- res_SURV$estimation$signal.estimate
signal_SE_SURV              <- res_SURV$estimation$SE
names(signal_SURV)          <- rownames(DATA_top_SURV)
names(signal_SE_SURV)       <- rownames(DATA_top_SURV)


signal_sorted_SURV    <- sort(signal_SURV, decreasing=TRUE)
ids <- match(names(signal_sorted_SURV),names(signal_SE_SURV))
signal_sorted_SE_SURV <- signal_SE_SURV[ids]

# Plots - Elbow Plot #####################
##########################################
## ggplot
#library(ggplot2)
#pdf("SURVIVED_ELBOW_quadratic.pdf")
#plot.elbow <- ggplot(data = data.frame(),aes(x = 1:length(signal_sorted_SURV),y = signal_sorted_SURV)) + geom_point() + geom_line() +
#        geom_errorbar(aes(ymin=signal_sorted_SURV-signal_sorted_SE_SURV, ymax=signal_sorted_SURV + signal_sorted_SE_SURV), width=.2,
#                 position=position_dodge(0.05),colour="black") + 
#        xlab("Genes (ordered)") +
#        ylab("Signal estimate (theta)") +
#        ylim(c(-2,4)) +
#	theme_light()
#plot.elbow #
#dev.off()



#### SUBFUNCTIONS ###############################################################
#################################################################################

readTCGA <- function(cancertype=c("kidney"), loc="~/TCGA_data/NAR Data/"){

TCGA_death <- list()
TCGA_surv  <- list()

 for(xx in 1:length(cancertype)){

        ctype <- cancertype[xx]

        # NAR TCGA  -------------------------
        path1 <- paste(loc, ctype, "/exp", sep="")
        path2 <- paste(loc, ctype, "/survival", sep="")

        tcga_exp      <- read.table(path1)
        tcga_survival <- read.table(path2)

        # Fix names 
        colnames(tcga_exp) <- gsub(".", "-", colnames(tcga_exp), fixed=TRUE) 

        # Get the Death and Survived Outcomes 
        death_ids <- which(tcga_survival[,3]==1)
        surv_ids  <- which(tcga_survival[,3]==0)

        death_pat <- as.character(tcga_survival[,1][death_ids])
        surv_pat  <- as.character(tcga_survival[,1][surv_ids])
        
        # death
        iii       <- match(death_pat, colnames(tcga_exp))
        iii       <- iii[!is.na(iii)]
 
        # survived
        jjj       <- match(surv_pat, colnames(tcga_exp))
        jjj       <- jjj[!is.na(jjj)]

        tcga_exp_death  <- tcga_exp[,iii]
        tcga_exp_surv   <- tcga_exp[,jjj]
        
        # fix the rownames
        rnames <- rownames(tcga_exp_death)
        rnames <- gsub(".", "|", rnames, fixed=TRUE)
        rownames(tcga_exp_death) <- rnames 

        rnames <- rownames(tcga_exp_surv)
        rnames <- gsub(".", "|", rnames, fixed=TRUE)
        rownames(tcga_exp_surv) <- rnames       
        
        TCGA_death[[xx]] <- tcga_exp_death
        TCGA_surv[[xx]]  <- tcga_exp_surv 
 }

# Concatenate
## Intersection of genes across cancer types
if(length(cancertype)==1){
return(list(FULL_TCGA_DEATH=TCGA_death[[1]], FULL_TCGA_SURV=TCGA_surv[[1]]))
}

genes_inter <- NULL
 for(xx in 1:(length(cancertype)-1)){
    genes_inter <- c(genes_inter,intersect(rownames(TCGA_death[[xx]]),
                                           rownames(TCGA_death[[xx+1]])))
 }

# Sample the same number of patients accross the cancer types
n.samples <- sapply(TCGA_death, function(x){dim(x)[2]})
n.samples <- min(n.samples)

# Concatenate cancer types (exp) DEATH
FULL_TCGA_DEATH <- TCGA_death[[1]][genes_inter,1:n.samples]
 for(xx in 2:length(cancertype)){
    FULL_TCGA_DEATH <- cbind(FULL_TCGA_DEATH, TCGA_death[[xx]][genes_inter,1:n.samples])        
 }
# Concatenate cancer types (exp) SURVIVED
FULL_TCGA_SURV <- TCGA_surv[[1]][genes_inter,1:n.samples]
 for(xx in 2:length(cancertype)){
    FULL_TCGA_SURV <- cbind(FULL_TCGA_SURV, TCGA_surv[[xx]][genes_inter,1:n.samples])   
 }

return(list(FULL_TCGA_DEATH=FULL_TCGA_DEATH,FULL_TCGA_SURV=FULL_TCGA_SURV))

}# end of function 

