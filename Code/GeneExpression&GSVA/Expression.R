#load required packages
library(BiocManager)
library(readxl)
library(ggplot2)
library(ggrepel)
library(gage)
library(xlsx)
library(fgsea)
library(KEGGprofile)
library(writexl)
library(tibble)
library(limma)
library(reshape2)
library(pheatmap)
library(grid)

#load experimental data
setwd("/Users/julia/Documents/Internship/Boolean Model/Version 3 (advanced)/Verification")
expdata <- read.csv2("exp_data.csv")

#load list of enzymeGenes, all genes and enzyme names from metabolic model
enzGene <- read.table("enzGenes.txt", quote="\"", comment.char="")
enzGene <- lapply(enzGene, as.character)
enzGene <- enzGene$V1
Enz <- read.table("POI.txt", quote="\"", comment.char="")
Enz <- lapply(Enz, as.character)
Enz <- Enz$V1

dict <- data.frame(enzGene, Enz)

#calculate mean of reference data (C-source:glucose, N-source:Ammonium, limited:Zinc)
ref <- apply(expdata[c("X53","X54","X55")],1, mean)

#calculate mean for glucose limitation
gluc_lim <- apply(expdata[c("X101", "X102", "X103")],1, mean)

#calculate mean for nitrogen limitation
nitr_lim <- apply(expdata[c("X44", "X46", "X56")],1, mean)

#merge data into datafame
gene <- expdata$yeast.systematic.name
gene <- apply(as.data.frame(gene), 1, function(genename){
  genename<-gsub("'", "", genename)
})

expdata2 <- data.frame(gene, expdata[c("X53", "X54", "X55", "X101", "X102", "X103", "X44", "X46", "X56")])

expdata <- data.frame(gene, ref, gluc_lim, nitr_lim)


#extract genes of interest (part of metabolic model)
bool<-lapply(gene, function(genename){
  bool <- is.element(genename,enzGene);
})

expdata <- expdata[unlist(bool),]
expdata2 <- expdata2[unlist(bool),]
expdata<- expdata[!duplicated(expdata$gene),]
expdata2<- expdata2[!duplicated(expdata2$gene),]

#change rownames to enzyme names
rownames(expdata) <- sapply(expdata[,1], function(name){
  dict[name == dict[,1],2]
})
rownames(expdata2) <- sapply(expdata2[,1], function(name){
  dict[name == dict[,1],2]
})

expdata$gene <- NULL
expdata2$gene <- NULL

#order rownames
expdata2 <- expdata2[order(row.names(expdata2)),]
expdata <- expdata[order(row.names(expdata)),]

#quality control: boxplot of expr data
plotdata <- melt(expdata2, id.var=0)
class <- data.frame(matrix(ncol=1, nrow=dim(plotdata)[1]))
class <- sapply(plotdata[,1], function(name){
  if (name == 'X53'| name == 'X54' | name == 'X55'){
    class = '1|1'
  } else if (name == 'X101'| name == 'X102' | name == 'X103'){
    class='0|1'
  } else if (name == 'X44'| name == 'X46' | name == 'X56'){
    class='1|0'
  }
})

plotdata <- data.frame(plotdata,class)

g<- ggplot(data=plotdata, aes(x=variable, y=value, fill=plotdata$class)) + geom_boxplot(outlier.shape=4, outlier.size=3) 
g <- g + theme_classic(base_size=16) + labs(y='log2(Expression)', x='Array ID') + scale_fill_manual(name='Nutrient available \n (glc|nitr)', values=c('grey67','gray35','gray93'))
g <- g + scale_x_discrete(labels= c('53', '54', '55', '101', '102', '103', '44', '46', '56')) 
g
ggsave('boxplot_exp.png', width=7.5, height=5)


################load simulated data################
gluc_lim <- read.delim("~/Documents/Internship/Boolean Model/Version 3 (advanced)/Verification/simulatedData/gluc0nitr1.txt")
nitr_lim <- read.delim("~/Documents/Internship/Boolean Model/Version 3 (advanced)/Verification/simulatedData/gluc1nitr0.txt")
ref <- read.delim("~/Documents/Internship/Boolean Model/Version 3 (advanced)/Verification/simulatedData/gluc1nitr1.txt")
simdata <- data.frame(ref, gluc_lim[,2], nitr_lim[,2])
colnames(simdata) <- c('gene','ref', 'gluc_lim','nitr_lim')
rownames(simdata) <- simdata$gene
simdata$gene <- NULL
#simdata <- 2^(simdata)

plotdata <- melt(simdata, id.var=0)
class <- data.frame(matrix(ncol=1, nrow=dim(simdata)[1]))
class <- sapply(plotdata[,1], function(name){
  if (name == 'ref'){
    class = '1|1'
  } else if (name == 'nitr_lim'){
    class='1|0'
  } else if (name == 'gluc_lim'){
    class='0|1'
  }
})

plotdata <- data.frame(plotdata,class)

g <- ggplot(data=plotdata, aes(x=variable, y=value, fill=plotdata$class)) + geom_boxplot(outlier.shape=4, outlier.size=3) 
g <- g + theme_classic(base_size=16) + labs(y='gene expression rank', x='Simulation') + scale_fill_manual(name='Nutrient available \n (glc|nitr)', values=c('grey67','gray35','gray93'))
g <- g + scale_x_discrete(labels= c('1', '2', '3')) 
g
ggsave('boxplot_simulation.png', width=5, height=5)

#filter for genes that are port of exp data and metabolic model
filtered_simdata <- simdata[rownames(simdata) %in% rownames(expdata),]
filtered_simdata <- filtered_simdata[order(row.names(filtered_simdata)),]

################heatmaps#############################
#normalize data to range between 0 and 1
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

norm_expdata <- normalize(expdata)
#norm_expdata <-sapply(expdata, function(x) {
#  return ((x - min(x)) / (max(x) - min(x)))
#})

norm_simdata <- normalize(filtered_simdata)
#norm_simdata <-apply(filtered_simdata,2, function(x) {
#  return ((x - min(x)) / (max(x) - min(x)))
#})

plotdata <- data.frame(norm_expdata[,1], norm_simdata[,1], norm_expdata[,2], norm_simdata[,2], norm_expdata[,3], norm_simdata[,3])
rownames(plotdata) <- rownames(norm_simdata)
collabels <- c('1|1', '1|1', '0|1', '0|1', '1|0', '1|0')
col_annotation <- data.frame(sample = rep(c("experimental", "simulated"), c(3)))
row.names(col_annotation) <- colnames(plotdata)
my_colour = list(sample = c(experimental = "gray75", simulated = "grey45"))


setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.95, name="vp", just=c("right","top"))), action="prepend")
p <- pheatmap(plotdata, fontsize_row=7.5,border_color=NA,cellwidth=50, color = colorRampPalette(c("darkblue", "lightgrey", "darkorange2"))(50), cluster_cols = FALSE, annotation_col = col_annotation, annotation_colors = my_colour, gaps_col=c(2,4), labels_col=collabels, angle_col=0)
setHook("grid.newpage", NULL, "replace")
grid.text("Nutrient available (glc|nitr)", y=-0.02, x=0.45, gp=gpar(fontsize=12))

###############gsva##############################
kegg <-gmtPathways("c2.cp.kegg.v7.0.symbols.gmt")
gobp <- gmtPathways("c5.bp.v7.0.symbols.gmt")
gomol <- gmtPathways("c5.mf.v7.0.symbols.gmt")
react <- gmtPathways("c2.cp.reactome.v7.0.symbols.gmt")
tohuman <- read.csv("yeastToHuman.csv", header=FALSE, sep=";")

humannames <- lapply(rownames(plotdata), function(name){
  name <- tohuman[which(name == tohuman[,2]),3]
})
humannames <- unlist(humannames)
rownames(plotdata) <- humannames

#calculate mean for duplicated genes (ENO3, AIFM2, ILVBL, PC, GAPDH, PGD, TKTL1, ALDH1A2)
ALDH1A2 <- apply(plotdata[c(9,62),], 2, mean)
AIFM2 <- apply(plotdata[c(74,75),], 2, mean)
ILVBL <- apply(plotdata[c(81,82,83),], 2, mean)
PC<- apply(plotdata[c(91,92),], 2, mean)
GAPDH <- apply(plotdata[c(110,111,112),], 2, mean)
PGD <- apply(plotdata[c(55,56,107),], 2, mean)
TKTL1 <- apply(plotdata[c(114,115),], 2, mean)
rownames(plotdata)[which(is.element(rownames(plotdata), 'ENO3-2'))] <- 'ENO3'

plotdata[c(9,62,74,75,81,82,83,91,92,110,111,112,55,56,114,115,107),] <- NA
plotdata <- na.omit(plotdata)
join <- rbind(ALDH1A2, AIFM2, ILVBL, PC, GAPDH, PGD, TKTL1)
plotdata <- rbind(plotdata, join)


#heatmap(gsva(as.matrix(plotdata[,c(1,3,5)]), gobp, min.sz=4), cluster_cols = FALSE,  color = colorRampPalette(c("darkblue", "lightgrey", "darkorange2"))(50))# cluster_cols = FALSE, fontsize_row=7.5,border_color=NA,cellwidth=50, angle_col=0, color = colorRampPalette(c("darkblue", "lightgrey", "darkorange2"))(50), gaps_col=c(2,4))
#pheatmap(gsva(as.matrix(plotdata), gobp, abs.ranking = FALSE, min.sz=4), cluster_cols = FALSE,  color = colorRampPalette(c("darkblue", "lightgrey", "darkorange2"))(50))# cluster_cols = FALSE, fontsize_row=7.5,border_color=NA,cellwidth=50, angle_col=0, color = colorRampPalette(c("darkblue", "lightgrey", "darkorange2"))(50), gaps_col=c(2,4))

#calculate enrichment scores (via gsva) for both sets(experimental and simulated, score indicates pathway activity)
exp <-gsva(as.matrix(plotdata[,c(1,3,5)]), gobp, min.sz=10, abs.ranking=FALSE)
sim <-gsva(as.matrix(plotdata[,c(2,4,6)]), gobp, min.sz=10, abs.ranking=FALSE)
ind <- which(rownames(sim) %in% rownames(exp))
exp <- exp[ind,]
gsva_res <- data.frame(exp[,1],sim[,1], exp[,2], sim[,2], exp[,3], sim[,3])

collabels <- c('1|1', '1|1', '0|1', '0|1', '1|0', '1|0')
col_annotation <- data.frame(sample = rep(c("experimental", "simulated"), c(3)))
row.names(col_annotation) <- colnames(gsva_res)
my_colour = list(sample = c(experimental = "gray75", simulated = "grey45"))

setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.95, name="vp", just=c("right","top"))), action="prepend")
pheatmap(gsva_res, color = colorRampPalette(c("darkblue", "lightgrey", "darkorange2"))(50),fontsize_row = 10, cluster_cols=FALSE, border_color=NA, cellwidth=50, annotation_col = col_annotation, annotation_colors = my_colour, gaps_col=c(2,4), labels_col=collabels, angle_col=0)
setHook("grid.newpage", NULL, "replace")
grid.text("Nutrient available (glc|nitr)", y=-0.02, x=0.26, gp=gpar(fontsize=12))

#which genes are contained in go sets?
sel_indx <- which(names(gobp) %in% rownames(gsva_res))
selgo <- gobp[sel_indx]

go_res <- list()
for (i in 1:(dim(gsva_res)[1])){
  go_res[rownames(gsva_res)[i]] <- list(rownames(plotdata)[is.element(rownames(plotdata), unlist(selgo[i]))])
}

for (i in 1:length(go_res)){
  for (j in 1:length(go_res[[i]])){
    if (go_res[[i]][j] %in% tohuman[,3]){
      go_res[[i]][j] = as.character(tohuman[which(tohuman[,3] %in% go_res[[i]][j]), 2])
    } else if (go_res[[i]][j] == 'ENO3'){
      go_res[[i]][j] = 'ERR1-3'
    } else if (go_res[[i]][j] == 'PGD'){
      go_res[[i]][j] ='PGD1-3'
    } else if (go_res[[i]][j] == 'AIFM2'){
      go_res[[i]][j] = 'NDE1-2'
    } else if (go_res[[i]][j] == 'ILVBL'){
      go_res[[i]][j] = 'PDC1,5,6'
    } else if (go_res[[i]][j] == 'PC'){
      go_res[[i]][j] = 'PYC1-2'
    } else if (go_res[[i]][j] == 'GAPDH'){
      go_res[[i]][j] = 'TDH1-3'
    } else if (go_res[[i]][j] == 'TKTL1'){
      go_res[[i]][j] = 'TKL1-2'
    } else if (go_res[[i]][j] == 'ALDH1A2'){
      go_res[[i]][j] = 'ALD5/HFD1'
    }
  }
}


#was not included, could be done in general
################limma analysis#####################
sample <- c('ref', 'ref', 'ref', 'gluc_lim', 'gluc_lim', 'gluc_lim', 'nitr_lim', 'nitr_lim', 'nitr_lim')
replicate <- c(1,2,3,1,2,3,1,2,3)
targets <- data.frame(sample, replicate)
f <- factor(targets$sample, levels=c('ref', 'gluc_lim', 'nitr_lim')) 
design <- model.matrix(~0+f)
colnames(design) <- c("ref","gluc_lim","nitr_lim")

fit <- lmFit(expdata, design)

contrast.matrix <- makeContrasts(gluc_lim-ref, nitr_lim-ref, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


