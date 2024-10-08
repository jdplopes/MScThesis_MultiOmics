##Master's Thesis
##João Lopes
set.seed(1)
setwd("C:\\Users\\jdpl2\\OneDrive\\Ambiente de Trabalho\\Mestrado\\2º Ano\\Both Omics")

##Load a R package
library(limma)
library(edgeR)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(seqinr)
library(data.table)
library(tidyverse)
library(tximport)
library(readxl)
library(tidyr)
library(multiGSEA)
windowsFonts(Calibri = windowsFont("Arial"))

dir.create("Plots")
dir.create("Plots/Dotplot")
dir.create("Files")

pathDotplot <- "Plots/Dotplot/"
pathFiles <- "Files/"

dotplot <- function (data,path) {
  o <- ggplot(data, aes(x = contrast, y = pathway)) +
    geom_point(aes(size = 1-padj, color = NES, fill = NES), shape = 21, stroke = 0.5) +  
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  
    scale_size_continuous(range = c(3, 15)) +
    labs(x = "Contrast", y = "Pathway", 
         size = "1-p.adjust", 
         color = "NES",
         fill = "NES") +  
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),  
          legend.position = "bottom",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          panel.grid.major = element_line(linewidth = 1.2, color = "grey90"),  
          panel.grid.minor = element_line(linewidth = 0.8, color = "grey90"),  
          panel.background = element_rect(fill = "grey", color = NA))  
  ggsave(filename = paste(path, "Dotplot.svg", sep = ""), plot = o, device = "svg", width = 15, height = 15)
  ggsave(filename = paste(path, "Dotplot.tiff", sep = ""), plot = o, device = "tiff", width = 15, height = 15)
}
 
#######################################Files/Data######################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

DEG_pathway <- read.csv("DEG_pathway.csv",sep =";")
DEP_pathway <- read.csv("DEP_pathway.csv",sep = ";")
genes_table <- read.csv("genes_ID.csv",sep =";")
DEG_acc <- read.csv("DEG.csv",sep =";")
DEP_acc <- read.csv("DEP.csv",sep=";")
#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#######################################Files/Data######################################

contrasts <- c("CTL_10vMHW2_10","CTL_25vMHW1_25","CTL_25vMHW2_25","MHW1_25vMHW2_25","MHW2_10vMHW2_25")

proteins <- unique(DEP_pathway$Accession)
genes <- unique(genes_table$Accession)
acc_common <- intersect(proteins,genes);length(acc_common)
de_acc_common <- intersect(unique(DEG_acc$ID),unique(DEP_acc$Accession));length(de_acc_common)

####################################Pathway analysis###################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

##Create data frames with Accession, logFC and Pvalue for the pathway analysis
allg_CTL_10vMHW2_10<- data.frame(DEG_pathway$ID,DEG_pathway$`logFC.CTL_10vMHW2_10`,DEG_pathway$`Pvalue.CTL_10vMHW2_10`)
names(allg_CTL_10vMHW2_10) <- c("Accession","logFC","Pvalue")
allp_CTL_10vMHW2_10<- data.frame(DEP_pathway$Accession,DEP_pathway$`logFC.CTL_10vMHW2_10`,DEP_pathway$`Pvalue.CTL_10vMHW2_10`)
names(allp_CTL_10vMHW2_10) <- c("Accession","logFC","Pvalue")

allg_CTL_25vMHW1_25<- data.frame(DEG_pathway$ID,DEG_pathway$`logFC.CTL_25vMHW1_25`,DEG_pathway$`Pvalue.CTL_25vMHW1_25`)
names(allg_CTL_25vMHW1_25) <- c("Accession","logFC","Pvalue")
allp_CTL_25vMHW1_25<- data.frame(DEP_pathway$Accession,DEP_pathway$`logFC.CTL_10vMHW2_10`,DEP_pathway$`Pvalue.CTL_10vMHW2_10`)
names(allp_CTL_25vMHW1_25) <- c("Accession","logFC","Pvalue")

allg_CTL_25vMHW2_25<- data.frame(DEG_pathway$ID,DEG_pathway$`logFC.CTL_25vMHW2_25`,DEG_pathway$`Pvalue.CTL_25vMHW2_25`)
names(allg_CTL_25vMHW2_25) <- c("Accession","logFC","Pvalue")
allp_CTL_25vMHW2_25<- data.frame(DEP_pathway$Accession,DEP_pathway$`logFC.CTL_25vMHW2_25`,DEP_pathway$`Pvalue.CTL_25vMHW2_25`)
names(allp_CTL_25vMHW2_25) <- c("Accession","logFC","Pvalue")

allg_MHW1_25vMHW2_25<- data.frame(DEG_pathway$ID,DEG_pathway$`logFC.MHW1_25vMHW2_25`,DEG_pathway$`Pvalue.MHW1_25vMHW2_25`)
names(allg_MHW1_25vMHW2_25) <- c("Accession","logFC","Pvalue")
allp_MHW1_25vMHW2_25<- data.frame(DEP_pathway$Accession,DEP_pathway$`logFC.MHW1_25vMHW2_25`,DEP_pathway$`Pvalue.MHW1_25vMHW2_25`)
names(allp_MHW1_25vMHW2_25) <- c("Accession","logFC","Pvalue")

allg_MHW2_10vMHW2_25<- data.frame(DEG_pathway$ID,DEG_pathway$`logFC.MHW2_10vMHW2_25`,DEG_pathway$`Pvalue.MHW2_10vMHW2_25`)
names(allg_MHW2_10vMHW2_25) <- c("Accession","logFC","Pvalue")
allp_MHW2_10vMHW2_25<- data.frame(DEP_pathway$Accession,DEP_pathway$`logFC.MHW2_10vMHW2_25`,DEP_pathway$`Pvalue.MHW2_10vMHW2_25`)
names(allp_MHW2_10vMHW2_25) <- c("Accession","logFC","Pvalue")
############################################################################

##Create a data structure
layers <- c("transcriptome","proteome")
odataCTL_10vMHW2_10 <- initOmicsDataStructure(layer=layers)
odataCTL_25vMHW1_25 <- initOmicsDataStructure(layer=layers)
odataCTL_25vMHW2_25 <- initOmicsDataStructure(layer=layers)
odataMHW1_25vMHW2_25 <- initOmicsDataStructure(layer=layers)
odataMHW2_10vMHW2_25 <- initOmicsDataStructure(layer=layers)
#########################

##Add transcriptome layer
odataCTL_10vMHW2_10$transcriptome <- rankFeatures(allg_CTL_10vMHW2_10$logFC,allg_CTL_10vMHW2_10$Pvalue)
names(odataCTL_10vMHW2_10$transcriptome) <- allg_CTL_10vMHW2_10$Accession
odataCTL_10vMHW2_10$transcriptome <- sort(odataCTL_10vMHW2_10$transcriptome)
head(odataCTL_10vMHW2_10$transcriptome)

odataCTL_25vMHW1_25$transcriptome <- rankFeatures(allg_CTL_25vMHW1_25$logFC,allg_CTL_25vMHW1_25$Pvalue)
names(odataCTL_25vMHW1_25$transcriptome) <- allg_CTL_25vMHW1_25$Accession
odataCTL_25vMHW1_25$transcriptome <- sort(odataCTL_25vMHW1_25$transcriptome)
head(odataCTL_25vMHW1_25$transcriptome)

odataCTL_25vMHW2_25$transcriptome <- rankFeatures(allg_CTL_25vMHW2_25$logFC,allg_CTL_25vMHW2_25$Pvalue)
names(odataCTL_25vMHW2_25$transcriptome) <- allg_CTL_25vMHW2_25$Accession
odataCTL_25vMHW2_25$transcriptome <- sort(odataCTL_25vMHW2_25$transcriptome)
head(odataCTL_25vMHW2_25$transcriptome)

odataMHW1_25vMHW2_25$transcriptome <- rankFeatures(allg_MHW1_25vMHW2_25$logFC,allg_MHW1_25vMHW2_25$Pvalue)
names(odataMHW1_25vMHW2_25$transcriptome) <- allg_MHW1_25vMHW2_25$Accession
odataMHW1_25vMHW2_25$transcriptome <- sort(odataMHW1_25vMHW2_25$transcriptome)
head(odataMHW1_25vMHW2_25$transcriptome)

odataMHW2_10vMHW2_25$transcriptome <- rankFeatures(allg_MHW2_10vMHW2_25$logFC,allg_MHW2_10vMHW2_25$Pvalue)
names(odataMHW2_10vMHW2_25$transcriptome) <- allg_MHW2_10vMHW2_25$Accession
odataMHW2_10vMHW2_25$transcriptome <- sort(odataMHW2_10vMHW2_25$transcriptome)
head(odataMHW2_10vMHW2_25$transcriptome)
####################

##Add proteome layer
odataCTL_10vMHW2_10$proteome <- rankFeatures(allp_CTL_10vMHW2_10$logFC,allp_CTL_10vMHW2_10$Pvalue)
names(odataCTL_10vMHW2_10$proteome) <- allp_CTL_10vMHW2_10$Accession
odataCTL_10vMHW2_10$proteome <- sort(odataCTL_10vMHW2_10$proteome)
head(odataCTL_10vMHW2_10$proteome)

odataCTL_25vMHW1_25$proteome <- rankFeatures(allp_CTL_25vMHW1_25$logFC,allp_CTL_25vMHW1_25$Pvalue)
names(odataCTL_25vMHW1_25$proteome) <- allp_CTL_25vMHW1_25$Accession
odataCTL_25vMHW1_25$proteome <- sort(odataCTL_25vMHW1_25$proteome)
head(odataCTL_25vMHW1_25$proteome)

odataCTL_25vMHW2_25$proteome <- rankFeatures(allp_CTL_25vMHW2_25$logFC,allp_CTL_25vMHW2_25$Pvalue)
names(odataCTL_25vMHW2_25$proteome) <- allp_CTL_25vMHW2_25$Accession
odataCTL_25vMHW2_25$proteome <- sort(odataCTL_25vMHW2_25$proteome)
head(odataCTL_25vMHW2_25$proteome)

odataMHW1_25vMHW2_25$proteome <- rankFeatures(allp_MHW1_25vMHW2_25$logFC,allp_MHW1_25vMHW2_25$Pvalue)
names(odataMHW1_25vMHW2_25$proteome) <- allp_MHW1_25vMHW2_25$Accession
odataMHW1_25vMHW2_25$proteome <- sort(odataMHW1_25vMHW2_25$proteome)
head(odataMHW1_25vMHW2_25$proteome)

odataMHW2_10vMHW2_25$proteome <- rankFeatures(allp_MHW2_10vMHW2_25$logFC,allp_MHW2_10vMHW2_25$Pvalue)
names(odataMHW2_10vMHW2_25$proteome) <- allp_MHW2_10vMHW2_25$Accession
odataMHW2_10vMHW2_25$proteome <- sort(odataMHW2_10vMHW2_25$proteome)
head(odataMHW2_10vMHW2_25$proteome)
####################

##Select the databases we want to query and download pathway definitions
databases <- c("kegg")
pathways <- getMultiOmicsFeatures(dbs = databases, layer = layers,
                                  returnTranscriptome = "UNIPROT",
                                  returnProteome = "UNIPROT",
                                  organism = "drerio",
                                  useLocal =  FALSE)
pathways_short <- lapply(names(pathways), function(name) {
  head(pathways[[name]], 2)
})
names(pathways_short) <- names(pathways)
pathways_short
pathways$transcriptome[8]
########################################################################

##Run the pathway enrichment
enrichment_scoresCTL_10vMHW2_10 <- multiGSEA(pathways,odataCTL_10vMHW2_10)
Tenrichment_scoresCTL_10vMHW2_10 <- as.data.frame(enrichment_scoresCTL_10vMHW2_10$transcriptome)
Tenrichment_scoresCTL_10vMHW2_10$leadingEdge <- sapply(Tenrichment_scoresCTL_10vMHW2_10$leadingEdge, function(x) paste(x, collapse = ";"))
Tenrichment_scoresCTL_10vMHW2_10 <- Tenrichment_scoresCTL_10vMHW2_10[order(Tenrichment_scoresCTL_10vMHW2_10$padj),]
Tenrichment_scoresCTL_10vMHW2_10$contrast <- rep("CTL_10vMHW2_10")
top_10_esCTL_10vMHW2_10 <- head(Tenrichment_scoresCTL_10vMHW2_10,10)
write.table(top_10_esCTL_10vMHW2_10,paste(pathFiles,"top10_enrichment_scoresCTL_10vMHW2_10.csv",sep=""),sep=";",row.names = FALSE)

enrichment_scoresCTL_25vMHW1_25 <- multiGSEA(pathways,odataCTL_25vMHW1_25)
Tenrichment_scoresCTL_25vMHW1_25 <- as.data.frame(enrichment_scoresCTL_25vMHW1_25$transcriptome)
Tenrichment_scoresCTL_25vMHW1_25$leadingEdge <- sapply(Tenrichment_scoresCTL_25vMHW1_25$leadingEdge, function(x) paste(x, collapse = ";"))
Tenrichment_scoresCTL_25vMHW1_25 <- Tenrichment_scoresCTL_25vMHW1_25[order(Tenrichment_scoresCTL_25vMHW1_25$padj),]
Tenrichment_scoresCTL_25vMHW1_25$contrast <- rep("CTL_25vMHW1_25")
top_10_esCTL_25vMHW1_25 <- head(Tenrichment_scoresCTL_25vMHW1_25,10)
write.table(top_10_esCTL_25vMHW1_25,paste(pathFiles,"top10_enrichment_scoresCTL_25vMHW1_25.csv",sep=""),sep=";",row.names = FALSE)

enrichment_scoresCTL_25vMHW2_25 <- multiGSEA(pathways,odataCTL_25vMHW2_25)
Tenrichment_scoresCTL_25vMHW2_25 <- as.data.frame(enrichment_scoresCTL_25vMHW2_25$transcriptome)
Tenrichment_scoresCTL_25vMHW2_25$leadingEdge <- sapply(Tenrichment_scoresCTL_25vMHW2_25$leadingEdge, function(x) paste(x, collapse = ";"))
Tenrichment_scoresCTL_25vMHW2_25 <- Tenrichment_scoresCTL_25vMHW2_25[order(Tenrichment_scoresCTL_25vMHW2_25$padj),]
Tenrichment_scoresCTL_25vMHW2_25$contrast <- rep("CTL_25vMHW2_25")
top_10_esCTL_25vMHW2_25 <- head(Tenrichment_scoresCTL_25vMHW2_25,10)
write.table(top_10_esCTL_25vMHW2_25,paste(pathFiles,"top10_enrichment_scoresCTL_25vMHW2_25.csv",sep=""),sep=";",row.names = FALSE)

enrichment_scoresMHW1_25vMHW2_25 <- multiGSEA(pathways,odataMHW1_25vMHW2_25)
Tenrichment_scoresMHW1_25vMHW2_25 <- as.data.frame(enrichment_scoresMHW1_25vMHW2_25$transcriptome)
Tenrichment_scoresMHW1_25vMHW2_25$leadingEdge <- sapply(Tenrichment_scoresMHW1_25vMHW2_25$leadingEdge, function(x) paste(x, collapse = ";"))
Tenrichment_scoresMHW1_25vMHW2_25 <- Tenrichment_scoresMHW1_25vMHW2_25[order(Tenrichment_scoresMHW1_25vMHW2_25$padj),]
Tenrichment_scoresMHW1_25vMHW2_25$contrast <- rep("MHW1_25vMHW2_25")
top_10_esMHW1_25vMHW2_25 <- head(Tenrichment_scoresMHW1_25vMHW2_25,10)
write.table(top_10_esMHW1_25vMHW2_25,paste(pathFiles,"top10_enrichment_scoresMHW1_25vMHW2_25.csv",sep=""),sep=";",row.names = FALSE)

enrichment_scoresMHW2_10vMHW2_25 <- multiGSEA(pathways,odataMHW2_10vMHW2_25)
Tenrichment_scoresMHW2_10vMHW2_25 <- as.data.frame(enrichment_scoresMHW2_10vMHW2_25$transcriptome)
Tenrichment_scoresMHW2_10vMHW2_25$leadingEdge <- sapply(Tenrichment_scoresMHW2_10vMHW2_25$leadingEdge, function(x) paste(x, collapse = ";"))
Tenrichment_scoresMHW2_10vMHW2_25 <- Tenrichment_scoresMHW2_10vMHW2_25[order(Tenrichment_scoresMHW2_10vMHW2_25$padj),]
Tenrichment_scoresMHW2_10vMHW2_25$contrast <- rep("MHW2_10vMHW2_25")
top_10_esMHW2_10vMHW2_25 <- head(Tenrichment_scoresMHW2_10vMHW2_25,10)
write.table(top_10_esMHW2_10vMHW2_25,paste(pathFiles,"top10_enrichment_scoresMHW2_10vMHW2_25.csv",sep=""),sep=";",row.names = FALSE)
############################

##Create dataset with all enrichment scores
combined_df <- rbind(top_10_esCTL_10vMHW2_10, top_10_esCTL_25vMHW1_25, top_10_esCTL_25vMHW2_25, top_10_esMHW1_25vMHW2_25, top_10_esMHW2_10vMHW2_25)
combined_df$pathway <- sub("^\\(KEGG\\) ", "", combined_df$pathway)
###########################################

#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
####################################Pathway analysis###################################



##########################################Plots########################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

#Dotplot
dotplot(combined_df,pathDotplot)
########

#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
##########################################Plots########################################