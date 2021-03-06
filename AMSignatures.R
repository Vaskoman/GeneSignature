## FUNCTIONS ####
uniqGeneTabUp <- function (tableUp, tableDown, celltype) {
  aup <- subset(tableUp, tableUp$cellType==celltype)
  adown <- subset(tableDown, tableDown$cellType==celltype)
  aup.u <- aup[!(aup$Gene.ID %in% adown$Gene.ID),]
  return (aup.u)
}

uniqGeneTabDown <- function (tableUp, tableDown, celltype) {
  aup <- subset(tableUp, tableUp$cellType==celltype)
  adown <- subset(tableDown, tableDown$cellType==celltype)
  adown.u <- adown[!(adown$Gene.ID %in% aup$Gene.ID),]
  return (adown.u)
}

averageFC <- function (tab, excols) {
  # Outputs a table with average FC and the first of 
  # other values from input table
  # remove P-values, dataset, species, treatmentTime, sex, GEOlink, papers, and GEO
  tab <- tab[,
             !(colnames(tab) %in% excols)]
  # Sort table based on GENE.ID:
  tab <- tab[order(tab$Gene.ID),]
  # extract gene names
  geneids <- names(table(tab$Gene.ID))
  # extract appearances of each gene name
  geneFC <- table(tab$Gene.ID)
  # Make empty data frame with same # of columns and rows = genes (NAs)
  averFC <- as.data.frame(matrix(data = NA, nrow = length(geneids),
                                 ncol = length(names(tab))+1))
  colnames(averFC) <- c(colnames(tab), "appearances")
  # Populate table with GENE.NAME and GENE.ID:
  averFC$Gene.ID <- unique(tab$Gene.ID)
  
  # loop over every single gene
  for (i in geneids){
    print (i)
    # subset table to a mini0-table (a) for each gene
    a <- tab[tab$Gene.ID==i,]
    # Take average FC and populate table
    averFC$FC[averFC$Gene.ID==i] <- mean(a$FC)
    # take first entry from rest of input 
    averFC[averFC$Gene.ID==i, c("Gene.Name","appearances")] <- c(
      as.character(a$Gene.Name[1]), dim(a)[1])
  }
  return(averFC)
}

## setwd, load libraries, & load data ####
setwd("/home/vassil/Documents/A.Martineau/GeneSignature")
# load tables with all genes
load("//home/vassil/Documents/A.Martineau/GeneSignature/masterTableDown") 
load("/home/vassil/Documents/A.Martineau/GeneSignature/masterTableUp")
# load tables with same-direction genes
load("//home/vassil/Documents/A.Martineau/GeneSignature/tabDown.u")
load("/home/vassil/Documents/A.Martineau/GeneSignature/tabUp.u")

## Datasets metadata ####

# load datasets table:
datasets <- read.delim("Datasets")
olddatasets <- read.delim("olddatasets")

# merge the 2 datasets tables:
metadata <- merge (datasets, olddatasets, by.x = "Old.name",
                   by.y = "olddataset", all.x = TRUE)
names(metadata)[2] <- "name"

# clean up:
rm(list = c("datasets","olddatasets"))

# save table:
write.table(metadata, file = "metadata.txt", sep = "\t",
            row.names = F)

## exclude columns ####
excol <- c("Old.name","cancer","primary","organ","type.y")

## PBMCs ####

# subset genes for PBMCs
tab.up.pbmc <-tabUp.u[tabUp.u$cellType == 'pbmc',]
tab.down.pbmc <- tabDown.u[tabDown.u$cellType == 'pbmc',]

# Get the UP genes that appear in 4/6 datasets:
pbmc.up <- names(table(tab.up.pbmc$Gene.ID)[table(tab.up.pbmc$Gene.ID)>=4])
pbmc.sign.up <- tab.up.pbmc[tab.up.pbmc$Gene.ID %in% pbmc.up,]
# Get the DOWN genes that appear in 4/6 datasets:
pbmc.down <- names(table(tab.down.pbmc$Gene.ID)[table(tab.down.pbmc$Gene.ID)>=4])
pbmc.sign.down <- tab.down.pbmc[tab.down.pbmc$Gene.ID %in% pbmc.down,]
# combine in one table (signature):
pbmc.sign.4.6 <- rbind(pbmc.sign.up, pbmc.sign.down)
# add metadata and order:
pbmc.sign.4.6 <- merge (pbmc.sign.4.6, metadata, by.x = "dataset",
                        by.y = "name", all.x = T)
pbmc.sign.4.6 <- pbmc.sign.4.6[,
                               !(colnames(pbmc.sign.4.6) %in% excol)]
pbmc.sign.4.6 <- pbmc.sign.4.6[order(pbmc.sign.4.6$Gene.Name),]
# clean up
rm(list=c("pbmc.up","pbmc.down","pbmc.sign.up","pbmc.sign.down"))

# Get the UP genes that appear in 5/6 datasets:
pbmc.up <- names(table(tab.up.pbmc$Gene.ID)[table(tab.up.pbmc$Gene.ID)>=5])
pbmc.sign.up <- tab.up.pbmc[tab.up.pbmc$Gene.ID %in% pbmc.up,]
# Get the DOWN genes that appear in 5/6 datasets:
pbmc.down <- names(table(tab.down.pbmc$Gene.ID)[table(tab.down.pbmc$Gene.ID)>=5])
pbmc.sign.down <- tab.down.pbmc[tab.down.pbmc$Gene.ID %in% pbmc.down,]
# combine in one table (signature)
pbmc.sign.5.6 <- rbind(pbmc.sign.up, pbmc.sign.down)
# add metadata and order:
pbmc.sign.5.6 <- merge (pbmc.sign.5.6, metadata, by.x = "dataset",
                        by.y = "name", all.x = T)
pbmc.sign.5.6 <- pbmc.sign.5.6[,
                               !(colnames(pbmc.sign.5.6) %in% excol)]
pbmc.sign.5.6 <- pbmc.sign.5.6[order(pbmc.sign.5.6$Gene.Name),]
# clean up
rm(list=c("pbmc.up","pbmc.down","pbmc.sign.up","pbmc.sign.down"))

# save as tables:
write.table(pbmc.sign.4.6, file="pbmc.sign.4.6.txt", sep="\t", row.names = F)
write.table(pbmc.sign.5.6, file="pbmc.sign.5.6.txt", sep="\t", row.names = F)

# clean up:
rm (list = c("tab.down.pbmc", "tab.up.pbmc", "tabUp.u", "tabDown.u",
             "pbmc.sign.down", "pbmc.down", "pbmc.sign.up", "pbmc.up",
             "pbmc.sign.4.6", "pbmc.sign.5.6"))


## Remove opposite direction genes per cell type in epithelial datasets####
# subset for epithelial only excluding pancreatic, liver, and corneal:
epith.data <- c("E.Br.20h.1", "E.Br.24h.1", "E.Br.24h.2", "E.Br.24h.3", "E.Br.24h.4",
                "E.Br.30d.1", "E.Br.50h.1", "E.Br.8h.1", "E.HN.24h.1", "E.Lu.24h.1",
                "E.Pr.24h.1", "E.Pr.24h.2", "E.Pr.24h.3", "E.Pr.48h.1", "E.Pr.48h.2",
                "E.Pr.48h.3", "E.Pr.6h.1")
tabUp.epith <- masterTableUp[masterTableUp$dataset %in% epith.data,]
tabDown.epith <- masterTableDown[masterTableDown$dataset %in% epith.data,]

# obtain uniquely UP- or DOWN-regulated genes:
tabUp.epith.u <- tabUp.epith[!(tabUp.epith$Gene.ID %in% tabDown.epith$Gene.ID),]
tabDown.epith.u <- tabDown.epith[!(tabDown.epith$Gene.ID %in% tabUp.epith$Gene.ID),]

# genes 
genes.up <- names(table(tabUp.epith.u$Gene.ID)[table(tabUp.epith.u$Gene.ID)>=8])
genes.down <- names(table(tabDown.epith.u$Gene.ID)[table(tabDown.epith.u$Gene.ID)>=4])

# get signature genes in tables:
epith.sig.up <- tabUp.epith.u[tabUp.epith.u$Gene.ID %in% genes.up,]
epith.sig.down <- tabDown.epith.u[tabDown.epith.u$Gene.ID %in% genes.down,]

# combine in one table (signature):
epith.sig <- rbind(epith.sig.up, epith.sig.down)

# add metadata and order:
epith.sig <- merge (epith.sig, metadata, by.x = "dataset",
                        by.y = "name", all.x = T)
epith.sig <- epith.sig[,!(colnames(epith.sig) %in% excol)]

# order by gene name:
epith.sig <- epith.sig[order(epith.sig$Gene.Name),]

# clean up
rm(list=c("epith.data","tabUp.epith","tabDown.epith","tabUp.epith.u",
          "tabDown.epith.u", "genes.down", "genes.up", "epith.sig.down",
          "epith.sig.up"))

# save as tables:
write.table(epith.sig, file="epith.sig.txt", sep="\t", row.names = F)


# Get only the 8 genes that also appear in lung:
lung <- epith.sig[
  epith.sig$Gene.ID %in% epith.sig$Gene.ID[epith.sig$dataset=="E.Lu.24h.1"],]
write.table (lung, file = "lung.txt", sep="\t", row.names=F)


## Make lists with averages ####

# load data:
pbmc.sign.4.6 <- read.delim("pbmc.sign.4.6.txt")
pbmc.sign.5.6 <- read.delim("pbmc.sign.5.6.txt")
epith.sig <- read.delim("epith.sig.txt")
pbmc.sign.5.6 <- read.delim("pbmc.sign.5.6.txt")
lung <- read.delim("lung.txt")

# define columns to be excluded:
excols <- c("dataset","P.Value","species","cellType",
            "treatmentTime","cellName","sex","GEOlink","papers",
            "GEO")

# get average tables
pbmc4out6 <- averageFC(tab = pbmc.sign.4.6, excols = excols)
pbmc5out6 <- averageFC(tab = pbmc.sign.5.6, excols = excols)
epith <- averageFC(tab = epith.sig, excols = excols)
lung.sig <- averageFC(tab = lung, excols = excols)

# sort tables based on FC:
pbmc4out6 <- pbmc4out6[order(pbmc4out6$FC,decreasing=T),]
pbmc5out6 <- pbmc5out6[order(pbmc5out6$FC,decreasing=T),]
epith <- epith[order(epith$FC,decreasing=T),]
epith.lung <- lung.sig[order(lung.sig$FC,decreasing=T),]

# save tables:
write.table(pbmc4out6, file="pbmc4out6.txt", sep="\t", row.names = F)
write.table(pbmc5out6, file="pbmc5out6.txt", sep="\t", row.names = F)
write.table(epith, file="epith.txt", sep="\t", row.names = F)
write.table(epith.lung, file="epith.lung.txt", sep="\t", row.names = F)

# average 2-FC and save
pbmc4out6.2 <- pbmc4out6[abs(pbmc4out6$FC)>=2,]
pbmc5out6.2 <- pbmc5out6[abs(pbmc5out6$FC)>=2,]
epith.2 <- epith[abs(epith$FC)>=2,]
epith.lung.2 <- lung.sig[abs(lung.sig$FC)>=2,]
write.table(pbmc4out6.2, file="pbmc4out6.2.txt", sep="\t", row.names = F)
write.table(pbmc5out6.2, file="pbmc5out6.2.txt", sep="\t", row.names = F)
write.table(epith.2, file="epith.2.txt", sep="\t", row.names = F)
write.table(epith.lung.2, file="epith.lung.2.txt", sep="\t", row.names = F)
