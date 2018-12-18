rm(list=ls())
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

## setwd, load libraries, & load data ####
setwd("/home/vassil/Documents/Bcells/Meta-analysis")
# source ("http://www.bioconductor.org/biocLite.R")
# biocLite("clusterProfiler") # also choose to update all packages
# biocLite("pathview")
library(clusterProfiler)
# install.packages("GOplot")
# library(GOplot)
library(pathview)

load("/home/vassil/Documents/Bcells/BioInfo/Meta-analysis/masterTableDown")
load("/home/vassil/Documents/Bcells/BioInfo/Meta-analysis/masterTableUp")
tabUp <- masterTableUp
tabDown <- masterTableDown
rm (list = c("masterTableUp", "masterTableDown"))

## Add cell type to tables ####
# up:
epith <- rep("epith", length(grep("^E.", tabUp$dataset)))
fibr <- rep("fibr", length(grep("^F.", tabUp$dataset)))
b <- rep("B", length(grep("^L.B.", tabUp$dataset)))
pbmc <- rep("pbmc", length(grep("^L.Pbmc.", tabUp$dataset)))
t <- rep("T", length(grep("^L.T.", tabUp$dataset)))
mono <- rep("mono", length(grep("^M.APC", tabUp$dataset)))
gran <- rep("gran", length(grep("^M.Gran", tabUp$dataset)))
gran2 <- rep("gran", length(grep("^M.Neutr", tabUp$dataset)))
musc <- rep("musc", length(grep("^Musc.", tabUp$dataset)))
hsc <- rep("hsc", length(grep("^O.HSC.", tabUp$dataset)))
test <- rep("test", length(grep("^O.Test.", tabUp$dataset)))
cellType <- c(epith,fibr,b,pbmc,t,mono,gran,gran2,musc,hsc,test)
tabUp$cellType <- cellType
# down:
epith <- rep("epith", length(grep("^E.", tabDown$dataset)))
fibr <- rep("fibr", length(grep("^F.", tabDown$dataset)))
b <- rep("B", length(grep("^L.B.", tabDown$dataset)))
pbmc <- rep("pbmc", length(grep("^L.Pbmc.", tabDown$dataset)))
t <- rep("T", length(grep("^L.T.", tabDown$dataset)))
mono <- rep("mono", length(grep("^M.APC", tabDown$dataset)))
gran <- rep("gran", length(grep("^M.Gran", tabDown$dataset)))
gran2 <- rep("gran", length(grep("^M.Neutr", tabDown$dataset)))
musc <- rep("musc", length(grep("^Musc.", tabDown$dataset)))
hsc <- rep("hsc", length(grep("^O.HSC.", tabDown$dataset)))
test <- rep("test", length(grep("^O.Test.", tabDown$dataset)))
cellType <- c(epith,fibr,b,pbmc,t,mono,gran,gran2,musc,hsc,test)
tabDown$cellType <- cellType

# clean up
rm (list = c("epith","fibr","gran","gran2","hsc","mono",
             "musc","pbmc","t","test","b","cellType"))

## Remove opposite direction genes per cell type ####

# up:
tabUp.u <- data.frame(matrix (ncol=length(names(tabUp)), # define holder 
                              nrow = 0)) # dataframe for unique upregulated genes
names(tabUp.u) <- names(tabUp)
cellsup <- unique(tabUp$cellType) # all cell types to go through
for (i in cellsup) {
  up <- uniqGeneTabUp (tableUp=tabUp,
                       tableDown=tabDown,
                       celltype=i)
  tabUp.u <- rbind(tabUp.u, up)
  rm(up)
}
rm(i)
rm(cellsup)

# down:
tabDown.u <- data.frame(matrix (ncol=length(names(tabDown)), # define holder 
                                nrow = 0)) # dataframe for unique upregulated genes
names(tabDown.u) <- names(tabDown)
cellsdown <- unique(tabDown$cellType) # all cell types to go through
for (i in cellsdown) {
  down <- uniqGeneTabDown (tableUp=tabUp,
                           tableDown=tabDown,
                           celltype=i)
  tabDown.u <- rbind(tabDown.u, down)
  rm(down)
}
rm(i)
rm(cellsdown)

save (tabUp.u, file = "tabUp.u")
save (tabDown.u, file = "tabDown.u")

