## ----options, include=FALSE----------------------------------------------
if(!require(knitr)) install.packages("knitr")
require(knitr)
opts_chunk$set(fig.path = 'images/grafic', tidy = FALSE, cache = FALSE, message = TRUE, echo = TRUE, print = TRUE) 
options(width=80, warn=0, digits=4)
options(scipen=100)
options(java.parameters = "-Xmx4g" )

## ----librerias, echo=FALSE, message=FALSE, warning=FALSE, results='hide'----
installifnot <- function(pckgName){
 if (!(require(pckgName, character.only = TRUE))) {
 install.packages(pckgName)
  }else{
    print(paste("Package", pckgName, "already installed", sep = " "))
  }
}
installBiocifnot <- function(pckgName){
 if (!(require(pckgName, character.only = TRUE))) {
    source("https://bioconductor.org/biocLite.R")
    biocLite(eval(pckgName), suppressUpdates = TRUE)
  }else{
    print(paste("Package", pckgName, "already installed", sep = " "))
  }
}

installBiocifnot("annotate")
installBiocifnot("GOstats")
installBiocifnot("org.Hs.eg.db")
installBiocifnot("hgu133a.db")
installBiocifnot("biomaRt")

## ----printr, eval=FALSE--------------------------------------------------
## if(!(require(printr))) {
##   install.packages(
##     'printr',
##     type = 'source',
##     repos = c('http://yihui.name/xran', 'http://cran.rstudio.com')
##   )
## }

## ----getGenes------------------------------------------------------------
topTab <- read.table("https://raw.githubusercontent.com/alexsanchezpla/scripts/master/Exemple_Analisis_BioC/results/ExpressAndTop_AvsB.csv2", head=TRUE, sep=";", dec=",", row.names = 1)
colnames(topTab)
head(topTab)

## ----setProxy------------------------------------------------------------
atVHIR <- FALSE
if (atVHIR){
    http_proxy="http://conf_www.ir.vhebron.net:8080/"
    https_proxy="http://conf_www.ir.vhebron.net:8080/"
    }

## ----annotateFromArrayPackage--------------------------------------------
probeIDsAll <- rownames(topTab)
probeIDsUp <- probeIDsAll [topTab$adj.P.Val<0.05 & topTab$logFC > 0] 
probeIDsDown <- probeIDsAll [topTab$adj.P.Val<0.05 & topTab$logFC < 0]

require(hgu133a.db)
keytypes(hgu133a.db)

geneEntrezsUp <- select(hgu133a.db, keys = probeIDsUp, columns=c("ENTREZID", "SYMBOL"))
geneEntrezsDown <- select(hgu133a.db, keys = probeIDsUp, columns=c("ENTREZID", "SYMBOL"))
geneEntrezsUniverse <- select(hgu133a.db, keys = probeIDsAll, columns=c("ENTREZID", "SYMBOL"))

head(geneEntrezsUp)

## ----annotateWithBiomart-------------------------------------------------
biodataset <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
listDatasets(biodataset)$dataset

filters<-listFilters(biodataset)
# We need to find the filter to link with Affymetrx arrays hgu133a
u133aFilters<- grep("u133a", filters[,1] )
u133aFilters <- filters[u133aFilters,]
myu133aFilter <- u133aFilters[3,1]
myu133aFilter

atributs<- listAttributes(biodataset)
entrezAtributs<- grep("entrez", atributs[,1])
entrezAtribut <- atributs[entrezAtributs,]
myentrezAtribut <- entrezAtribut[2,1]
myentrezAtribut

# Now we can do the search
entrezfromProbesUp <- getBM(filters= myu133aFilter, 
                          attributes= c(myentrezAtribut, myu133aFilter),
                          values= probeIDsUp,
                          mart= biodataset,uniqueRows=TRUE)
head(entrezfromProbesUp)

## ----geneLists-----------------------------------------------------------
geneListUp <- topTab$EntrezsA [topTab$adj.P.Val<0.05 & topTab$logFC > 0]
head(geneListUp)
geneListDown <- topTab$EntrezsA [topTab$adj.P.Val<0.05 & topTab$logFC < 0]
length(geneListDown)
geneUniverse <- topTab$EntrezsA
length(geneUniverse)
write.csv(geneListUp, file="selectedAvsB.up.csv")
write.csv(geneListDown, file="selectedAvsB.down.csv")
write.csv(geneUniverse, file="geneUniverse.csv")

## ----prepareEntrezs------------------------------------------------------
# Remove potential NA's values
geneEntrezsUp <- unique(geneListUp[!is.na(geneListUp)])
geneEntrezsDown <- unique(geneListDown[!is.na(geneListDown)])
geneEntrezsUniverse <- unique(geneUniverse[!is.na(geneUniverse)])

## ----GOAnalysis1---------------------------------------------------------
require(GOstats)
## Creamos los "hiperparametros" en que se basa el analisis
GOparams = new("GOHyperGParams",
               geneIds=geneEntrezsUp, universeGeneIds=geneEntrezsUniverse,
               annotation="org.Hs.eg.db", # might have use hgu133a.db instead
               ontology="BP",
               pvalueCutoff=0.001, conditional=FALSE,
               testDirection="over")
KEGGparams = new("KEGGHyperGParams",
                 geneIds=geneEntrezsUp, universeGeneIds=geneEntrezsUniverse,
                 annotation="org.Hs.eg.db", # might have use hgu133a.db instead 
                 pvalueCutoff=0.01, testDirection="over")

## ----GOAnalysis2---------------------------------------------------------
GOhyper = hyperGTest(GOparams)
KEGGhyper = hyperGTest(KEGGparams)
cat("GO\n")
print(head(summary(GOhyper)))
cat("KEGG\n")
print(head(summary(KEGGhyper)))

## ----GOAnalysis3---------------------------------------------------------
# Creamos un informe html con los resultados
GOfilename =file.path(paste("GOResults.AvsB.up",".html", sep=""))
KEGGfilename =file.path(paste("KEGGResults.AvsB.up",".html", sep=""))
htmlReport(GOhyper, file = GOfilename, summary.args=list("htmlLinks"=TRUE))
htmlReport(KEGGhyper, file = KEGGfilename, summary.args=list("htmlLinks"=TRUE))

