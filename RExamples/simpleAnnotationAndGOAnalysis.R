## ----options, include=FALSE----------------------------------------------
if(!require(knitr)) install.packages("knitr")
require(knitr)
opts_chunk$set(fig.path = 'images/grafic', tidy = FALSE, cache = FALSE)
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

## ----getGenes------------------------------------------------------------
topTab <- read.table("https://raw.githubusercontent.com/alexsanchezpla/Ejemplo_de_MDA_con_Bioconductor/master/results/ExpressAndTop_AvsB.csv2", head=TRUE, sep=";", dec=",", row.names = 1)
colnames(topTab)
head(topTab)

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
writeGeneLists<- FALSE
if(writeGeneLists){
  write.csv(geneListUp, file="selectedAvsB.up.csv")
  write.csv(geneListDown, file="selectedAvsB.down.csv")
  write.csv(geneUniverse, file="geneUniverse.csv")
}

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

## ----goProfiles1---------------------------------------------------------
require(goProfiles)
BPprofile1<- basicProfile(genelist=geneListUp, onto="BP", orgPackage="org.Hs.eg.db", empty.cats=FALSE, level=2)[[1]]
head(BPprofile1)

## ----GOANots-------------------------------------------------------------
require(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
entrezsUp2GO <- select(org.Hs.eg.db, keys = as.character(geneListUp), columns=c("SYMBOL", "GOALL"))
head(entrezsUp2GO)
entrezsUp2GOBP<- entrezsUp2GO[entrezsUp2GO$ONTOLOGY=="BP",]
BPprofileWithGenes<- cbind(BPprofile1, genes=rep("", nrow(BPprofile1)))
BPprofileWithGenes$genes<- as.character(BPprofileWithGenes$genes)
for (i in 1:nrow(BPprofile1)){
  GOIDi<- BPprofile1[i,"GOID"]
  genesi <-unique(entrezsUp2GOBP[entrezsUp2GOBP$GO==GOIDi,"ENTREZID"])
  genesi <- paste(genesi[!is.na(genesi)], collapse = " ")
  BPprofileWithGenes[i,"genes"]=genesi
}
head(BPprofileWithGenes)

