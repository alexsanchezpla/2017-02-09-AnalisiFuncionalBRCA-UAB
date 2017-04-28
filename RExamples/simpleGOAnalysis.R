# Rpackages
source("http://bioconductor.org/biocLite.R")
if (!(require(annotate))) biocLite("annotate")
if (!(require(GOstats))) biocLite("GOstats")
if (!(require(org.Hs.eg.db))) biocLite("org.Hs.eg.db")
# Data
## Read data
topTab <- read.table("https://raw.githubusercontent.com/alexsanchezpla/scripts/master/Exemple_Analisis_BioC/results/ExpressAndTop_AvsB.csv2",
                     head=TRUE, sep=";", dec=",")
colnames(topTab)
head(topTab)
geneListUp <- unique(topTab$EntrezsA [topTab$adj.P.Val<0.05 & topTab$logFC > 0] )
length(geneListUp)
geneListDown <- unique(topTab$EntrezsA [topTab$adj.P.Val<0.05 & topTab$logFC < 0] )
length(geneListDown)
geneUniverse <- unique(topTab$EntrezsA)
length(geneUniverse)
write.csv(geneListUp, file="selectedAvsB.up.csv")
write.csv(geneListDown, file="selectedAvsB.down.csv")
write.csv(geneUniverse, file="geneUniverse.csv")

### Symbols to Entrez
# If we hadn't had 'Entrez' Identifiers, but only 'Gene Symbols' we would have done as follows:
require(org.Hs.eg.db)
geneSymbolsUp <- unique(topTab$SymbolsA [topTab$adj.P.Val<0.05 & topTab$logFC > 0] )
geneSymbolsDown <- unique(topTab$SymbolsA [topTab$adj.P.Val<0.05 & topTab$logFC < 0] )
geneSymbolsUniverse <- unique(topTab$SymbolsA)

require(org.Hs.eg.db)
#keytypes(org.Hs.eg.db)
geneEntrezsUp <- select(org.Hs.eg.db, keys = as.character(geneSymbolsUp), columns=c("ENTREZID"), keytype="SYMBOL")$ENTREZID
geneEntrezsDown <- select(org.Hs.eg.db, keys = as.character(geneSymbolsDown), columns=c("ENTREZID"), keytype="SYMBOL")$ENTREZID
geneEntrezsUniverse <- select(org.Hs.eg.db, keys = as.character(geneSymbolsUniverse), columns=c("ENTREZID"), keytype="SYMBOL")$ENTREZID

# Remove potential NA's values
geneEntrezsUp <- geneEntrezsUp[!is.na(geneEntrezsUp)]
geneEntrezsDown <- geneEntrezsDown[!is.na(geneEntrezsDown)]
geneEntrezsUniverse <- geneEntrezsUniverse[!is.na(geneEntrezsUniverse)]

# GOAnalysis for "up regulated" genes

require(GOstats)
entrezUniverse <-  geneEntrezsUniverse
geneIds <- geneEntrezsUp

## Creamos los "hiperparametros" en que se basa el analisis
GOparams = new("GOHyperGParams",
               geneIds=geneIds, universeGeneIds=entrezUniverse,
               annotation="org.Hs.eg.db", ontology="BP",
               pvalueCutoff=0.001, conditional=FALSE,
               testDirection="over")
KEGGparams = new("KEGGHyperGParams",
                 geneIds=geneIds, universeGeneIds=entrezUniverse,
                 annotation="org.Hs.eg.db",  
                 pvalueCutoff=0.01, testDirection="over")
  
## Ejecutamos los analisis
  
GOhyper = hyperGTest(GOparams)
KEGGhyper = hyperGTest(KEGGparams)
cat("GO\n")
print(head(summary(GOhyper)))
cat("KEGG\n")
print(head(summary(KEGGhyper)))
  
# Creamos un informe html con los resultados
GOfilename =file.path(paste("GOResults.AvsB.up",".html", sep=""))
KEGGfilename =file.path(paste("KEGGResults.AvsB.up",".html", sep=""))
htmlReport(GOhyper, file = GOfilename, summary.args=list("htmlLinks"=TRUE))
htmlReport(KEGGhyper, file = KEGGfilename, summary.args=list("htmlLinks"=TRUE))

