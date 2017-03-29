# Rpackages
source("http://bioconductor.org/biocLite.R")
if (!(require(annotate))) biocLite("annotate")
if (!(require(GOstats))) biocLite("GOstats")
if (!(require(org.Mm.eg.db))) biocLite("org.Mm.eg.db")
# Data
## Read data
topTab <- read.table("http://ueb.vhir.org/tiki-download_file.php?fileId=2689", 
                     head=TRUE, sep=",", dec=".")
expres <- read.table("http://ueb.vhir.org/tiki-download_file.php?fileId=2690", 
                     head=TRUE, sep="\t", dec=".")
## Unique identifiers
geneList <- unique(as.character(topTab$SYMBOL)); length(geneList)
backgrd <- unique(as.character(expres$SYMBOL)); length(backgrd)
## Symbols2Entrezs
require(org.Mm.eg.db)
geneListEntrezs <- unlist(mget(geneList, org.Mm.egSYMBOL2EG, ifnotfound=NA))
length(!is.na(geneListEntrezs))
backgrdEntrezs <-  unlist(mget(backgrd, org.Mm.egSYMBOL2EG, ifnotfound=NA)) 
length(!is.na(backgrdEntrezs))
geneListEntrezs <- geneListEntrezs[!is.na(geneListEntrezs)]
backgrdEntrezs <- backgrdEntrezs[!is.na(backgrdEntrezs)]

# GOAnalysis
require(GOstats)
entrezUniverse <-  backgrdEntrezs
geneIds <- geneListEntrezs

## Creamos los "hiperparametros" en que se basa el analisis
GOparams = new("GOHyperGParams",
               geneIds=geneIds, universeGeneIds=entrezUniverse,
               annotation="org.Mm.eg.db", ontology="BP",
               pvalueCutoff=0.001, conditional=FALSE,
               testDirection="over")
KEGGparams = new("KEGGHyperGParams",
                 geneIds=geneIds, universeGeneIds=entrezUniverse,
                 annotation="org.Mm.eg.db",  
                 pvalueCutoff=0.01, testDirection="over")
  
## Ejecutamos los analisis
  
GOhyper = hyperGTest(GOparams)
KEGGhyper = hyperGTest(KEGGparams)
cat("GO\n")
print(head(summary(GOhyper)))
cat("KEGG\n")
print(head(summary(KEGGhyper)))
  
# Creamos un informe html con los resultados
GOfilename =file.path(paste("GOResults.",".html", sep=""))
KEGGfilename =file.path(paste("KEGGResults.",".html", sep=""))
htmlReport(GOhyper, file = GOfilename, summary.args=list("htmlLinks"=TRUE))
htmlReport(KEGGhyper, file = KEGGfilename, summary.args=list("htmlLinks"=TRUE))

