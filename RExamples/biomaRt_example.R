library(biomaRt)

atVHIR <- FALSE
if (atVHIR){
  http_proxy="http://conf_www.ir.vhebron.net:8080/"
  https_proxy="http://conf_www.ir.vhebron.net:8080/"
}

biodataset <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#antes de utilizar la funcion biomaRt hay que definir: filtros, atributos y valores.
filtros<-listFilters(biodataset)
attributos<- listAttributes(biodataset)
#fem la cerca
entrezfromgimsin <- getBM(filters= "affy_ragene_2_1_st_v1", 
                          attributes= c("entrezgene","affy_ragene_2_1_st_v1"),
                          values= affyID,
                          mart= biodataset,uniqueRows=TRUE)


###ESTO VIENE DE UN EJERCICIO DE RNASEQ PARA NORMALIZAR POR LONGITUD DEL GEN Y %GC
#buscamos el inicio y el final de cada gen y el %GC utilizando biomart
biodataset <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

entrezID <- as.character(rownames(counts.f))
annot <-getBM(filters = "entrezgene",
              attributes =c("hgnc_symbol","start_position",
                            "end_position","entrezgene",
                            "percentage_gc_content"),
              values=entrezID,
              mart=biodataset,
              uniqueRows=TRUE,
              checkFilters = TRUE)

### OTRO RECURSO:
# http://davetang.org/muse/2012/04/27/learning-to-use-biomart/
