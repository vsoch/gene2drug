# This script will parse GSEA Drug results,
inputfile = "/home/vanessa/Documents/Work/GENE_EXPRESSION/cmap/gsea/DRUGS_GSEA_FDRpt05_RMA_Filter.Rda" 
load(inputfile)

# IDEA 0 - we want a matrix of TERM_gene (rows) and medication (columns), and the value inside each to be the
# enrichment score

make list of drugs, for each drug, include list of core genes (how to deal with overlap?)
use genes as "features" of drug - put 1 if gene is sig. enriched for the drug, 0 if not
then do some correlation... 

# IDEA 1 - For all significant gene sets, download NES for ALL significant terms for subset of enriched drugs, do GEPHI image
# First let's create a matrix of medications (rows), and gene sets (columns), and
# we want to see if we can cluster.
tmp = as.character(reportFinal$NAME)
tmp= sapply(tmp,strsplit,"_")

medications = c()
for (t in tmp){
  medications = c(medications,t[[1]][1])
}  

result = cbind(medications,reportFinal)
colnames(result)[2] = "FOLDER"
colnames(result)[3] = "GENESET"
save(result,file=paste("/home/vanessa/Documents/Work/GENE_EXPRESSION/cmap/gsea/DRUGS_GSEA_FDRpt05_RMA_Filter_Final.Rda",sep=""))

# Matrix - drugs in rows, gene sets in columns
geneset = as.character(unique(result$GENESET))

df = array(data=0,dim=c(length(unique(result$medications)),length(geneset)))
rownames(df) = as.character(sort(unique(result$medications)))
colnames(df) = as.character(geneset)

# Fill in the NES score
for (r in 1:nrow(result)) {
  cat(r,"of",nrow(result),"\n")
  df[as.character(result$medications[r]),as.character(result$GENESET[r])] = result$NES[r]
}

# Try maude's method...
# IDEA 2 - create matrix of term by drug, put NES score in matrix (is it too sparse?)
# IDEA 3 -  find drugs with common gene sets - calculate correlation between
# enrichment scores, GEPHI to visualize similarity
