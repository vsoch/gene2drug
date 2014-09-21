# This script will compare the similarity matrix from pubChem (based on structural similarity) to my clustering based
# on common enriched genes in gsea.

# SIMILARITY BASED ON STRUCTURE ---------------------------------------------------------------------------------
setwd("/home/vanessa/Documents/Dropbox/Code/R/gene2drug/data/cluster/2d")
load("/home/vanessa/Documents/Dropbox/Code/R/gene2drug/data/DrugGeneLists74Final.Rda")
tanimoto = read.table("tanimoto_cluster_reqid_1069925969951886533.csv",sep=",",skip=1)
rownames(tanimoto) = tanimoto$V1
tanimoto = tanimoto[,-c(1,73)]
colnames(tanimoto) = rownames(tanimoto)
save(tanimoto,file="tanimoto_cluster_reqid_1069925969951886533.Rda")

# Load lookup table
protdb = read.csv("/home/vanessa/Documents/Dropbox/Code/R/gene2drug/pubChem/protDBLookup.csv",sep="\t",head=TRUE)

# Now let's recreate tree
pubchem.dist = as.dist(tanimoto)
my.dist = as.dist(result$tanimotos)
pubchem.hc = hclust(pubchem.dist)
my.hc = hclust(my.dist)

# Fix labels to be drugs
pubchem.hc$labels = protdb$MEDICATION[match(pubchem.hc$labels,protdb$PUBCHEM)]
plot(pubchem.hc,main="PubChem Tanimoto Structure Similarity Clustering")
plot(my.hc,main="gene2Drug Tanimoto Genes Similarity Clustering")
# I'm not sure how to interpret - assess this - need to meet with someone who knows about drugs

# SIMILARITY BASED ON PROPERTIES ---------------------------------------------------------------------------------
load("../../drugProperties74.Rda")
numeric.properties = properties[,c(3,9,seq(10,40))]
numeric.properties.matrix  = array(0,dim(numeric.properties))
rownames(numeric.properties.matrix) = rownames(numeric.properties)
for (n in 1:nrow(numeric.properties)){
  numeric.properties.matrix[n,] = as.numeric(numeric.properties[n,])  
}

disty = dist(numeric.properties.matrix)
hc = hclust(disty)
plot(hc,main="Clustering of Drugs based on Chemical 33 Properties")