# This script will determine if associations between a group of genes associated with a drug have been 
# found in the literature, by way of bioassay results.  We will do the following:

# parse the bioassay results for each drug - I want a matrix of bioassay proteins by drugs, and some score? in the cell to indicate the association. OR I could just make a list of bioassays for each drug.
# associate each bioassay with a gene.  I should take the unique list of the above, assign IDs, genes, and then recreate the above using those unique identifiers.
# for each drug, rank order the list of GSEA genes (use NES score?)
# rank order the proteins for each bioassay (each one will have assigned to it a gene)
# figure out the GSEA algorithm, and implement it for this case
# do the walk to calculate the enrichment score!
# Permutations!
# The enrichment score will give a significance for the whole set, but I can also just do intersect to see if I have found any known interactions (and calculate a percentage?)

library(rjson)
setwd("/home/vanessa/Documents/Dropbox/Code/R/gene2drug/data/json/bioassay")

# Load our drug data
load("/home/vanessa/Documents/Dropbox/Code/R/gene2drug/data/DrugGeneLists74Final.Rda")

# First we will create a matrix of unique bioassays - the unique id is the protein (target) id - the GID
# We will then parse the results again and catalog which targets are associated with each drug

# Let's get list of AIDs
ba = list.files(pattern="*_bioassay.json")
gids = c()
for (f in 1:length(ba)){
  file = ba[f]
  cat("Starting",f,"of",length(ba),"\n")
  json_data <- fromJSON(paste(readLines(file), collapse=""))
  for (row in 1:length(json_data$Table$Row)) {
    gid = json_data$Table$Row[[row]][1]$Cell[8]
    gids = c(gids,gid)
  }
}
gids = unique(gids)
cat("Found",length(gids),"associated with",length(ba),"drugs.","\n")

# Prepare matrix to hold assay info
drug.names = gsub("_bioassay.json","",ba)

# This will be a binary matrix of associations -
# 1: active
#-1: inactive
# 0 : anything else
associations.binary = array(0,dim=c(length(gids),length(ba)))
colnames(associations.binary) = sort(drug.names)
rownames(associations.binary) = sort(gids)

# This will be a matrix of Activity Value [uM]

# Here is the info that we have 
# [1] AID
# [2] AID Version
# [3] AID Revision
# [4] Panel Member ID
# [5] SID
# [6] CID
# [7] Bioactivity Outcome
# [8] Target GI
# [9] Activity Value [uM]
# [10] Activity Name
# [11] Assay Name
# [12] Bioassay Type
# [13] PubMed ID
# [14] RNAi
# [15] Gene Target if RNAi

for (file in ba){
  json_data = fromJSON(paste(readLines(file), collapse=""))
  drug = gsub("_bioassay.json","",file)
  for (row in 1:length(json_data$Table$Row)) {
    outcome = json_data$Table$Row[[row]]$Cell[7]
    gid = json_data$Table$Row[[row]]$Cell[8]  
    if (!is.null(gid)){
       if (gid != ""){
          if (outcome == "Inactive"){
            associations.binary[gid,drug] = -1
          }
          if (outcome == "Active"){
            associations.binary[gid,drug] = 1
          }
          # If inactive, code as -1
          # If active, code as 1
          # If not conclusive, already coded as 0
       }
    }
  }
}

associations.binary = associations.binary[-1,]
save(associations.binary,file="../../associated_gidToDrug72.Rda")

# Vanessa: if this method turns out good, turn this into an R package!
# Write a function to match gid to gene from our data
giToGene <- function(data){

  require(XML)
  
  # Here we will keep the ids to gene mapping
  gidlookup = list()

  for (i in 1:length(data)){
    for (j in 1:length(data[i]$GBSeq["GBSeq_other-seqids"][[1]])) {
      contender = as.character(data[i]$GBSeq["GBSeq_other-seqids"][[1]][j])
      contender = strsplit(contender,"[|]")[[1]]
      if (contender[1]=="gi"){
        gi = contender[2]
      }
    }
    
    # Now find the gene name
    for (j in 1:length(data[i]$GBSeq["GBSeq_feature-table"][[1]])){
      contender = data[i]$GBSeq["GBSeq_feature-table"][[1]][j]
      idx = which(names(contender[[1]])=="GBFeature_quals")
      contender[[1]][idx]["GBFeature_quals"]
      for (k in 1:length(contender[[1]][idx]["GBFeature_quals"][[1]])){
        tmp = as.character(contender[[1]][idx]["GBFeature_quals"][[1]][k][[1]])
        if (tmp[1] == "gene"){
          gene = tmp[2]
        }
      }
    }
    gidlookup[gi] = gene
    cat(gi,":",gene,"\n")
  }
  return(gidlookup)

}

# Now we want to get a gene for each target (gid) using our function
load("/home/vanessa/Documents/Dropbox/Code/R/gene2drug/data/associated_gidToDrug72.Rda")
gidsearch = paste(rownames(associations.binary)[1:899],collapse=",")
query = paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&retmode=xml&id=",gidsearch,sep="")
result = xmlParse(query)
data <- xmlToList(result)
gidlist1 = giToGene(data)
# Do in two parts so query is <1000 each time
gidsearch = paste(rownames(associations.binary)[900:length(rownames(associations.binary))],collapse=",")
query = paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&retmode=xml&id=",gidsearch,sep="")
result = xmlParse(query)
data <- xmlToList(result)
gidlist2 = giToGene(data)

# Combine the two lists
gidlookup = c(gidlist1,gidlist2)
setwd("/home/vanessa/Documents/Dropbox/Code/R/gene2drug/data")
save(gidlookup,file="giToGeneLookup1550.Rda")

# ENRICHMENT ANALYSIS ------------------------------------------------------------------------------

# Load our drug result again
load("/home/vanessa/Documents/Dropbox/Code/R/gene2drug/data/DrugGeneLists74Final.Rda")

# Now for each medication, let's make a list of active and inactive genes
# We should award points for active genes, take off for inactive, do nothing for neither
scores = list()
posscore = list()
negscore = list()
nocalls = list()
overlap.genes = list()
for (i in 1:ncol(associations.binary)){
  drug = colnames(associations.binary)[i]
  active = names(which(associations.binary[,i] == 1))
  inactive = names(which(associations.binary[,i] == -1))
  nocall = names(which(associations.binary[,i] == 0))
  
  # Now get gene names for all
  genes.active = as.character(gidlookup[active])
  genes.inactive = as.character(gidlookup[inactive])
  genes.nocall = as.character(gidlookup[nocall])
  
  # Now compare to our list
  mygenes = strsplit(result$geneLists[which(result$meds == drug)][[1]]," ")[[1]]
  pospoints = length(which(mygenes %in% genes.active))
  negpoints = length(which(mygenes %in% genes.inactive))
  zeropoints = length(which(mygenes %in% genes.nocall))
  overlap.genes[[drug]] = mygenes[which(mygenes %in% genes.active)]
  scores[[drug]] = pospoints - negpoints
  posscore[[drug]] = pospoints
  negscore[[drug]] = negpoints 
  nocalls[[drug]] = zeropoints
}

# Conclusion from above - most of these genes have no evidence either way :(
# NOTE: The drugs with active genes found above, eg:

# chrysin: "CALM1"
# novobiocin: "ALB"
# sulconazole: "NR1H4"
# tanespimycin: "HSP90B1"

# are exactly the same for this analysis and the analysis below when
# I map from genes --> pathways --> more genes first!
# ATTEMPT # 2 --- looking for pathway enrichment ---------------------------------------------------------
# It could be the case that the differences in RNA expression are due to indirect changes in the pathway, and
# not necessarily the exact gene.  To test this, I will map the genes enriched in my drugs to proteins, and then 
# the proteins to pathways.

# IDEA 2:
# 1) Map drug --> proteins --> pathways: enrichment of pathways (look for common pathways)
# 2) Map GSEA genes --> pathways 
# 3) Find some way to assess if they are related (overlap)

library("KEGGREST")
library(RCurl)
options(RCurlOptions = list(cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")))

setwd("/home/vanessa/Documents/Dropbox/Code/R/gene2drug/data/json/bioassay")

# Load our drug data
load("/home/vanessa/Documents/Dropbox/Code/R/gene2drug/data/DrugGeneLists74Final.Rda")
names(result$geneLists) = result$meds

# Just get all hsa pathways in one go
kegg = keggList("hsa")

# LOOKUP by hsa ID
hsa2genes = list()
hsa2path = list()

# Lookup by gene
gene2hsa = list()
gene2path = list()

# Parse kegg result - create list of genes with hsaid
for (k in 1:length(kegg)){
  cat("Processing",k,"of",length(kegg),"\n")
  hsa = names(kegg[k])
  genes = strsplit(kegg[k],";")[[1]][1]
  pathway = gsub("^ ","",strsplit(kegg[k],";")[[1]][2])
  hsa2genes[[hsa]] = genes
  hsa2path[[hsa]] = pathway
  genes = gsub(" ","",strsplit(genes,",")[[1]])
  for (g in genes){
    if (g %in% names(gene2hsa)) {
      holder = gene2hsa[[g]]
      gene2hsa[[g]] = c(holder,hsa)
    } else {
      gene2hsa[[g]] = hsa
    }
    if (g %in% names(gene2path)) {
      holder = gene2path[[g]]
      gene2path[[g]] = c(holder,pathway)
      } else {
      gene2path[[g]] = pathway
      }
  }
}

# Save our kegg objects to file
path.hsa2hsa = keggLink("pathway", "hsa")
hsa2path.hsa = keggLink("hsa", "pathway")
readme = "KEGG database with 30739 entries, downloaded 9/22/2014/n gene2hsa: maps gene names to hsa identifiers \n gene2path: maps gene symbols to pathway descriptions \n hsa2gene: maps hsa pathway identifiers to gene symbols \n hsa2path: maps hsa identifiers to pathways. \n questions: email vsochat@stanford.edu"
kegg = list(gene2hsa = gene2hsa,gene2path=gene2path,hsa2gene=hsa2genes,hsa2path=hsa2path,readme=readme,pathhsa2hsa=path.hsa2hsa,hsa2pathhsa=hsa2path.hsa)
save(kegg,file="/home/vanessa/Documents/Dropbox/Code/R/gene2drug/data/kegg.Rda")

# Now for each drug, create list of hsa identifiers
hsas = list()
for (g in 1:length(result$geneLists)){
  med = names(result$geneLists[g])
  # Here are the genes
  tmp = result$geneLists[[g]]
  tmp = strsplit(tmp," ")[[1]]
  hsa = list()
  for (t in tmp){
    hsa[[t]] = kegg$gene2hsa[[t]]
  }
  hsas[[med]] = hsa
}

# Now for each med, get a list of pathways
pathwaylookup = keggLink("pathway", "hsa")
pathways = list()
for (h in 1:length(hsas)){
  med = names(hsas[h])
  tmp = hsas[[med]]
  allhsa = c()
  for (t in 1:length(tmp)){
    if (length(tmp[[t]]) == 1) {
      allhsa = c(allhsa,as.character(tmp[t]))
    } else {
      for (tt in tmp[[t]]){
      allhsa = c(allhsa,as.character(tmp[[t]][tt]))
     }
    }
  }
  allhsa = allhsa[-which(is.na(allhsa))]
  pathways[[med]] = allhsa
}
  
# Save both to our result object
result$pathways.kegg = pathways
result$hsa = hsas
save(result,file="DrugGeneLists74Final.Rda")

# Now we have complete pathways for enriched genes for each med
# Now we need to get ALL the possible genes based on the pathways!

expandedgenes = list()
for (e in 1:length(result$pathways.kegg)){
  med = names(result$pathways.kegg[e])
  hsa = result$pathways.kegg[[med]]
  tmp = as.character(kegg$hsa2gene[hsa])
  tmp2 = c()
  for (t in tmp){
    tmp2 = c(tmp2,gsub(" ","",strsplit(t,",")[[1]]))  
  }
  tmp = unique(c(tmp2,strsplit(result$geneLists[[med]]," ")[[1]]))
  expandedgenes[[med]] = tmp
}

# Save to result object
result$expanded.pathway.genes = expandedgenes
save(result,file="DrugGeneLists74Final.Rda")

# NOW we should look for overlapping genes with association.
# Now for each medication, let's make a list of active and inactive genes
# We should award points for active genes, take off for inactive, do nothing for neither

# Load the associations from file, associations.binary
load("associated_gidToDrug72.Rda")
load("giToGeneLookup1550.Rda")

activeGenes = list()
for (i in 1:ncol(associations.binary)){
  drug = colnames(associations.binary)[i]
  active = names(which(associations.binary[,i] == 1))
  
  # Now get gene names for all
  genes.active = as.character(gidlookup[active])
  
  # Now compare to our list
  mygenes = result$expanded.pathway.genes[[drug]]
  ag = mygenes[which(mygenes %in% genes.active)]
  activeGenes[[drug]] = ag
  cat("Found",length(ag),"active genes for",drug,"\n")
}
activeGenes = activeGenes[which(as.character(activeGenes)!="character(0)")]
result$activeGenes = activeGenes
save(result,file="DrugGeneLists74Final.Rda")

# The result is equivalent to above:

# chrysin: "CALM1"
# novobiocin: "ALB"
# sulconazole: "NR1H4"
# tanespimycin: "HSP90B1"

# Next I will look up these drug/gene pairs, try to understand the pathway/how interact, and create gene lists for those pathways to possibly test against data.
# TODO: GET MESH HIERARCHY FOR EACH DRUG