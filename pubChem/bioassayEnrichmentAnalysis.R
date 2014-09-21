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
  scores[[drug]] = pospoints - negpoints
  posscore[[drug]] = pospoints
  negscore[[drug]] = negpoints 
  nocalls[[drug]] = zeropoints
}

# Conclusion from above - most of these genes have no evidence either way :(

# GET MESH HIERARCHY FOR EACH DRUG