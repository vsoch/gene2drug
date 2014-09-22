setwd("/home/vanessa/Documents/Work/GENE_EXPRESSION/drugs")
list.files()

# Here is the drug result, FDR .05, for 74 drugs
# The gene lists are only for CORE genes of each set, combined across brainterms
load("DrugGeneLists74.Rda")

# Create a matrix of genes by drugs
uniquegenes = c()
for (c in 1:length(result$geneLists)){
  uniquegenes = c(uniquegenes,strsplit(result$geneLists[[c]]," ")[[1]])
}
uniquegenes = unique(uniquegenes)
genematrix = array(0,dim=c(length(result$geneLists),length(uniquegenes)))
rownames(genematrix) = sort(result$meds)
colnames(genematrix) = sort(uniquegenes)

# Now go through data, make a matrix of core genes by medications
for (g in 1:length(result$geneLists)){
  genes = strsplit(result$geneLists[[g]]," ")[[1]]
  genematrix[result$meds[g],genes] = 1
}

result$genematrix = genematrix
save(result,file="DrugGeneLists74Final.Rda")

# Let's also calculate tanimoto scores for the gene sets
tanimotos = array(0,dim=c(length(result$geneLists),length(result$geneLists)))
rownames(tanimotos) = sort(result$meds)
colnames(tanimotos) = sort(result$meds)

# Now go through data, make a matrix of core genes by medications
for (g in 1:length(result$meds)){
  med1 = result$meds[g] 
  med1genes = strsplit(result$geneLists[[g]]," ")[[1]]
  for (c in 1:length(result$meds)){
    med2 = result$meds[c]
    med2genes = strsplit(result$geneLists[[c]]," ")[[1]]
    tanimotos[med1,med2] = length(intersect(med1genes,med2genes))/length(union(med1genes,med2genes))
  }
}

result$tanimotos = tanimotos
save(result,file="DrugGeneLists74Final.Rda")

disty = as.dist(tanimotos)
hc = hclust(disty)
plot(hc,main="Drug Similarity Based on Tanimoto Scores of Core Gene Sets")

# Now let's try loading the SOM, and coloring map for drug gene sets
load("brainMap.Rda")
colorscale = brewer.pal(9,"YlOrRd")
colorscale = colorRampPalette(brewer.pal(8,"YlOrRd"))(100)

# This is match scores for one compobent to all images - the range is the max match score for all images
# For each drug, get the unique list of terms
medications = c()
for (m in 1:nrow(result$report)){
  medications = c(medications, strsplit(as.character(result$report$NAME[m]),"_")[[1]][1])
}
result$report = result$report[-14]
result$report = cbind(result$report,medications)
save(result,file="DrugGeneLists74Final.Rda")

upterms = list()
downterms = list()
allterms = list()
for (r in result$meds){
  idx = which(result$report$medications %in% r)
  tmp = (tolower(as.character(result$report$NAME.1[idx])))
  up = c()
  down = c()
  all = c()
  for (t in 1:length(tmp)){
    if (strsplit(tmp[t],"_")[[1]][2] == "down"){
      down = c(down,strsplit(tmp[t],"_")[[1]][1])
    } else {
      up = c(up,strsplit(tmp[t],"_")[[1]][1])
    }
    all = c(all,strsplit(tmp[t],"_")[[1]][1])
  }
  upterms[[r]] = unique(up)
  downterms[[r]] = unique(down)
  allterms[[r]] = unique(all)
}

terms = list(up=upterms,down=downterms,all=allterms)
result$terms = terms
save(result,file="DrugGeneLists74Final.Rda")

# Now let's map to som, for each of up, down, and all
# First create a lookup for index in the SOM
labelidx = list()
for (l in 1:length(brainMap$labels)){
  tmp = strsplit(brainMap$labels[l],"\n")[[1]]
  for (t in tmp){
    labelidx[[t]] = l
  }
}

brainMap$idx = labelidx
save(brainMap,file="brainMap.Rda")

# Now create a map for each drug!
for (t in 1:length(terms$all)){
  med = names(terms$all[t])
  ups = terms$up[[med]]
  downs = terms$down[[med]]
  # If both aren't null
  if (!is.null(ups) && !is.null(downs)){
    both = intersect(ups,downs)
    ups = setdiff(both,ups)
    downs = setdiff(both,downs)
  } else {
    both = ""
  }
  
  # Now find the index of the nodes
  upidx = unique(as.numeric(brainMap$idx[which(names(brainMap$idx) %in% ups)]))
  downidx = unique(as.numeric(brainMap$idx[which(names(brainMap$idx) %in% downs)]))
  bothidx = unique(as.numeric(brainMap$idx[which(names(brainMap$idx) %in% both)]))
  
  # Create a vector of colors
  colors = array("#FFFFFF",dim=506)
  colors[upidx] = "#00FF00"
  colors[downidx] = "#FF0040"
  colors[bothidx] ="#FE9A2E"

  png(paste("img/",med,".png",sep=""),width=1000,height=800)
  plot(brainMap$som$grid$pts,main=med,col=colors,xlab="Nodes",ylab="Nodes",pch=15,cex=7.5)
  text(brainMap$som$grid$pts,brainMap$labels,cex=.6)
  dev.off()  
}

# Count number of up and down genes
tmp = tolower(as.character(unique(result$report$NAME.1))
upcount = 0
downcount = 0
for (t in 1:length(tmp)) {
  if(strsplit(tmp[t],"_")[[1]][2]=="down"){
    downcount = downcount +1
  }else{
    upcount = upcount+1
  }
}

# Finally, calculate tanimoto scores based on overlap of gene sets
# Get unique terms
uniqueterms = as.character(unique(result$report$NAME.1))
tanimotos.terms = array(dim=c(length(result$meds),length(result$meds)))
rownames(tanimotos.terms) = sort(result$meds)
colnames(tanimotos.terms) = sort(result$meds)

for (t in 1:length(result$meds)){
  med1 = result$meds[t]
  sets1 = as.character(result$report$NAME.1[which(result$report$medications ==med1)])
  for (u in 1:length(result$meds)) {
    med2 = result$meds[u]
    sets2 = as.character(result$report$NAME.1[which(result$report$medications ==med2)])
    tanimotos.terms[med1,med2] = length(intersect(sets1,sets2)) / length(union(sets1,sets2))
  }
}

result$tanimotos.terms = tanimotos.terms
save(result,file="DrugGeneLists74Final.Rda")

# Cluster based on gene sets
load("drugNESMatrix.Rda")
disty = dist(df)
hc = hclust(disty)
plot(hc,main="Clustering based on Gene Sets")