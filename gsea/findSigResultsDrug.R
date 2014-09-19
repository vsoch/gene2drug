# Get list of files in the output directory
outdir = "/scratch/PI/dpwall/DATA/DRUG/CONNECTIVITY_MAP/gseaRMA"
folders = list.files(outdir)
#folders = folders[-which(folders == "ERROR")]
#folders = folders[-which(folders=="nocodazole_HGU133A.Gsea.1402258840744")]

# Threshold results at FDR q value:
threshold = .05
# This is the report directory to create and write results to
writedir = "/scratch/PI/dpwall/DATA/DRUG/CONNECTIVITY_MAP/report/" 

reportFinal = c()
for (f in folders){
  reportdir = paste(outdir,"/",f,"/",sep="")
  report = list.files(reportdir,pattern="gsea_report_for_DRUG")
  idx = grep("*.xls",report)
  report = report[idx]
  report = read.csv(paste(reportdir,"/",report,sep=""),sep="\t",head=TRUE)
  report = report[report$FDR.q.val < threshold,]
  tmp = cbind(rep(f,dim(report)[1]),report)
  reportFinal = rbind(reportFinal,tmp)
}
colnames(reportFinal)[1] = "NAME"

# Save to report folder
dir.create(writedir, showWarnings = FALSE)
save(reportFinal,file=paste(writedir,"DRUGS_GSEA_FDRpt05_RMA_Filter.Rda",sep=""))
write.table(reportFinal,file=paste(writedir,"DRUGS_GSEA_FDRpt05_RMA_Filter.txt",sep=""),sep="\t",row.names=FALSE)

# Extract the unique medications
result = reportFinal

# Threshold results at FDR q value:
threshold = .05
result = result[which(result$FDR.q.val <= .05),]
tmp = as.character(result$NAME)
tmp= sapply(tmp,strsplit,"_")

# Now find medications
medications = c()
for (t in tmp){
  medications = c(medications,t[[1]][1])
}

result = cbind(medications,result)
colnames(result)[2] = "FOLDER"
colnames(result)[3] = "GENESET"
meds = sort(unique(medications))
# Next, let's download the genes for each of our core results
# This will be the same order as the medications "meds"
geneLists = list()

# This is the top directory for the results
topdir = "/scratch/PI/dpwall/DATA/DRUG/CONNECTIVITY_MAP/gseaRMA/"
outfile = paste(writedir,"drugGeneLists.txt",sep="")

for (mm in 1:length(meds)) {
   cat("Processing",mm,"of",length(meds),"\n")
   m = meds[mm]
   tmplist = c()
   # Filter data to that medication
   subset = result[which(result$medications==as.character(m)),]
   # Now for each significant term, get the genes!
   for (ss in 1:length(subset$FOLDER)){
     cat("  folder",ss,"of",length(subset$FOLDER),"\n")
     s = subset$FOLDER[ss]
     folder = as.character(s)
     tmp = subset[which(subset$FOLDER %in% folder),]
     geneset = as.character(tmp$GENESET)
     for (g in geneset){
       filepath = paste(topdir,folder,sep="")
       filepath = paste(filepath,"/",list.files(filepath,pattern=paste("^",g,"*.xls",sep="")),sep="")
       if (file.exists(filepath)){
         filey = read.csv(filepath,sep="\t",head=TRUE)
         tmplist = unique(c(tmplist,as.character(filey$PROBE[which(filey$CORE.ENRICHMENT == "Yes")]))) 
       } else {
        cat(filepath,"does not exist\n")
       }
     }
   }
   tmplist = unique(tmplist)
   # Write to output file
   tmplist = paste(sort(tmplist),collapse=" ")
   cat(m,"\t",tmplist,"\n",file=outfile,append=TRUE)
   geneLists = c(geneLists,list(tmplist))
}

result = list(geneLists=geneLists,meds=meds,report=reportFinal)
save(result,file="DrugGeneLists74.Rda")
# Now let's explore these results!!
#result = read.table(file=paste(writedir,"DRUGS_GSEA_FDRpt05_RMA_Filter.txt",sep=""),head=TRUE)
#subset = reportFinal[which(reportFinal$NAME == unique(reportFinal$NAME)[18]),]

load('DRUGS_GSEA_FDRpt05_RMA_Filter.Rda')
result = reportFinal

# First let's create a matrix of medications (rows), and gene sets (columns), and
# we want to see if we can cluster.
geneset = unique(sort(result$NAME.1))

result = result[,-c(4,5,14)]
result = cbind(medications,result)
colnames(result)[2] = "FOLDER"
colnames(result)[3] = "GENESET"
save(result,file=paste("DRUGS_GSEA_FDRpt05_RMA_Filter_Final.Rda",sep=""))

load(paste(writedir,"DRUGS_GSEA_FDRpt05.Rda",sep=""))
# Matrix - drugs in rows, gene sets in columns

geneset = as.character(unique(result$GENESET))

df = array(data=0,dim=c(length(unique(result$medications)),length(geneset)))
rownames(df) = as.character(sort(unique(result$medications)))
colnames(df) = as.character(geneset)

for (g in geneset){
  gene = as.character(g)
  med = as.character(result$medications[which(result$GENESET %in% gene)])
  df[which(rownames(df) %in% as.character(med)),g] = 1
}

disty = dist(df)
hc = hclust(disty)
plot(hc,main="Clustering drugs by BrainTerms Gene Sets")

# RXNORM WORK IS BELOW - WE DON'T ACTUALLY NEED IT
library(RCurl)
library(XML)

# We will save a final list of medications
actions = array(dim=74)   # list of actions
foundmeds = array(dim=74) # list of meds we found
missmeds = c()  # missing medications

medications = unique(medications)

# Now we need to look up medication names in RxNorm:
idx=1
for (i in 1:length(medications)){
  m = medications[i]
  med = tolower(as.character(m))
  xml = getURL(paste('http://rxnav.nlm.nih.gov/REST/approximateTerm?term=',med,"&maxEntries=20&option=1",sep=""))
  xml = xmlParse(xml)
  xml = xmlToList(xml)
  if (length(xml$approximateGroup$candidate$rxcui) == 0){
    missmeds = c(missmeds,med)
  } else { # Look up actions with cuid
    cuid = xml$approximateGroup$candidate$rxcui
    xml = getURL(paste('http://rxnav.nlm.nih.gov/REST/rxcui/',cuid,"/hierarchy?src=MESH&oneLevel=1",sep=""))
    xml = xmlParse(xml)
    xml = xmlToList(xml)
    if (xml[[1]]$title != "***  No MeSH identifier found ***"){
      if (!is.null(xml[[1]][3]$node$nodeName)) {
        actions[idx] = xml[[1]][3]$node$nodeName
      } else {
        actions[idx] = "NA" 
      }
    }
    foundmeds[idx] = med
    idx = idx+1
  }
}

medlookup = c(actions,others)
allmeds= c(foundmeds,missmeds)

medtable = data.frame(actions=medlookup,medication=allmeds)
save(medtable,file=paste(writedir,"DRUGS_GSEA_FDRpt05_mediction.Rda",sep=""))
load(paste(writedir,"DRUGS_GSEA_FDRpt05_medaction.Rda",sep=""))
# we find 72 of the set! Look up others manually

# IDEA - find drugs with common gene sets - calculate correlation between
# enrichment scores, GEPHI to visualize similarity

# Look at larger sets
subset = reportFinal[which(reportFinal$SIZE>50),]

# Look up drug names
labels = c()
for (i in 1:length(result$medications)){
  idx = which(as.character(medtable$medication) %in% as.character(result$medications[i]))
  idx = idx[1]
  labels = c(labels,as.character(medtable$actions[idx]))
}
result$actions = labels


# Now redo  clustering with drug indications
load(paste(writedir,"DRUGS_GSEA_FDRpt05.Rda",sep=""))
geneset = unique(sort(result$GENESET))

# Matrix - drugs in rows, gene sets in columns
df = array(data=0,dim=c(length(unique(result$medications)),length(geneset)))
rownames(df) = as.character(sort(unique(result$medications)))
colnames(df) = as.character(geneset)

# Make new row labels for unique drugs
newrows = c()
for (m in rownames(df)){
  newrows = c(newrows,result$actions[which(result$medications==m)[1]])  
}

for (g in geneset){
  gene = as.character(g)
  med = as.character(result$medications[which(result$GENESET %in% gene)])
  df[which(rownames(df) %in% as.character(med)),g] = result$NES[which(result$GENESET %in% gene)]
}
colors = as.numeric(as.factor(rownames(df)))
rownames(df) = newrows
disty = dist(df)
hc = hclust(disty)
plot(hc,main="Clustering drugs by BrainTerms Gene Sets")

# Output file for gephi
tmp =as.matrix(dist(df))
tmp2 =as.matrix(dist(df2))

# Make drug by term matrix with NES scores
seuil<-0.4
index<-which(abs(test$co[,1])>seuil)
selection<-test$co[index,]
s.arrow(selection)

write.table(tmp2,file="/scratch/PI/dpwall/SCRIPT/R/gsea/gseaR/gephiDrug.csv",sep=";")
)
