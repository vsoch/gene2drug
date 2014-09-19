# cmapGSEA.R Will read in CMAP instances, normalize, and output table for gsea, intended to be run on the Sherlocj cluster
# by run_cmapGSEA.R.  The input data are compressed .CEL instance files, and a lookup table is used to determine perturbation and
# vehicle (control) files to look for over expression of a particular gene set (with GSEA).

# May 2014 Vanessa Sochat
# Wall Lab
 
library('oligo')
library('affy')
library('R.utils')

args <- commandArgs(TRUE)
name = args[1]
celfiles = args[2]
datadir = args[3]
chip = args[4]

cat(name,":drug\n")
cat(celfiles,":cels\n")
cat(datadir,":output\n")

# Split the cel paths by comma
celfiles = strsplit(celfiles,",")[[1]]
setwd(datadir)

# Unzip cel files
files = c()
for (c in celfiles){  
  if (!file.exists(gsub(".bz2","",c))) {
    bunzip2(c) 
  } 
  file = gsub(".bz2","",c)
  files = c(files,file)    
}

# Summarize and normalize with quantile normalize
affy.data = ReadAffy(filenames = files)
# Summarize and normalize with MAS5
eset.mas5 = rma(affy.data)
exprset.nologs = as.data.frame(exprs(eset.mas5))

# Save the probes to file - will need to look up
# and make chip file for each
probes = rownames(exprset.nologs)
#write.table(probes,file=paste("/scratch/PI/dpwall/DATA/DRUG/CONNECTIVITY_MAP/probes/probes",chip,".dat",sep=""),row.names=FALSE,col.names=FALSE,sep="\n",quote=FALSE)

# Write data to file (do we need a class file?)
#classfile = paste(topdir,"/cls/",gsub(".gct",".cls",o),sep="")

# Add probe names, description(with na)
DESCRIPTION = rep('na',dim(exprset.nologs)[1])
NAME = rownames(exprset.nologs)
filtered = cbind(NAME,DESCRIPTION,exprset.nologs)

# Now print new file
outfile = paste(datadir,"/",name,"_",chip,".gct",sep="")
dimprint = paste(dim(filtered)[1],dim(filtered)[2]-2,sep="\t")
fileConn = file(outfile)
writeLines(c("#1.2",dimprint),sep="\n",fileConn)
close(fileConn)  
write.table(filtered,file=outfile,append=TRUE,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")   
