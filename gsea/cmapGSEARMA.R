# This script will read in un-normalized cmap data, and normalize with rma (quantile),
# and run gsea

# June 2014 Vanessa Sochat
# Wall Lab

library('oligo')
library('preprocessCore')
library('R.utils')

datadir = "/scratch/PI/dpwall/DATA/DRUG/CONNECTIVITY_MAP/inputRMA"
indir = "/scratch/PI/dpwall/DATA/DRUG/CONNECTIVITY_MAP/input"
files = list.files(indir)
setwd('/scratch/PI/dpwall/SCRIPT/R/gsea')

# First prepare data
for (i in 1: length(files)){
  f=files[i]
  # read in file
  data = read.csv(paste(indir,"/",f,sep=""),skip=2,sep="\t")
  meta = data[,c(1,2)]
  data = data[,-c(1,2)]
  norm = normalize.quantiles(as.matrix(data))
  colnames(norm) = colnames(data)
  # Add probe names, description(with na)
  DESCRIPTION = rep('na',dim(norm)[1])
  NAME = meta$NAME
  filtered = cbind(as.character(NAME),DESCRIPTION,norm)
  colnames(filtered)[1] = "NAME"
  
  # Now print new file
  outfile = paste(datadir,"/",f,sep="")
  dimprint = paste(dim(filtered)[1],dim(filtered)[2]-2,sep="\t")
  fileConn = file(outfile)
  writeLines(c("#1.2",dimprint),sep="\n",fileConn)
  close(fileConn)  
  write.table(filtered,file=outfile,append=TRUE,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")   
}


files = list.files(datadir)
# Now run GSEA when the top above is done
topdir = "/scratch/PI/dpwall/DATA/DRUG/CONNECTIVITY_MAP/"
for (f in files){
  name = gsub(".gct","",f)
  chip = strsplit(name,"_")[[1]]
  chip = paste(chip[2:length(chip)],collapse="")
  chip = paste(topdir,"chip/",chip,".chip",sep="")
  input = paste(topdir,"inputRMA/",name,".gct",sep="")
  cls = paste(topdir,"cls/",name,".cls",sep="")
  
  jobby = paste(f,"_gseaN.job",sep="")
  sink(file=paste(".job/",jobby,sep=""))
  cat("#!/bin/bash\n")
  cat("#SBATCH --job-name=",jobby,"\n",sep="")  
  cat("#SBATCH --output=.out/",jobby,".out\n",sep="")  
  cat("#SBATCH --error=.out/",jobby,".err\n",sep="")  
  cat("#SBATCH --time=2-00:00\n",sep="")
  cat("#SBATCH --mem=8000\n",sep="")
  cat("java -cp",gseadir,"xtools.gsea.Gsea -res",input,"-cls",as.character(cls),"-gmx",inputdb,"-chip",as.character(chip),"-collapse true -mode Max_probe -norm None -nperm 1000 -permute genes -rnd_type no_balance -scoring_scheme weighted -rpt_label",f,"-metric Signal2Noise -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out",outdir,"-gui false\n")
  sink()
  
  system(paste("sbatch",paste(".job/",jobby,sep="")))
}

