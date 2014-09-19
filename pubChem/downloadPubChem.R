# This script will read in a list of drugs and associated CID, and download target tables
# from pubchem.  We will then read in genes enriched in drugs from our gsea, map those
# genes to proteins, and calculate a validation score.

library(RJSONIO)
library(RCurl)
options(RCurlOptions = list(cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")))

setwd("/home/vanessa/Documents/Dropbox/Code/R/gene2drug/pubChem")

# Load our drug data
load("/home/vanessa/Documents/Dropbox/Code/R/gene2drug/data/DrugGeneLists74Final.Rda")

# Load our CID mappings of drugs to CID identifiers
protdb = read.csv("protDBLookup.csv",sep="\t",head=TRUE)

# Get rid of the NA
protdb = protdb[-which(is.na(protdb$PUBCHEM)),]

# Here we will store structures, properties
property.names = "CID,MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChI,InChIKey,IUPACName,XLogP,ExactMass,MonoisotopicMass,TPSA,Complexity,Charge,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,HeavyAtomCount,IsotopeAtomCount,AtomStereoCount,DefinedAtomStereoCount,UndefinedAtomStereoCount,BondStereoCount,DefinedBondStereoCount,UndefinedBondStereoCount,CovalentUnitCount,Volume3D,XStericQuadrupole3D,YStericQuadrupole3D,ZStericQuadrupole3D,FeatureCount3D,FeatureAcceptorCount3D,FeatureDonorCount3D,FeatureAnionCount3D,FeatureCationCount3D,FeatureRingCount3D,FeatureHydrophobeCount3D,ConformerModelRMSD3D,EffectiveRotorCount3D,ConformerCount3D,Fingerprint2D"
property.names = strsplit(property.names,",")[[1]]
properties = array(dim=c(nrow(protdb),length(property.names)))
rownames(properties) = protdb$MEDICATION
colnames(properties) = property.names

# For each drug, download from pubChem
for (d in 1:nrow(protdb)){
  drug = as.character(protdb$MEDICATION[d])
  cid = as.character(protdb$PUBCHEM[d])

  cat("Processing",as.character(drug),d,"of",nrow(protdb),"\n")
  
  # Get properties
  query = paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",cid,"/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChI,InChIKey,IUPACName,XLogP,ExactMass,MonoisotopicMass,TPSA,Complexity,Charge,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,HeavyAtomCount,IsotopeAtomCount,AtomStereoCount,DefinedAtomStereoCount,UndefinedAtomStereoCount,BondStereoCount,DefinedBondStereoCount,UndefinedBondStereoCount,CovalentUnitCount,Volume3D,XStericQuadrupole3D,YStericQuadrupole3D,ZStericQuadrupole3D,FeatureCount3D,FeatureAcceptorCount3D,FeatureDonorCount3D,FeatureAnionCount3D,FeatureCationCount3D,FeatureRingCount3D,FeatureHydrophobeCount3D,ConformerModelRMSD3D,EffectiveRotorCount3D,ConformerCount3D,Fingerprint2D/JSON",sep="")
  result = getURL(query)
  data <- fromJSON(result)
  # Here are the names
  # names(data$PropertyTable$Properties[[1]])

  # Find fields that are not null
  notnull = names(data$PropertyTable$Properties[[1]])
  for (n in notnull){
    properties[drug,n] = as.character(data$PropertyTable$Properties[[1]][which(names(data$PropertyTable$Properties[[1]])==n)])
  }
    
  # This is structure
  query = paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",cid,"/JSON",sep="")
  result = getURL(query)
  sink(paste("/home/vanessa/Documents/Dropbox/Code/R/gene2drug/data/json/structure/",drug,"_structure.json",sep=""))
  cat(result)
  sink()
  
  # Get picture - not sure how to do this
  #outfile = paste("/home/vanessa/Documents/Dropbox/Code/R/gene2drug/data/img/",med,".png",sep="")
  #query = paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",cid,"/PNG",sep="")
  # download.file(query, outfile)

  # Get associated bioassays
  query = paste("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",cid,"/assaysummary/JSON",sep="")
  result = getURL(query)
  sink(paste("/home/vanessa/Documents/Dropbox/Code/R/gene2drug/data/json/bioassay/",drug,"_bioassay.json",sep=""))
  cat(result)
  sink()
  }  

# Save to file
save(properties,file="../data/drugProperties74.Rda")
