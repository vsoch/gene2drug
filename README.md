# gene2drug

UNDER DEVELOPMENT

This toolbox will take a set of genes, and allow a researcher to determine if any of the genes are meaningfully related to drugs. The algorithm works as follows:

## Gene Set Enrichment Analysis
1) Perform GSEA (gene set enrichment analysis) to test gene sets against a database of thousands of drug compounds
  
  GSEA/
  cmapGSEA.R: convert expression objects (.bz2 files) to expression matrices
  cmapGSEARMA.R: runs gene set enrichment analysis with quantile normalization of above
  findSigResultsDrug.R: parses GSEA results and returns significant genes, drugs, etc. for a specified FDR threshold
  parseDrugResult.R: format results into object that includes:
     geneLists: a list of core enriched genes for each drug
     meds: a list of drugs with significantly enriched genes
     report: The full report table from GSEA
     genematrix: a binary matrix of drugs by genes, with 1 indicated significant enrichment
     tanimotos: tanimoto scores to assess similarity of drugs based on enriched genes
     terms: gene set names that are significantly enriched for each drug

## PubChem
2) Extract structure and properties of drugs from Pubchem

  pubChem/
  downloadPubChem.R: create table of drug properties, and download structure and bioassay json files

## gene2drug Enrichment Scoring
3) 

TODO:
  - map genes to proteins
  - calculate enrichment score of proteins in my drugs
  - generate final report

## Visualize & Summarize Results
4) Visualize

  gephi/
  export data for gephi (not useful!)
  exportd3: todo: will export drug to visualize in d3
