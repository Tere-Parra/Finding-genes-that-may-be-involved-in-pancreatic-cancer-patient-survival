# Differential-Expression-Analysis-of-Pancreatic-Cancer-Genes
The Genomic Data Commons (GDC) is a portal that contains a lot of projects related to cancer diseases. The GDC has transcriptomic, proteomic, genomic, and metabolomic data from Cancer patients. Here, I developed a pipeline to analyze a pancreatic cancer project.  [Oficcial GDC page: ]([https://github.com/CdeCMx-org/templates_paginaweb/edit/main/README.md](https://portal.gdc.cancer.gov/])  

The R program has a library named "TCGA_biolinks" with free access to GDC data. 

First, I introduced the project (ductal pancreatic cancer subtype) and downloaded the data to create a counts data frame that we will use for the Differential Analysis Expression. Then, I perform an Enrichment Analysis over the principal DEG genes.

## Methods 

About the project

Here, I use the TCGA-PADD project. This project has 185 cases with 19,556 genes. We have different metadata available. We are going to focus only on dead and live patients. We are conducting a Differential Expression Analysis of pancreatic cancer to find differential dead and live genes.

1. Downloading data
   
TCGA-bio links have their codes and structure. Here, I combined edgeR and limma to perform a DEG analysis.
We downloaded the data with query.raw. 

The function query.raw has different options as data.category, data.type, and workflow.  We choose "Transcriptome Profiling", but this project has "Simple Nucleotide Variation", "mi RNA profile", and more. To consult the available data: ([https://portal.gdc.cancer.gov/projects/TCGA-PAAD]) 

 ``` R
#Project Name
project <- "TCGA-PAAD"

#Category
Category <- "Transcriptome Profiling"

#We choose the STAR counts 
query.raw <- GDCquery(project = project,
                      data.category = category,
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts"
)

# We download all the data 
GDCdownload(query.raw)

# We organize the data and save it in the RAM memory 
SKCM.counts <- GDCprepare(query = query.raw,
                          summarizedExperiment = TRUE)

#Detele it from the RAM 
rm(query.raw)

 ```





REFERENCES

Colaprico, Antonio, et al. “TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data.” Nucleic acids research 44.8 (2015): e71-e71.

Silva, Tiago C., et al. “TCGA Workflow: Analyze cancer genomics and epigenomics data using Bioconductor packages.” F1000Research 5 (2016). (https://f1000research.com/articles/5-1542/v2)

Mounir, Mohamed, et al. “New functionalities in the TCGAbiolinks package for the study and integration of cancer data from GDC and GTEx.” PLoS computational biology 15.3 (2019): e1006701. (https://doi.org/10.1371/journal.pcbi.1006701)
Other useful links


