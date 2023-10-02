# Differential-Expression-Analysis-of-Pancreatic-Cancer-Genes
The Genomic Data Commons (GDC) is a portal that contains a lot of projects related to cancer diseases. The GDC has transcriptomic, proteomic, genomic, and metabolomic data from Cancer patients. Here, I developed a pipeline to analyze a pancreatic cancer project.  [Oficcial GDC page: ]([https://github.com/CdeCMx-org/templates_paginaweb/edit/main/README.md](https://portal.gdc.cancer.gov/)  R has a library name "TCGA_biolinks" with free access to all the data. 

First, I introduced the project (ductal pancreatic cancer subtype) and downloaded the data to create a counts data frame that we will use for the Differential Analysis Expression. Then I perform an Enrichment Analysis over the principal DEG genes.

## Methods 

1. About the project
Here we use the TCGA-PADD project. This project has 185 cases with 19,556 genes. We have different metadata available, but we focus only on dead and live cases. For that reason, we are conducting a Differential Expression Analysis of pancreatic cancer to find differential dead and live genes.

2. Downloading data
TCGA-biolinks have their own codes.

REFERENCES

TCGAbiolinks. (2023). Bioconductor. https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html


