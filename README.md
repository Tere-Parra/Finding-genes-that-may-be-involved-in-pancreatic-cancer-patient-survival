# Differential-Expression-Analysis-of-Pancreatic-Cancer-Genes
The Genomic Data Commons (GDC) is a portal that contains a lot of projects related to cancer diseases. The GDC has transcriptomic, proteomic, genomic, and metabolomic data from Cancer patients. Here, I developed a pipeline to analyze a pancreatic cancer project.  [Oficcial GDC page: ]([https://github.com/CdeCMx-org/templates_paginaweb/edit/main/README.md](https://portal.gdc.cancer.gov/)  

The R program has a library named "TCGA_biolinks" with free access to GDC data. 

First, I introduced the project (ductal pancreatic cancer subtype) and downloaded the data to create a counts data frame that we will use for the Differential Analysis Expression. Then, I perform an Enrichment Analysis over the principal DEG genes.

## Methods 

About the project
Here, I use the TCGA-PADD project. This project has 185 cases with 19,556 genes. We have different metadata available. We are going to focus only on dead and live patients. We are conducting a Differential Expression Analysis of pancreatic cancer to find differential dead and live genes.

1. Downloading data
TCGA-bio links have their codes and structure. Here, I combined edgeR and limma to perform a DEG analysis.
We downloaded the data with query.raw  

 ``` R

proyecto1 <- "TCGA-PAAD"
categorias_deseadas <- "Transcriptome Profiling"
query.raw <- GDCquery(project = proyecto1,
                      data.category = categorias_deseadas,
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts"
)

 ```

REFERENCES

TCGAbiolinks. (2023). Bioconductor. https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html


