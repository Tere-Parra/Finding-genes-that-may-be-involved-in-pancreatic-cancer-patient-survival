# Differential-Expression-Analysis-of-Pancreatic-Cancer-Genes
The Genomic Data Commons (GDC) is a portal that contains a lot of projects related to cancer diseases. The GDC has transcriptomic, proteomic, genomic, and metabolomic data from Cancer patients. Here, I developed a pipeline to analyze a pancreatic cancer project.  [Oficcial GDC page: ]([https://github.com/CdeCMx-org/templates_paginaweb/edit/main/README.md](https://portal.gdc.cancer.gov/])  

The R program has a library named "TCGA_biolinks" with free access to GDC data. 

First, I introduced the project (ductal pancreatic cancer subtype) and downloaded the data to create a counts data frame that we will use for the Differential Analysis Expression. Then, I perform an Enrichment Analysis over the principal DEG genes.

## Methods 

About the project

Here, I use the TCGA-PADD project. This project has 185 cases with 19,556 genes. We have different metadata available. We are going to focus only on dead and live patients. We are conducting a Differential Expression Analysis of pancreatic cancer to find differential dead and live genes.

**1. Downloading data**
   
TCGA-bio links have their codes and structure. Here, I combined edgeR and limma to perform a DEG analysis.
We downloaded the data with query.raw. 

The function query.raw has different options as data.category, data.type, and worflow.type (counts data).  We choose "Transcriptome Profiling", but this project has different data.category such as: "Simple Nucleotide Variation", "mi RNA profile", and more. To consult the available data: ([https://portal.gdc.cancer.gov/projects/TCGA-PAAD]) . Data.Type is the quantification form. workflow.type has the available software counts data. 


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

Next, we create the count's data frame with the assay() function. 

 ``` R
# We generate the Expression Matrix (TPM counts)
data2<-assay(SKCM.counts)

#save the data frame
write.table(data2, "conteos_TCGA.csv")

#Save the R data session
save(data2, SKCM.counts, file="conteos.RData")

 ```

**2. Pre-process**

In these steps, we eliminate outliers and perform the data normalization.  

First, we use the TCGAanalyze_Preprocessing() function to find outliers. Then, with the TCGAanalyze_Normalization() normalized the data (samples) and observed the different distribution with boxplots. Finally, TCGAanalyze_Filtering() filter genes, and with tmm() we applied the TMM normalization.

We use TMM normalization. 

``` R

# 
# 1- Function TCGAanalyze_Preprocessing
dataPrep<-TCGAanalyze_Preprocessing(object=SKCM.counts,
                                    cor.cut = 0.6
                                    )

#No outliers !! 

# 2- Function TCGAanalyze_Normalization
dataNorm<-TCGAanalyze_Normalization(tabDF=data2,
                                    geneInfo = TCGAbiolinks::geneInfoHT,
                                    method="gcContent")


#Revisamos que todo vaya bien con la normalizacion y metodo

boxplot(dataPrep, outline = FALSE)
boxplot(dataNorm, outline = FALSE)

# 3- Function TCGAanalyze_Filtering 

dataFilt<-TCGAanalyze_Filtering(tabDF=dataNorm,
                                method="quantile", 
                                qnt.cut = 0.25)

# 4- Normalized Method -- >  TMM
dataTMMnorm<- tmm(dataFilt)

#We save the normalized data frame 
write.table(dataTMMnorm, "Tabla_Normalizada.csv")

#Save the project
save(dataPrep, dataNorm, dataFilt, dataTMMnorm, file="Preprocesamiento.RData")
``` 

REFERENCES

Colaprico, Antonio, et al. “TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data.” Nucleic acids research 44.8 (2015): e71-e71.

Silva, Tiago C., et al. “TCGA Workflow: Analyze cancer genomics and epigenomics data using Bioconductor packages.” F1000Research 5 (2016). (https://f1000research.com/articles/5-1542/v2)

Mounir, Mohamed, et al. “New functionalities in the TCGAbiolinks package for the study and integration of cancer data from GDC and GTEx.” PLoS computational biology 15.3 (2019): e1006701. (https://doi.org/10.1371/journal.pcbi.1006701)
Other useful links


