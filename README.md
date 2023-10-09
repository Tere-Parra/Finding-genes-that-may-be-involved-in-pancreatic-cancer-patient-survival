# Finding differential genes that may be involved in the overall survival of pancreatic cancer patients

The Genomic Data Commons (GDC) is a portal that contains a lot of projects related to cancer diseases. The GDC has transcriptomic, proteomic, genomic, and metabolomic data from Cancer patients. Here, I developed a pipeline to analyze a pancreatic cancer project.  [Oficcial GDC page: ]([https://github.com/CdeCMx-org/templates_paginaweb/edit/main/README.md](https://portal.gdc.cancer.gov/])  

The R program has a library named "TCGA_biolinks" with free access to GDC data. 

### About the project

Here, I use the TCGA-PADD project. This project has 185 cases with 19,556 genes. We have different metadata available. We are going to focus only on dead and live patients. We are conducting a Differential Expression Analysis of pancreatic cancer to find differential dead and live genes.

First, I introduced the project (ductal pancreatic cancer subtype) and downloaded the data to create a counts data frame that we will use for the Differential Analysis Expression. Then, I perform an Enrichment Analysis over the principal DEG genes.

## Methods 


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

![Before Normalization](Antes_Norm.png)

![After Normalization](After_Norm.png)


**3. Expression analysis**

We are going to see the effects of the normalization with the PCA graph. 

``` R

# We take the first 1500 genes.
varianza <-apply(dataTMMnorm, 1, var)
varianza <-sort(varianza, decreasing=TRUE) #decreasing a TRUE to pick up the more variant genes
milquinientosgenes <-varianza[1:1500]
genes <-names(milquinientosgenes)
milquinientosgenesdata <-dataTMMnorm[genes,]

# PCA
library(factoextra)
pca <- prcomp(milquinientosgenesdata[1:50,1:20])
fviz_eig(pca)

#Genes distribution
library(scales)
fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
 
```
![](PCA_After_Normalization_genes.png)


**4.DEG**

First, we need to verify the available metadata. then, we choose the Dead and Alive projects. 

``` R

#1 We see the metadata 
sample.info<-colData(SKCM.counts)
print(sample.info)
head(sample.info)

#We need the Dead and Alive patients. 
TMdata<-dataTMMnorm[,which(sample.info@listData[["paper_Follow up vital status"]]=="Dead")]

PSTdata<-dataTMMnorm[,which(sample.info@listData[["paper_Follow up vital status"]]=="Alive")]

#Save the R session
save(TMdata, PSTdata, file="Conditions.RData")

```
Now we have the samples and we can perform the DEG with the condition (DEAD over ALIVE patients). The function is TCGAanalyze_DEA()

``` R
load("Conditions.RData")

#Perform DEG
dataDEGs <- TCGAanalyze_DEA(mat1 = TMdata,
                            mat2 = PSTdata,
                            Cond1type = "Dead",
                            Cond2type = "Alive",
                            fdr.cut = 0.01 ,  #False Discovery Rate
                            logFC.cut = 1,   #Fold Change
                            method = "glmLRT") 


# 79 samples in Cond1type Dead
# 71 samples in Cond2type Alive
# 45,362 features as miRNA or genes

summary(dataDEGs)

``` 

By defect, we have miRNA, mRNA, and other transcriptomic data. We are going to separate the protein-coding genes because we want to perform an enrichment analysis.

``` R
#Coding genes
DEG <- filter(dataDEGs, gene_type=="protein_coding")
#861 coding genes 

#save the data frames
write.table(dataDEGs, "Genes_Diferenciados_TCGA.csv",
            row.names=T)
write.table(DEG, "Genes_Dif_ProteinCoding.csv",
            row.names=T)

#save the session
save(DEG, dataDEGs, file="DEG.RData")

summary(DEG)

``` 
Let's see the DEG results with the volcano plot and the heatmap 

``` R
######################   Volcano plot and heatmap ##############

### A Basic volcano plot
par(mfrow=c(1,1))

# volcano plot
with(DEG, plot(logFC, -log10(PValue), pch=20, main="Volcano plot", xlim=c(-7,7)))

#We select these values, but it is important to adjust these values according to your data and perform an analysis of the values #beforehand.

# blue padj<0.1, red si log2FC>1 y padj<0.1)
with(subset(DEG, PValue<0.05 ), points(logFC, -log10(PValue), pch=20, col="blue"))
with(subset(DEG, PValue<0.05 & abs(logFC)>2), 
     points(logFC, -log10(PValue), pch=20, col="red"))
  
#heatmap
pheatmap(DEG, main="Heatmap", color = heat.colors, cluster_rows = T,
         show_rownames=F, border_color=NA, scale="row",
         fontsize_row = 8, fontsize_col = 12, angle_col = "45")

```
![](Volcano_Plot2.png)



## Enrichment Analysis 

Now, we will perform an Enrichment Analysis. First and foremost, we will make the GO function and graphs 

``` R

#Libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library (DOSE)
library(GSEABase)
library(pathview)
library(enrichplot)
library(ggridges)

# We need to select the LogFC genes 
original_gene_list <- DegGenes$logFC

# Now, the genes list
names(original_gene_list) <- DegGenes$gene_name

# We omit NA values 
gene_list<-na.omit(original_gene_list)

# We sort the list in descending order (requirement of the clusterProfiler package)
gene_list = sort(gene_list, decreasing = TRUE)
print(gene_list)

#3. Let's make the GSEA 

#Humans package
keytypes(org.Hs.eg.db)

# GSEA Parameters
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", #We want CC, FM, MM 
             keyType = "SYMBOL", #ID type of the Genes
             nPerm = 10000, #Permutation number
             minGSSize = 3, #Minimium genes 
             maxGSSize = 200, #Max number of genes
             pvalueCutoff = 0.05, #Cut 
             verbose = TRUE, #Print the process
             OrgDb = "org.Hs.eg.db", #Organism
             pAdjustMethod = "none") #Ajust method

#283 GO Terms

#Save the results
enrichment <- as.data.frame(gse)
write.table(enrichment, "tabla_enriquecimiento.csv")


#4. We built graphics

head(gse)

#Enrich Map
emapplot(gse, showCategory = 10)

#(see the code attached for more)
```
![](GO_Terms.png)


## Acknowledgment

Tere Parra (2022).  GitHub. [https://github.com/Tere-Parra/Finding-genes-that-may-be-involved-in-pancreatic-cancer-patient-survival.git]

## REFERENCES

Colaprico, Antonio, et al. “TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data.” Nucleic acids research 44.8 (2015): e71-e71.

Silva, Tiago C., et al. “TCGA Workflow: Analyze cancer genomics and epigenomics data using Bioconductor packages.” F1000Research 5 (2016). (https://f1000research.com/articles/5-1542/v2)

Mounir, Mohamed, et al. “New functionalities in the TCGAbiolinks package for the study and integration of cancer data from GDC and GTEx.” PLoS computational biology 15.3 (2019): e1006701. (https://doi.org/10.1371/journal.pcbi.1006701)
Other useful links

Carmona Pedo & Cano Carlos, Analisis bioinformatico para un problems de ómicas. Google Colaboratory. (2019). Google.com. https://colab.research.google.com/drive/1sgkBRHqzUOxUHh4j_li3daHWqyNjEDOk?usp=sharing#scrollTo=Oe0PL76cRV23

‌


