

#Paquetes a instalar

install.packages(c("scales", "pheatmap", "factoextra", "BiocManager", "EDASeq", "tidyverse"))
BiocManager::install(c ("NOISeq", "ComplexHeatmap", "TCGAbiolinks", "limma"))
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "DOSE", "enrichplot", "edgeR", "SummarizedExperiment"))
BiocManager::install("TCGAbiolinks")


#librerias 
library(TCGAbiolinks)
library(SummarizedExperiment)
library(pheatmap)
library(limma)
library(scales)
library(BiocManager)
library(SummarizedExperiment)
library(NOISeq)
library(ComplexHeatmap)
library(EDASeq)
library(factoextra)
library(tidyverse)
library(edgeR)



########################################################################
############################  STEP 1 #########################
######################  Descarga de datos del TCGA ###########
#################################################################

# 1. Descarga de datos del GDC para el proyecto llamado TCGA-PAAD

proyecto1 <- "TCGA-PAAD"
categorias_deseadas <- "Transcriptome Profiling"
query.raw <- GDCquery(project = proyecto1,
                      data.category = categorias_deseadas,
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts"
)

#2- Descarga los datos
GDCdownload(query.raw)

# 3- Organiza los datos y los guarda en una variable en la memoria RAM de tu ordenador
SKCM.counts <- GDCprepare(query = query.raw,
                          summarizedExperiment = TRUE)

#Eliminar del RAM el query 
rm(query.raw)

# 4. Generar la Matriz de expresión (conteos en TPM)
data2<-assay(SKCM.counts)

#guardar tabla de conteos 
write.table(data2, "conteos_TCGA.csv")

save(data2, SKCM.counts, file="conteos.RData")

###########################################################################
######################### STEP 2 #########################

load("conteos.RData")

######################## Pre-procesamiento ################

# Pre-procesamiento
# 1- Función TCGAanalyze_Preprocessing
dataPrep<-TCGAanalyze_Preprocessing(object=SKCM.counts,
                                    cor.cut = 0.6
                                    )

#No detecto outliers !! 

# 2- Función TCGAanalyze_Normalization
dataNorm<-TCGAanalyze_Normalization(tabDF=data2,
                                    geneInfo = TCGAbiolinks::geneInfoHT,
                                    method="gcContent")


#Revisamos que todo vaya bien con la normalizacion y metodo

boxplot(dataPrep, outline = FALSE)
boxplot(dataNorm, outline = FALSE)

# 3- Función TCGAanalyze_Filtering (se filtran genes)
dataFilt<-TCGAanalyze_Filtering(tabDF=dataNorm,
                                method="quantile", 
                                qnt.cut = 0.25)

# 4- Método de normalización -- >  TMM
dataTMMnorm<- tmm(dataFilt)


write.table(dataTMMnorm, "Tabla_Normalizada.csv")

save(dataPrep, dataNorm, dataFilt, dataTMMnorm, file="Preprocesamiento.RData")

###############################################################
####################### STEP 3  ##################
############### graficas ############################

load("Preprocesamiento.RData")

# Boxplots

boxplot(data2[,1:50], outline=FALSE, main="Antes de la normalización", xaxt="n")
boxplot(dataTMMnorm[,1:50], outline=FALSE, main="Después de la normalización", xaxt="n")

# Tomar los 1500 genes expresados con mas variabilidad.
varianza<-apply(dataTMMnorm, 1, var)
varianza<-sort(varianza, decreasing=TRUE) #decreasing a TRUE para coger los 1500 genes de más varianza
milquinientosgenes<-varianza[1:1500]
genes<-names(milquinientosgenes)
milquinientosgenesdata<-dataTMMnorm[genes,]

# Aplicar PCA aplicado sobre los genes
library(factoextra)
pca <- prcomp(milquinientosgenesdata[1:50,1:20])
fviz_eig(pca)

#Distribucion de los genes
library(scales)
fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

 

############################## STEP 4 ##########################
################## Analisis de Expresion diferencial ###############


#Información de muestras (metadatos disponibles del proyecto)
sample.info<-colData(SKCM.counts)
print(sample.info)
head(sample.info)

#En este caso queremos evaluar genes DEG entre pacientes vivos y muertos
# Separar los datos en pacientes vivos y muertos.
TMdata<-dataTMMnorm[,which(sample.info@listData[["paper_Follow up vital status"]]=="Dead")]

PSTdata<-dataTMMnorm[,which(sample.info@listData[["paper_Follow up vital status"]]=="Alive")]

save(TMdata, PSTdata, file="Conditions.RData")


###############################################################
################################################################
############# Expresión Diferencial ####################

load("Conditions.RData")

#Aplicar las condiciones
dataDEGs <- TCGAanalyze_DEA(mat1 = TMdata,
                            mat2 = PSTdata,
                            Cond1type = "Dead",
                            Cond2type = "Alive",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT") 


# 79 samples in Cond1type Dead
# 71 samples in Cond2type Alive
# 45,362 features as miRNA or genes

summary(dataDEGs)

#Filtramos solo aquellos genes codificadores
DEG <- filter(dataDEGs, gene_type=="protein_coding")
#861 genes codificantes 

#guardamos las tablas
write.table(dataDEGs, "Genes_Diferenciados_TCGA.csv",
            row.names=T)
write.table(DEG, "Genes_Dif_ProteinCoding.csv",
            row.names=T)

save(DEG, dataDEGs, file="DEG.RData")

summary(DEG)

genesexpresadosdif <- as.character(rownames(DEG))
genesexpresadosdif[1:10]


######################   Volcano plot and heatmap ##############

### A Basic volcano plot
par(mfrow=c(1,1))

# volcano plot
with(DEG, plot(logFC, -log10(PValue), pch=20, main="Volcano plot", xlim=c(-7,7)))

#Seleccionamos estos valores, pero es importante adecuar estos valores de acuerdo a tus datos y realizar antes un analisis de los valores
# azul si padj<0.1, rojo si log2FC>1 y padj<0.1)
with(subset(DEG, PValue<0.05 ), points(logFC, -log10(PValue), pch=20, col="blue"))
with(subset(DEG, PValue<0.05 & abs(logFC)>2), 
     points(logFC, -log10(PValue), pch=20, col="red"))
  

#heatmap
pheatmap(DEG, main="Heatmap", color = heat.colors, cluster_rows = T,
         show_rownames=F, border_color=NA, scale="row",
         fontsize_row = 8, fontsize_col = 12, angle_col = "45")

                                                

####################################################33
      ### Enriquecimiento funcional ###

#1. Seleccionar solo los genes que están involucrados en codificacion proteica

DegGenes <- dataDEGs[dataDEGs$gene_type == "protein_coding", ]

#2. Esto dió como resultado 861 genes en total expresados diferencialmente. 
#Cargamos las librerias para el Enriquecimiento funcional

BiocManager::install(c("org.Hs.eg.db", "DOSE", "GSEABase","ggridges"))
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("clusterProfiler")
BiocManager::install("HDO.db")

library(clusterProfiler)
library(org.Hs.eg.db)
library (DOSE)
library(GSEABase)
library(pathview)
library(enrichplot)
library(ggridges)

# Seleccionamos el valor log2 fold change 
original_gene_list <- DegGenes$logFC

# Ahora la lista de genes 
names(original_gene_list) <- DegGenes$gene_name

# omitimos cualquier valor NA  
gene_list<-na.omit(original_gene_list)

# Ordenamos la lista de manera decresiente (requisito del paquete clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
print(gene_list)

#3. Realizamos el GSEA con el organimo "homo sapiens"

keytypes(org.Hs.eg.db)

# Parametros para el GSEA
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", #Queremos CC, FM, MM (leer sobre GO TERMS, en caso de no saber)
             keyType = "SYMBOL", #ID tipo de entrada de genes
             nPerm = 10000, #Numero de permutaciones
             minGSSize = 3, #Numero minimo de genes
             maxGSSize = 800, #Maximo numero de genes
             pvalueCutoff = 0.05, #Corte 
             verbose = TRUE, #Mostrar lo que esta pasando
             OrgDb = "org.Hs.eg.db", #Organimo
             pAdjustMethod = "none") #Metodo de ajuste

#283 terminos GO como resultado

#guardamos el enriquecimiento 
enrichment <- as.data.frame(gse)
write.table(enrichment, "tabla_enriquecimiento.csv")


#4. Realizamos graficas para ver los resultados
head(gse)

#Matriz de puntos 
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

#Mapa de enriquecimiento
emapplot(gse, showCategory = 10)

#Diagrama de cresta
ridgeplot(gse) + labs(x = "enrichment distribution")

#Diagrama de red
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 10)

###############################################
#Enriqueciminento funcional KEGG ####

#Convertimos los ID's a Entrez ID's

ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", 
          OrgDb="org.Hs.eg.db")
head(gene_list)

#eliminar duplicados 
dedup_ids = ids[!duplicated(ids[c("ENTREZID")]),]

#Creamos un data frame solo con los genes mapeados
df2 = DegGenes[DegGenes$gene_name %in% dedup_ids$SYMBOL,]

#cremos una columnas con los identificadores
colnames(df2)[colnames(df2) == "gene_name"] <- "SYMBOL"

df2 <- merge(dedup_ids, df2, by="SYMBOL")

#Con el fold change
kegg_gene_list <- df2$logFC

#Los nomnbramos con los IDs
names(kegg_gene_list) <- df2$ENTREZID

#Omite el valor NA
kegg_gene_list<-na.omit(kegg_gene_list)

#ordenalo en orden decreciente 
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

print(kegg_gene_list)

### Analisis KEEG
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = "hsa",
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

#Guardamos la tabla KEGG
KEGG <- as.data.frame(kk2)
write.table(KEGG, "KEGG_Table.csv")

#Grafica de puntos
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + 
  facet_grid(.~.sign)

# native KEGG plot (PNG)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04613", 
                species = "hsa")

knitr::include_graphics("hsa04613.pathview.png")



