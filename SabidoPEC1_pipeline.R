#Definimos directorios
workingDir <-getwd()
dataDir <-file.path(workingDir, "data")
resultsDir <- file.path(workingDir, "results")

#Cargamos librerias
require(Biobase)
require(affy)
require(GEOquery)
library(arrayQualityMetrics)
library(oligo)
library(knitr)
library(ggplot2)
library(ggrepel)
library(genefilter)
library(limma)
library(hgu133plus2.db)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(magrittr)


#Identificación de grupos
sampleInfo <- read.AnnotatedDataFrame(file.path(dataDir,"targets.txt"),header=TRUE,
                                      row.names=1,sep="\t")
targets <- read.table(file.path(dataDir,"targets.txt"), header = TRUE, sep = "\t")
fileNames <- pData(sampleInfo)$FileName
rawData <- read.affybatch(filenames=file.path(dataDir,fileNames), phenoData=sampleInfo)

targets_table <- read.table("data/targets.txt", header = TRUE)
sampleNames <- targets_table$SampleName
Group <- targets_table$Group
kable(targets_table, caption = "Resumen de los individuos y sus condiciones experimentales")


#Control de calidad

plotPCA3 <- function (datos, labels, factor, title, scale,colores, size = 1.5, glineas = 0.25) {
  data <- prcomp(t(datos),scale=scale)
  # plot adjustments
  dataDf <- data.frame(data$x)
  Group <- factor
  loads <- round(data$sdev^2/sum(data$sdev^2)*100,1)
  # main plot
  p1 <- ggplot(dataDf,aes(x=PC1, y=PC2)) +
    theme_classic() +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_point(aes(color = Group), alpha = 0.55, size = 3) +
    coord_cartesian(xlim = c(min(data$x[,1])-5,max(data$x[,1])+5)) +
    scale_fill_discrete(name = "Group")
  # avoiding labels superposition
  p1 + geom_text_repel(aes(y = PC2 + 0.25, label = labels),segment.size = 0.25, size = size) + 
    labs(x = c(paste("PC1",loads[1],"%")),y=c(paste("PC2",loads[2],"%"))) +  
    ggtitle(paste("Principal Component Analysis for: ",title,sep=" "))+ 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=colores)
}


plotPCA3(exprs(rawData), labels = targets$SampleName, factor = targets$Group, 
         title="Raw data", scale = FALSE, size = 3,
         colores = c("green", "blue", "red", "yellow"))

boxplot(rawData, cex.axis=0.5, las=2,  which="both", 
        col = c(rep("green", 4), rep("blue", 4), rep("red", 4), rep("yellow", 4)),
        main="Distribution of raw intensity values")

#Normalización
eset_rma <- affy::rma(rawData)

#Control de calidad de los datos normalizados
plotPCA3(exprs(eset_rma), labels = targets$SampleName, factor = targets$Group, 
         title="Normalized data", scale = FALSE, size = 3,
         colores = c("green", "blue", "red", "yellow"))

boxplot(eset_rma, cex.axis=0.5, las=2,  which="all", 
        col = c(rep("green", 4), rep("blue", 4), rep("red", 4), rep("yellow", 4)),
        main="Boxplot for arrays intensity: Normalized Data")

#Filtraje
annotation(eset_rma) <- "hgu133plus2.db"
filtered <- nsFilter(eset_rma, 
                     require.entrez = TRUE, remove.dupEntrez = TRUE,
                     var.filter=TRUE, var.func=IQR, var.cutoff=0.75,filterByQuantile=TRUE, feature.exclude = "^AFFX")
print(filtered$filter.log)

eset_filtered <-filtered$eset

### Matriz de diseño y contrastes
designMat<- model.matrix(~0+Group, pData(eset_filtered))
colnames(designMat) <- c("GFP_tr", "GFP_un", "p53_tr", "p53_un")
print(designMat)

cont.matrix <- makeContrasts (GFP_trvsGFP_un = GFP_tr-GFP_un,
                              p53_trvsp53_un = p53_tr-p53_un,
                              INT = (p53_tr-p53_un)-(GFP_tr-GFP_un),
                              levels=designMat)
print(cont.matrix)

#Genes diferencialmente expresados
fit<-lmFit(eset_filtered, designMat)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)
class(fit.main)

topTab_GFP_trvsGFP_un <- topTable (fit.main, number=nrow(fit.main), coef="GFP_trvsGFP_un", adjust="fdr") 
head(topTab_GFP_trvsGFP_un)

topTab_p53_trvsp53_un <- topTable (fit.main, number=nrow(fit.main), coef="p53_trvsp53_un", adjust="fdr") 
head(topTab_p53_trvsp53_un)

topTab_INT  <- topTable (fit.main, number=nrow(fit.main), coef="INT", adjust="fdr") 
head(topTab_INT, n = 5)

#Anotación de resultados
annotatedTopTable <- function(topTab, anotPackage)
{
  topTab <- cbind(PROBEID=rownames(topTab), topTab)
  myProbes <- rownames(topTab)
  thePackage <- eval(parse(text = anotPackage))
  geneAnots <- select(thePackage, myProbes, c("SYMBOL", "ENTREZID", "GENENAME"))
  annotatedTopTab<- merge(x=geneAnots, y=topTab, by.x="PROBEID", by.y="PROBEID")
  return(annotatedTopTab)
}

topAnnotated_GFP_trvsGFP_un <- annotatedTopTable(topTab_GFP_trvsGFP_un,
                                                 anotPackage="hgu133plus2.db")
topAnnotated_p53_trvsp53_un <- annotatedTopTable(topTab_p53_trvsp53_un,
                                                 anotPackage="hgu133plus2.db")
topAnnotated_INT <- annotatedTopTable(topTab_INT,
                                      anotPackage="hgu133plus2.db")
write.csv(topAnnotated_GFP_trvsGFP_un, file="./results/topAnnotated_GFP_trvsGFP_un.csv")
write.csv(topAnnotated_p53_trvsp53_un, file="./results/topAnnotated_p53_trvsp53_un.csv")
write.csv(topAnnotated_INT, file="./results/topAnnotated_INT.csv")

shortGFP<- head(topAnnotated_GFP_trvsGFP_un[1:5,1:3])
show(shortGFP)

shortp53<- head(topAnnotated_p53_trvsp53_un[1:5,1:3])
show(shortp53)

shortINT<- head(topAnnotated_INT[1:5,1:3])
show(shortINT)

#Comparaciones múltiples
res<-decideTests(fit.main, method="separate", adjust.method="fdr", p.value=0.1, lfc=1)

sum.res.rows<-apply(abs(res),1,sum)
res.selected<-res[sum.res.rows!=0,] 
print(summary(res))

vennDiagram (res.selected[,1:3], cex=0.9)

#Análisis de significación biológica
datasetGFP <- topAnnotated_GFP_trvsGFP_un$logFC
names(datasetGFP) <- as.character(topAnnotated_GFP_trvsGFP_un$ENTREZID)
datasetGFP <- sort(datasetGFP, decreasing = TRUE)

datasetp53 <- topAnnotated_p53_trvsp53_un$logFC
names(datasetp53) <- as.character(topAnnotated_p53_trvsp53_un$ENTREZID)
datasetp53 <- sort(datasetp53, decreasing = TRUE)

geneGFP <- names(datasetGFP)[abs(datasetGFP) > 1]
genep53 <- names(datasetp53)[abs(datasetp53) > 1]

edoGFP <- enrichDGN(geneGFP)
barplot(edoGFP, showCategory=10)

edop53 <- enrichDGN(genep53)
barplot(edop53, showCategory=10)

edoxGFP <- setReadable(edoGFP, 'hgu133plus2.db', 'ENTREZID')
heatplot(edoxGFP, foldChange=datasetGFP, showCategory = 10)

edoxp53 <- setReadable(edop53, 'hgu133plus2.db', 'ENTREZID')
heatplot(edoxp53, foldChange=datasetp53, showCategory = 10)

wpgmtfile <- system.file("extdata/wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
wp2gene <- read.gmt(wpgmtfile)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

ewpGFP <- enricher(geneGFP, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
head(ewpGFP[1:5,2:5])

ewpp53 <- enricher(genep53, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
head(ewpp53[1:5,2:5])

