library(NOISeq)
library(GO.db)
library(AnnotationHub)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(gplots)
library(plotly)
library(topGO)
library(Rgraphviz)
library(limma)
library(biomaRt)

# Se leen los datos y se crea un data frame
datospac <- read.csv('Paci.tabular', header = TRUE, sep="\t",stringsAsFactors = FALSE)
datospac <- as.data.frame(datospac[1:7])


rownames(datospac)<- datospac$tracking_id



datospac<- data.frame(datospac[,2:7])

datos_reducidos <- datospac


factores = data.frame(Muestras= c("SANO","SANO","SANO","CANCER","CANCER","CANCER"))
misdatos <- NOISeq::readData(data = datos_reducidos,  factors = factores)
head(assayData(misdatos)$exprs)

#COMPARACIÓN sano, cancer
comparacion=c("SANO","CANCER")

mynoiseq1 = noiseq(misdatos, conditions=comparacion,k = 0.5, norm="n", factor = "Muestras",nss = 0, lc = 1, replicates = "no")
    #con esto nos genera la matriz de la diapo 66 del tema.

#######LO QUE VIENE A PARTIR DE AHORA SON FILTROS SOBRE ESA MATRIZ GENERADA

# Genes expresados diferencialmente 
mynoiseq1.deg = degenes(mynoiseq1, q = 0.9, M = NULL) #Aui valen los expresados up y down regulated


# Genes upregulated
mynoiseq1.deg1 = degenes(mynoiseq1, q = 0.9, M = "up")

# Genes downregulated
mynoiseq1.deg2 = degenes(mynoiseq1, q = 0.9, M = "down")

#Expression plot
png("expresionplot1.png",width = 600, height = 600)
DE.plot(mynoiseq1, q = 0.9, graphic = "expr", log.scale = TRUE)
dev.off() #en engro no tienen diferencial y rojo si, los negros suelen estar en la diagonal

#MD plot.
png("MD.png",width = 600, height = 600)
DE.plot(mynoiseq1, q = 0.9, graphic = "MD")
dev.off() #Este es similar down a la izq y up regulated a la derecha, cuanto mas alto mas dif en valor absoluto (de izq a derecah en valor relativo)



geneList<-rownames(mynoiseq1.deg)
# Convertir los identificadores
gene.df <- bitr(geneList, fromType = "REFSEQ", 
                toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)


# Agrupar los t??rminos go al nivel de profundidad deseado, cuanto m??s grande t??rminos m??s espec??ficos aparecer??n
# De momento no estamos enriqueciendo

ggo <- groupGO(gene     =  geneList,
               OrgDb    = org.Hs.eg.db,
               ont      = "MF",
               keyType = "REFSEQ",
               level    = 4,
               readable = TRUE)
barplot(ggo, drop=TRUE, showCategory=10)



# Enriquecimiento GO
ego2 <- enrichGO(gene         = geneList,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'REFSEQ',
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)

# Mapeamos ids a SYMBOL
ego2 <- setReadable(ego2, OrgDb = org.Hs.eg.db)


#Gr??ficos del Enriquecimiento
barplot(ego2)
dotplot(ego2)

gsea_genes=data.frame(rownames(mynoiseq1.deg),mynoiseq1.deg[3:3])

geneList = gsea_genes[,2]
## feature 2: named vector
names(geneList) = as.character(gsea_genes[,1])
## feature 3: decreasing orde
geneList = sort(geneList, decreasing = TRUE)

ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              keyType       = 'REFSEQ',
              ont          = "MF",
        #      nPerm        = 1000,
              minGSSize    = 10,
              maxGSSize    = 500,
              pvalueCutoff = 0.1,
              verbose      = FALSE)

ego3 <- setReadable(ego3, OrgDb = org.Hs.eg.db)
dotplot(ego3)






########### HEAT MAP ########################

x<-as.matrix(datos_reducidos[1:6])
datos_heatmap<-apply(x, 2, as.double)
rownames(datos_heatmap)<-rownames(datos_reducidos)
group <- as.factor(c("SAN","SAN","SAN", "CAN", "CAN", "CAN" ))
# creates a own color palette from red to green
mycol <- colorpanel(9,"red","green")
png("heatmap.png",width = 1024, height = 1024)
#rowside_colors=c(rep("black", 29), rep("blue", 38))
col_breaks = c(seq(-1,0,length=5),  # for red
               seq(0.01,1,length=5) )

heatmap.2(datos_heatmap,
          col = mycol,
          xlab = "Experiment",
          ylab = "Genes",
          main = "Heatmap by expression",
          labRow = rownames(datos_heatmap),
          labCol = group,
          dendrogram = c("both"),
          margins=c(12,8),
          trace="column",
          tracecol="black",
          #  breaks=col_breaks,
          density.info="none",
          scale = "row",
          notecol="black",
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          distfun=function(x) as.dist(1-cor(t(x))),
          #  distfun = dist(datos_fusionados_FC, method = "manhattan"),
          hclustfun=function(x) hclust(x, method="complete")
          #     RowSideColors = rowside_colors ,   # line width
          
)

dev.off()

