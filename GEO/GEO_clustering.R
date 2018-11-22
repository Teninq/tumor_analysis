setwd("D:/GEO/COAD/")
library(GEOquery)
# library(arrayQualityMetrics)
library(pcaMethods)
library(ggplot2)
library(sva)
library(ConsensusClusterPlus)

lncRNA.probe <- read.table("lncRNA_Probe.txt",sep = "\t",header = F)
table(row.names(expr.GSE14333) %in% lncRNA.probe$V1)
expr.lnc14333 <- expr.GSE14333[row.names(expr.GSE14333) %in% lncRNA.probe$V1,]
#GSE5851 <- getGEO(filename = "D:/GEO/COAD/GSE5851/GSE5851_series_matrix.txt.gz",GSEMatrix = TRUE)
load("GSE5851.RData")
phenotype.GSE5851 <- pData(GSE5851@phenoData)
expr.GSE5851 <- GSE5851@assayData$exprs
expr.GSE5851 <- log2(expr.GSE5851)
expr.GSE5851 <- expr.GSE5851[apply(expr.GSE5851,1,function(x){sum(is.na(x))<(length(x)*0.9)}),]

#GSE14333 <- getGEO(filename = "D:/GEO/COAD/GSE14333/GSE14333_series_matrix.txt.gz",GSEMatrix = TRUE)
load("GSE14333.RData")
expr.GSE14333 <- GSE14333@assayData$exprs
expr.GSE14333 <- expr.GSE14333[apply(expr.GSE14333,1,function(x){sum(is.na(x))<(length(x)*0.9)}),]

#GSE17538 <- getGEO("GSE17538",destdir = "D:/GEO/COAD/GSE17538/",GSEMatrix = TRUE)
load("GSE17538.RData")
expr.GSE17538 <- GSE17538$`GSE17538-GPL570_series_matrix.txt.gz`@assayData$exprs
expr.GSE17538 <- expr.GSE17538[apply(expr.GSE17538,1,function(x){sum(is.na(x))<(length(x)*0.9)}),]

#GSE13294 <- getGEO("GSE13294",destdir = "D:/GEO/COAD/GSE13294/",GSEMatrix = TRUE)
load("GSE13294.RData")
expr.GSE13294 <- GSE13294$GSE13294_series_matrix.txt.gz@assayData$exprs
expr.GSE13294 <- expr.GSE13294[apply(expr.GSE13294,1,function(x){sum(is.na(x))<(length(x)*0.9)}),]

#############################batch effect remove
expr.GSE13294 <- expr.GSE13294[match(intersect(row.names(expr.GSE5851),row.names(expr.GSE13294)),
                                     row.names(expr.GSE13294)),]
expr.GSE14333 <- expr.GSE14333[match(row.names(expr.GSE13294),row.names(expr.GSE14333)),]
expr.GSE5851 <- expr.GSE5851[match(row.names(expr.GSE14333),row.names(expr.GSE5851)),]
expr.GSE17538 <- expr.GSE17538[match(row.names(expr.GSE14333),row.names(expr.GSE17538)),]

expr.coad <- cbind(expr.GSE13294,expr.GSE14333,expr.GSE17538)
expr.coad.zscore <- scale(expr.coad)

data.label <- c(rep("GSE13294",ncol(expr.GSE13294)),rep("GSE14333",ncol(expr.GSE14333)),
                rep("GSE17538",ncol(expr.GSE17538)))
names(data.label) <- c(colnames(expr.GSE13294),colnames(expr.GSE14333),
                       colnames(expr.GSE17538))
#coad.result <- pca(t(expr.coad),method = "bpca",nPcs = 2)
#coad.result1 <- pca(expr.coad,method = "bpca",nPcs = 2)
#coad.result2 <- pca(coad.result1@completeObs,method = "bpca",nPcs = 2)

pca.zscore <- pca(expr.coad.zscore,method = "bpca",nPcs = 2)
pca.coad <- cbind(pca.zscore@loadings,data.label)
ggplot(data.frame(pca.coad), aes(PC1, PC2, shape=data.label, color=data.label)) +
  geom_point(size = 4) +
  labs(x=paste("PC1", pca.zscore@R2[1] * 100, "% of variance"),
       y=paste("PC2", pca.zscore@R2[2] * 100, "% of variance"),
       title="Before remove batch effect") +
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 20),
        axis.text.x= element_blank(),
        axis.text.y = element_blank())

expr.combat <- ComBat(pca.zscore@completeObs,data.label)
combat.pca <- pca(expr.combat,method = "bpca",nPcs = 2)
pca.coad.combat <- cbind(combat.pca@loadings,data.label)
ggplot(data.frame(pca.coad.combat), aes(PC1, PC2, shape=data.label, color=data.label)) +
  geom_point(size = 4) +
  labs(x=paste("PC1", combat.pca@R2[1] * 100, "% of variance"),
       y=paste("PC2", combat.pca@R2[2] * 100, "% of variance"),
       title="After remove batch effect") +
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 20),
        axis.text.x= element_blank(),
        axis.text.y = element_blank())

#################################

############
expr.combat.mad <- sort(apply(expr.combat,1,mad),decreasing = TRUE,index.return = TRUE)
cluster.matrix <- expr.combat[expr.combat.mad$x>0.5,]
load("clusterresult.RData")
# coad.res <- ConsensusClusterPlus(cluster.matrix,maxK=8,pItem = 0.9,pFeature = 0.9,
#                                  reps = 1000,distance = "euclidean",clusterAlg = "km")
# calcICL(coad.res,title="untitled_consensus_cluster",plot=NULL,writeTable=FALSE)

# ###
# index <- expr.combat.mad$ix[1:5000]
# probe <- row.names(expr.combat)[index]
# expr.cluster.matrix <- expr.combat[probe,]
# coad.res.1 <- ConsensusClusterPlus(expr.cluster.matrix,maxK=8,pItem = 0.9,pFeature = 0.9,
#                                    reps = 1000,distance = "euclidean",clusterAlg = "km")
# calcICL(coad.res.1,title="untitled_consensus_cluster",plot=NULL,writeTable=FALSE)

subcluster.list <- list(one=names(coad.res[[5]]$consensusClass)[coad.res[[5]]$consensusClass==1],
                        two=names(coad.res[[5]]$consensusClass)[coad.res[[5]]$consensusClass==2],
                        three=names(coad.res[[5]]$consensusClass)[coad.res[[5]]$consensusClass==3],
                        four=names(coad.res[[5]]$consensusClass)[coad.res[[5]]$consensusClass==4],
                        five=names(coad.res[[5]]$consensusClass)[coad.res[[5]]$consensusClass==5])
subcluster.expression.list <- sapply(subcluster.list,function(x){expr.combat[,match(x,colnames(expr.combat))]})

x <-c(sapply(subcluster.list,length))
label <-c(names(subcluster.list))
piepercent<-round(100*x/sum(x), 1)
piepercent <-paste(piepercent, "%", sep = "")
pie(x,labels=piepercent, main="Subtype pie chart",col= terrain.colors(length(x)))
legend("topright",label, cex=0.8, fill=terrain.colors(length(x)))
#legend(locator(1),label, cex=0.8, fill=terrain.colors(length(x)))


##############################################
# GSE5851.meta <- Meta(GSE5851)
# GSE5851.gsmlist <- GSMList(GSE5851)
# GSE5851.gpllist <- GPLList(GSE5851)
# sample1 <- Table(GSE5851.gsmlist[[1]])
# sample1.annotation <- Table(GSE5851.gpllist[[1]])
# probesets<-Table(GSE5851.gpllist[[1]])$ID
# GSEMerge <- do.call("cbind",lapply(GSE5851.gsmlist,function(x){
#   tab <- Table(x)
#   #mymatch <- match(probeID,tab$ID_REF)
#   return(tab$VALUE)
# }))