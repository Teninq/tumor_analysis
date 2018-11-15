##generate a classifier
setwd("D:/TCGA/TCGA_COAD_exp_HiSeqV2-2015-02-24")
coad.tcga <- read.table("genomicMatrix",sep="\t",header = TRUE,row.names = 1)
coad.tcga <- coad.tcga[which(apply(coad.tcga,1,sum)>0),]
coad.tcga <- log2(coad.tcga+0.01)

coad.tcga.tumor <- coad.tcga[,grep(".0[0-9]$",colnames(coad.tcga))]
coad.tumor.zscore <- scale(coad.tcga.tumor)

coad.tumor.var <- sort(apply(coad.tcga.tumor,1,var),decreasing = TRUE,index.return = TRUE)
coad.tumor.mad <- sort(apply(coad.tcga.tumor,1,mad),decreasing = TRUE,index.return = TRUE)
index <- tumor.var$ix[1:5000]
gene <- row.names(tumor)[index]
index1 <- tumor.mad$ix[1:5000]
gene1 <- row.names(tumor)[index1]
tumor.classify <- coad.tumor.zscore[match(intersect(gene,gene1),colnames(coad.tumor.zscore)),]
res <- ConsensusClusterPlus(tumor.classify,maxK=8,pItem = 0.8,pFeature = 0.8,
                            reps = 1000,distance = "euclidean",clusterAlg = "km")
calcICL(res,title="untitled_consensus_cluster",plot=NULL,writeTable=FALSE)

library(siggenes)
sam.coad.out<-sam(coad.tumor.zscore,res[[4]]$consensusClass,rand=123,gene.names=rownames(coad.tumor.zscore))












#clinical data
clinic <- read.table("clinical_data",sep = "\t",header = TRUE)
data.sample <- colnames(data)
clc.sample <- chartr("-",".",clinic$sampleID)
clinic <- clinic[match(data.sample,clc.sample),]
tumor.clinic <- clinic[grep(".0[0-9]$",clinic$sampleID),]
