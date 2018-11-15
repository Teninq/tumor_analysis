setwd("D:/TCGA/TCGA_COAD_exp_HiSeqV2-2015-02-24")
######################library#########################
library(pipeR)
library(ConsensusClusterPlus)
library(limma)
library(pheatmap)
library(survival)
library(randomForest)
library(gplots)
library(scales)
library(data.table)
library(SNFtool)

##############################Data Preparation###############################
data <- read.table("genomicMatrix",sep="\t",header = TRUE,row.names = 1)
data <- data[which(apply(data,1,sum)>0),]
data <- log2(data+0.01)

#clinical data
clinic <- read.table("clinical_data",sep = "\t",header = TRUE)
data.sample <- colnames(data)
clc.sample <- chartr("-",".",clinic$sampleID)
clinic <- clinic[match(data.sample,clc.sample),]
tumor.clinic <- clinic[grep(".0[0-9]$",clinic$sampleID),]

#clinical label
CIMP=tumor.clinic$CIMP
mirsta <- tumor.clinic$CDE_ID_3226963
names(mirsta) <- tumor.clinic$sampleID
Anatomic=tumor.clinic$anatomic_neoplasm_subdivision
location=as.character(Anatomic)
location[which(location=="Cecum")]="RCC"
location[which(location=="Ascending Colon")]="RCC"
location[which(location=="Transverse Colon")]="RCC"
location[which(location=="Hepatic Flexure")]="RCC"
location[which(location=="Descending Colon")]="LCC"
location[which(location=="Sigmoid Colon")]="LCC"
location[which(location=="Splenic Flexure")]="LCC"
location[which(location=="Rectosigmoid Junction")]="LCC"
location[which(location=="")]="NA"

#select active genes
tumor <- data[,grep(".0[0-9]$",colnames(data))]
normal <- data[,grep(".1[0-9]$",colnames(data))]
tumor.var <- sort(apply(tumor,1,var),decreasing = TRUE,index.return = TRUE)
tumor.mad <- sort(apply(tumor,1,mad),decreasing = TRUE,index.return = TRUE)
index <- tumor.var$ix[1:5000]
gene <- row.names(tumor)[index]
index1 <- tumor.mad$ix[1:5000]
gene1 <- row.names(tumor)[index1]
#length(intersect(gene,gene1))

tumor.zscore <- scale(tumor)
# row.names(tumor.zscore) <- row.names(tumor)
# colnames(tumor.zscore) <- colnames(tumor)

tumor.matrix.zscore <- tumor.zscore[match(intersect(gene,gene1),row.names(data)),]
# normal.matrix.zscore <- normal.zscore[match(intersect(gene,gene1),row.names(data)),]
# COAD.matrix <- data[match(intersect(gene,gene1),row.names(data)),]

################################classification##############################
##ConsensusCluster
#load("cluster.RData")
#load("coad_cluster(log).RData")
res <- ConsensusClusterPlus(tumor.matrix.zscore,maxK=8,pItem = 0.8,pFeature = 0.8,
                            reps = 1000,distance = "euclidean",clusterAlg = "km")
calcICL(res,title="untitled_consensus_cluster",plot=NULL,writeTable=FALSE)
#save(res,file = "coad_cluster(log).RData")

###Survival Analysis
clusterID <- data.frame(res[[4]]$consensusClass)
clusterID$sampleID <- rownames(clusterID)
tumor.clinic$sampleID <- chartr("-",".",tumor.clinic$sampleID)
tumor.clinic1 <- merge(tumor.clinic,clusterID,by="sampleID",all.y = TRUE)

event <- tumor.clinic1$X_EVENT
days <- tumor.clinic1$X_OS
sdf <- survdiff(Surv(days, event) ~ tumor.clinic1$res..4...consensusClass)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
kmsurvival <- survfit(Surv(days, event) ~ tumor.clinic1$res..4...consensusClass)

plot(kmsurvival,lwd=2, col=rainbow(4),xlab="Time",ylab="Survival Probability",
     main=paste("Survival Curse",'\n',"P=",format(pvalue,digits=2)) ,mark=3)
legend("topright",legend=names(table(tumor.clinic1$res..4...consensusClass)),col=rainbow(4),lty=1,bty="n")

#Feature Selection
setwd("D:/TCGA/TCGA_COAD_miRNA_HiSeq-2015-02-24")
matrix <-  read.table("genomicMatrix",sep="\t",header = TRUE,row.names = 1)
matrix <- matrix[which(apply(matrix,1,sum)>0),]
#matrix <- log2(as.matrix(matrix)+0.01)
#data_qt=normalizeQuantiles(data)
matrix.zscore <- t(apply(matrix,1,scale))
row.names(matrix.zscore) <- row.names(matrix)
colnames(matrix.zscore) <- colnames(matrix)

##############################Differential Expression##############################
###unpaired differential expression
COAD.matrix <- matrix.zscore[,na.omit(match(colnames(data),colnames(matrix.zscore)))]
tumor.rna <- COAD.matrix[,grep(".0[0-9]$",colnames(COAD.matrix))]
normal.rna <- COAD.matrix[,grep(".1[0-9]$",colnames(COAD.matrix))]
# paired.tumor <- tumor.rna[,na.omit(match(sub(".1[0-9]$","",colnames(normal.rna)),
#                                          sub(".0[0-9]$","",colnames(tumor.rna))))]
# paired.normal <- normal.rna[,na.omit(match(sub(".0[0-9]$","",colnames(paired.tumor)),
#                                            sub(".1[0-9]$","",colnames(normal.rna))))]
rna.COAD <- cbind(normal.rna,tumor.rna[row.names(tumor.rna),])
p.value <- apply(rna.COAD,1,function(x){t.test(x[1:8],x[9:254])$p.value})
p.adjust <- p.adjust(p.value,method="fdr",n=length(p.value))
feature.matrix <- tumor.rna[p.adjust<0.05,]

#############################################

###limma
# design <- rep("Tumor",length(colnames(COAD.matrix)))
# design[grep("1[0-9]$",colnames(matrix.zscore))] = "Normal"
# design.mat <- factor(design, levels = c("Tumor","Normal"))
# design.matrix <- model.matrix(~0 + design.mat)
# colnames(design.matrix) <- levels(design.mat)
# data.fit <- lmFit(COAD.matrix,design.matrix)
# cont.matrix <- makeContrasts(data.TvsN=Tumor-Normal,levels = design.matrix)
# data.fit2 <- contrasts.fit(data.fit, cont.matrix)
# data.fit2 <- eBayes(data.fit2)
# diff.result <- topTable(data.fit2,coef=1,n=nrow(data.fit2),lfc=log2(1))
# diff.result.final <- diff.result[diff.result[,"P.Value"]<0.05,]
# nrow(diff.result.final) #470
# up.sample=diff.result.final[diff.result.final[,"logFC"]>0,]
# #up.gene <- unique(up.sample[,3])
# nrow(up.sample)
# down.sample=diff.result.final[diff.result.final[,"logFC"]<0,]
# #down.gene <- unique(down.sample[,3])
# nrow(down.sample)

clinical=read.table("clinical_data",header=T,sep="\t")
clc_sample = chartr("-",".",clinical$sampleID)
sample_id=intersect(colnames(feature.matrix),clc_sample)
clinic = clinical[match(sample_id,clc_sample),]
c_time=clinic$X_OS
c_event=clinic$X_EVENT
cox_pvalue <- apply(feature.matrix,1,function(x){y=coxph(Surv(c_time,c_event)~x);y=summary(y)$sctest[3]})
choose_cox_matrix <- feature.matrix[cox_pvalue<0.05,]
cox2 <- coxph(Surv(c_time, c_event)~t(choose_cox_matrix))
coefficient=summary(cox2)$coef[,1]
risk_score=t(choose_cox_matrix)%*%coefficient
low_score=order(risk_score)[1:(length(risk_score)/2)]
hight_score=setdiff(order(risk_score),low_score)
lable=c()
lable[low_score]="L"
lable[hight_score]="H"
names(lable) <- colnames(choose_cox_matrix)

##survival curse
sdf1 <- survdiff(Surv(c_time, c_event) ~ lable)
pvalue <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
kmsurvival1 <- survfit(Surv(c_time, c_event) ~ lable)

plot(kmsurvival1,lwd=2, col=rainbow(2),xlab="Time",ylab="Survival Probability",
     main=paste("Survival Curse",'\n',"P=",format(pvalue,digits=2)) ,mark=3)
legend("topright",legend=names(table(lable)),col=rainbow(2),lty=1,bty="n")

##map the label to mRNA
label <- rep(NA,288)
names(label) <- colnames(tumor.matrix.zscore)
label[names(lable)] <- lable
label <- as.matrix(label)

#Heatmap
# forest.data <- rbind(tumor.matrix.zscore[,colnames(tumor.matrix.zscore)],res[[4]]$consensusClass)
# row.names(forest.data)=c(rownames(tumor.matrix.zscore),"cluster")
# forest.data <- t(forest.data)
# forest.data.m <- randomForest(cluster ~.,data = forest.data,importance = TRUE)
# for(i in 1:100){
#   forest.data.m <- randomForest(cluster ~.,data = forest.data,mtry=i,importance = TRUE)
#   i
# }

design.cluster <- sub("1","One",res[[4]]$consensusClass)
design.cluster <- sub("2","Two",design.cluster)
design.cluster <- sub("3","Three",design.cluster)
design.cluster <- sub("4","Four",design.cluster)
cluster.mat <- factor(design.cluster, levels = c("One","Two","Three","Four"))
cluster.matrix <- model.matrix(~0 + cluster.mat)
colnames(cluster.matrix) <- levels(cluster.mat)
cluster.fit <- lmFit(tumor.zscore,cluster.matrix)
cont.cluster <- makeContrasts(One-(Two+Three+Four),
                              Two-(One+Three+Four),
                              Three-(One+Two+Four),
                              Four-(One+Two+Three),
                              levels = cluster.matrix)
cluster.fit2 <- contrasts.fit(cluster.fit, cont.cluster)
cluster.fit2 <- eBayes(cluster.fit2)

##decidetest
diff.res <- decideTests(cluster.fit2,method = "global",adjust.method = "BH",
                        p.value = 0.01,lfc = log2(2))
summary(diff.res)

###Cluster One
One.diff <- topTable(cluster.fit2,coef=1,n=nrow(cluster.fit2),lfc=log2(1.5))
One.diff.final <- One.diff[One.diff[,"adj.P.Val"]<0.01,]
nrow(One.diff.final)
up.sample.One <- One.diff.final[One.diff.final[,"logFC"]>0,]
nrow(up.sample.One)

###Cluster Two
Two.diff <- topTable(cluster.fit2,coef=2,n=nrow(cluster.fit2),lfc=log2(1.5))
Two.diff.final <- Two.diff[Two.diff[,"adj.P.Val"]<0.01,]
nrow(Two.diff.final)
up.sample.Two <- Two.diff.final[Two.diff.final[,"logFC"]>0,]
nrow(up.sample.Two)

###Cluster Three
Three.diff <- topTable(cluster.fit2,coef=3,n=nrow(cluster.fit2),lfc=log2(1.5))
Three.diff.final <- Three.diff[Three.diff[,"adj.P.Val"]<0.01,]
nrow(Three.diff.final)
up.sample.Three=Three.diff.final[Three.diff.final[,"logFC"]>0,]
nrow(up.sample.Three)

###Cluster Four
Four.diff <- topTable(cluster.fit2,coef=4,n=nrow(cluster.fit2),lfc=log2(1.5))
Four.diff.final <- Four.diff[Four.diff[,"adj.P.Val"]<0.01,]
nrow(Four.diff.final)
up.sample.Four <-Four.diff.final[Four.diff.final[,"logFC"]>0,]
nrow(up.sample.Four)

up.genes <- list(one=row.names(up.sample.One),
                 two=row.names(up.sample.Two),
                 three=row.names(up.sample.Three),
                 four=row.names(up.sample.Four))
gplots::venn(up.genes)

one.specific.gene <- setdiff(setdiff(setdiff(row.names(up.sample.One),
                                             row.names(up.sample.Two)),
                                     row.names(up.sample.Three)),
                             row.names(up.sample.Four))
two.specific.gene <- setdiff(setdiff(setdiff(row.names(up.sample.Two),
                                             row.names(up.sample.One)),
                                     row.names(up.sample.Three)),
                             row.names(up.sample.Four))
three.specific.gene <- setdiff(setdiff(setdiff(row.names(up.sample.Three),
                                               row.names(up.sample.One)),
                                       row.names(up.sample.Two)),
                               row.names(up.sample.Four))
four.specific.gene <- setdiff(setdiff(setdiff(row.names(up.sample.Four),
                                              row.names(up.sample.One)),
                                      row.names(up.sample.Two)),
                              row.names(up.sample.Three))
all.specific.gene <- c(one.specific.gene,two.specific.gene,three.specific.gene,four.specific.gene)
specific.gene.list <- list(one=one.specific.gene,
                           two=two.specific.gene,
                           three=three.specific.gene,
                           four=four.specific.gene)
tumor.all.specific <- tumor.matrix.zscore[match(all.specific.gene,row.names(tumor.matrix.zscore)),]

all.diff <- union(union(union(row.names(up.sample.Four),
                              row.names(up.sample.One)),
                        row.names(up.sample.Three)),
                  row.names(up.sample.Two))
all.diff.gene <- union(union(union(row.names(Four.diff.final),
                                   row.names(One.diff.final)),
                             row.names(Three.diff.final)),
                       row.names(Two.diff.final))
tumor.all.diff <- tumor.matrix.zscore[match(all.diff,row.names(tumor.matrix.zscore)),]

cluster.result <- data.frame(cbind(names(res[[4]]$consensusClass),
                                   res[[4]]$consensusClass,
                                   res[[2]]$consensusClass,
                                   label[row.names(label),]))
cluster.result.final <- dplyr::arrange(cluster.result,X2)
annotation.col <- data.frame(cluster=factor(cluster.result.final$X2,labels = c("1","2","3","4")),
                             cluster1=factor(cluster.result.final$X3,labels = c("1","2")),
                             label1=factor(cluster.result$X4,labels = c("H","L")),
                             location_label=factor(location,labels = c("LCC","NA","RCC")),
                             microsatellite=factor(mirsta,labels = c("","Indeterminate","MSI-H","MSI-L","MSS")))
cluster <- c("red","green","yellow","blue")
names(cluster) <- c("1","2","3","4")
cluster1 <- c("green","yellow")
names(cluster1) <- c("1","2")
label1 <- c("red","blue")
names(label1) <- c("H","L")
location_label <- c("yellow","green","black")
names(location_label) <- c("RCC","LCC","NA")
microsatellite <- c("white","black","red","blue","yellow")
names(microsatellite) <- c("","Indeterminate","MSI-H","MSI-L","MSS")
ann.colors <- list(cluster=cluster,cluster1=cluster1,label1=label1,
                   location_label=location_label,microsatellite=microsatellite)
rownames(annotation.col) <- cluster.result.final$X1

###plot heatmap
pheatmap(tumor.all.specific[,as.character(cluster.result.final$X1)],
         annotation_col=annotation.col,
         annotation_colors=ann.colors,
         cluster_cols = FALSE,
         cluster_rows =FALSE,
         show_rownames = F, 
         show_colnames = F,
         color = redblue(100)[100:1])
#colorRampPalette(c("blue","white","red"))(25)
pheatmap(t(apply(tumor.all.specific[,as.character(cluster.result.final$X1)],1,function(x){rescale(x,to=c(-3,3))})),
         annotation_col=annotation.col,
         annotation_colors=ann.colors,
         cluster_cols = FALSE,
         cluster_rows =FALSE,
         show_rownames = F, 
         show_colnames = F,
         color =redblue(100)[100:1])

pheatmap(t(coad),show_colnames = F,color = brewer.pal(9,"YlOrRd"),
         cluster_cols = T,clustering_method = "complete",
         main="immune",fontsize = 10, fontsize_row = 8, fontsize_col = 9)

###pairwise
cont.matrix <- makeContrasts(One-Two,One-Three,One-Four,Two-Three,Two-Four,Three-Four,levels = cluster.matrix)
matrix.fit <- contrasts.fit(cluster.fit,cont.matrix)
matrix.fit2 <- eBayes(matrix.fit)


###########################SNF Analysis######################
setwd("D:/TCGA/methylation")
# coad.methylation <- fread("COAD.txt")
# row.names(coad.methylation) <- coad.methylation$V1
# coad.methylation <- coad.methylation[,-1]
load("coad_methy.RData")
coad.methylation <- as.matrix(coad.methylation)
coad.methy <- coad.methylation[,match(intersect(colnames(tumor),colnames(coad.methylation)),colnames(coad.methylation))]
coad.rna <- tumor[,match(intersect(colnames(tumor),colnames(coad.methylation)),colnames(tumor))]

coad.rna.norm <- standardNormalization(t(coad.rna))
coad.methy.norm <- standardNormalization(t(coad.methy))
Dist.rna <- dist2(as.matrix(coad.rna.norm),as.matrix(coad.rna.norm))
Dist.methy <- dist2(as.matrix(coad.methy.norm),as.matrix(coad.methy.norm))

K = 35
alpha = 0.7
T = 35
C = 5
W1 <- affinityMatrix(Dist.rna,K,alpha)
W2 <- affinityMatrix(Dist.methy,K,alpha)
W = SNF(list(W1,W2), K, T)
SNF.labels = spectralClustering(W, C)
names(SNF.labels) <- row.names(Dist.methy)
table(SNF.labels)
displayClusters(W, SNF.labels)
displayClusters(W1,SNF.labels)
displayClusters(W2,SNF.labels)

plot(coad.rna.norm, col=SNF.labels, main='RNA')
plot(coad.methy.norm, col=SNF.labels, main='Methylation')


cluster.result.SNF <- data.frame(cbind(names(res[[4]]$consensusClass)[na.omit(match(names(SNF.labels),
                                                                                    names(res[[4]]$consensusClass)))],
                                   res[[4]]$consensusClass[na.omit(match(names(SNF.labels),names(res[[4]]$consensusClass)))],
                                   SNF.labels[names(SNF.labels)]))
cluster.result.final.SNF <- dplyr::arrange(cluster.result.SNF,X2)
annotation.col.SNF <- data.frame(cluster.SNF=factor(cluster.result.final.SNF$X2[na.omit(match(names(SNF.labels),
                                                                                              names(res[[4]]$consensusClass)))],
                                                    labels = c("1","2","3","4")),
                             cluster1.SNF=factor(SNF.labels,labels = c("1","2","3","4","5")))
cluster.SNF <- c("red","green","yellow","blue")
names(cluster.SNF) <- c("1","2","3","4")
cluster1.SNF <- c("green","yellow","blue","red","pink")
names(cluster1.SNF) <- c("1","2","3","4","5")

ann.colors.SNF <- list(cluster.SNF=cluster.SNF,cluster1.SNF=cluster1.SNF)
rownames(annotation.col.SNF) <- cluster.result.final.SNF$X1

pheatmap(t(apply(tumor.all.specific[,as.character(cluster.result.final.SNF$X1)],1,function(x){rescale(x,to=c(-3,3))})),
         annotation_col=annotation.col.SNF,
         annotation_colors=ann.colors.SNF,
         cluster_cols = FALSE,
         cluster_rows =FALSE,
         show_rownames = F, 
         show_colnames = F,
         color = redblue(10)[10:1])


#################################drug response###################################
library(affy)
library(affyPLM)
fns <-list.celfiles(path="D:/cell line/coad_ic50",full.names=TRUE)
data.affy <- ReadAffy(filenames=fns)
data.rma <- rma(data.affy)
expr.rma <- exprs(data.rma)


library(hgu219.db)
#library(hgu219hsrefseq.db)
# Annot <- data.frame(REFSEQ=sapply(contents(hgu219REFSEQ), paste, collapse=", "),
#                     SYMBOL=sapply(contents(hgu219SYMBOL), paste, collapse=", "),
#                     DESC=sapply(contents(hgu219GENENAME), paste, collapse=", "))
probeset <- rownames(expr.rma)
Symbol <- as.character(as.list(hgu219SYMBOL[probeset]))
expr.symbol <- cbind.data.frame(Symbol,expr.rma)
rmDupID <-function(a){
  exprSet <- a[,-1]
  rowMeans <- apply(exprSet,1,function(x) mean(as.numeric(x),na.rm=T))
  a <- a[order(rowMeans,decreasing=T),]
  exprSet <- a[!duplicated(a[,1]),]
  #exprSet=apply(exprSet,2,as.numeric)
  exprSet <- exprSet[-grep("^NA$",exprSet[,1]),]
  rownames(exprSet) <- exprSet[,1]
  exprSet <- exprSet[,-1]
  return(exprSet)
}
exprSet <- rmDupID(expr.symbol)

gene.label <- c(rep("one",length(na.omit(match(one.specific.gene,row.names(exprSet))))),
                rep("two",length(na.omit(match(two.specific.gene,row.names(exprSet))))),
                rep("three",length(na.omit(match(three.specific.gene,row.names(exprSet))))),
                rep("four",length(na.omit(match(four.specific.gene,row.names(exprSet))))))
gene.label.mat <- factor(gene.label, levels = c("one","two","three","four"))
gene.label.matrix <- model.matrix(~0 + gene.label.mat)
colnames(gene.label.matrix) <- levels(gene.label.mat)
row.names(gene.label.matrix) <- c(one.specific.gene[na.omit(match(row.names(exprSet),one.specific.gene))],
                                  two.specific.gene[na.omit(match(row.names(exprSet),two.specific.gene))],
                                  three.specific.gene[na.omit(match(row.names(exprSet),three.specific.gene))],
                                  four.specific.gene[na.omit(match(row.names(exprSet),four.specific.gene))])
specific.gene.expr <- exprSet[na.omit(match(all.specific.gene,row.names(exprSet))),]

###NTP Classification Predict
NTP_Fuc <- function(ref.exp,template,nresmpl=1000,dist.selection = "cosine"){
  num.samples <- ncol(ref.exp)
  num.features.extract <- nrow(template)
  num.cls <- ncol(template)
  predict.label<-vector(length=num.samples,mode="numeric")
  dist.to.template<-vector(length=num.samples,mode="numeric")
  dist.to.cls1<-vector(length=num.samples,mode="numeric")
  
  rnd.feature.matrix<-matrix(0,nrow=num.features.extract,ncol=nresmpl)
  sample.names<-as.vector(as.matrix(colnames(ref.exp)))
  expression <- t(apply(ref.exp,1,scale))
  row.names(expression) <- row.names(ref.exp)
  colnames(expression) <- colnames(ref.exp)
  
  perm.dist.vector<-vector(length=nresmpl*num.cls,mode="numeric")
  nominal.p<-vector(length=num.samples,mode="numeric")
  BH.FDR<-vector(length=num.samples,mode="numeric")
  Bonferroni.p<-vector(length=num.samples,mode="numeric")
  
  for (i in 1:num.cls){
    temp.temp<-as.numeric(as.vector(template[,i]))
    eval(parse(text=paste("temp.",i,"<-temp.temp",sep="")))
  }
  
  for (i in 1:num.samples){
    
    print(paste("sample # ",i,sep=""))
    
    current.sample <- as.vector(expression[,i])
    
    # compute original distance
    
    orig.dist.to.all.temp <- vector(length=num.cls,mode="numeric")
    
    if (dist.selection=="cosine"){
      for (o in 1:num.cls){      # compute distance to all templates
        eval(parse(text=paste("current.temp <- temp.",o,sep="")))
        orig.dist.to.all.temp[o]<-sum(current.temp*current.sample)/
          (sqrt(sum(current.temp^2))*sqrt(sum(current.sample^2)))
      }
    }
    if (dist.selection=="correlation"){
      for (o in 1:num.cls){      # compute distance to all templates
        eval(parse(text=paste("current.temp <- temp.",o,sep="")))
        orig.dist.to.all.temp[o] <- cor(current.temp,current.sample,method="pearson",use="complete.obs")
      }
    }
    for (o in 1:num.cls){       # find nearest neighbor (>2 classes)
      if (is.na(orig.dist.to.all.temp[o])!=T){
        if (orig.dist.to.all.temp[o]==max(orig.dist.to.all.temp,na.rm=T)){
          predict.label[i]<-o
          dist.to.template[i]<-1-orig.dist.to.all.temp[o]
          dist.to.cls1[i]<-(1-orig.dist.to.all.temp[o])+o
        }
      }
    }
    
    for (p in 1:nresmpl){
      rnd.feature.matrix[,p]<-sample(expression[,i],num.features.extract,replace=F)
    }
    
    if (dist.selection=="cosine"){          # cosine
      for (res in 1:num.cls){
        eval(parse(text=paste("temp.resmpl<-temp.",res,sep="")))
        
        prod.sum<-apply(t(t(rnd.feature.matrix)*temp.resmpl),2,sum)
        
        data.sq.sum<-apply(rnd.feature.matrix^2,2,sum)
        temp.sq.sum<-sum(temp.resmpl^2)
        
        perm.dist.vector[(1+(nresmpl*(res-1))):(nresmpl*res)]<-
          (1-(prod.sum/(sqrt(data.sq.sum)*sqrt(temp.sq.sum))))
      }
    }
    
    if (dist.selection=="correlation"){          # correlation
      for (res in 1:num.cls){
        eval(parse(text=paste("temp.resmpl<-temp.",res,sep="")))
        perm.dist.vector[(1+(nresmpl*(res-1))):(nresmpl*res)]<-
          (1-as.vector(cor(rnd.feature.matrix,temp.resmpl,method="pearson",use="complete.obs")))
      }
    }
    
    # compute nominal p-value
    
    combined.stats.rank<-rank(c(dist.to.template[i],perm.dist.vector))
    nominal.p[i]<-combined.stats.rank[1]/length(combined.stats.rank)
    
    
  }
  BH.FDR<-nominal.p*num.samples/rank(nominal.p)
  Bonferroni.p<-nominal.p*num.samples
  
  BH.FDR[BH.FDR>1]<-1
  Bonferroni.p[Bonferroni.p>1]<-1
  
  dist.to.cls1.rank <- rank(dist.to.cls1)
  pred.summary <- cbind(sample.names,predict.label,dist.to.template,dist.to.cls1.rank,
                        nominal.p,BH.FDR,Bonferroni.p)
  
  return(pred.summary)
}
res.pre <- NTP_Fuc(specific.gene.expr,gene.label.matrix)
celid <- as.character(res.pre[,1])
res.label <- as.matrix(res.pre[,2])
row.names(res.label) <- sapply(celid,function(x){unlist(strsplit(x,"[.]"))[1]})

###drug discovery
setwd("D:/TCGA/TCGA_COAD_exp_HiSeqV2-2015-02-24")
ic50 <- read.table("IC50.txt",sep = "\t",header = T)
ic <- cbind(ic50,res.label[row.names(ic50),])
kw.pvalue <- apply(ic50,2,function(x){kruskal.test(x,ic$`res.label[row.names(ic50),]`)$p.value})




#############################Mutation Analysis##############################
setwd("D:/TCGA/TCGA_COAD_mutation_bcm_gene-2015-02-24")
#coad.mutation <- fread("genomicMatrix",sep = "\t")
#row.names(coad.mutation) <- coad.mutation$sample
load("COAD_mutation.Rdata")
MutualExclusivity <- function(Part.mutual.exclusity.mat,cc){
  E <- t(Part.mutual.exclusity.mat)
  g <- intersect(colnames(E),as.character(specific.gene.list[[1]]))
  #s <- intersect(colnames(E),as.character(specific.gene.list[[2]]))
  resu <- NULL
  for(i in 1:length(g)){
    for(j in i:length(g)){
      if(i!=j){
        f <- fisher.test(E[,i], E[,j]) 
        resu <- rbind(resu,cbind(geneA=i,geneB=j,oddsRatio=f$estimate,pvalue=f$p.value))
      }
    }
  }
  # some formatting
  res.new <- as.data.frame(resu)
  #res.new$geneA <- factor(res.new$geneA,levels=colnames(E))
  #res.new$geneB <- factor(res.new$geneB,levels=colnames(E))
  res.new$oddsRatio <- as.numeric(as.character(res.new$oddsRatio))
  res.new$pvalue <- as.numeric(as.character(res.new$pvalue))
  row.names(res.new) <- NULL
  # use p.adjust to correct for multi testing using a FDR
  res2 <- cbind(res.new,fdr=p.adjust(res.new$pvalue,"fdr"))
  # change the FDR in labels for plotting
  res2$stars <- cut(res2$fdr, 
                    breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), 
                    label=c("***", "**", "*", ""))
  res2$events <- cut(res2$oddsRatio,
                     breaks=c(0, 0.1, 0.5, 2, 10, Inf), 
                     label=c("SE", "TE", "NO", "TC","SC"),
                     include.lowest = T)
  # plot with ggplot 2
  require(ggplot2)
  require(cowplot) # not necessary but the plot is nicer
  p <- ggplot(res2, aes(geneA, geneB)) + 
    geom_tile(aes(fill = events),colour = "white") + 
    geom_text(aes(label=stars), color="black", size=5) + 
    xlab("Gene A") + 
    ylab("Gene B") +
    scale_fill_brewer(palette = "Set3") + 
    ggtitle(paste(cc,"Mutual exclusivity",sep = "_")) + 
    theme(legend.key.size = unit(0.4, "cm"),
          axis.text.y=element_text(size=9,face = "bold"),
          axis.text.x=element_text(angle=70,size=9,hjust=1,face = "bold",
                                   margin = margin(0,0,0,20)))+
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, 
                                    margin = margin(l=100,r=50,t=10,b=10),
                                    face = "bold", colour = "black",size=10))
  p
}
MutualExclusivity(COAD,"coad")
colnames(COAD) <- paste(colnames(COAD),".01",sep = "")
COAD.sub1 <- COAD[,na.omit(match(colnames(tumor)))]

p <- ggplot(res2, aes(geneA, geneB)) + 
  geom_tile(aes(fill = events),colour = "white") + 
  geom_text(aes(label=stars), color="black", size=5) + 
  xlab("")+ylab("")+
  scale_fill_brewer(palette = "Set3")+
  theme(legend.key.size = unit(0.4, "cm"),
        axis.text.y=element_text(size=7,face = "bold"),
        axis.text.x=element_text(angle=70,size=9,hjust=1,face = "bold",margin = margin(0,0,0,20)))
p


subcluster.list <- list(one=names(res[[4]]$consensusClass)[res[[4]]$consensusClass==1],
                        two=names(res[[4]]$consensusClass)[res[[4]]$consensusClass==2],
                        three=names(res[[4]]$consensusClass)[res[[4]]$consensusClass==3],
                        four=names(res[[4]]$consensusClass)[res[[4]]$consensusClass==4])
subcluster.expression.list <- sapply(subcluster.list,function(x){tumor[,match(x,colnames(tumor))]})


One.gseaset <- topTable(cluster.fit2,coef=1,n=nrow(cluster.fit2),lfc=0)




###########################Immume Infiltration###############################
setwd("D:/TCGA/TCGA_COAD_exp_HiSeqV2-2015-02-24")
load("COAD_Infil.Rdata")
row.names(COAD) <- chartr("-",".",row.names(COAD))
COAD.infil <- COAD[na.omit(match(paste(colnames(tumor),"A",sep = ""),row.names(COAD))),]

cluster.result.infil <- data.frame(cbind(paste(names(res[[4]]$consensusClass),"A",sep = "")
                                         [na.omit(match(row.names(COAD.infil),
                                                        paste(names(res[[4]]$consensusClass),"A",sep = "")))],
                                       res[[4]]$consensusClass[na.omit(match(row.names(COAD.infil),
                                                                             paste(names(res[[4]]$consensusClass),"A",sep = "")))]))
cluster.result.final.infil <- dplyr::arrange(cluster.result.infil,X2)
annotation.col.infil <- data.frame(cluster.imm=factor(cluster.result.final.infil$X2,labels = c("1","2","3","4")))
cluster.infil <- c("red","green","yellow","blue")
names(cluster.infil) <- c("1","2","3","4")


ann.colors.infil <- list(cluster.infil=cluster.infil)
rownames(annotation.col.infil) <- cluster.result.final.infil$X1

pheatmap(t(apply(t(COAD.infil[as.character(cluster.result.final.infil$X1),]),1,function(x){rescale(x,to=c(-3,3))})),
         annotation_col=annotation.col.infil,
         annotation_colors=ann.colors.infil,
         cluster_cols = FALSE,
         cluster_rows =FALSE,
         show_rownames = T, 
         show_colnames = F,
         color = redblue(100)[100:1])



######################################RPPA##################################
setwd("D:/TCGA/gdac.broadinstitute.org_COAD.RPPA_AnnotateWithGene.Level_3.2016012800.0.0/")
coad.rppa <- read.table("COAD.rppa.txt",header = T,sep = "\t",row.names = 1)
coad.rppa <- coad.rppa[is.na(coad.rppa)]
rppa.res <- ConsensusClusterPlus(as.matrix(coad.rppa),maxK = 8,pItem = 0.8,pFeature = 0.8,
                                 reps = 1000,distance = "euclidean",clusterAlg = "km")
calcICL(rppa.res,title="untitled_consensus_cluster",plot=NULL,writeTable=FALSE)


