setwd("D:/GEO/COAD/")
library(GEOquery)
load("GSE5851.RData")
phenotype.GSE5851 <- pData(GSE5851@phenoData)
expr.GSE5851 <- GSE5851@assayData$exprs
expr.GSE5851 <- log2(expr.GSE5851)
expr.GSE5851 <- expr.GSE5851[apply(expr.GSE5851,1,function(x){sum(is.na(x))<(length(x)*0.9)}),]
liver.metas <- expr.GSE5851[,grep("Liver$",phenotype.GSE5851$characteristics_ch1)]

load("GSE14333.RData")
phenotype.GSE14333 <- pData(GSE14333@phenoData)
expr.GSE14333 <- GSE14333@assayData$exprs
expr.GSE14333 <- expr.GSE14333[apply(expr.GSE14333,1,function(x){sum(is.na(x))<(length(x)*0.9)}),]
primary.expr <- expr.GSE14333[match(row.names(liver.metas),row.names(expr.GSE14333)),]

label <- c(rep("metastase",ncol(liver.metas)),rep("primary",ncol(primary.expr)))
names(label) <- c(colnames(liver.metas),colnames(primary.expr))
tumor.expr <- cbind(liver.metas,primary.expr)
tumor.zscore <- scale(tumor.expr)

library(pcaMethods)
pca.zscore <- pca(tumor.zscore,method = "bpca",nPcs = 2)
pca.coad <- cbind(pca.zscore@loadings,label)
library(ggplot2)
ggplot(data.frame(pca.coad), aes(PC1, PC2, shape=label, color=label)) +
  geom_point(size = 4) +
  labs(x=paste("PC1", pca.zscore@R2[1] * 100, "% of variance"),
       y=paste("PC2", pca.zscore@R2[2] * 100, "% of variance"),
       title="Before remove batch effect") +
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 20),
        axis.text.x= element_blank(),
        axis.text.y = element_blank())

library(sva)
expr.combat <- ComBat(pca.zscore@completeObs,label)
combat.pca <- pca(expr.combat,method = "bpca",nPcs = 2)
pca.coad.combat <- cbind(combat.pca@loadings,as.factor(label))
COLORS=c("darkred","orchid1","darkorange4","darkorange","darkolivegreen4",
         "darkolivegreen1","darkmagenta","darkkhaki","darkgreen","darkgrey",
         "darkgoldenrod1","cyan4","cyan","cornflowerblue","chartreuse",
         "burlywood4","brown1","blueviolet","rosybrown4","aquamarine4",
         "darkslateblue","darksalmon","forestgreen","firebrick3","blue2",
         "dodgerblue",'deeppink',"gray65","yellow","springgreen","mediumseagreen")
ggplot(data.frame(pca.coad.combat), aes(PC1, PC2, shape=label, color=label)) +
  geom_point(size = 4) +
  labs(x=paste("PC1", combat.pca@R2[1] * 100, "% of variance"),
       y=paste("PC2", combat.pca@R2[2] * 100, "% of variance"),
       title="After remove batch effect") +
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 20),
        axis.text.x= element_blank(),
        axis.text.y = element_blank())

combat.pca1 <- pca(expr.combat,method = "bpca",nPcs = 3)
pca.coad.combat1 <- cbind(combat.pca1@loadings,as.factor(label))
library(threejs)
scatterplot3js(pca.coad.combat1[,1:3],
               color = COLORS[factor(pca.coad.combat[,"label"])],
               size = 0.3,pch=10,renderer = "canvas")



library(siggenes)
sam.out <- sam(tumor,label,rand=123,gene.names = row.names(tumor))
sam.out.sum <- summary(sam.out,1)
tumor.sam <- tumor[sam.out.sum@row.sig.genes,]

######
tumor.sam.pvalue <- tumor[sam.out@p.value<0.01,]
######

library(pamr)
mydata <- list(x=tumor.sam, y=label,geneid=rownames(tumor.sam),
               genenames=rownames(tumor.sam))
mytrain <-pamr.train(mydata)
mycv <- pamr.cv(mytrain,mydata, nfold=10)
#pamr.plotcv(mycv)
pamr.confusion(mycv, threshold=1.1)
pam.out <- pamr.listgenes(mytrain, mydata, threshold=1.1)
tumor.sam.pam <- tumor.sam[match(pam.out[,1],row.names(tumor.sam)),]

######
mydata <- list(x=tumor.sam.pvalue, y=label,geneid=rownames(tumor.sam.pvalue),
               genenames=rownames(tumor.sam.pvalue))
mytrain <-pamr.train(mydata)
mycv <- pamr.cv(mytrain,mydata, nfold=10)
#pamr.plotcv(mycv)
pamr.confusion(mycv, threshold=1.5)
pam.out.pvalue <- pamr.listgenes(mytrain, mydata, threshold=1.5)
tumor.pvalue.pam <- tumor.sam.pvalue[match(pam.out.pvalue[,1.5],row.names(tumor.sam.pvalue)),]
######

library(pROC)
roc.auc <- apply(tumor.sam.pam,1,function(x){y <- roc(label,x);y$auc})
tumor.roc <- tumor.sam.pam[roc.auc>0.7,]

######
roc.auc.pvalue <- apply(tumor.pvalue.pam,1,function(x){y <- roc(label,x);y$auc})
tumor.roc.pvalue <- tumor.pvalue.pam[roc.auc.pvalue>0.7,]
######

library(hgu133plus2.db)
probeset <- rownames(tumor.roc)
Symbol <- as.character(as.list(hgu133plus2SYMBOL[probeset]))
expr.symbol <- cbind.data.frame(Symbol,tumor.roc)
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

######
probeset.pvalue <- rownames(tumor.roc.pvalue)
Symbol.pvalue <- as.character(as.list(hgu133plus2SYMBOL[probeset.pvalue]))
expr.symbol.pvalue <- cbind.data.frame(Symbol.pvalue,tumor.roc.pvalue)
exprSet.pvalue <- rmDupID(expr.symbol.pvalue)
length(intersect(row.names(exprSet),row.names(exprSet.pvalue)))
######

library(randomForest)
tumor.rf <- cbind(t(exprSet),samp.cli[,"First.site.of.recurrence_liver.0......1....",drop=FALSE])
colnames(tumor.rf) <- c(row.names(exprSet),"metastatic")
tumor.rf$metastatic[tumor.rf$metastatic==1]="Trans"
tumor.rf$metastatic[tumor.rf$metastatic==0]="Primary"
tumor.rf$metastatic <- factor(tumor.rf$metastatic)
n<-length(names(tumor.rf))
Error<-NULL
set.seed(123)
for(i in 1:(n-1)){
  rf<-randomForest(metastatic ~ .,data = tumor.rf,mtry=i,importance=TRUE,proximity=TRUE)
  err<-mean(rf$err.rate)
  Error[i]<-err
}
m=which.min(Error)
rf.res <- randomForest(metastatic ~ .,data = tumor.rf,mtry=m,importance=TRUE,proximity=TRUE)
plot(rf.res)
rf.res <- randomForest(metastatic ~ .,data = tumor.rf,mtry=m,ntree=100,importance=TRUE,proximity=TRUE)

##SMOTE
library(unbalanced)
n <- ncol(tumor.rf)
output<-tumor.rf$metastatic
input<-tumor.rf[ ,-n]
smote.res <- ubBalance(X=input,Y=output,positive = "Trans",type="ubSMOTE", percOver=300, percUnder=150, verbose=TRUE)
balanced.data <- cbind(smote.res$X,smote.res$Y)
colnames(balanced.data) <- c(row.names(exprSet),"metastatic")
n1<-length(names(balanced.data))
Error<-NULL
set.seed(123)
for(i in 1:(n1-1)){
  rf<-randomForest(metastatic ~ .,data = balanced.data,mtry=i,importance=TRUE,proximity=TRUE)
  err<-mean(rf$err.rate)
  Error[i]<-err
}
m1=which.min(Error)
rf.res1 <- randomForest(metastatic ~ .,data = balanced.data,mtry=m1,importance=TRUE,proximity=TRUE)
plot(rf.res1)
rf.res1 <- randomForest(metastatic ~ .,data = balanced.data,mtry=m1,ntree=250,importance=TRUE,proximity=TRUE)

##remove batch effect
TCGA.expr <- TCGA.data[na.omit(match(row.names(exprSet),row.names(TCGA.data))),]
expr.data <- exprSet[na.omit(match(row.names(TCGA.expr),row.names(exprSet))),]
data.expr <- cbind(TCGA.expr,expr.data)
data.label <- c(rep("TCGA",ncol(TCGA.expr)),rep("GSE62254",ncol(expr.data)))
library(pcaMethods)
library(ggplot2)
library(sva)
pca.zscore <- pca(data.expr,method = "bpca",nPcs = 2)
pca.res <- cbind(pca.zscore@loadings,data.label)
ggplot(data.frame(pca.res), aes(PC1, PC2, shape=data.label, color=data.label)) +
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

TCGA.rf <- rbind(combat.pca@completeObs[,1:ncol(TCGA.expr)],TCGA.label)
row.names(TCGA.rf) <- c(row.names(TCGA.expr),"metastatic")
TCGA.rf <- data.frame(t(TCGA.rf))
TCGA.rf$metastatic[TCGA.rf$metastatic==1]="Trans"
TCGA.rf$metastatic[TCGA.rf$metastatic==0]="Primary"
# GSM.label <- label
# GSM.label[GSM.label==1]="Trans"
# GSM.label[GSM.label==0]="Primary"
GSE62254.rf <- rbind(combat.pca@completeObs[,(ncol(TCGA.expr)+1):ncol(data.expr)],label)
row.names(GSE62254.rf) <- c(row.names(TCGA.expr),"metastatic")
GSE62254.rf <- data.frame(t(GSE62254.rf))
GSE62254.rf$metastatic[GSE62254.rf$metastatic==1]="Trans"
GSE62254.rf$metastatic[GSE62254.rf$metastatic==0]="Primary"
data.rf <- rbind(TCGA.rf,GSE62254.rf)
index <- sample(2,nrow(data.rf),replace = TRUE,prob=c(0.6,0.4))
train.rf <- data.rf[index==1,]
test.rf <- data.rf[index==2,]

nbal <- ncol(train.rf)
output.combat<-as.factor(train.rf$metastatic)
input.combat<-train.rf[ ,-nbal]
smote.combat <- ubBalance(X=input.combat,Y=output.combat,positive = "Trans",type="ubSMOTE", percOver=360, percUnder=55, verbose=TRUE)
balanced.combat <- cbind(smote.combat$X,smote.combat$Y)
colnames(balanced.combat) <- c(row.names(expr.data),"metastatic")
# n2<-length(names(balanced.combat))
# Error<-NULL
# set.seed(2017)
# for(i in 1:(n2-1)){
#   rf<-randomForest(metastatic ~ .,data = balanced.combat,mtry=i,importance=TRUE,proximity=TRUE)
#   err<-mean(rf$err.rate)
#   Error[i]<-err
# }
# m2=which.min(Error)
rf.res.combat <- randomForest(metastatic ~ .,data = balanced.combat,importance=TRUE,proximity=TRUE)
plot(rf.res.combat)
rf.res1.combat <- randomForest(metastatic ~ .,data = balanced.combat,ntree=400,importance=TRUE,proximity=TRUE)
rf.pred <- predict(rf.res1.combat,test.rf)
table(rf.pred,test.rf$metastatic)

######################
###NTP
##differential expression
library(limma)
design.lab <- sub(1,"Trans",label)
design.lab <- sub(0,"Primary",design.lab)
design.mat <- factor(design.lab,levels = c("Trans","Primary"))
design.matrix <- model.matrix(~0 + design.mat)
colnames(design.matrix) <- levels(design.mat)
design.fit <- lmFit(tumor,design.matrix)
cont.cluster <- makeContrasts(Trans-Primary,levels = design.matrix)
design.fit2 <- contrasts.fit(design.fit, cont.cluster)
design.fit2 <- eBayes(design.fit2)
diff.res <- decideTests(design.fit2,method = "global",adjust.method = "BH",
                        p.value = 0.05,lfc = log2(1))
summary(diff.res)
tumor.diff <- topTable(design.fit2,coef=1,n=nrow(design.fit2),lfc=log2(1))
tumor.diff.final <- tumor.diff[tumor.diff[,"adj.P.Val"]<0.05,]
nrow(tumor.diff.final)
up.sample.tumor <- tumor.diff.final[tumor.diff.final[,"logFC"]>0,]
nrow(up.sample.tumor)
#NTP
gene.label <- c(rep("Trans",nrow(up.sample.tumor)),
                rep("Primary",nrow(tumor.diff.final)-nrow(up.sample.tumor)))
gene.label.mat <- factor(gene.label, levels = c("Trans","Primary"))
gene.label.matrix <- model.matrix(~0 + gene.label.mat)
colnames(gene.label.matrix) <- levels(gene.label.mat)
row.names(gene.label.matrix) <- row.names(tumor.diff.final)
specific.gene.expr <- tumor[match(row.names(gene.label.matrix),row.names(tumor)),]
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

##################################
##balanced
smote.res1 <- ubBalance(X=t(tumor),Y=design.mat,positive = "Trans",type="ubSMOTE", percOver=500, percUnder=50, verbose=TRUE)
balanced.data1 <- cbind(smote.res1$X,smote.res1$Y)
design.lab1 <- sub(1,"Trans",balanced.data$metastatic)
design.lab1 <- sub(0,"Primary",design.lab1)
design.mat1 <- factor(design.lab1,levels = c("Trans","Primary"))
design.matrix1 <- model.matrix(~0 + design.mat1)
colnames(design.matrix1) <- levels(design.mat1)
design.fit1 <- lmFit(t(balanced.data[,-ncol(balanced.data)]),design.matrix1)
cont.cluster1 <- makeContrasts(Trans-Primary,levels = design.matrix1)
design.fit2.1 <- contrasts.fit(design.fit1, cont.cluster1)
design.fit2.1 <- eBayes(design.fit2.1)
diff.res1 <- decideTests(design.fit2.1,method = "global",adjust.method = "BH",
                         p.value = 0.05,lfc = log2(1))
summary(diff.res1)
tumor.diff1 <- topTable(design.fit2.1,coef=1,n=nrow(design.fit2.1),lfc=log2(1))
tumor.diff.final1 <- tumor.diff1[tumor.diff1[,"adj.P.Val"]<0.05,]
nrow(tumor.diff.final1)
up.sample.tumor1 <- tumor.diff.final1[tumor.diff.final1[,"logFC"]>0,]
nrow(up.sample.tumor1)
#NTP
gene.label <- c(rep("Trans",nrow(up.sample.tumor)),
                rep("Primary",nrow(tumor.diff.final)-nrow(up.sample.tumor)))
gene.label.mat <- factor(gene.label, levels = c("Trans","Primary"))
gene.label.matrix <- model.matrix(~0 + gene.label.mat)
colnames(gene.label.matrix) <- levels(gene.label.mat)
row.names(gene.label.matrix) <- row.names(tumor.diff.final)
specific.gene.expr <- tumor[match(row.names(gene.label.matrix),row.names(tumor)),]
res.pre <- NTP_Fuc(specific.gene.expr,gene.label.matrix)