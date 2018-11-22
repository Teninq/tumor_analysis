setwd("D:/GEO/COAD/")
load("diffneed.RData")
library(limma)
library(pheatmap)
library(gplots)
library(scales)
##########################specific genes
clusterlabel <- coad.res[[5]]$consensusClass

################subtype one
one.design <- sub("1","one",clusterlabel)
one.design <- sub("2","other",one.design)
one.design <- sub("3","other",one.design)
one.design <- sub("4","other",one.design)
one.design <- sub("5","other",one.design)
one.mat <- factor(one.design, levels = c("one","other"))
one.matrix <- model.matrix(~0 + one.mat)
colnames(one.matrix) <- levels(one.mat)
one.fit <- lmFit(expr.combat,one.matrix)
cont.one <- makeContrasts(one-other,
                          levels = one.matrix)
one.fit2 <- contrasts.fit(one.fit, cont.one)
one.fit2 <- eBayes(one.fit2)
##decidetest
one.res <- decideTests(one.fit2,method = "global",adjust.method = "BH",
                       p.value = 0.01,lfc = log2(1.5))
summary(one.res)
##diff expr
one.diff <- topTable(one.fit2,coef=1,resort.by = "logFC",n=nrow(one.fit2),lfc=log2(1.15))
one.diff.final <- one.diff[one.diff[,"adj.P.Val"]<0.01,]
nrow(one.diff.final)
up.sample.one <- one.diff.final[one.diff.final[,"logFC"]>0,]
nrow(up.sample.one)
# one.up <- as.matrix(up.sample.one[,1])
# row.names(one.up) <- row.names(up.sample.one)
# write.table(one.up,file = "one_up.rnk",quote = F,sep = "\t",row.names = TRUE,col.names = FALSE)
# one.gene <- as.matrix(one.diff.final[,1])
# row.names(one.gene) <- row.names(one.diff.final)
# write.table(one.gene,file = "one_gene.rnk",quote = F,sep = "\t",row.names = TRUE,col.names = FALSE)

################subtype two
two.design <- sub("2","two",clusterlabel)
two.design <- sub("1","other",two.design)
two.design <- sub("3","other",two.design)
two.design <- sub("4","other",two.design)
two.design <- sub("5","other",two.design)
two.mat <- factor(two.design, levels = c("two","other"))
two.matrix <- model.matrix(~0 + two.mat)
colnames(two.matrix) <- levels(two.mat)
two.fit <- lmFit(expr.combat,two.matrix)
cont.two <- makeContrasts(two-other,
                          levels = two.matrix)
two.fit2 <- contrasts.fit(two.fit, cont.two)
two.fit2 <- eBayes(two.fit2)
##decidetest
two.res <- decideTests(two.fit2,method = "global",adjust.method = "BH",
                       p.value = 0.01,lfc = log2(1))
summary(two.res)
##diff expr
two.diff <- topTable(two.fit2,coef=1,resort.by = "logFC",n=nrow(two.fit2),lfc=log2(1.15))
two.diff.final <- two.diff[two.diff[,"adj.P.Val"]<0.01,]
nrow(two.diff.final)
up.sample.two <- two.diff.final[two.diff.final[,"logFC"]>0,]
nrow(up.sample.two)
# two.up <- as.matrix(up.sample.two[,1])
# row.names(two.up) <- row.names(up.sample.two)
# write.table(two.up,file = "two_up.rnk",quote = F,sep = "\t",row.names = TRUE,col.names = FALSE)
# two.gene <- as.matrix(two.diff.final[,1])
# row.names(two.gene) <- row.names(two.diff.final)
# write.table(two.gene,file = "two_gene.rnk",quote = F,sep = "\t",row.names = TRUE,col.names = FALSE)


################subtype three
three.design <- sub("3","three",clusterlabel)
three.design <- sub("1","other",three.design)
three.design <- sub("2","other",three.design)
three.design <- sub("4","other",three.design)
three.design <- sub("5","other",three.design)
three.mat <- factor(three.design, levels = c("three","other"))
three.matrix <- model.matrix(~0 + three.mat)
colnames(three.matrix) <- levels(three.mat)
three.fit <- lmFit(expr.combat,three.matrix)
cont.three <- makeContrasts(three-other,
                            levels = three.matrix)
three.fit2 <- contrasts.fit(three.fit, cont.three)
three.fit2 <- eBayes(three.fit2)
##decidetest
three.res <- decideTests(three.fit2,method = "global",adjust.method = "BH",
                         p.value = 0.01,lfc = log2(1))
summary(three.res)
##diff expr
three.diff <- topTable(three.fit2,coef=1,resort.by = "logFC",n=nrow(three.fit2),lfc=log2(1.15))
three.diff.final <- three.diff[three.diff[,"adj.P.Val"]<0.01,]
nrow(three.diff.final)
up.sample.three <- three.diff.final[three.diff.final[,"logFC"]>0,]
nrow(up.sample.three)
# three.up <- as.matrix(up.sample.three[,1])
# row.names(three.up) <- row.names(up.sample.three)
# write.table(three.up,file = "three_up.rnk",quote = F,sep = "\t",row.names = TRUE,col.names = FALSE)
# three.gene <- as.matrix(three.diff.final[,1])
# row.names(three.gene) <- row.names(three.diff.final)
# write.table(three.gene,file = "three_gene.rnk",quote = F,sep = "\t",row.names = TRUE,col.names = FALSE)


################subtype four
four.design <- sub("4","four",clusterlabel)
four.design <- sub("1","other",four.design)
four.design <- sub("2","other",four.design)
four.design <- sub("3","other",four.design)
four.design <- sub("5","other",four.design)
four.mat <- factor(four.design, levels = c("four","other"))
four.matrix <- model.matrix(~0 + four.mat)
colnames(four.matrix) <- levels(four.mat)
four.fit <- lmFit(expr.combat,four.matrix)
cont.four <- makeContrasts(four-other,
                           levels = four.matrix)
four.fit2 <- contrasts.fit(four.fit, cont.four)
four.fit2 <- eBayes(four.fit2)
##decidetest
four.res <- decideTests(four.fit2,method = "global",adjust.method = "BH",
                        p.value = 0.01,lfc = log2(1))
summary(four.res)
##diff expr
four.diff <- topTable(four.fit2,coef=1,resort.by = "logFC",n=nrow(four.fit2),lfc=log2(1.15))
four.diff.final <- four.diff[four.diff[,"adj.P.Val"]<0.01,]
nrow(four.diff.final)
up.sample.four <- four.diff.final[four.diff.final[,"logFC"]>0,]
nrow(up.sample.four)
# four.up <- as.matrix(up.sample.four[,1])
# row.names(four.up) <- row.names(up.sample.four)
# write.table(four.up,file = "four_up.rnk",quote = F,sep = "\t",row.names = TRUE,col.names = FALSE)
# four.gene <- as.matrix(four.diff.final[,1])
# row.names(four.gene) <- row.names(four.diff.final)
# write.table(four.gene,file = "four_gene.rnk",quote = F,sep = "\t",row.names = TRUE,col.names = FALSE)


################subtype five
five.design <- sub("5","five",clusterlabel)
five.design <- sub("1","other",five.design)
five.design <- sub("2","other",five.design)
five.design <- sub("3","other",five.design)
five.design <- sub("4","other",five.design)
five.mat <- factor(five.design, levels = c("five","other"))
five.matrix <- model.matrix(~0 + five.mat)
colnames(five.matrix) <- levels(five.mat)
five.fit <- lmFit(expr.combat,five.matrix)
cont.four <- makeContrasts(five-other,
                           levels = five.matrix)
five.fit2 <- contrasts.fit(five.fit, cont.four)
five.fit2 <- eBayes(five.fit2)
##decidetest
five.res <- decideTests(five.fit2,method = "global",adjust.method = "BH",
                        p.value = 0.01,lfc = log2(1))
summary(five.res)
##diff expr
five.diff <- topTable(five.fit2,coef=1,resort.by = "logFC",n=nrow(five.fit2),lfc=log2(1.15))
five.diff.final <- five.diff[five.diff[,"adj.P.Val"]<0.01,]
nrow(five.diff.final)
up.sample.five <- five.diff.final[five.diff.final[,"logFC"]>0,]
nrow(up.sample.five)
# five.up <- as.matrix(up.sample.five[,1])
# row.names(five.up) <- row.names(up.sample.five)
# write.table(five.up,file = "five_up.rnk",quote = F,sep = "\t",row.names = TRUE,col.names = FALSE)
# five.gene <- as.matrix(five.diff.final[,1])
# row.names(five.gene) <- row.names(five.diff.final)
# write.table(five.gene,file = "five_gene.rnk",quote = F,sep = "\t",row.names = TRUE,col.names = FALSE)

###specific genes
up.genes <- list(one=row.names(up.sample.one),
                 two=row.names(up.sample.two),
                 three=row.names(up.sample.three),
                 four=row.names(up.sample.four),
                 five=row.names(up.sample.five))
gplots::venn(up.genes)

one.specific <- setdiff(setdiff(setdiff(setdiff(up.genes[[1]],up.genes[[2]]),up.genes[[3]]),up.genes[[4]]),up.genes[[5]])
two.specific <- setdiff(setdiff(setdiff(setdiff(up.genes[[2]],up.genes[[1]]),up.genes[[3]]),up.genes[[4]]),up.genes[[5]])
three.specific <- setdiff(setdiff(setdiff(setdiff(up.genes[[3]],up.genes[[1]]),up.genes[[2]]),up.genes[[4]]),up.genes[[5]])
four.specific <- setdiff(setdiff(setdiff(setdiff(up.genes[[4]],up.genes[[2]]),up.genes[[3]]),up.genes[[1]]),up.genes[[5]])
five.specific <- setdiff(setdiff(setdiff(setdiff(up.genes[[5]],up.genes[[2]]),up.genes[[3]]),up.genes[[4]]),up.genes[[1]])

specific.genes <- c(one.specific,two.specific,three.specific,four.specific,five.specific)
tumor.all.specific <- expr.combat[match(specific.genes,row.names(expr.combat)),]
cluster.result <- data.frame(cbind(names(coad.res[[5]]$consensusClass),
                                   coad.res[[5]]$consensusClass))
cluster.result.final <- dplyr::arrange(cluster.result,X2)
annotation.col <- data.frame(cluster=factor(cluster.result.final$X2,labels = c("1","2","3","4","5")))
cluster <- c("red","green","yellow","blue","pink")
names(cluster) <- c("1","2","3","4","5")
ann.colors <- list(cluster=cluster)
rownames(annotation.col) <- cluster.result.final$X1

###plot heatmap
pheatmap(t(apply(tumor.all.specific[,as.character(cluster.result.final$X1)],1,function(x){rescale(x,to=c(-3,3))})),
         annotation_col=annotation.col,
         annotation_colors=ann.colors,
         cluster_cols = FALSE,
         cluster_rows =FALSE,
         show_rownames = F, 
         show_colnames = F,
         color = redblue(100)[100:1])

ntp.specific.genes <- c(one.specific[1:200],two.specific[1:200],three.specific[1:200],four.specific[1:200],five.specific[1:200])
test.expr <- expr.combat[match(ntp.specific.genes,row.names(expr.combat)),]

library(hgu133plus2.db)
probeset <- rownames(test.expr)
Symbol <- as.character(as.list(hgu133plus2SYMBOL[probeset]))
expr.symbol <- cbind.data.frame(Symbol,test.expr)
rmDupID <-function(a){
  exprSet <- a[,-1]
  rowMeans <- apply(exprSet,1,function(x){mean(as.numeric(x),na.rm=T)})
  gene.symbol <- row.names(a[!duplicated(a[,1]),])
  b <- a[order(rowMeans,decreasing=T),]
  exprSet <- b[match(gene.symbol,row.names(b)),]
  exprSet <- exprSet[-grep("^NA$",exprSet[,1]),]
  #exprSet=apply(exprSet,2,as.numeric)
  rownames(exprSet) <- exprSet[,1]
  exprSet <- exprSet[,-1]
  return(exprSet)
}
exprSet <- rmDupID(expr.symbol)
up.symbol <- list(one=unique(as.character(as.list(hgu133plus2SYMBOL[one.specific[1:200]]))),
                  two=unique(as.character(as.list(hgu133plus2SYMBOL[two.specific[1:200]]))),
                  three=unique(as.character(as.list(hgu133plus2SYMBOL[three.specific[1:200]]))),
                  four=unique(as.character(as.list(hgu133plus2SYMBOL[four.specific[1:200]]))),
                  five=unique(as.character(as.list(hgu133plus2SYMBOL[five.specific[1:200]]))))
gplots::venn(up.symbol)
ntp.label <- c(rep("one",163),
               rep("two",139),
               rep("three",151),
               rep("four",165),
               rep("five",163))
ntp.label.mat <- factor(ntp.label, levels = c("one","two","three","four","five"))
ntp.label.matrix <- model.matrix(~0 + ntp.label.mat)
colnames(ntp.label.matrix) <- levels(ntp.label.mat)
row.names(ntp.label.matrix) <- row.names(exprSet)

source("D:/TCGA/Tumor Analysis/NTP.R")
res.test <- NTP_Fuc(exprSet,ntp.label.matrix)
t <- table(res.test[,2],coad.res[[5]]$consensusClass)
(t[1,1]+t[2,2]+t[3,3]+t[4,4]+t[5,5])/length(coad.res[[5]]$consensusClass)
save(ntp.label.matrix,file = "NTPneed.RData")
