setwd("D:/GEO/COAD/")
load("diffneed.RData")
library(cluster)
sil <- silhouette(coad.res[[5]]$consensusClass,dist(t(expr.combat)))
sample.filter <- coad.res[[5]]$consensusClass[sil[,3] > 0]
expr.filter <- expr.combat[,match(names(sample.filter),names(coad.res[[5]]$consensusClass))]

library(siggenes)
sam.out<-sam(expr.filter,sample.filter,rand=123,gene.names=row.names(expr.filter))
sum.sam.out <- summary(sam.out,30)
expr.filter.sam <- expr.filter[sum.sam.out@row.sig.genes,]

library(pamr)
mydata <- list(x=expr.filter.sam, y=sample.filter,geneid=rownames(expr.filter.sam),
               genenames=rownames(expr.filter.sam))
mytrain <-pamr.train(mydata)
mycv <- pamr.cv(mytrain,mydata, nfold=10)
#pamr.plotcv(mycv)
pamr.confusion(mycv, threshold=5)
pam.out <- pamr.listgenes(mytrain, mydata, threshold=5)
expr.filter.sam.pam <- expr.filter.sam[match(pam.out[,1],row.names(expr.filter.sam)),]


##########################################
######specific gene selection
library(limma)
clusterlabel <- coad.res[[5]]$consensusClass
expr.matrix.filter <- expr.combat[match(row.names(expr.filter.sam.pam),row.names(expr.combat)),]
################subtype one
one.design <- sub("1","one",clusterlabel)
one.design <- sub("2","other",one.design)
one.design <- sub("3","other",one.design)
one.design <- sub("4","other",one.design)
one.design <- sub("5","other",one.design)
one.mat <- factor(one.design, levels = c("one","other"))
one.matrix <- model.matrix(~0 + one.mat)
colnames(one.matrix) <- levels(one.mat)
one.fit <- lmFit(expr.matrix.filter,one.matrix)
cont.one <- makeContrasts(one-other,
                          levels = one.matrix)
one.fit2 <- contrasts.fit(one.fit, cont.one)
one.fit2 <- eBayes(one.fit2)
##decidetest
one.res <- decideTests(one.fit2,method = "global",adjust.method = "BH",
                       p.value = 0.01,lfc = log2(1))
summary(one.res)
##diff expr
one.diff <- topTable(one.fit2,coef=1,resort.by = "logFC",n=nrow(one.fit2),lfc=log2(1.1))
one.diff.final <- one.diff[one.diff[,"adj.P.Val"]<0.01,]
nrow(one.diff.final)
up.sample.one <- one.diff.final[one.diff.final[,"logFC"]>0,]
nrow(up.sample.one)


################subtype two
two.design <- sub("2","two",clusterlabel)
two.design <- sub("1","other",two.design)
two.design <- sub("3","other",two.design)
two.design <- sub("4","other",two.design)
two.design <- sub("5","other",two.design)
two.mat <- factor(two.design, levels = c("two","other"))
two.matrix <- model.matrix(~0 + two.mat)
colnames(two.matrix) <- levels(two.mat)
two.fit <- lmFit(expr.matrix.filter,two.matrix)
cont.two <- makeContrasts(two-other,
                          levels = two.matrix)
two.fit2 <- contrasts.fit(two.fit, cont.two)
two.fit2 <- eBayes(two.fit2)
##decidetest
two.res <- decideTests(two.fit2,method = "global",adjust.method = "BH",
                       p.value = 0.01,lfc = log2(1))
summary(two.res)
##diff expr
two.diff <- topTable(two.fit2,coef=1,resort.by = "logFC",n=nrow(two.fit2),lfc=log2(1.1))
two.diff.final <- two.diff[two.diff[,"adj.P.Val"]<0.01,]
nrow(two.diff.final)
up.sample.two <- two.diff.final[two.diff.final[,"logFC"]>0,]
nrow(up.sample.two)


################subtype three
three.design <- sub("3","three",clusterlabel)
three.design <- sub("1","other",three.design)
three.design <- sub("2","other",three.design)
three.design <- sub("4","other",three.design)
three.design <- sub("5","other",three.design)
three.mat <- factor(three.design, levels = c("three","other"))
three.matrix <- model.matrix(~0 + three.mat)
colnames(three.matrix) <- levels(three.mat)
three.fit <- lmFit(expr.matrix.filter,three.matrix)
cont.three <- makeContrasts(three-other,
                            levels = three.matrix)
three.fit2 <- contrasts.fit(three.fit, cont.three)
three.fit2 <- eBayes(three.fit2)
##decidetest
three.res <- decideTests(three.fit2,method = "global",adjust.method = "BH",
                         p.value = 0.01,lfc = log2(1))
summary(three.res)
##diff expr
three.diff <- topTable(three.fit2,coef=1,resort.by = "logFC",n=nrow(three.fit2),lfc=log2(1.1))
three.diff.final <- three.diff[three.diff[,"adj.P.Val"]<0.01,]
nrow(three.diff.final)
up.sample.three <- three.diff.final[three.diff.final[,"logFC"]>0,]
nrow(up.sample.three)


################subtype four
four.design <- sub("4","four",clusterlabel)
four.design <- sub("1","other",four.design)
four.design <- sub("2","other",four.design)
four.design <- sub("3","other",four.design)
four.design <- sub("5","other",four.design)
four.mat <- factor(four.design, levels = c("four","other"))
four.matrix <- model.matrix(~0 + four.mat)
colnames(four.matrix) <- levels(four.mat)
four.fit <- lmFit(expr.matrix.filter,four.matrix) 
cont.four <- makeContrasts(four-other,
                           levels = four.matrix)
four.fit2 <- contrasts.fit(four.fit, cont.four)
four.fit2 <- eBayes(four.fit2)
##decidetest
four.res <- decideTests(four.fit2,method = "global",adjust.method = "BH",
                        p.value = 0.01,lfc = log2(1))
summary(four.res)
##diff expr
four.diff <- topTable(four.fit2,coef=1,resort.by = "logFC",n=nrow(four.fit2),lfc=log2(1.1))
four.diff.final <- four.diff[four.diff[,"adj.P.Val"]<0.01,]
nrow(four.diff.final)
up.sample.four <- four.diff.final[four.diff.final[,"logFC"]>0,]
nrow(up.sample.four)


################subtype five
five.design <- sub("5","five",clusterlabel)
five.design <- sub("1","other",five.design)
five.design <- sub("2","other",five.design)
five.design <- sub("3","other",five.design)
five.design <- sub("4","other",five.design)
five.mat <- factor(five.design, levels = c("five","other"))
five.matrix <- model.matrix(~0 + five.mat)
colnames(five.matrix) <- levels(five.mat)
five.fit <- lmFit(expr.matrix.filter,five.matrix)
cont.four <- makeContrasts(five-other,
                           levels = five.matrix)
five.fit2 <- contrasts.fit(five.fit, cont.four)
five.fit2 <- eBayes(five.fit2)
##decidetest
five.res <- decideTests(five.fit2,method = "global",adjust.method = "BH",
                        p.value = 0.01,lfc = log2(1))
summary(five.res)
##diff expr
five.diff <- topTable(five.fit2,coef=1,resort.by = "logFC",n=nrow(five.fit2),lfc=log2(1.1))
five.diff.final <- five.diff[five.diff[,"adj.P.Val"]<0.01,]
nrow(five.diff.final)
up.sample.five <- five.diff.final[five.diff.final[,"logFC"]>0,]
nrow(up.sample.five)

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
gene.label <- c(rep("one",length(one.specific)),
                rep("two",length(two.specific)),
                rep("three",length(three.specific)),
                rep("four",length(four.specific)),
                rep("five",length(five.specific)))
gene.label.mat <- factor(gene.label, levels = c("one","two","three","four","five"))
gene.label.matrix <- model.matrix(~0 + gene.label.mat)
colnames(gene.label.matrix) <- levels(gene.label.mat)

source("D:/TCGA/Tumor Analysis/NTP.R")
test.expr <- expr.combat[match(specific.genes,row.names(expr.combat)),]
res.test <- NTP_Fuc(test.expr,gene.label.matrix)
(165+179+91+25+87)/length(coad.res[[5]]$consensusClass)


##################################
library(ROCR)
pred <- prediction(expr.filter.sam.pam[1,],sample.filter)

#roc(factor(sample.filter),expr.filter.sam.pam[1,],levels = c(1,2,3,4,5))
