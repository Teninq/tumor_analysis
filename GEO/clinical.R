setwd("D:/GEO/COAD/")
library(GEOquery)
load("diffneed.RData")

#age,time,gender,ethnicity,dss,dfs,overall,
load("GSE17538.RData")
library(survival)
library(ggplot2)
library(ggpubr)
library(GGally)
phenotype.GSE17538 <- pData(GSE17538$`GSE17538-GPL570_series_matrix.txt.gz`@phenoData)
GSE17538.classlabel <- coad.res[[5]]$consensusClass[row.names(phenotype.GSE17538)]
Subtype <- sub("1","S-I",GSE17538.classlabel)
Subtype <- sub("2","S-II",Subtype)
Subtype <- sub("3","S-III",Subtype)
Subtype <- sub("4","S-IV",Subtype)
Subtype <- sub("5","S-V",Subtype)
names(Subtype) <- names(GSE17538.classlabel)
###Gender
gender.subtype <- data.frame(table(phenotype.GSE17538$characteristics_ch1.1,Subtype))
total.gender <- rep(as.numeric(table(Subtype)),each=2)
gender.subtype$propotion <- gender.subtype$Freq/total.gender
ggplot(gender.subtype, aes(x=Subtype, y=propotion, fill=Var1)) +
  geom_bar(stat="identity")+
  theme_pubr()+
  labs(x="Subtype",y="Gender")+
  theme(axis.text.x = element_text(angle=0,size=8,face = "bold",hjust = 1))

male.subtype <- gender.subtype[gender.subtype$Var1=="gender: male",]
female.subtype <- gender.subtype[gender.subtype$Var1=="gender: female",]
ggplot(male.subtype, aes(Subtype, propotion,fill=Subtype)) + 
  geom_bar(stat="identity", position="dodge")+
  labs(x="Subtype",y="Frequency",title="Male Propotion")+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 20),
        axis.text.x=element_text(face="bold",size=12,angle = 45),
        axis.text.y=element_text(face="bold",size=12),
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15),
        legend.text = element_text(size = 12, hjust = 3, vjust = 3))

###Age
age <- as.numeric(unlist(lapply(phenotype.GSE17538$characteristics_ch1,function(x){a=strsplit(as.character(x),":")[[1]];a[2]})))
age.subtype <- data.frame(cbind(age,Subtype))
row.names(age.subtype) <- names(Subtype)
ggplot(age.subtype, aes(x=factor(Subtype),y=as.numeric(age),fill=Subtype))+ 
  geom_boxplot()+
  labs(x="Subtype",y="Ages",title="Boxplot")+
  theme(plot.title = element_text(face="bold",hjust = 0.5,size = 20),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12,angle=45),
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15),
        plot.margin = unit(c(2, 2, 2, 2), "lines"))

###Survival Analysis
event <- unlist(lapply(phenotype.GSE17538$characteristics_ch1.5,function(x){a=strsplit(as.character(x),":")[[1]];a[2]}))
event.dds <- unlist(lapply(phenotype.GSE17538$characteristics_ch1.6,function(x){a=strsplit(as.character(x),":")[[1]];a[2]}))
overall.event <- sub(" no death",0,event)
overall.event <- sub(" no recurrence",0,overall.event)
overall.event <- sub(" recurrence",1,overall.event)
overall.event <- as.numeric(sub(" death",1,overall.event))
overall.time <- as.numeric(unlist(lapply(phenotype.GSE17538$characteristics_ch1.8,function(x){a=strsplit(as.character(x),":")[[1]];a[2]})))
overall.subtype <- survfit(Surv(overall.time, overall.event) ~ Subtype) 
ggsurv(overall.subtype,size.est = 1.2,
       surv.col=c("orchid1","brown1","blueviolet","yellow","aquamarine4"),
       cens.size = 3,
       cens.shape = 3,
       order.legend = F) +
  labs(x="Time (days)",y="Survival",title="Overall Survival")+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 15),
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))

#################dds
event.dds <- unlist(lapply(phenotype.GSE17538$characteristics_ch1.6,function(x){a=strsplit(as.character(x),":")[[1]];a[2]}))
dds.event <- sub(" no death",0,event.dds)
dds.event <- sub(" no recurrence",0,dds.event)
dds.event <- sub(" recurrence",1,dds.event)
dds.event <- as.numeric(sub(" death",1,dds.event))
dds.time <- as.numeric(unlist(lapply(phenotype.GSE17538$characteristics_ch1.9,function(x){a=strsplit(as.character(x),":")[[1]];a[2]})))
dds.subtype <- survfit(Surv(dds.time, dds.event) ~ Subtype) 
ggsurv(dds.subtype,size.est = 1.2,
       surv.col=c("orchid1","brown1","blueviolet","yellow","aquamarine4"),
       cens.size = 3,
       cens.shape = 3,
       order.legend = F) +
  labs(x="Time (days)",y="Survival",title="Disease Specific Survival")+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 15),
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))

###################dfs
event.dfs <- unlist(lapply(phenotype.GSE17538$characteristics_ch1.7,function(x){a=strsplit(as.character(x),":")[[1]];a[2]}))
dfs.event <- sub(" no death",0,event.dfs)
dfs.event <- sub(" no recurrence",0,dfs.event)
dfs.event <- sub(" recurrence",1,dfs.event)
dfs.event <- as.numeric(sub(" death",1,dfs.event))
dfs.time <- as.numeric(unlist(lapply(phenotype.GSE17538$characteristics_ch1.10,function(x){a=strsplit(as.character(x),":")[[1]];a[2]})))
dds.subtype <- survfit(Surv(dfs.time, dfs.event) ~ Subtype) 
ggsurv(dds.subtype,size.est = 1.2,
       surv.col=c("orchid1","brown1","blueviolet","yellow","aquamarine4"),
       cens.size = 3,
       cens.shape = 3,
       order.legend = F) +
  labs(x="Time (days)",y="Survival",title="Disease Free Survival")+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 15),
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))


#Microsatellite
load("GSE13294.RData")
phenotype.GSE13294 <- pData(GSE13294$GSE13294_series_matrix.txt.gz@phenoData)
GSE13294.classlabel <-  coad.res[[5]]$consensusClass[row.names(phenotype.GSE13294)]
#msi.fisher <- fisher.test(phenotype.GSE13294$characteristics_ch1,GSE13294.classlabel)
Subtype1 <- sub("1","S-I",GSE13294.classlabel)
Subtype1 <- sub("2","S-II",Subtype1)
Subtype1 <- sub("3","S-III",Subtype1)
Subtype1 <- sub("4","S-IV",Subtype1)
Subtype1 <- sub("5","S-V",Subtype1)
names(Subtype1) <- names(GSE13294.classlabel)
msi.mat <- table(phenotype.GSE13294$characteristics_ch1,Subtype1)
msi.pval <- pairwise.prop.test(msi.mat[1,],apply(msi.mat,2,sum))$p.value
msi.subtype <- data.frame(table(phenotype.GSE13294$characteristics_ch1,Subtype1))
total.msi <- as.numeric(rep(table(Subtype1),each=2))
msi.subtype$propotion <- msi.subtype$Freq/total.msi
ggplot(msi.subtype, aes(x=Subtype1, y=propotion, fill=Var1)) +
  geom_bar(stat="identity")+
  theme_pubr()+
  labs(x="Subtype",y="MSI/MSS")+
  theme(axis.text.x = element_text(angle=0,size=8,face = "bold",hjust = 1))

MSI.subtype <- msi.subtype[msi.subtype$Var1=="MSI",]
MSS.subtype <- msi.subtype[msi.subtype$Var1=="MSS",]
my_comparisons <- list(c("S-I","S-II"),c("S-I","S-III"),c("S-II","S-III"),
                       c("S-III","S-IV"),c("S-III","S-V"))
ggplot(msi.subtype[msi.subtype$Var1=="MSI",], aes(Subtype1, propotion,fill=Subtype1)) + 
  geom_bar(stat="identity", position="dodge")+
  geom_signif(comparisons=my_comparisons, annotations=c("3.97e-3","6.36e-11","2.00e-3","6.81e-6","9.44e-7"),
              y_position = c(0.65,1.1,1.2,1.3,1.4), tip_length = 0, vjust=0.4) +
  labs(x="Subtype",y="Frequency",title="MSI Propotion")+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 20),
        axis.text.x=element_text(face="bold",size=12,angle = 45),
        axis.text.y=element_text(face="bold",size=12),
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15),
        legend.text = element_text(size = 12, hjust = 3, vjust = 3))

#location,age,gender,DFS
load("GSE14333.RData")
phenotype.GSE14333 <- pData(GSE14333@phenoData)
GSE14333.classlabel <- coad.res[[5]]$consensusClass[row.names(phenotype.GSE14333)]
Subtype2 <- sub("1","S-I",as.character(GSE14333.classlabel))
Subtype2 <- sub("2","S-II",Subtype2)
Subtype2 <- sub("3","S-III",Subtype2)
Subtype2 <- sub("4","S-IV",Subtype2)
Subtype2 <- sub("5","S-V",Subtype2)
names(Subtype2) <- names(GSE14333.classlabel)
location <- as.character(unlist(lapply(phenotype.GSE14333$characteristics_ch1,
                                       function(x){
                                         a=strsplit(as.character(x),";")[[1]];
                                         b=strsplit(a[1],":")[[1]];
                                         b[2]})))
Age <- as.numeric(unlist(lapply(phenotype.GSE14333$characteristics_ch1,
                                function(x){
                                  a=strsplit(as.character(x),";")[[1]];
                                  b=strsplit(a[3],":")[[1]];
                                  b[2]})))
Gender <- as.character(unlist(lapply(phenotype.GSE14333$characteristics_ch1,
                                     function(x){
                                       a=strsplit(as.character(x),";")[[1]];
                                       b=strsplit(a[4],":")[[1]];
                                       b[2]})))
DFS_Time <- as.numeric(unlist(lapply(phenotype.GSE14333$characteristics_ch1,
                                     function(x){
                                       a=strsplit(as.character(x),";")[[1]];
                                       b=strsplit(a[5],":")[[1]];
                                       b[2]})))
DFS_Event <- as.numeric(unlist(lapply(phenotype.GSE14333$characteristics_ch1,
                                      function(x){
                                        a=strsplit(as.character(x),";")[[1]];
                                        b=strsplit(a[6],":")[[1]];
                                        b[2]})))
GSE14333.clinic <- data.frame(cbind(location,Age,Gender,DFS_Time,DFS_Event,Subtype2))
row.names(GSE14333.clinic) <- names(Subtype2)
###location
location.subtype <- data.frame(table(location,Subtype2))
total.location <- rep(as.numeric(table(Subtype2)),each=5)
location.subtype$propotion <- location.subtype$Freq/total.location
ggplot(location.subtype, aes(x=Subtype2, y=propotion, fill=location)) +
  geom_bar(stat="identity")+
  theme_pubr()+
  labs(x="Subtype",y="MSI/MSS")+
  theme(axis.text.x = element_text(angle=0,size=8,face = "bold",hjust = 1))

NA.subtype <- location.subtype[location.subtype$location==" ",]
Rectum.subtype <- location.subtype[location.subtype$location==" Rectum",]
Colon.subtype <- location.subtype[location.subtype$location==" Colon",]
Left.subtype <- location.subtype[location.subtype$location==" Left",]
Right.subtype <- location.subtype[location.subtype$location==" Right",]
######
#library(plyr)
# library(dplyr)
# data <- GSE14333.clinic %>% group_by(location,Subtype2) %>% summarise(Count=n())
# data1 <- GSE14333.clinic %>% group_by(Subtype2) %>% summarise(Count=length(Subtype2))

###Gender
Gender.subtype <- data.frame(table(Gender,Subtype2))
total.Gender <- rep(as.numeric(table(Subtype2)),each=2)
Gender.subtype$propotion <- Gender.subtype$Freq/total.Gender
ggplot(Gender.subtype, aes(x=Subtype2, y=propotion, fill=Gender)) +
  geom_bar(stat="identity")+
  theme_pubr()+
  labs(x="Subtype",y="MSI/MSS")+
  theme(axis.text.x = element_text(angle=0,size=8,face = "bold",hjust = 1))

###Age
Age.subtype <- data.frame(cbind(Age,Subtype2))
ggplot(Age.subtype, aes(x=factor(Subtype2),y=as.numeric(Age),fill=Subtype2))+ 
  geom_boxplot()+
  labs(x="Subtype",y="Ages",title="Boxplot")+
  theme(plot.title = element_text(face="bold",hjust = 0.5,size = 20),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12,angle=45),
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15),
        plot.margin = unit(c(2, 2, 2, 2), "lines"))

###DFS
DFS.subtype <- survfit(Surv(DFS_Time, DFS_Event) ~ Subtype2) 
ggsurv(DFS.subtype,size.est = 1.2,
       surv.col=c("orchid1","brown1","blueviolet","yellow","aquamarine4"),
       cens.size = 3,
       cens.shape = 3,
       order.legend = F) +
  labs(x="Time (days)",y="Survival",title="Disease Free Survival")+
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 15),
        axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12))

#######################

###All Gender,Age
colnames(gender.subtype) <- colnames(Gender.subtype)
Gender.subtype$Gender <- gender.subtype$Gender
All.Gender <- rbind(Gender.subtype,gender.subtype)
ggplot(All.Gender, aes(x=Subtype2, y=propotion, fill=Gender)) +
  geom_bar(stat="identity")+
  theme_pubr()+
  labs(x="Subtype",y="Gender")+
  theme(axis.text.x = element_text(angle=0,size=8,face = "bold",hjust = 1))

colnames(age.subtype) <- colnames(Age.subtype)
All.Age <- rbind(age.subtype,Age.subtype)
ggplot(All.Age, aes(x=factor(Subtype2),y=as.numeric(Age),fill=Subtype2))+ 
  geom_boxplot()+
  labs(x="Subtype",y="Ages",title="Boxplot")+
  theme(plot.title = element_text(face="bold",hjust = 0.5,size = 20),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12,angle=45),
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15),
        plot.margin = unit(c(2, 2, 2, 2), "lines"))



#pheatmap
GSE13294.combat <- expr.combat[specific.genes,names(sort(GSE13294.classlabel))]
library(gplots)
library(pheatmap)
annotation.col <- data.frame(cluster=factor(sort(GSE13294.classlabel),labels = c("1","2","3","4","5")))
annotation.row <- data.frame(spec.gene=factor(c(rep("one",length(one.specific)),rep("two",length(two.specific)),
                                                rep("three",length(three.specific)),rep("four",length(four.specific)),
                                                rep("five",length(five.specific))),labels = c("one","two","three","four","five")))
cluster <- c("red","green","yellow","blue","black")
names(cluster) <- c("1","2","3","4","5")
spec.gene <- c("red","green","yellow","blue","black")
names(spec.gene) <- c("one","two","three","four","five")
ann.colors <- list(cluster=cluster,spec.gene=spec.gene)
pheatmap(GSE13294.combat,cluster_cols = FALSE,cluster_rows =FALSE,
         annotation_col=annotation.col,annotation_row = annotation.row,annotation_colors=ann.colors,
         show_rownames = F,show_colnames = F,color = redblue(100)[100:1])

