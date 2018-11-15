files1=list.files(path="/media/_EXTEND2/zhangli/out_work/SRP069212",full.names=T,recursive=T,pattern="circular")
files2=list.files(path="../../outwork",full.names=T,recursive=T,pattern="circular")
files=c(files1,files2)
data=sapply(c(files),function(x){y=as.matrix(read.table(x,sep='\t'))})
names(data)=sapply(strsplit(names(data),"\\/"),function(x){x[grep("out_dir",x)-1]})
labels=sapply(data,function(x){y=gsub(" ","",paste(x[,1],x[,2],x[,3],x[,6],sep="_"))})
uniq_labels=unique(unlist(labels))
val_circRNA_mat=matrix(0,nrow=length(uniq_labels),ncol=length(data))
for(i in 1:length(data))
{
val_circRNA_mat[match(labels[[i]],uniq_labels),i]=as.numeric(data[[i]][,5])
}
rownames(val_circRNA_mat)=uniq_labels
colnames(val_circRNA_mat)=names(data)
val_circRNA_mat=val_circRNA_mat[,sort(colnames(val_circRNA_mat))]
save(val_circRNA_mat,file="val_circRNA_mat.Rdata")

