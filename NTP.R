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