genefilter <- function(exprdata, Cutoff) {
  # Restrict analysis to genes with sample var above threshold
  nLastCol <- dim(exprdata)[2]
  print("Removing genes with low sample variance...")
  
  sampleVar <- apply(exprdata,1,var, na.rm=TRUE)
  # compute sample variance
  idxHighVar <- sampleVar > Cutoff
  data.highVar <- exprdata[idxHighVar,]
  sampleVar.highVar <- sampleVar[idxHighVar]
  print(paste(" Found ",dim(data.highVar)[1]," genes exceeding SD threshold of
              ",Cutoff,".",sep=""))
  nGenesHighVar <- dim(data.highVar)[1]
  # Cannot have replicate gene names for clustering; choose gene with largest
  #sample variance
  genesUnique <- as.vector(unique(data.highVar[,1]))
  nGenesUnique <- length(genesUnique)
  nSamples <- dim(data.highVar)[2]-1; # first col is still gene name
  data <- array(dim=c(nGenesUnique,nLastCol))
  data[,1] <- genesUnique
  colnames(data) <- colnames(data.highVar)
  print("Removing duplicate genes (selecting for max standard deviation)...")
  
  for (gene in genesUnique) {
    # index/indices of all genes matching gene (detect duplicates)
    idxGenes <- seq(along=1:nGenesHighVar)[data.highVar[,1]==gene]
    
    data.slice <- data.highVar[idxGenes,2:nLastCol]
    if (length(idxGenes)>1) {
      idxMaxVar <- which.max(sampleVar.highVar[idxGenes])
      # find dupls with max var
      data.slice <- data.slice[idxMaxVar,]
    }
    data[data[,1]==gene,2:nLastCol] <- as.matrix(data.slice)
  }
  print(paste(" ",nGenesUnique," unique genes IDs remain.",sep=""))
  outFile <- paste(strsplit(exprFile,".txt"),"_sd",sdCutoff,".txt",sep="")
  write.table(file=outFile, data, quote=FALSE, row.names=FALSE, sep="\t")
  print(paste("Screened data written to file: ", outFile, sep=""))
  print("")
  rm(data.all,data.mapped,data.highVar,data)
}