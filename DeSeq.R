library(DESeq2)
library(EnhancedVolcano)
theme_set(theme_minimal())
countMatrix<-readRDS("processed_data/countMatrix.RDS")

counts<-countMatrix$counts
samplenames<-matrix(nrow=2,unlist(strsplit(colnames(counts),"\\.")))[1,]
colnames(counts)<-samplenames



coldata = data.frame(condition=paste0(c(rep("Fulvestrant.",8),rep("Vehicle.",8)),rep(c("Hypoxia","Normoxia"),8)),rep=gsub("[^0-9.-]", "", samplenames))
rownames(coldata)<-samplenames
coldata

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)
dds<- DESeq(dds)

#PCA

se<-SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                     colData=colData(dds))
plotPCA( DESeqTransform( se ) )

######
#Rep4 looks like an Outlier, removing for now.
#############


removeSamples<-grep("4",colnames(counts), value=TRUE)

filteredCounts<-counts[,colnames(counts)[!colnames(counts) %in% removeSamples]]



coldata2 = data.frame(condition=paste0(c(rep("Fulvestrant.",6),rep("Vehicle.",6)),rep(c("Hypoxia","Normoxia"),6)),rep=gsub("[^0-9.-]", "", colnames(filteredCounts)))
rownames(coldata2)<-colnames(filteredCounts)
coldata2

dds2 <- DESeqDataSetFromMatrix(countData = filteredCounts,
                              colData = coldata2,
                              design = ~ condition)
dds2 <- DESeq(dds2)


#PCA


se2<-SummarizedExperiment(log2(counts(dds2, normalized=TRUE) + 1),
                         colData=colData(dds2))
plotPCA( DESeqTransform( se2) )
plotPCA( DESeqTransform( se2),  intgroup=c("rep")  )



#Function to convert entrez ids to gene symbols.
library(org.Hs.eg.db)
eg2sym<-function(x){
  
  
  eg2symmap<-as.list( org.Hs.egSYMBOL[mappedkeys( org.Hs.egSYMBOL)])
  x<-as.character(x)
  out<-eg2symmap[x]
  names(out)<-x
  out2<-sapply(out,function(x){
    if(is.null(x)){
      return(NA)
    } else {
      return(x[1])
    }
  })
  out3<-unlist(out2)
  out<-setNames(out3,names(out))
  return(out)
}


resVNvFN<- results(dds2, contrast = c("condition","Vehicle.Normoxia", "Fulvestrant.Normoxia"))
resVNvFN
EnhancedVolcano(resVNvFN,
                lab = eg2sym(rownames(resVNvFN)),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "Vehicle.Normoxia - Fulvestrant.Normoxia" 
                )


resVHvFH<- results(dds2, contrast = c("condition","Vehicle.Hypoxia", "Fulvestrant.Hypoxia"))
resVHvFH
EnhancedVolcano(resVHvFH,
                lab = eg2sym(rownames(resVHvFH)),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "Vehicle.Hypoxia - Fulvestrant.Hypoxia" 
                )


resFHvFN<- results(dds2, contrast = c("condition", "Fulvestrant.Hypoxia","Fulvestrant.Normoxia"))
resFHvFN
EnhancedVolcano(resFHvFN,
                lab = eg2sym(rownames(resFHvFN)),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "Fulvestrant.Hypoxia - Fulvestrant.Normoxia"
                )

resVHvVN<- results(dds2, contrast = c("condition", "Vehicle.Hypoxia","Vehicle.Normoxia"))
resVHvVN
EnhancedVolcano(resVHvVN,
                lab = eg2sym(rownames(resVHvVN)),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "Vehicle.Hypoxia - Vehicle.Normoxia"
                )


resVHvVN<- results(dds2, contrast = c("condition", "Vehicle.Hypoxia","Vehicle.Normoxia"))
resVHvVN
EnhancedVolcano(resVHvVN,
                lab = eg2sym(rownames(resVHvVN)),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "Vehicle.Hypoxia - Vehicle.Normoxia"
)


#What are the signifant gene sets
FulvGenes<-resVNvFN[!is.na(resVNvFN$padj),]
FulvGenes<-FulvGenes[FulvGenes$padj<0.05,]
FulvGenes

HypoxGenes<-resVHvVN[!is.na(resVHvVN$padj),]
HypoxGenes<-HypoxGenes[HypoxGenes$padj<0.05,]
HypoxGenes


