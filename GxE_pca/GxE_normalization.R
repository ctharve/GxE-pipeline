### Step 1
##The first step is to get the annotation for the genes, (i.e. transcript length) you may already have that in an 
## object called anno, but if not look at the following code:

geneModelFile <- "/wsu/home/groups/piquelab/data/RefTranscriptome/ensGene.hg19.faCount.gz"
geneAnno <- read.table(geneModelFile,sep="\t",as.is=T)
geneAnno <- geneAnno[,-c(9:12)]
rownames(geneAnno) <- geneAnno$V4
colnames(geneAnno) <-  c("chr","start","stop","t.id","score","strand","c.start","c.stop","ensg","g.id")

## Annotating trnascript lenght in bp and GC content.  
gcContentFile <- "/wsu/home/groups/piquelab/data/RefTranscriptome/ensGene.hg19.faCount.gz";
anno2 <- read.table(gcContentFile,as.is=T,sep="\t",header=T,comment="")
anno2 <- anno2[-191892, ]
rownames(anno2) <- gsub("hg19_ensGene_","",anno2$X.seq)
colnames(anno2) <- c('t.id', 'len', 'A', 'C', 'G', 'T', 'N', 'cpg', 'avg.cg')
anno2$avg.cg <- (anno2$C+anno2$G)/anno2$len
## coding length             
anno2$codLen <- anno2[anno2$t.id,"len"]
## transcript length  
# geneAnno$txLen <- (geneAnno$c.start-geneAnno$start)+(geneAnno$stop-geneAnno$c.stop)+geneAnno$codLen
#Avg. CG on the coding part.                                                                                                   
# geneAnno$avg.cg <- anno2[geneAnno$t.id,"avg.cg"]
rm(anno2)

##===============================
## STEP 2
## d is my matrix with the counts, rows transcripts, columns samples
cs <- colSums(d)/1E6

d2 <- mysapply(1:ncol(d),function(jj){d[,jj]/cs[jj]*1000/geneAnno$codLen})
cs <- total number of reads after filtering for a single transcript / 1E6 across all samples 
## d2 has now reads per kb of transcript, you can divide additional by the total number of reads mapped in millions and you will have rpkm, 
## but I don't do it because the I column normalize which makes dividing by a constant kind of useless. 

## ======================
## Optional step 3
## gc correction and transcript length correction. 

d3 <- sapply(1:ncol(d2),function(ii){
  fit.lm <- lm( log10(d2[,ii]+0.001) ~ geneAnno2$avg.cg + geneAnno2$codLen + geneAnno2$avg.cg * geneAnno2$codLen + geneAnno2$utrLen)
  x <- residuals(fit.lm)
})


##Column normalization                                                                                                              
d4 <- apply(d3,2,function(x){qqnorm(rank(x, ties.method = "random"), plot = F)$x})
d4 <- as.matrix(d4)