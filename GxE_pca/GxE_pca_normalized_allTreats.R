##################################################################
##
## Chris Harvey
## 07182014
## identify control only columns and conduct PCA on only controls 
##
##################################################################
require(parallel)
require(ggplot2)
cores <- as.integer(Sys.getenv("NCPUS"))
if(cores<1 || is.na(cores)){cores <- 1}
## need to get this working for lapply
ParallelSapply <- function(...,mc.cores=cores){
  simplify2array(mclapply(...,mc.cores=mc.cores))
}

##################################################################
x11(display="localhost:10.0" ,type="Xlib")
#################
plate_seq <- 1:12
plates <- paste('P', plate_seq[-c(3, 6, 9)], sep='')
cellLine <- c('LCL', 'PBMC', 'HUVEC', 'PBMC', 'LCL', 'HUVEC', 'SMC', 'MEL', 'MEL')

#################
## all_data.Rd is from /wsu/home/groups/piquelab/charvey/GxE/derived_data.old/scripts/GxE_assemble_tables.R
load('all_data_centered.Rd')
#load('all_data.Rd')
annotations <- data_centered$anno
#all_data$anno <- NULL
data_centered$anno <- NULL

####################################################################
## Annotating trnascript length in bp and GC content.  
gcContentFile <- "/wsu/home/groups/piquelab/data/RefTranscriptome/ensGene.hg19.faCount.gz";
anno2 <- read.table(gcContentFile,as.is=T,sep="\t",header=T,comment="")
anno2 <- anno2[-191892, ]
rownames(anno2) <- gsub("hg19_ensGene_","",anno2$X.seq)
anno2$avg.cg <- (anno2$C+anno2$G)/anno2$len
colnames(anno2) <- c('t.id', 'len', 'A', 'C', 'G', 'T', 'N', 'cpg', 'avg.cg')
## coding length             
anno2$codLen <- anno2[,"len"]

## combine raw data into a single matrix
full_matrix <- do.call(cbind, data_centered[1:9])   ## remove annotations; hence [1:9]

####################################################################
## combine all of the covariates files into a single file
####################################################################
all_cvs <- lapply(plates, function(plate){
  #plate <- plates[1]
  filename <- paste('/wsu/home/groups/piquelab/charvey/GxE/derived_data.old/cvs/GxE_', plate, '_covariates.txt', sep='')
  cov <- read.table(filename, sep='\t', as.is=TRUE, header=TRUE)
  plate <- cov$Plate.ID
  barcodes <- as.integer(cov$Barcode.ID)
  cellLine <- cov$CellLine
  treat <-  cov$Treatment.ID
  ctrl <- cov$Control.ID
  list(plate, barcodes, cellLine, treat, ctrl)
})
names(all_cvs) <- plates

## unlist the above results
plate <- unlist(sapply(plates, function(plate){all_cvs[[plate]][1]}))
barcodes <- unlist(sapply(plates, function(plate){all_cvs[[plate]][2]}))
cellLine <- unlist(data.frame(sapply(plates, function(plate){all_cvs[[plate]][3]}), stringsAsFactors=FALSE))
treat <- unlist(data.frame(sapply(plates, function(plate){all_cvs[[plate]][4]}), stringsAsFactors=FALSE))
ctrl <- unlist(data.frame(sapply(plates, function(plate){all_cvs[[plate]][5]}), stringsAsFactors=FALSE))
names(plate) <- NULL
names(barcodes) <- NULL
names(cellLine) <- NULL
names(treat) <- NULL
names(ctrl) <- NULL
cvs <- data.frame(plate, barcodes, cellLine, treat, ctrl, stringsAsFactors=FALSE)
rm(plate, barcodes, cellLine, treat, ctrl)
####################################################################
# normalization to transcript length, and quantile normalization
####################################################################
## 1.) use the indicator to select columns which correspond to controls
## LOOK HERE: in this case we are not subsetting controls so ctrl_matrix is still the full_matrix
ctrl_matrix <- full_matrix[, ]
##ctrl_centered <- t(apply(ctrl_matrix, 1, function(row){row-mean(row)}))

## 2.) normalize the ctrl_matrix to transcript length
sample_sums <- colSums(ctrl_matrix)/1E6

ctrl_normal <- sapply(1:ncol(ctrl_matrix),function(jj){
  ctrl_matrix[,jj]/sample_sums[jj]*1000/anno2$codLen
})

## 2.a.)

d2 <- ctrl_normal

rs <- rowMeans(d2>10)
ind <- rs>1/5
sum(ind)

d2 <- d2[ind,]
anno3 <- anno2[ind,]

##d3 <- sapply(1:ncol(d2),function(ii){
##    fit.lm <- lm( log10(d2[,ii]+0.01) ~ anno3$avg.cg + anno3$codLen + anno3$avg.cg * anno3$codLen)
##      x <- residuals(fit.lm)
##  })
d4 <- apply(d2, 2, function(x){qqnorm(rank(x, ties.method = "random"), plot = F)$x})
d4 <- as.matrix(d4)

## 3.) normalize each feature/transcript
ctrl_centered <- t(apply(d4, 1, function(row){row-mean(row)}))
colnames(ctrl_centered) <- colnames(ctrl_matrix)
####################################################################
# create a list of averaged controls for each plate, cellLine, & control
# ctrl_centered ~ normalized data
# cvs ~ all covariates files
####################################################################
str(cvs)
plates <- unique(cvs$plate)
avg_ctrl <- list()
for(pl in plates){
  #pl <- 'P1'
  lines <- unlist(unique(subset(cvs, plate==pl, cellLine)))
  
  for(line in lines){
    #line <- "19239"
    ctrls <- unlist(unique(subset(cvs, plate==pl & cellLine==line, ctrl)))

    for(cl in ctrls){
      #cl <- 'CO1'
      out_tag <- paste(pl, line, cl, sep='_')
      print(out_tag)
      bcs <- subset(cvs, plate==pl & cellLine==line & treat==cl, barcodes)
      samples <- sapply(bcs, function(bc){paste(pl, '-HT', bc, sep='')})
      avg_ctrl[[out_tag]] <- apply(ctrl_centered[, samples], 1, mean)
    }
  }
}
####################################################################
# create a list of averaged controls for each plate, cellLine, & control
# ctrl_centered ~ normalized data
# cvs ~ all covariates files
# avg_ctrl ~ average control for each plate, cellLine, and control
####################################################################
all_samples <- colnames(ctrl_centered)
normed_dat <- ParallelSapply(all_samples, function(this.sam){
  #this.sam <- all_samples[1]

  aux <- strsplit(this.sam, '-')
  pl <- aux[[1]][1]
  bc <- gsub('HT', '', aux[[1]][2])
  fields <- subset(cvs, plate==pl & barcodes==bc)
  this.ctrl <- paste(fields$plate, fields$cellLine, fields$ctrl, sep='_')

  ctrl_centered[, this.sam] - avg_ctrl[[this.ctrl]] 

})

cor.mat <- cor(normed_dat)

mycol.plate <- c("red","blue","green","pink","purple","gold","darkgreen","darkred","lightblue")
names(mycol.plate) <- c("P1","P2","P4","P5","P7","P8","P10","P11","P12")


mycol <- c("red","blue","green","pink","purple")
names(mycol) <- c("LCL","PBMC","HUVEC","SMC","MEL")
#########
colvec <- mycol[cellLine[ctrl_plates]]
colrow <- mycol.plate[ctrl_plates]

#pdf('heatmap_corr_mat.pdf')
hc <- hclust(dist(ctrl_centered), method = "average", members = NULL)
##hc <- hclust(as.dist(1-(1+cor.mat)/2), method = "average", members = NULL)

hd <- as.dendrogram(hc)

pdf('heatmap_corr_plates_cellTypes_noDendro.pdf')

heatmap.2(cor.mat,Rowv=hd,Colv=hd,trace="none",ColSideColor=colvec,RowSideColor=colrow)

dev.off()
#########
## 4.) create the covariance matrix for the t.id normalized, and transcript normalized controls then do SVD
ctrl_cov <- t(normed_dat) %*% normed_dat
ctrl_svd <- svd(ctrl_cov)
n.samples <- length(ctrl_svd$d)
save(ctrl_svd, file='data_cov_allDat_contrlNormalized_svd.Rd')

## percentage of variance explained by PCs
## 07182014
(ctrl_svd$d[1]^2)/sum(ctrl_svd$d^2) ## squared 0.893188 0.3979114
(ctrl_svd$d[2]^2)/sum(ctrl_svd$d^2) ## squared 0.101896 0.182582

## 5.) get plates, barcodes, and treatments
indic_mat <- data.frame(sapply(plates, function(plate){ctrl_indic[[plate]][1]}))
numb_controls <- colSums(indic_mat)
ctrl_plates <- unlist(sapply(plates, function(plate){n <- numb_controls[plate]; rep(plate, n)}))
ctrl_barcodes <- barcodes[indic]
ctrl_treatments <-  controls[indic]

##plot_dat <- data.frame(cellLine=cellLine[ctrl_plates],plates=ctrl_plates, bcs=ctrl_barcodes, treats=ctrl_treatments, PC1=ctrl_svd$v[, 1], PC2=ctrl_svd$v[, 2], PC3=ctrl_svd$v[, 3], PC4=ctrl_svd$v[, 4], PC5=ctrl_svd$v[, 5], stringsAsFactors=FALSE)

v <- ctrl_svd$v
colnames(v) <- paste("PC",1:ncol(v),sep="")

plot_dat <- data.frame(cellLine=cellLine[ctrl_plates],plates=ctrl_plates, bcs=ctrl_barcodes, treats=ctrl_treatments, v, stringsAsFactors=FALSE)


pdf('PCA_cellLines_normalized_5_6.pdf', width=15, height=15)
my_plot <- ggplot(plot_dat, aes(x=PC5, y=PC6)) + geom_point(aes(color=cellLine)) + ggtitle('Control only PCA')
my_plot + geom_text(aes(label=treats, color=cellLine))
dev.off()

pdf('variance_explained.pdf')
plot(data_cov_svd$d^2/sum(data_cov_svd$d^2), xlim = c(0, 15), type = "b", pch = 16, xlab = "principal components", 
    ylab = "variance explained", main='Variance Explained')
dev.off()


##
## The End
##
