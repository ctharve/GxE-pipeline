##################################################################
##
## Chris Harvey
## June 9th 2014
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
load('all_data_centered.Rd')
load('all_data.Rd')
annotations <- data_centered$anno
all_data$anno <- NULL
data_centered$anno <- NULL

####################################################################
## Annotating trnascript lenght in bp and GC content.  
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
## need to identify the controls for control_only PCA
####################################################################
# extract a list of vectors containinng the indicator and controls for each plate_barcode
## TODO extract cell type & cell line
ctrl_indic <- lapply(plates, function(plate){
  #plate <- plates[1]
  filename <- paste('/wsu/home/groups/piquelab/charvey/GxE/derived_data/cvs/GxE_', plate, '_covariates.txt', sep='')
  cov <- read.table(filename, sep='\t', as.is=TRUE, header=TRUE)
  barcodes <- as.integer(cov$Barcode.ID)
  treat <-  cov$Treatment.ID
  indic <- treat %in% c('CO1', 'CO2', 'CO3')
  list(indic, treat, barcodes)
})
names(ctrl_indic) <- plates

## unlist the above results
indic <- unlist(sapply(plates, function(plate){ctrl_indic[[plate]][1]}))
controls <- unlist(data.frame(sapply(plates, function(plate){ctrl_indic[[plate]][2]}), stringsAsFactors=FALSE))
barcodes <-  unlist(sapply(plates, function(plate){ctrl_indic[[plate]][3]}))

## 1.) use the indicator to select columns which correspond to controls
ctrl_matrix <- full_matrix[, indic]
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

cor.mat <- cor(ctrl_centered)

mycol.plate <- c("red","blue","green","pink","purple","gold","darkgreen","darkred","lightblue")
names(mycol.plate) <- c("P1","P2","P4","P5","P7","P8","P10","P11","P12")


mycol <- c("red","blue","green","pink","purple")
names(mycol) <- c("LCL","PBMC","HUVEC","SMC","MEL")

colvec <- mycol[cellLine[ctrl_plates]]
colrow <- mycol.plate[ctrl_plates]

#pdf('heatmap_corr_mat.pdf')
hc <- hclust(dist(ctrl_centered), method = "average", members = NULL)
##hc <- hclust(as.dist(1-(1+cor.mat)/2), method = "average", members = NULL)

hd <- as.dendrogram(hc)

pdf('heatmap_corr_plates_cellTypes_noDendro.pdf')

heatmap.2(cor.mat,Rowv=hd,Colv=hd,trace="none",ColSideColor=colvec,RowSideColor=colrow)

dev.off()

## 4.) create the covariance matrix for the t.id normalized, and transcript normalized controls then do SVD
ctrl_cov <- t(ctrl_centered) %*% ctrl_centered
ctrl_svd <- svd(ctrl_cov)
n.samples <- length(ctrl_svd$d)
save(ctrl_svd, file='data_cov_controls_normalized_svd.Rd')

## percentage of variance explained by PCs
(ctrl_svd$d[1]^2)/sum(ctrl_svd$d^2) ## squared 0.893188
(ctrl_svd$d[2]^2)/sum(ctrl_svd$d^2) ## squared 0.101896

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
