##################################################################
## Chris Harvey
##
## Assemble a Complete matrix of GxE data
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

plate_seq <- 1:12
plates <- paste('P', plate_seq[-c(3, 6, 9)], sep='')

## aggregate all count objects and output annotations and count data into .gz
all_data <- invisible(
            	lapply(plates, function(ii){       ## across all samples within 
            		#ii <-'P1'	
            		cat(paste('Plate: ', ii, sep=''))
            		cat('\n')
            	
            		load(paste('./', ii, '.data.Rd', sep=''))

            		cat(paste('Colnames & Transcript.IDs identical? ', identical(rownames(data), anno$t.id), sep=''))
            		cat('\n')
            		cat(paste('Annotation transcript count: ', dim(anno)[1], sep=''))
            		cat('\n')

                    #annoFile <- paste(ii,"_anno.gz",sep="")
                    #write.table(anno,gzfile(annoFile),sep="\t",quote=FALSE, row.names=FALSE)
                    #dataFile <- paste(ii,"_data.gz",sep="")
                    #write.table(data,gzfile(dataFile),sep="\t",quote=FALSE, row.names=FALSE)

            		data
            	})
)
names(all_data) <- plates 
load(paste('./P1.data.Rd', sep=''))
all_data$anno <- anno
all_data_file <- paste("all_data.Rd",sep="")
save(all_data, file=all_data_file)
all_data$anno <- NULL

## row normalize all data
data_centered <- invisible(

lapply(plates, function(ii){       ## across all samples within 
                    #ii <- plates[1]  
                    cat(paste('Row centering Plate: ', ii, sep=''))
                    cat('\n')
                    t(apply(all_data[[ii]], 1, function(row){row-mean(row)}))
                })
)
names(data_centered) <- plates
data_centered$anno <- anno
data_centered_file <- paste("all_data_centered.Rd",sep="")
save(data_centered ,file=data_centered_file)
data_centered$anno <- NULL

#################
load('all_data_centered.Rd')
load('all_data.Rd')
annotations <- data_centered$anno
#all_data$anno <- NULL
data_centered$anno <- NULL
potatoes <- do.call(cbind, data_centered)

cvs.directory <- '/wsu/home/fl/fl97/fl9788/piquelab/charvey/GxE/derived_data/cvs'

idx_keep <- sapply(plates, function(plate){
                                        #plate <- plates[1]
                                        table <- paste('../../cvs/', 'GxE_', plate, '_covariates.txt', sep='')
                                        fields <- read.table(table, as.is=TRUE, header=TRUE, sep='\t')
                                        treatment <- fields[, c('Treatment.ID')]
                                        keep <- treatment %in% c('CO1', 'CO2', 'CO3') 
})

####################################################################
## Normalize the data to transcript length ##
## Annotating trnascript lenght in bp and GC content.  
gcContentFile <- "/wsu/home/groups/piquelab/data/RefTranscriptome/ensGene.hg19.faCount.gz";
anno2 <- read.table(gcContentFile,as.is=T,sep="\t",header=T,comment="")
anno2 <- anno2[-191892, ]
rownames(anno2) <- gsub("hg19_ensGene_","",anno2$X.seq)
colnames(anno2) <- c('t.id', 'len', 'A', 'C', 'G', 'T', 'N', 'cpg', 'avg.cg')
anno2$avg.cg <- (anno2$C+anno2$G)/anno2$len
## coding length             
anno2$codLen <- anno2[,"len"]

## combine raw data into a single matrix
full_matrix <- do.call(cbind, data_centered[1:9])   ## remove annotations; hence [1:9]
sample_sums <- colSums(full_matrix)/1E6

full_normal <- sapply(1:ncol(full_matrix),function(jj){
    full_matrix[,jj]/sample_sums[jj]*1000/anno2$codLen
    })
## d2 has now reads per kb of transcript, you can divide additional by the total number of reads mapped in millions and you will have rpkm, 
## but I don't do it because the I column normalize which makes dividing by a constant kind of useless.
full_cov_norm <- t(full_normal) %*% full_normal
full_norm_svd <- svd(full_cov_norm)
save(full_norm_svd, file='data_cov_normalized_svd.Rd')

####
potatoes <- do.call(cbind, data_centered)
data_svd <- svd(potatoes)
save(data_svd, file='data_svd.Rd')

####
data_cov <- t(potatoes) %*% potatoes
data_cov_svd <- svd(data_cov)
save(data_cov_svd, file='data_cov_svd.Rd')

#######################################################################
# data_svd$d^2
# data_cov_svd$d
#######################################################################
load('data_cov_svd.Rd')
load('all_data.Rd')
n.samples <- length(data_cov_svd$d)
(data_cov_svd$d[1]^2)/sum(data_cov_svd$d^2) ## .4675 and squared 0.8154672
(data_cov_svd$d[2]^2)/sum(data_cov_svd$d^2) ## .1966 and squared 0.1443409

## extract group labels from 
group_labels <- c(
    sapply(plates, function(this_plate){
        colnames(all_data[[this_plate]])
        }
    )
)

plate <- sapply(1:n.samples, function(ii){
    #ii <- 1
    temp <- strsplit(group_labels[ii], '-')
    temp[[1]][1]
})

bc <- sapply(1:n.samples, function(ii){
    #ii <- 1
    temp <- strsplit(group_labels[ii], '-')
    temp[[1]][2]
})

cell_lines <- c(rep('LCL.P1', 96), rep('PBMC.P2', 96), rep('HUVEC.P4', 96), rep('PBMC.P5', 96), rep('LCL.P7', 96), rep('HUVEC.P8', 96), rep('SMC.P10', 96), rep('Melanocytes.P11', 96), rep('Melanocytes.P12', 96))

initial_dat <- data.frame(IDS=group_labels, Plate=plate, Lines=cell_lines, PC1=data_cov_svd$v[, 1], PC2=data_cov_svd$v[, 2], stringsAsFactors=FALSE)

## extract the treatment names from the covariates file. 
treatments_full <- c(sapply(plates, function(this_plate){
    #this_plate <- plates[1]

    barcodes <- as.integer(
                    gsub(paste(this_plate, '-HT', sep=''), '', initial_dat[initial_dat$Plate==this_plate, 'IDS'])
                )
    covs <- read.table(paste('/wsu/home/groups/piquelab/charvey/GxE/derived_data/cvs/GxE_', this_plate, '_covariates.txt', sep=''), header=TRUE, sep='\t', as.is=TRUE)
    t(covs[, c('Treatment.ID')])

}))

treats_split <- strsplit(treatments_full, 'C')
aux <- sapply(treats_split, function(this_treat){
    cat(this_treat)
})

n.treats <- length(treats_split)
treat_con <- cbind(rep(NA, n.treats), rep(NA, n.treats))
colnames(treat_con) <- c('treatment', 'concentration')

for(ii in 1:n.treats){
    #ii <- 1
    temp_treat <- treats_split[[ii]][1]
    temp_control <- treats_split[[ii]][2]
    if(nchar(temp_treat)==0){
        treat_con[ii, 'treatment'] <- paste('C', temp_control, sep='')
        treat_con[ii, 'concentration'] <- NAll

    }else{
        treat_con[ii, 'treatment'] <- temp_treat
        treat_con[ii, 'concentration'] <- temp_control
    }
}


plot_dat <- data.frame(treatments_full, treat_con, initial_dat, stringsAsFactors=FALSE)


pdf('PCA_cellLines_annotated.pdf', width=15, height=15)
my_plot <- ggplot(plot_dat, aes(x=PC1, y=PC2)) + geom_point(aes(color=Lines)) + ggtitle('Cell Line Specific PCA')
my_plot + geom_text(aes(label=treatment, color=Lines))
dev.off()

pdf('PCA_plates.pdf')
my_plot <- ggplot(plot_dat, aes(x=PC1, y=PC2)) + geom_point(aes(color=Plate)) + ggtitle('Plate Specific PCA')
my_plot
dev.off()


pdf('variance_explained.pdf')
plot(data_cov_svd$d^2/sum(data_cov_svd$d^2), xlim = c(0, 15), type = "b", pch = 16, xlab = "principal components", 
    ylab = "variance explained", main='Variance Explained')
dev.off()



## for each plate, and sample within a plate, take the median count for each transcript within a gene
gene_data <- lapply(all_data, function(this_data){       # all plates
                #this_data <- all_data[[1]]
                ## TODO: '-97' is dirty dirty, fix it, but not now, becuase it works, for now.
                samples <- colnames(this_data)[-97]
                unique_genes <- unique(this_data$ensg)

                potato <- sapply(samples, function(this_sample){      # all samples within a plate
                    #this_sample <- samples[1]
                    ParallelSapply(unique_genes, function(this_gene){       # across all genes within a sample within a plate
                        #this_gene <- unique_genes[1]
                        cat(this_gene)
                        median(
                            subset(this_data[, c(this_sample, 'ensg')], ensg==this_gene, select=this_sample)[, 1]
                        )

                    })  
                })
})


format(object.size(all_data), units='GB')

##
## The End
##
