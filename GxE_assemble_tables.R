##################################################################
## Chris Harvey
##
## Assemble a Complete matrix of GxE data
##
##################################################################
require(parallel)
require(ggplot2)
cores <- as.integer(Sys.getenv("NCPUS"))
if(cores<1){cores <- 1}
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
all_data$anno <- NULL
data_centered$anno <- NULL

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
data_cov_svd$d[1]/sum(data_cov_svd$d) ## .4675
data_cov_svd$d[2]/sum(data_cov_svd$d) ## .1966

percent_var <- data_cov_svd$d


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

cell_lines <- c(rep('LCL', 96), rep('PBMC.1', 96), rep('HUVEC.1', 96), rep('PBMC.2', 96), rep('LCL', 96), rep('HUVEC.2', 96), rep('SMC', 96), rep('Melanocytes', 96), rep('Melanocytes', 96))

plot_dat <- data.frame(IDS=group_labels, Plate=plate, Lines=cell_lines, PC1=data_cov_svd$v[, 1], PC2=data_cov_svd$v[, 2])

pdf('PCA_cellLines.pdf')
my_plot <- ggplot(plot_dat, aes(x=PC1, y=PC2)) + geom_point(aes(color=Lines)) + ggtitle('Cell Line Specific PCA')
my_plot
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