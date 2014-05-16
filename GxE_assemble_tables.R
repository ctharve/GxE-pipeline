##################################################################
## Chris Harvey
##
## Assemble a Complete matrix of GxE data
##
##################################################################
require(parallel)
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

                    annoFile <- paste(ii,"_anno.gz",sep="")
                    write.table(anno,gzfile(annoFile),sep="\t",quote=FALSE, row.names=FALSE)
                    dataFile <- paste(ii,"_data.gz",sep="")
                    write.table(data,gzfile(dataFile),sep="\t",quote=FALSE, row.names=FALSE)

            		data
            	})
)
names(all_data) <- plates 

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

potatoes <- do.call(cbind, data_centered)

data_svd <- svd(potatoes)

## 5:52pm
data_cov <- t(potatoes) %*% potatoes

data_cov_svd <- svd(data_cov)


# data_cov_svd should be equal to  data_svd$d

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



dataFile=paste(platePrefix,".data.gz",sep="")
write.table(data,gzfile(dataFile),sep="\t",quote=F)
