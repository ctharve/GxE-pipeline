##################################################################
## Chris Harvey
##
## Assemble a Complete matrix of GxE data
##
##################################################################
require(parallel)
cores <- as.integer(Sys.getenv("NCPUS"))
if(cores<1){cores <- 1}
ParallelSapply <- function(...,mc.cores=cores){
  simplify2array(mclapply(...,mc.cores=mc.cores))
}

##################################################################
x11(display="localhost:10.0" ,type="Xlib")

plate_seq <- 1:12
plates <- paste('P', plate_seq[-c(3, 6, 9)], sep='')

all_data <- invisible(
            	lapply(plates, function(ii){       ## across all samples within 
            		#ii <-'P1'	
            		cat(paste('Plate: ', ii, sep=''))
            		cat('\n')
            	
            		load(paste('./', ii, '/counts/', ii, '.data.Rd', sep=''))

            		cat(paste('Colnames & Transcript.IDs identical? ', identical(rownames(data), anno$t.id), sep=''))
            		cat('\n')
            		cat(paste('Annotation transcript count: ', dim(anno)[1], sep=''))
            		cat('\n')

            		data <- data.frame(data)
            		data$ensg <- anno$ensg

            		data
            	})
)
names(all_data) <- plates 

## for each plate, and sample within a plate, take the median count for each transcript within a gene
gene_data <- lapply(all_data, function(this_data){       # all plates
                #this_data <- all_data[[1]]
                ## TODO: '-97' is dirty dirty, fix it, but not now, becuase it works, for now.
                samples <- colnames(this_data)[-97]
                unique_genes <- unique(this_data$ensg)

                sapply(samples, function(this_sample){      # all samples within a plate
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
