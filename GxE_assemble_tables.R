##################################################################
## Chris Harvey
##
## Assemble a Complete matric of GxE data
##
##################################################################
plate_seq <- 1:12
plates <- paste('P', plate_seq[-c(3, 6, 9)], sep='')

all_data <- invisible(
            	lapply(plates, function(ii){
            		#ii <-1	
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

gene_data <- invisible(
                lapply(all_data, function(this_data){
                    this_data <- all_data[[1]]

                    unique_genes <- unique(this_data$ensg)
                    


                })
)


unique_genes <- unique(data.P1$ensg)



format(object.size(all_data), units='MB')
