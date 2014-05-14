##################################################################
## Chris Harvey
##
## May 14, 2013
##
##################################################################
plates <- 1:12
plates <- plates[-c(3, 6, 9)]

invisible(
	sapply(plates, function(ii){
		#ii <-1	
		cat(paste('Plate ', ii, sep=''))
		cat('\n')
	
		load(paste('./P', ii, '/counts/P', ii, '.data.Rd', sep=''))

		cat(paste('Colnames & Transcript.IDs identical? ', identical(rownames(data), anno$t.id), sep=''))
		cat('\n')
		cat(paste('Annotation transcript count: ', dim(anno)[1], sep=''))
		cat('\n')

		data <- data.frame(data)
		data$ensg <- anno$ensg

		assign(paste('data.P', ii, sep=''), data, envir=globalenv())

		assign('anno', anno, envir=globalenv())
		rm(data)

	})
)







unique_genes <- unique(data.P1$ensg)



format(object.size(data.P1), units='MB')
