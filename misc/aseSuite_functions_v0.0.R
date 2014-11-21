##################################################################
## Combines pileups accross samples by first creating the union 
## of all loci, then combining reads.
##
## Created: 08/06/2013
##
## Author: RPR
##
## Version0.1: use to determine if samples are from the same individual
##################################################################  

##################################################################  
## Find the union of and extract relevant fields from a set of samples
UnionExtractFields <- function(fileList, combine=FALSE){			
	#browser()
	tmpFile <- scan(pipe("mktemp -t"),character(0))
	#cmd <- paste("less ", paste(fileList, collapse=" "), " | grep -v -w '^chr\\|^NA' | cut -f 1-3,6-7 | bedSort stdin stdout | uniq | gzip > ", tmpFile)
        ## this command in includes the reference and alternate allele in the annotations
        cmd <- paste("less ", paste(fileList, collapse=" "), " | grep -v -w '^chr\\|^NA' | cut -f 1-7 | bedSort stdin stdout | uniq | gzip > ", tmpFile)
	cat(cmd,"\n");
	system(cmd)
	##TODO a better way to do this is to have the directory and the file list as separate arguments.
	#sNames <- strsplit(gsub("/nfs/hpnfs/groups/piquelab/charvey/ASE/prod2/", "", gsub(".pileup.clean.bed.gz", "", fileList)), "/")
	#sNames <- gsub("_NEB_Rep1.q10.p10", "", sapply(sNames, `[`, 1))	
	
	anno <- read.table(gzfile(tmpFile),sep="\t",as.is=T)	

        ## returns a 2D (3Xn.files) list, a list of lists, where rows are ref, alt, err, and each ref list is the ref from a sample
	aux <- sapply(fileList,function(fn){			
				#fn <- fileList[1]
				cat("Processing:",fn,"\n")
				#command=paste("intersectBed -a ",tmpFile," -b ",fn," -wao | cut -f 1-3,13-15 ",sep="")
                                command=paste("intersectBed -a ",tmpFile," -b ",fn," -wao | cut -f 1-3,15-17",sep="")
                                cat(command,"\n")
				aa <- read.table(pipe(command),sep="\t",as.is=T,na.strings=".")
				aa[is.na(aa)] <- 0
				stopifnot(identical(aa[,1:3],anno[,1:3]))
				aa[,-(1:3)]
			})
	
	#colnames(anno) = c("chr","pos0","pos", "rsID", "af")
        colnames(anno) = c("chr","pos0","pos", "ref", "alt", "rsID", "af")
	
	Ref <- as.matrix(do.call(cbind,aux[1,]))
	#colnames(Ref) <- sNames
	Alt <- as.matrix(do.call(cbind,aux[2,]))
	#colnames(Alt) <- sNames
	Err <- as.matrix(do.call(cbind,aux[3,]))
	#colnames(Err) <- sNames
	
	returnlist <- list(ref=Ref,alt=Alt,err=Err,anno=anno);
	
	if(combine==TRUE){
		allRef<-apply(Ref, MARGIN=1, sum)
		allAlt<-apply(Alt, MARGIN=1, sum)
		allErr<-apply(Err, MARGIN=1, sum)
		
		returnlist$all<-as.matrix(cbind(allRef, allAlt, allErr))
		
	}
	return(returnlist)
}

##################################################################  
## TODO, set .pileup.clean.bed.gz as a default but add an option for
## a different stem
## TODO, create an option for no file.name
GetSamples <- function(file.name, directory.name){
	target <- paste(directory.name, file.name, sep="/")
	system(paste("cd ", directory.name,"; ls *pileup.clean.bed.gz > ", target, "; cd -",sep=""))
	pileups.raw <- as.character(scan(target, what = character(), strip.white=TRUE))
	pileups <- paste(rep(directory.name, length(pileups.raw)), pileups.raw, sep="/")
	names <- gsub(".pileup.clean.bed.gz", "", pileups.raw)
	return(list(pileups=pileups, names=names))
}

##################################################################  
## takes a list of inference results and recursively merges them  
## across common variables
##################################################################  
MergeList <-  function(index, list){
	stopifnot(index >= 2)
	stopifnot(is.list(list))
	l.len <- length(list)
	stopifnot(l.len >= 2)
	cat(paste("Merging sample ", l.len, " \n",sep=""))
	if(l.len==2){
		return(merge(list[[2]][, c(1, index)], list[[1]][, c(1, index)], by = "annotations.rsID", all=TRUE))
	}else{
		temp <- list[[l.len]]
		list[[l.len]] <- NULL
		merge(temp[, c(1, index)], Recall(index, list), by = "annotations.rsID", all=TRUE)
	}
}

##################################################################  
## given a cell.line and plate, identify barcodes in a covariate
## table
################################################################## 
GetBarcodes <- function(cell.line, plate, cov.table){
	stopifnot(cell.line %in% levels(cov.table$CellLine))
	stopifnot(plate %in% levels(cov.table$Plate.ID))
	
	cov.table[which((cov.table$Plate.ID == plate) & (cov.table$CellLine == cell.line)), "Barcode.ID"]
}


##################################################################  
## select samples using target barcodes
################################################################## 
SelectByBarcode <- function(target.barcodes, sample.list){
	stopifnot(is.numeric(target.barcodes))
	#target.barcodes <- barcodes
	#sample.list <- all.samples
	
	## get samples names and extract the barcode number
	names <- sample.list$names
	start.position <- sapply(seq_along(1:length(names)), function(ii){gregexpr("H", names[ii])[[1]][1]})
	end.position <- sapply(seq_along(1:length(names)), function(ii){gregexpr("_", names[ii])[[1]][1]})
	barcode.numbs <- as.numeric(substr(names, start.position+2, end.position-1))
	
	barcode.positions <- which(barcode.numbs %in% target.barcodes)
	
	sample.list$pileups[barcode.positions]
}

##################################################################  
## prepare the ASE data for genotyping
################################################################## 
PrepForGenotyping <- function(ase.dat, min.coverage, dampen.priors=TRUE){
	ref.all <- ase.dat[[1]]
	alt.all <- ase.dat[[2]]
	phi.all <- ase.dat[[4]]$af 
	n.samples <-dim(ref.all)[2]
	
	##################################################################  
	## reads.floor ~ minumum coverage across all samples
	## collapsed.indicator ~ loci with minimum coverage
	## dat.collapssed.final ~ dat.collapsed with >reads.floor minimum coverage
	## ref.all.final ~ sample wise reference counts for loci with sufficient covergae
	## alt.all.final ~ sample wise alternate counts for loci with sufficient covergae
	## annotations.included.samples ~ annotations for loci with sufficient coverage
	## gmat ~ genotype priors from 1k genomes for loci with combined coverage >15
	## gmat.collapsed ~ genotype priors from 1k genomes for all loci
	reads.floor <- min.coverage
	collapsed.indicator <- (rowSums(ref.all + alt.all) > reads.floor)
	ref.all.final <- ref.all[collapsed.indicator, ]		
	alt.all.final <- alt.all[collapsed.indicator, ]
	phi.all.final <- phi.all[collapsed.indicator]
	annotations.included.samples <- ase.dat$anno[collapsed.indicator, ]
	gmat <- cbind(g0=(1-phi.all.final)^2, g1=2*phi.all.final*(1-phi.all.final), g2=phi.all.final^2)
	
	if(dampen.priors){
		gmat <- (gmat+0.0001)/1.0003
	}
	
	list(ref=ref.all.final, alt=alt.all.final, gmat=gmat, annotations=annotations.included.samples)
}
	
