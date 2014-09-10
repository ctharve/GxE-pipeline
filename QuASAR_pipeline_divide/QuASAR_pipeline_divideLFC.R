##################################################################
## Software suite for joint genotyping and ASE inference on multiple
## experiments and samples
##
## Created 09.10.2014 
## select samples with logFC cutoff
##
## CTH 090214
## made changes to include 
##
## derived from aseSuite_v0.0_D1P6_meshprep2.R
## Author: CTH
##
## Version 0.0: Preliminary
## Arguments: plate, cell.line, covariate.file, pileup.dir, out.dir
## Return Values: Genotypes, model.convergence, inference, metaData
##################################################################
## automate this in a Python or shell script
LPG <- '/wsu/home/groups/piquelab'
## make sure to include a directory to download and source all external packages
## in stampede
#myRlib <- paste(LPG, '/charvey/tools/Rlib', sep='')
## helper functions for data processing
source(paste(LPG, '/charvey/GxE/jointGenotyping/scripts/aseSuite_functions_v0.0.R', sep=''))
## ASE model fitting functions funtions 
source(paste(LPG, '/charvey/source/ASE/fitAseModels.v4.R', sep=''))
## qqplot functions
source(paste(LPG, '/gmb/AI/results/qqman.r', sep=''))

##################################################################
#library('ggplot2', lib.loc=myRlib)
##x11(display="localhost:11.0" ,type="Xlib")

##################################################################
## use the covariate table and use the barcodes to select samples
##################################################################    
min.cov <- 15

cargs <- commandArgs(trail=TRUE);
if(length(cargs)>=1)
  plate <- cargs[1]
if(length(cargs)>=2)
  cell.line <- cargs[2]
if(length(cargs)>=3)
  cutOff <- cargs[3]

#plate <- "DP1" 
## Plate DP6 cell lines KP39334 KP39346 KP39351
#cell.line <- "18507"
#cutOff <- 0.5


output.folder <- paste('./output',sep='')
## extract covariates table
cov.file <- paste('~/piquelab/scratch/charvey/GxE/derived_data/covariates/GxE_', plate, '_covariates.txt', sep='')
cv <- read.table(file=cov.file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
cv$Treatment <- gsub(' ',  '_', cv$Treatment)
cov.file <- cv[cv$Plate.ID==plate & cv$CellLine==cell.line, ]

dat.name <- paste(output.folder, "/",plate, "_",cell.line, "_cov", min.cov, '_inference.RData', sep='')
load(dat.name)

## REQUIRES: load(dat.name)
n.sams <- length(inference.data)
groupings <- c('low', 'high')
ids <- unique(cov.file$Treatment.ID[!(cov.file$Treatment.ID %in% c('CO1', 'CO2', 'CO3'))])

for(grouping in groupings){

  load(dat.name)

  for(id in ids){
    #id <- ids[1]
    cat('This id: ', id, '\n')
    idx <- which(names(inference.data)==id)
    cat('Sample ID: ', names(inference.data)[idx], '\n')
    q_dat <- read.table(paste0('./output/', paste(plate, cell.line, id, 'allOutput.txt.tid', sep='_')), stringsAsFactors=FALSE)
    names(q_dat) <- c('chr', 'pos0', 'pos', 'rsID', 'beta', 'beta.se',
                    'pval', 'qval', 't.id', 'g.id', 'ensg')
    str(q_dat)

    d_dat <- read.table(paste0(LPG, '/charvey/GxE/differential_expression/DEseq2_results/out_data_',
           plate, '/stats/', plate, '_DEG_stats', '_', id, '.txt'), stringsAsFactors=FALSE, header=TRUE)
    rownames(d_dat) <- d_dat$t.id
    str(d_dat)

    dq_dat <- q_dat
    dq_dat$logFC <- d_dat[dq_dat$t.id, 'logFC']
    dq_dat$high <- dq_dat$logFC > cutOff 

    if(grouping=='high'){
      out_dat <- dq_dat[dq_dat$high, ]
    } else {
      out_dat <- dq_dat[!dq_dat$high, ]
    }
    str(out_dat)
  
    ## obtain inference.dat
    temp <- inference.data[[id]][[1]]
    temp$annotations.rsID <- as.character(temp$annotations.rsID)
    str(temp)
    idx_keep <- temp$annotations.rsID %in% out_dat$rsID  
    temp <- temp[idx_keep, ]
    temp$annotations.rsID <- as.factor(temp$annotations.rsID)
    inference.data[[id]][[1]] <- temp
  }

  relevant.samples <- which(names(inference.data) %in% ids)

  for(ii in relevant.samples){
	  #ii <- 3
	  treat <- names(inference.data[ii])
          control <- unlist(subset(cov.file, Treatment.ID==treat, Control.ID))

	  	  treatdata <- cbind(as.character(inference.data[[treat]]$dat$annotations.rsID), treat, inference.data[[treat]]$dat$betas, inference.data[[treat]]$dat$betas.se)
		  controldata <- cbind(as.character(inference.data[[control]]$dat$annotations.rsID), control, inference.data[[control]]$dat$betas, inference.data[[control]]$dat$betas.se)
		  colnames(treatdata) <- c('SNP_ID', 'label', 'beta', 'SE')
		  colnames(controldata) <- c('SNP_ID', 'label', 'beta', 'SE')

		  totaldata <- merge(treatdata, controldata, by='SNP_ID')

		  totaldata$SNP_ID <- as.character(totaldata$SNP_ID)
		  totaldata$label.x <- as.character(totaldata$label.x)
		  totaldata$beta.x <- as.numeric(levels(totaldata$beta.x)[as.integer(totaldata$beta.x)]) 
		  totaldata$SE.x <- as.numeric(levels(totaldata$SE.x)[as.integer(totaldata$SE.x)]) 
		  totaldata$label.y <- as.character(totaldata$label.y)
		  totaldata$beta.y <- as.numeric(levels(totaldata$beta.y)[as.integer(totaldata$beta.y)]) 
		  totaldata$SE.y <- as.numeric(levels(totaldata$SE.y)[as.integer(totaldata$SE.y)]) 

		  finaldata <- lapply(seq_along(1:dim(totaldata)[1]), FUN=function(ii){
			  #ii <- 1
			  control <- c(totaldata[ii, 'SNP_ID'], totaldata[ii, 'label.x'], totaldata[ii, 'beta.x'], totaldata[ii, 'SE.x'])
			  treatment <- c(totaldata[ii, 'SNP_ID'], totaldata[ii, 'label.y'], totaldata[ii, 'beta.y'], totaldata[ii, 'SE.y'])

			  returnval <- rbind(control, treatment)
			  rownames(returnval) <- NULL
			  returnval
			  })

		  for(ii in 1:length(finaldata)){
			  #ii <- 1
			  if(ii==1){
			  	  finallong <- finaldata[[ii]]
			  }else{
				  finallong <- rbind(finallong, finaldata[[ii]])
			  }
		  }
		
		  colnames(finallong) <- c('SNP_ID', 'label', 'beta', 'SE')

                  filename <- paste('./mesh/analysis_data/', plate, '_', cell.line, '_', treat, '.txt', sep='')
		  filename <- paste('./mesh/analysis_data/', plate, '_', cell.line, '_', treat, '_', grouping, '.txt', sep='')
		  write.table(finallong, file=filename, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t") 	
  }
}

##
cat("###### THE END ######")
##
