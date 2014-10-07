##################################################################
## Software suite for joint genotyping and ASE inference on multiple
## experiments and samples
##
## Created 10.07.2014 
## select controls
##
## CTH 100714
## made changes to include 
##
## derived from aseSuite_v0.0_D1P6_meshprep2.R
## Author: CTH
##
## Version 0.0: Preliminary
## Arguments: plate, cell.line, covariate.file, pileup.dir, out.dir
## Return Values: Genotypes, model.convergence, inference, metaData
##################################################################
#library('ggplot2', lib.loc=myRlib)
## x11(display="localhost:12.0" ,type="Xlib")
##################################################################    
cargs <- commandArgs(trail=TRUE);
if(length(cargs)>=1)
  plate <- cargs[1]
if(length(cargs)>=2)
  cell.line <- cargs[2]
## plate <- "DP9"
## cell.line <- "SM046-A"

## extract covariates table
cov.name <- paste('~/piquelab/scratch/charvey/GxE/derived_data/covariates/GxE_', plate, '_covariates.txt', sep='')
cv <- read.table(file=cov.name, sep="\t", header=TRUE, stringsAsFactors=FALSE)
cv <- cv[cv$Plate.ID==plate & cv$CellLine==cell.line, ]
ids <- unique(cv$Treatment.ID)[(unique(cv$Treatment.ID) %in% c('CO1', 'CO2', 'CO3'))] ## control only ids

out_dat <- Reduce(rbind, lapply(ids, function(id){
  ## id <- ids[1]
  ## cat('This id: ', id, '\n')
  aux <- read.table(paste0('./output/', paste(plate, cell.line, id, 'allOutput.txt', sep='_')), stringsAsFactors=FALSE, header=TRUE)
  aux$plate <- plate
  aux$cell.line <- cell.line
  aux$treatment <- id
  #head(aux)
  aux <- aux[,c('plate', 'cell.line', 'treatment', 'chr', 'pos0', 'pos', 'rsID', 'beta', 'beta.se', 'pval', 'qval')]
  aux
}))

fileName <- paste0('./output/', plate, '_', cell.line, '_masterControl.txt')
write.table(out_dat, fileName, col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')
                                
##
cat("###### THE END ######")
##
