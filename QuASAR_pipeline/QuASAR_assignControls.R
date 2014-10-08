##################################################################
## Software suite for joint genotyping and ASE inference on multiple
## experiments and samples
##
## Created 07.14.2014 
## changes from *meshprep: all controls are co
##
## derived from aseSuite_v0.0_D1P6_meshprep2.R
## Author: CTH
##
## Version 0.0: Preliminary
## Arguments: plate, cell.line, covariate.file, pileup.dir, out.dir
## Return Values: Genotypes, model.convergence, inference, metaData
##################################################################
#library('ggplot2', lib.loc=myRlib)
library('ggplot2')
##x11(display="localhost:12.0" ,type="Xlib")
top_dir <- '/wsu/home/groups/piquelab/charvey/GxE/jointGenotyping/'
##################################################################
## use the covariate table and use the barcodes to select samples
##################################################################    
cargs <- commandArgs(trail=TRUE);
if(length(cargs)>=1)
  plate <- cargs[1]
if(length(cargs)>=2)
  cell.line <- cargs[2]
## plate <- "DP1" 
## cell.line <- "18507"

##################################################################
## specify covariate file and pileup directories
##################################################################    

cov_data <- function(plate, cell.line, treat, fields='all'){
  cov.file <- paste('~/piquelab/scratch/charvey/GxE/derived_data/covariates/GxE_', plate, '_covariates.txt', sep='')
  cv <- read.table(file=cov.file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  cv$Treatment <- gsub(' ',  '_', cv$Treatment)
  if(fields=='all'){
    cv[cv$Plate.ID==plate & cv$CellLine==cell.line & cv$Treatment.ID==treat, ]
  } else {
    cv[cv$Plate.ID==plate & cv$CellLine==cell.line & cv$Treatment.ID==treat, fields]
  }
}

controls <- read.table('./data/MASTER_masterControl.txt', stringsAsFactors=FALSE, header=FALSE)
controls <- controls[complete.cases(controls), ]
names(controls) <- c('plate', 'cell.line', 'treatment', 'chr', 'pos0', 'pos', 'rsID', 'beta', 'beta.se', 'pval', 'qval')
treats <- read.table('./data/MASTER_masterTable.txt', stringsAsFactors=FALSE, header=FALSE)
treats <- treats[complete.cases(treats), ]
names(treats) <- c('plate', 'cell.line', 'treatment', 'chr', 'pos0', 'pos', 'rsID', 't.id', 'g.id', 'ensg', 'beta', 'beta.se', 'pval', 'qval', 'logFC')

ub <- 0.7
lb <- -ub
treat_ub <- treats[treats$logFC>ub,]
treat_lb <- treats[treats$logFC<lb, ]

my_data <- list()
my_data[[1]] <- treat_ub
my_data[[2]] <- treat_lb
names(my_data) <- c('treat_ub', 'treat_lb')


## get unique plates, cell.lines, treatments
## get the control data
## merge on rsID
## reformat for MESH

## for each partitioned data set
out_dat <- lapply(my_data, function(data){  
  #data <- my_data[[1]]
  
  plates <- unique(data$plate)

  ## for each plate
  aux3 <- Reduce(rbind, lapply(plates, function(plate){  
    #plate <- plates[1]

    lines <- unique(data[which(data$plate==plate), 'cell.line'])

    ## for each cell line
    aux2 <- Reduce(rbind, lapply(lines, function(line){
      #line <- lines[1]
      
      treats <- unique(data[which(data$plate==plate & data$cell.line==line), 'treatment'])

      ## for each treatment
      aux1 <- Reduce(rbind, lapply(treats, function(treat){
        #treat <- treats[1]

        # get the control ID
        cont <- cov_data(plate=plate, cell.line=line, treat=treat, fields='Control.ID')
        stopifnot(length(cont)==1)

        # get the treatment data
        treat_dat <- data[which(data$plate==plate & data$cell.line==line & data$treatment==treat), c('rsID', 'treatment', 'beta', 'beta.se')]
        #str(treat_dat)
        
        # get the control data
        target <- paste0(top_dir, 'QuASAR_results_', plate, '/output/', plate, '_', line, '_cov15_', 'inference.RData')
        load(target)
        cont_dat <- inference.data[[which(names(inference.data)==cont)]][[1]][, c('annotations.rsID', 'betas', 'betas.se')]
        names(cont_dat) <- c('rsID', 'beta', 'beta.se')
        cont_dat$rsID <- as.character(cont_dat$rsID) 
        cont_dat$treatment <- cont
        cont_dat <- cont_dat[, c('rsID', 'treatment', 'beta', 'beta.se')]
        #str(cont_dat)
 
        ## merge 
        totaldata <- merge(cont_dat, treat_dat, by='rsID')

        ## format for MESH
        finaldata <- Reduce(rbind, lapply(seq_along(1:dim(totaldata)[1]), FUN=function(ii){
          #ii <- 1
	  control <- data.frame(rsID=totaldata[ii, 'rsID'], treat=totaldata[ii, 'treatment.x'], beta=totaldata[ii, 'beta.x'], se=totaldata[ii, 'beta.se.x'])
	  treatment <- data.frame(rsID=totaldata[ii, 'rsID'], treat=totaldata[ii, 'treatment.y'], beta=totaldata[ii, 'beta.y'], se=totaldata[ii, 'beta.se.y'])

	  returnval <- rbind(control, treatment)
	  rownames(returnval) <- NULL
	  returnval
        }))

        finaldata                          
        
      }))

      aux1                         
                               
    }))

    aux2              
                  
  }))

  aux3
  
}); names(out_dat) <- names(my_data)

hold <- out_dat
out_dat$treat_ub$rsID <- as.character(out_dat$treat_ub$rsID)
out_dat$treat_ub$treat <- as.character(out_dat$treat_ub$treat)
out_dat$treat_lb$rsID <- as.character(out_dat$treat_lb$rsID)
out_dat$treat_lb$treat <- as.character(out_dat$treat_lb$treat)

write.table(out_dat$treat_ub, file=paste0('./mesh/allData_gtpos.7.txt'), row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(out_dat$treat_lb, file=paste0('./mesh/allData_ltneg.7.txt'), row.names=FALSE, col.names=FALSE, quote=FALSE)

##
cat("###### THE END ######")
##
