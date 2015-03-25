###########################################################################################
## System diagnostics
## CTH 4.8.2014
###########################################################################################
## all files made with commands such as
## ls *quality.bam | while read f; do echo $f `samtools view -c $f`; done > counts_quality.txt
require('stringr')
library(stringr)
require('ggplot2')
library(ggplot2)

cargs <- commandArgs(trail=TRUE);
if(length(cargs)>=1)
  platePrefix <- cargs[1]
if(length(cargs)>=2)
  bamDir <- cargs[2]

## for debugging
#platePrefix <- "D1P4"
bamDir <- "../../bams"

## Extract the covariates file to match cell-lines and treaments with barcodes.
cov.file <- paste('../../../covariates/GxE_', platePrefix, '_covariates.txt', sep='')
covariates <- read.table(cov.file, header=TRUE, sep='\t')
cov_short <- covariates[, c('Plate.ID', 'Barcode.ID', 'CellLine', 'Treatment.ID')]

## Function takes a bam directory, a bam prefix(which filtered bams), and a plate prefix to return a matrix of counts
CountBams <- function(bamDir, bamPrefix, platePrefix){
    inCommand <- paste("ls ", bamDir, "/*", bamPrefix, ".bam", sep='')
    outFile <- paste(platePrefix, "_", bamPrefix, "_reads.txt", sep='')
    outCommand <- paste(inCommand, " | while read f; do echo $f `samtools view -c $f`; done", " > ", outFile, sep='')
    cat(outCommand)
    system(outCommand)
    return(outFile)
}

## Bams filters: Raw reads merged across lanes, filtered for -q10, and filtered -q10 with duplicates removed
filter_options <- c('merged', 'quality', 'clean')

## apply the CountBams function across all filter levels
counts <- lapply(filter_options, function(filter_level){
#        filter_level <- filter_options[1]
        outFile <- CountBams(bamDir, filter_level, platePrefix)

        temp <- read.table(paste('./', outFile, sep=''), as.is=TRUE)
        ## nice reg exp example
        Barcode.ID <- as.numeric(sub('_.*','',sub(".*HT","",temp$V1)))
        counts <- temp[, 2]

        data.frame(Barcode.ID, counts, stringsAsFactors=FALSE)
        })

## Reduce the output from the 
## barcode X filterLevel table
counts_merged <- Reduce(function(...){merge(..., by='Barcode.ID')}, counts)
#counts_merged <- read.table('DP2_QC_counts.txt', header=TRUE, sep='\t') ## debugging
#counts_merged <- counts_merged[, -c(2,3)] ## debugging
names(counts_merged) <- c('Barcode.ID', filter_options) 
counts_final <- merge(cov_short, counts_merged, by='Barcode.ID')
counts_final <- counts_final[, c('Plate.ID', 'Barcode.ID', 'CellLine', 'Treatment.ID', 'merged', 'quality', 'clean')]

write.table(counts_final, file=paste(platePrefix, '_QC_counts.txt', sep=''), quote=FALSE, row.names=FALSE, sep='\t')

## prepare plot data
counts <- c(counts_merged$merged, counts_merged$quality, counts_merged$clean)
barcodes <- rep(counts_merged$Barcode.ID, 3)
n.barcodes <- length(counts_merged$Barcode.ID)
filter <- c(rep('initial', n.barcodes), rep('quality', n.barcodes), rep('clean', n.barcodes))
dat <- data.frame(barcodes=as.factor(barcodes), filter=as.factor(filter), counts)
dat$filter <- factor(dat$filter, levels=c('initial', 'quality', 'clean'))

## plot
plotName <- paste(platePrefix, '_QC_counts.pdf', sep='')
pdf(plotName, width=17);
my.plot <- ggplot(dat, aes(x=barcodes, y=counts, fill=filter))
my.plot + geom_bar(position='dodge', stat='identity') + ggtitle(paste(platePrefix, ' QC counts',sep='')) + theme_bw()
dev.off()

##
# THE END
##
