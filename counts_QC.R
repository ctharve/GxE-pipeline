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


#platePrefix <- "D1P1"
#bamDir <- "../../bams"

CountBams <- function(bamDir, bamPrefix, platePrefix){
    inCommand <- paste("ls ", bamDir, "/*", bamPrefix, ".bam", sep='')
    outFile <- paste(platePrefix, "_", bamPrefix, "_reads.txt", sep='')
    outCommand <- paste(inCommand, " | while read f; do echo $f `samtools view -c $f`; done", " > ", outFile, sep='')
    cat(outCommand)
    system(outCommand)

    return(outFile)
}

filter_options <- c('merged', 'quality', 'clean')

counts <- lapply(filter_options, function(filter_level){
#        filter_level <- filter_options[1]
        outFile <- CountBams(bamDir, filter_level, platePrefix)

        temp <- read.table(paste('./', outFile, sep=''), as.is=TRUE)
        ## nice reg exp example
        barcodes <- as.numeric(sub('_.*','',sub(".*HT","",temp$V1)))
        counts <- temp[, 2]

        data.frame(barcodes, counts, stringsAsFactors=FALSE)
        })

## barcode X filterLevel table
counts_merged <- Reduce(function(...){merge(..., by='barcodes')}, counts)
names(counts_merged) <- c('barcodes', filter_options)
write.table(counts_merged, file=paste(platePrefix, '_QC_counts.txt', sep=''), quote=FALSE, row.names=FALSE, sep='\t')

## prepare plot data
counts <- c(counts_merged$merged, counts_merged$quality, counts_merged$clean)
barcodes <- rep(counts_merged$barcodes, 3)
n.barcodes <- length(counts_merged$barcodes)
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
