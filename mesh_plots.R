##################################################################
## Plotting for QuASAR -> Mesh Output
##
## Created 05.03.2014 
##
## Author: CTH
##
##################################################################
library(ggplot2)

##x11(display="localhost:10.0" ,type="Xlib")

aux <- read.table("Mesh_Complete_results.txt", header=TRUE)

PlotMyBetas <- function(some.treat, some.snp, vertical=TRUE){

	
	temp.dat <- subset(aux, (Treatment==some.treat & annotations.rsID==some.snp))

	temp.clean <- temp.dat[, -c(4:6)]

	temp.treat <- temp.clean[, 1:5]
	temp.treat$condition <- some.treat
	colnames(temp.treat) <- c('Individual', 'Treatment', 'annotations.rsID', 'Beta', 'SE', 'Condition')

	temp.control <- temp.clean[, c(1:3, 6:7)]
	temp.control$condition <- 'Control'
	colnames(temp.control) <- c('Individual', 'Treatment', 'annotations.rsID', 'Beta', 'SE', 'Condition')

        plate='P1'
        some.snp <- "rs1059307"
        some.treat <- "Retinoic Acid"
        indv <- "19239"
        
        temp.control <- data.frame(Individual=indv, Treatment=some.treat, annotations.rsID=some.snp, Beta=-0.121324197616614, SE=0.284806652757819, Condition="Control", stringsAsFactors=FALSE)
        temp.treatment <- data.frame(Individual=indv, Treatment=some.treat, annotations.rsID=some.snp, Beta=-1.11861874100994, SE=0.303580696300933, Condition=some.treat, stringsAsFactors=FALSE)



        
	plot.dat <- rbind(temp.control, temp.treatment)
	plot.dat$Individual <- as.factor(plot.dat$Individual) 
	plot.dat$lb <- plot.dat$Beta - 2*plot.dat$SE
	plot.dat$ub <- plot.dat$Beta + 2*plot.dat$SE

	title <- paste(some.snp, ': ', some.treat, ' specific ASE', sep='')

	p <- ggplot(plot.dat, aes(x=Individual, y=Beta, colour=Condition))
	pd <- position_dodge(width=0.2,height=NULL)
pdf(paste(plate, indv, some.treat, some.snp, '.pdf', sep='_'))
	if(FALSE){
		p + geom_point(aes(shape=Condition), size=4, position=pd)  + scale_color_manual(name="Condition",values=c("lightgreen","steelblue")) + scale_shape_manual(name="Condition",values=c(17,19)) + theme_bw() + geom_errorbar(aes(ymin=lb,ymax=ub), width=0.1, position=pd) +geom_abline(intercept=0, slope=0, colour='red',linetype=4) + scale_y_continuous("Heterozygote Beta Estimate") + ggtitle(title)
	}else{	
		p + geom_point(aes(shape=Condition), size=4, position=pd)  + scale_color_manual(name="Condition",values=c("lightgreen","steelblue")) + scale_shape_manual(name="Condition",values=c(17,19)) + theme_bw() + geom_errorbar(aes(ymin=lb,ymax=ub), width=0.1, position=pd) +geom_abline(intercept=0, slope=0, colour='red',linetype=4) + scale_y_continuous("Heterozygote Beta Estimate") + ggtitle(title) + coord_flip()
	}
dev.off()
}
##some.treat <- 'Copper'
##some.snp <- 'rs2303883'

PlotMyBetas('Copper', 'rs2303883', vertical=FALSE)


