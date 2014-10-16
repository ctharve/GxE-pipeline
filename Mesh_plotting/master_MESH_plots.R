###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
library('ggplot2')
library('reshape')
x11(display="localhost:12.0" ,type="Xlib")
###########################################################################
###########################################################################
source('~/piquelab/charvey/GxE/jointGenotyping/scripts/mesh_analysis_functions.R')
####
## 0.) capture and append count data
####
dat <- read.table('test.txt',  header=FALSE, stringsAsFactors=FALSE)
names(dat) <- c('direction', 'bin', 'pi0', 'pi0.l', 'pi0.u', 'treatment', 'treatment.l', 'treatment.u', 'control', 'control.l', 'control.u', 'both', 'both.l', 'both.u')
head(dat)
countFiles <- scan(pipe("ls allPlates_*.txt"), character(0))
count_dat <- Reduce(rbind, lapply(countFiles, function(file){
  ## file <- countFiles[1]
  target <- paste0('less ', file, ' | wc -l | awk \'{print $1/2}\'')
  myDat <- as.numeric(scan(pipe(target), character(0)))
  myLabels <- strsplit(gsub('.txt', '',  gsub('allPlates_', '', file)), '_')[[1]]
  out_dat <- data.frame(direction=myLabels[1], bin=as.numeric(myLabels[2]), snp.count=myDat, stringsAsFactors=FALSE)
  out_dat
}))
dat <- merge(dat, count_dat, by=c('direction', 'bin'))

####
## 1.) plots of all posterior parameters from MESH
## plot all posteriors probabilities 
#### 
dat_long <- melt(dat[, c('direction', 'bin', 'pi0', 'treatment', 'control', 'both')], id=c("direction", "bin"))
dat_long$direction <- as.factor(dat_long$direction)
pdf('./plots/Mesh_posteriors_allParameters.pdf')
my_plot <- ggplot(dat_long, aes(x=bin, y=value))
my_plot + geom_point() + facet_grid(direction ~ variable)
dev.off()

####
## 2.) plots of treatment and condition only ASE number of analyzed SNPs as point size
## treatment only ASE with # SNPs as point size
####
plot_dat <- dat[which(dat$direction=="gtpos"), c('direction', 'bin', 'treatment', 'treatment.u', 'treatment.l', 'snp.count')]
p <- ggplot(plot_dat, aes(x=bin, y=treatment)) + theme_bw() + geom_errorbar(aes(ymin=treatment.l,ymax=treatment.u), width=0.02)
pdf('./plots/Mesh_posteriors_treatmentOnly_countSize.pdf')
p + geom_point(aes(size=log10(snp.count), alpha=.75), colour="steelblue") + ggtitle("Posterior Proportion of SNPs With Treatment Only ASE") + xlab("Fold change bin") + ylab("Posterior distribution") + theme_bw() + theme(legend.position="none") + scale_x_continuous(limits=c(1.0, 3.5))
dev.off()
## control only ASE with # SNPs as point size 
plot_dat <- dat[which(dat$direction=="gtpos"), c('direction', 'bin', 'control', 'control.u', 'control.l', 'snp.count')]
p <- ggplot(plot_dat, aes(x=bin, y=control)) + theme_bw() + geom_errorbar(aes(ymin=control.l,ymax=control.u), width=0.02)
pdf('./plots/Mesh_posteriors_controlOnly_countSize.pdf')
p + geom_point(aes(size=log10(snp.count), alpha=.75), colour="steelblue") + ggtitle("Posterior Proportion of SNPs With Control Only ASE") + xlab("Fold change bin") + ylab("Posterior distribution") + theme_bw() + theme(legend.position="none") + scale_x_continuous(limits=c(1.0, 3.5))
dev.off()

####
## 3.) plots of treatment and condition onle ASE number of analyzed SNPs as point size
## final plot: treatment only ASE
####
plot_dat <- dat[which(dat$direction=="gtpos"), c('direction', 'bin', 'treatment', 'treatment.u', 'treatment.l')]
p <- ggplot(plot_dat, aes(x=bin, y=treatment)) + theme_bw() + geom_errorbar(aes(ymin=treatment.l,ymax=treatment.u), width=0.05)
pdf('./plots/Mesh_posteriors_treatmentOnly.pdf')
p + geom_point(size=3.2, colour="steelblue") + ggtitle("Posterior Proportion of SNPs With Treatment Only ASE") + xlab("Fold change bin") + ylab("Posterior distribution") + theme_bw() + theme(legend.position="none") + scale_x_continuous(limits=c(1.0, 3.5))
dev.off()
## final plot: control only ASE
plot_dat <- dat[which(dat$direction=="gtpos"), c('direction', 'bin', 'control', 'control.u', 'control.l')]
p <- ggplot(plot_dat, aes(x=bin, y=control)) + theme_bw() + geom_errorbar(aes(ymin=control.l,ymax=control.u), width=0.05)
pdf('./plots/Mesh_posteriors_controlOnly.pdf')
p + geom_point(size=3.2, colour="steelblue") + ggtitle("Posterior Proportion of SNPs With Control Only ASE") + xlab("Fold change bin") + ylab("Posterior distribution") + theme_bw() + theme(legend.position="none") + scale_x_continuous(limits=c(1.0, 3.5))
dev.off()

####
## 4.) some violin plots of the density of treatment on ly ASE configurations
## files with draws from the posterior distribution 
####
myFiles <- scan(pipe("ls allPlates_*.hm_config.out"), character(0))
## gene expression direction and bin
plot_dat_treat <- Reduce(rbind, lapply(myFiles, function(file){
  ## file <- myFiles[1]
  target <- paste0('less ', file, ' | tail -10000 | awk \'{print $4}\'')
  myDat <- as.numeric(scan(pipe(target), character(0)))
  n.samples <- length(myDat)
  myLabels <- strsplit(gsub('.hm_config.out', '',  gsub('allPlates_', '', file)), '_')[[1]]
  out_dat <- data.frame(direction=rep(myLabels[1], n.samples), bin=rep(myLabels[2], n.samples), samples=myDat, stringsAsFactors=FALSE)
  out_dat$direction <- as.factor(out_dat$direction)
  out_dat$bin <- as.factor(out_dat$bin)
  out_dat
}))
plot_dat_treat$treat <- 'treatment'
plot_dat_treat$direction <- as.character(plot_dat_treat$direction) 
plot_dat_treat$bin <- as.character(plot_dat_treat$bin) 

plot_dat_control <- Reduce(rbind, lapply(myFiles, function(file){
  ## file <- myFiles[1]
  target <- paste0('less ', file, ' | tail -10000 | awk \'{print $3}\'')
  myDat <- as.numeric(scan(pipe(target), character(0)))
  n.samples <- length(myDat)
  myLabels <- strsplit(gsub('.hm_config.out', '',  gsub('allPlates_', '', file)), '_')[[1]]
  out_dat <- data.frame(direction=rep(myLabels[1], n.samples), bin=rep(myLabels[2], n.samples), samples=myDat, stringsAsFactors=FALSE)
  out_dat$direction <- as.factor(out_dat$direction)
  out_dat$bin <- as.factor(out_dat$bin)
  out_dat
}))
plot_dat_control$treat <- 'control'
plot_dat_control$direction <- as.character(plot_dat_control$direction) 
plot_dat_control$bin <- as.character(plot_dat_control$bin) 

plot_dat <- rbind(plot_dat_treat, plot_dat_control)
plot_dat <- plot_dat[which(plot_dat$direction=='gtpos'), ]

pdf('./plots/violin_plots_posterior_density_treatmentControl.pdf', height=5, width=7)
p <- ggplot(plot_dat, aes(x=bin, y=samples))
p + geom_violin(scale='width', aes(fill=treat)) + facet_grid(direction ~ .) + ggtitle("Posterior Distribution of Treatment Only ASE") + xlab("Fold change bin")
dev.off()

####
## 5.) plots of QuASAR beta estimtates from analysis of "significant" MESH BFs
##
####
## data for configurations that are Treatment only
betaTreat_dat <- read.table(pipe("less ./posteriors_bfs/allPlates_*bayesFactors.txt | awk \'($3>$4*2) && ($3>20) && ($2<0.5){print $0\"\t\"$3/$4}\' | sort -k6 -n | cut -f1 | tr \"_\" \"\t\" | cut -f1,2,3 | tr \"\t\" \"_\"| while read f; do grep $f allPlates*.txt | tr \" \" \"\t\"; done"), stringsAsFactors=FALSE)
names(betaTreat_dat) <- c('id', 'treat', 'beta', 'beta.se')
## parse ids
betaTreat_info <- get_info(betaTreat_dat$id)
betaTreat_info$treatment <- NA  ## this bad: for each SNP ensure the first row is the control and the second treatment
betaTreat_info$treatment[(as.numeric(row.names(betaTreat_dat)) %% 2) == 1] <- 'control'
betaTreat_info$treatment[(as.numeric(row.names(betaTreat_dat)) %% 2) == 0] <- 'treatment'

plot_dat <- cbind(betaTreat_info, betaTreat_dat[, c('beta', 'beta.se')]/log(2))
plot_dat$lb <- plot_dat$beta - 1.96*plot_dat$beta.se
plot_dat$ub <- plot_dat$beta + 1.96*plot_dat$beta.se
plot_dat[c(47, 48), 'id'] <- 'rs1059307_RetAc_II'
str(plot_dat)
plot_dat$direction <- as.factor(plot_dat$direction) 
plot_dat$bin <- as.factor(plot_dat$bin)
plot_dat$id <- as.factor(plot_dat$id)
plot_dat$treatment <- as.factor(plot_dat$treatment)
plot_dat <- plot_dat[which(plot_dat$direction=="gtpos"), ]

pd <- position_dodge(width=0.5, height=NULL)
p <- ggplot(plot_dat, aes(x=id, y=beta, colour=treatment)) + geom_abline(intercept=0, slope=0, colour='red',linetype=4) 
pdf('./plots/Mesh_allBetas_treatment_config.pdf')
p + geom_point(position=pd, size=3.2) + geom_errorbar(aes(ymin=lb , ymax=ub), width=0.25, position=pd) + theme_bw() + coord_flip() + ggtitle('Parameter Estimates of Treatment Only ASE Configurations')
dev.off()

## data for configurations that are Control only
betaControl_dat <- read.table(pipe("less ./posteriors_bfs/allPlates_*bayesFactors.txt | awk \'($2>$4*2) && ($2>20) && ($3<0.5){print $0\"\t\"$2/$4}\' | sort -k6 -n | cut -f1 | tr \"_\" \"\t\" | cut -f1,2,3 | tr \"\t\" \"_\"| while read f; do grep $f allPlates*.txt | tr \" \" \"\t\"; done"), stringsAsFactors=FALSE)
names(betaControl_dat) <- c('id', 'treat', 'beta', 'beta.se')
## parse ids
betaControl_info <- get_info(betaControl_dat$id)
betaControl_info$treatment <- NA  ## this bad: for each SNP ensure the first row is the control and the second treatment
betaControl_info$treatment[(as.numeric(row.names(betaControl_dat)) %% 2) == 1] <- 'control'
betaControl_info$treatment[(as.numeric(row.names(betaControl_dat)) %% 2) == 0] <- 'treatment'

plot_dat <- cbind(betaControl_info, betaControl_dat[, c('beta', 'beta.se')]/log(2))
plot_dat$lb <- plot_dat$beta - 1.96*plot_dat$beta.se
plot_dat$ub <- plot_dat$beta + 1.96*plot_dat$beta.se
str(plot_dat)
plot_dat$direction <- as.factor(plot_dat$direction) 
plot_dat$bin <- as.factor(plot_dat$bin)
plot_dat$id <- as.factor(plot_dat$id)
plot_dat$treatment <- as.factor(plot_dat$treatment)
plot_dat <- plot_dat[which(plot_dat$direction=="gtpos"), ]

pd <- position_dodge(width=0.5, height=NULL)
p <- ggplot(plot_dat, aes(x=id, y=beta, colour=treatment)) + geom_abline(intercept=0, slope=0, colour='red',linetype=4) 
pdf('./plots/Mesh_allBetas_control_config.pdf')
p + geom_point(position=pd, size=3.2) + geom_errorbar(aes(ymin=lb , ymax=ub), width=0.25, position=pd) + theme_bw() + coord_flip() + ggtitle('Parameter Estimates of Control Only ASE Configurations')
dev.off()


###########################################################################
###########################################################################
## dot plot for treatment only ASE
dat_long <- melt(subset(dat[, c('direction', 'bin', 'treatment', 'treatment.u', 'treatment.l')], direction=='gtpos'), id=c("direction", "bin"))
dat_long$variable <- as.character(dat_long$variable)
dat_long[which(dat_long$variable %in% c("treatment.l", "treatment.u")), "variable"] <- "Credible Interval"
dat_long[which(dat_long$variable %in% c("treatment")), "variable"] <- "Posterior Mean"
dat_long$direction <- as.factor(dat_long$direction)
dat_long$variable <- as.factor(dat_long$variable)
head(dat_long)
my_plot <- ggplot(dat_long, aes(x=bin, y=value, colour=variable))
#pdf('mesh_posteriors_treatmentOnly.pdf')
my_plot + geom_point(aes(size=2.5)) + ggtitle("Posterior Proportion of SNPs With Treatment Only ASE") + xlab("gene expression fold change") + ylab("posterior mean") + theme_bw()
+ theme(legend.position="none")
#dev.off()

## dot plot for control only ASE
dat_long <- melt(subset(dat[, c('direction', 'bin', 'control', 'control.u', 'control.l')], direction=='gtpos'), id=c("direction", "bin"))
dat_long$variable <- as.character(dat_long$variable)
dat_long[which(dat_long$variable %in% c("control.l", "control.u")), "variable"] <- "Credible Interval"
dat_long[which(dat_long$variable %in% c("control")), "variable"] <- "Posterior Mean"
dat_long$direction <- as.factor(dat_long$direction)
dat_long$variable <- as.factor(dat_long$variable)
head(dat_long)
my_plot <- ggplot(dat_long, aes(x=bin, y=value, colour=variable))
#pdf('mesh_posteriors_controlOnly.pdf')
my_plot + geom_point(aes(size=2.5)) + ggtitle("Posterior Porportion of SNPs With Control Only ASE") + xlab("logFC bin") + ylab("posterior mean") + theme_bw() + theme(legend.position="none")
#dev.off()
