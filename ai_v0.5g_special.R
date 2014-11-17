##################################################################
## Filters a pilup file from a single sample
## 
## Created: 08/06/2013
##
## Edited: 10/10/13 gmb + rpr: output new parameters (rho, rho vs phi)
##
## Edited: 10/17/13 gmb: added duplicate removal via x8b files
## Edited: 11/05/13 gmb: added outputs during x8b files, changed
##					source for ase functions, modified fitAseNull() params
## Edited: 11/08/13 gmb: added pval3 to the final outputs
## 
## Author: RPR, CH
##
## Version0.3: hacked from 'ai_v3.6_noImpute' to be run interactively
##################################################################
# setwd("/wsu/home/fl/fl97/fl9788/piquelab/charvey/ASE/prod2")


##################################################################
qual <- function(x) { strtoi(charToRaw(x),16L)-33 }
filt <- function(r,q,t) {
	ret=''
	for (i in 1:nchar(r)) {
		if (substr(r,i,i) != '$' & substr(r,i,i) != '^') {
			if (qual(substr(q,i,i)) >= t) { ret <- paste(ret,substr(r,i,i),sep='') }
		}
	}
	ret
}

##################################################################
thresh <- 13  # Base call quality threshold  NOT used
mincov <- 4   # Minimum coverage
maxcov <- 200000 # Max coverage

pileupFile <- "pileup/test5.OpenChromDnaseFibroblgm03348.pileup.bed.gz"
outDir <- "res5/OpenChromDnaseFibroblgm03348/"


##################################################################
cargs<-commandArgs(trail=TRUE);
if(length(cargs)>=1)
  pileupFile <-cargs[1];
if(length(cargs)>=2)
  outDir <-cargs[2];
if(length(cargs)>=3)
  outLociFile <- cargs[3];

##################################################################
## Output to make sure everything is as expected
#cat("#Pileup File:",pileupFile,"\n")
#cat("#Threshold:",thresh,"\n")
###cat("#fileName:",fileName,"\n")
sName <- gsub(".*/","",gsub(".pileup.bed.gz","",pileupFile));
cat("#sName:",sName,"\n")
cat("outDir:",outDir,"\n")

##################################################################
## Create an object to save data in
saveData <- data.frame(pileup.file=c(pileupFile)) #,impute.file=c(imputeFile))

##################################################################
## QC STEPS
# Remove all loci whose read coverage is below 4 (our picked threshold)
# cov = read coverage per allele to define heterozygote
# RPR: We could pipe the data in and apply these filters, to make all this go faster
## pileup <- pileup[pileup$num.reads >= mincov,]
## pileup <- pileup[pileup$num.reads <= maxcov,]

command=paste("less ",pileupFile,
		"| awk ' $5 >=",mincov," && $5 <=",maxcov,"'")
#cat("#pipe: ",command,"\n")

## RPR: pipe speeds up, as.is=T speeds up and avouids having to use character conversion later...
pileup <- read.table(file=pipe(command),header=F,quote="",comment.char="",as.is=T,sep="\t") 
# names(pileup) <- c("chr","pos-1","pos","ref","num.reads","read.alleles","read.quality","xchr","xpos-1","xpos","rsID","TKG.Ref","alt","xcov","xoverlap")
# names(pileup) <- c("chr","pos-1","pos","ref","num.reads","read.alleles","read.quality","xchr","xpos-1","xpos","rsID","TKG.Ref","alt","d1","af","d2")

## Restructured pileup output for res7
names(pileup) <- c("chr","pos-1","pos","ref","num.reads","read.alleles","read.quality","rsID","TKG.Ref","alt","af")



# See if the ref allels match, then discard uncesessary columns
##saveData$ref.match <- length(which(toupper(as.character(pileup$ref)) == as.character(pileup$TKG.Ref))) == dim(pileup)[1]
indMatch <- (toupper(pileup$ref) == pileup$TKG.Ref)
saveData$ref.match <- sum(!indMatch)==0;
saveData$ref.missmatches <- sum(!indMatch);
pileup <- pileup[indMatch,c("chr","pos-1","pos","ref","alt","rsID","num.reads","read.alleles","read.quality","af")]
## I'm also discarding mismatching rows.
stopifnot(mean(indMatch)>0.8) ##RPR stop if too many errors?
rm(indMatch)

# Duplicates arise from the impute data for 3 (so far) determined reasons: tri+ alleleic SNPs,
# indels (should already be filtered), and incongruencies between illumina and 1KG data
d1 <- duplicated(paste(pileup$chr,pileup$pos,sep=":"))
d2 <- duplicated(paste(pileup$chr,pileup$pos,sep=":"),fromLast=T)
saveData$num.duplicates <- length(which(d1)) + length(which(d2))
pileup <- pileup[!(d1 | d2),]
rm(d1,d2)
## RPR: Position did not inlude chr

## Filter the alleles to remove those at the beginning and end of a mapped read
## Unnecessary if data loaded with: as.is=T
## pileup$read.alleles <- as.character(pileup$read.alleles)
## pileup$read.quality <- as.character(pileup$read.quality)
pileup$read.alleles.filt <- mapply(gsub,'[a-zA-Z.,]\\$','$',pileup$read.alleles) ## Remove start base 
pileup$read.alleles.filt <- mapply(gsub,'\\^[[:punct:][:alnum:]][a-zA-Z.,]','^',pileup$read.alleles.filt)  ## Remove end base

##head((strsplit(pileup$read.quality,'')))
## RPR alternative way to extract base quality numbers. 
qual <- sapply(1:nrow(pileup),function(ii){qual(pileup$read.quality[ii])})
qual <- unlist(qual)
qtr <- quantile(qual,seq(0,1,0.1))
qual.table <- table(qual);
qtr
qual.table
##RPR: We could decide the the threshold based on the quantile...?
##RPR: For the moment I remove the function...

##pileup$read.alleles.filt <- mapply(filt,pileup$read.alleles.filt,pileup$read.quality,MoreArgs=list(t=thresh))
pileup <- pileup[nchar(pileup$read.alleles.filt)>0,]
## RPR: we could do some of the filters as a pipe. 

## Should be no more odd chars - no Ns, no +s, -s, ^s, $s, etc (still... assert that there are only ACGTs left?)
##pileup$read.alleles.clean <- gsub('[^AGCTagct.,]',"",pileup$read.alleles)
pileup$read.alleles.clean <- mapply(gsub,'[\\.\\,]',pileup$ref,pileup$read.alleles.filt)
pileup$read.alleles.clean <- toupper(pileup$read.alleles.clean)
pileup$ref <- toupper(pileup$ref) # for downstream analysis

# For each read at each location, calculate the number of ref and alt allele matches (alt defined by 1KG)
pileup$ref.matches <- as.integer(nchar(mapply(gsub,paste('[^',as.character(pileup$ref),']',sep=""),'',pileup$read.alleles.clean)))
pileup$alt.matches <- as.integer(nchar(mapply(gsub,paste('[^',as.character(pileup$alt),']',sep=""),'',pileup$read.alleles.clean)))

## RPR: I'm trying to see if we can calculate errors by counting the alleles not matching the ref,
## but it is not possible if the mapper does not allow errors. We could hash all alternate bases at the SNP location, this would allow
## to estimate the error rate, by counting the number of bases not matching REF/ALT alleles.
## This is relatively easy to do and could help to estimate accuracy in determining Hets.
pileup$errors <- as.integer(nchar(mapply(gsub,paste('[^ACGT]',sep=""),'',pileup$read.alleles.clean))) - (pileup$alt.matches + pileup$ref.matches)
#sum(pileup$errors)
eps0 <- sum(pileup$errors)/sum(pileup$ref.matches+pileup$alt.matches)
#eps0

## Reorder by position so chr names appear right

chrList <- unique(pileup$chr)
chrList <- paste("chr",sort(as.numeric(gsub("chr","",chrList))),sep="")
factor(chrList,levels=chrList)
pileup$chr <- factor(pileup$chr,levels=chrList)
pileup <- pileup[order(pileup$chr,pileup$pos),];

#####################################################################
##outFile=gsub("pileup.bed.gz", "pileup.clean.bed.gz", pileupFile)
outDir <- gsub(".pileup.clean.bed.gz","",outDir)
system(paste("mkdir -p",outDir))
#print(outDir)
#outDir <- "/wsu/home/fl/fl97/fl9788/piquelab/charvey/ASE/cleanBed1/beds"
fileName <- paste(outDir,"/",sName,"",sep="");
outFile <- paste(fileName,".pileup.clean.bed.gz",sep="");
#print(fileName)
#print(outFile)

## Write the original data just in case we want to view later
## chr, pos-1, pos, ref, alt, rsID, af, ref.matches, alt.matches, errors
# write.table(pileup[,c(1:6,10,13:15)],file=gzfile(paste(outFile,sep='')),row.names=F,col.names=F,sep='\t',quote=F)
##write.table(file=paste(fileName,'.stats.txt',sep=''),saveData,row.names=T,col.names=T,quote=F,sep='\t') # Must be data frame

#####################################################################
#####################################################################
#####################################################################
#####################################################################
# source("/wsu/home/fl/fl97/fl9788/piquelab/charvey/ASE/prod2/fitAseModels.v4.R")
# source("/nfs/hpnfs/groups/piquelab/charvey/ASE/prod2/fitAseModels.v4.R")

## Updated 11/05/13 to source the most recent fitASEModels.v4.R
source("/wsu/home/groups/piquelab/charvey/source/ASE/fitAseModels.v4.R")
source('/wsu/home/groups/piquelab/gmb/AI/results/qqman.r')
require(qvalue)

####################################
oraf <- pileup$ref.matches/(pileup$ref.matches + pileup$alt.matches)
omac <- pmin(pileup$alt.matches,pileup$ref.matches)
otac <- (pileup$ref.matches + pileup$alt.matches)

tr.l <- 15
tr.u <- 10000
ind <- otac>tr.l & otac<tr.u;
sum(ind)
names(ind) <- pileup$rsID
pileup <- pileup[ind,]
dim(pileup)

####################################################################
## Added 10/17/13 by GMB to test new duplicate removal method     ##
## Remove the duplicates using the x8b file                       ##
####################################################################

# expName='OpenChromDnaseFibroblgm03348'
x8Base='/wsu/home/groups/piquelab/allCentipede/x8bFiles.Dnase/'
# expName <- gsub(".*\\.","",sName)
expName <- gsub("test5.","",sName)
pileupTmp <- pileup[c("chr","pos","pos","ref","alt")]
pileupTmp$strand <- '+'

## Write a temporary pileup file to filter the x8b reads to retrieve
tmpFile=paste(tempdir(),'/tmp.pileup',sep='')
write.table(pileupTmp,file=tmpFile,row.names=F,col.names=F,sep='\t',quote=F)

## Get the x8b file
system(paste('cp ', x8Base,'/',expName,'.x8b.bz2 ',tempdir(),sep=''))
system(paste('pbzip2 -dvb1m2000 ', tempdir(), '/', expName, '.x8b.bz2', sep=''))

fileX8b <- paste(tempdir(),'/',expName,'.x8b',sep="")
con <- pipe(paste("bedXbFileXtract ",fileX8b, tmpFile," stdout -window=20"))
cutsite <- read.table(con,sep='\t',as.is=T,fill=T)

## Remove the x8b file when finished
system(paste('rm ',fileX8b, sep=''))

cutsite <- as.matrix(cutsite[,-1])
qFilt = quantile(unlist(cutsite),p=1-1E-4,na.rm=T)
cat('## qFilt:',qFilt,'\n')
m <- apply(cutsite, 1, max)
ind <- m < qFilt

## First 20 are forward strand, second 20 are reverse strand.
## Both sets of 20 end with the SNP position
cutsite <- cutsite[ind, c((1:20),(62:81))]
pileup <- pileup[ind,]
rownames(cutsite) <- pileup$rsID

## Apply a filter on the proprotion of reads coming from each position
## No position should account for more than 50% of the reads covering a SNP
r <- rowSums(cutsite)
m <- apply(cutsite, 1, max)
prop <- m/r

otac <- (pileup$ref.matches + pileup$alt.matches)
sum(r > otac)
mean(r > otac)

## What to do with Na's (36 in test case)? for now will make as 1 and remove
prop[is.na(prop)] <- 1
ind <- prop < 0.5

## redefine based on new pileup
oraf <- pileup$ref.matches/(pileup$ref.matches + pileup$alt.matches)
omac <- pmin(pileup$alt.matches,pileup$ref.matches)
otac <- (pileup$ref.matches + pileup$alt.matches)


####################################################################
####################################################################
####################################################################

ref=pileup$ref.matches[ind];
alt=pileup$alt.matches[ind];
phi <- pileup$af[ind]

##	create output data to asses ASE between samples
out.dat <- pileup[ind, c(1, 3, 6, 13:14)]
# out.dat <- pileup[ind, c(1, 3, 6, 8, 9)]
out.dat$logRA <- round(log(out.dat[, 4]/out.dat[, 5]), digits=4)
colnames(out.dat) <- c("chr", "pos", "rsID", "ref", "alt", "logRA") 
####################################

##bed <- pileup[otac>65,1:3]
##bed$chr <- sub("chr","",bed$chr)
##tmpFile <- paste(Sys.getenv("TMPDIR","/tmp"),"/",sName,'.bed.gz',sep='')
##write.table(bed,file=gzfile(tmpFile),row.names=F,col.names=F,sep='\t',quote=F)

##source('/wsu/home/groups/piquelab/testCentipede/AI/results/fitAseModels.v2.R')
##aux <- fitAse(ref,alt);

asenull <- fitAseNull(ref,alt);
##asenull <- fitAseNull(ref,alt,max.it=1,eps=0.1);

## Use allel frequencies for genotype prior
gmat <- cbind(g0=(1-phi)^2,g1=2*phi*(1-phi),g2=phi^2)
ase2 <- fitAseNull(ref,alt,log.gmat=log((gmat+0.0001)/1.0003))

## attach genotype data to output data
out.dat$gt0 <- round(ase2$gt[,1], digits=4)
out.dat$gt1 <- round(ase2$gt[,2], digits=4)
out.dat$gt2 <- round(ase2$gt[,3], digits=4)

eps <- ase2$eps;
##ase3 <- fitAseNull(ref,alt,log.gmat=log(gmat),fixGprior=F)
##ase2 <- fitAse(ref,alt,log.gmat=asenull$log.gt,eps=0.1,fixEps=TRUE,fixT1=TRUE,max.it=300,t1=1E-6);
##ase2 <- fitAse(ref,alt,eps=0.1,fixEps=TRUE,fixT1=TRUE,max.it=300,t1=1E-6);

##eps <- 0.2

## Calculate rho using eps
rho <- comp.rho(ref,alt,eps)

## Simple rho estimate
rho0 <- plogis(log(ref)-log(alt));

##logliktop = ref * log((1 - eps) * rho + eps * (1 - rho)) + alt * log((1 - eps) * (1 - rho) + eps * rho);
##logRatio <- logliktop-asenull$loglik
##pval <- (1-pchisq(2*logRatio,df=1))

lrt <- lrtEpsRhoBinom(ref,alt,eps,rho)
pval <- (1-pchisq(2*lrt,df=1))

qv <- qvalue(pval,pi0.method="bootstrap")
sum(qv$qv<0.1)

log.het <- pmin(ase2$log.gt[,2],asenull$log.gt[,2])
het <- exp(log.het)

hetInd <- het > 0.99
sum(hetInd)

## ## Find MLE for Beta-Binomial dispersion parameter
## aux <- optim(1,fn=fnLogLikBetaBinomialNull,gr=grLogLikBetaBinomialNull,R=ref[hetInd],A=alt[hetInd],method="L-BFGS-B",lower=1,upper=10,hessian=T)
## se <- sqrt(1/aux$hessian)
## Dmax <- exp(aux$par-0*se)
## Dmax
## str(aux)

############################################

D <- 1:100;
D <- exp((0:500)/50)
##D <- seq(1,40,0.01)
aux <- sapply(D,function(D){
			sum(logLikBetaBinomialRhoEps(0.5,eps,D,ref[hetInd],alt[hetInd]))
			##fnLogLikBetaBinomialNull(log(D),ref[hetInd],alt[hetInd])
			##grLogLikBetaBinomialNull(log(D),ref[hetInd],alt[hetInd])
		})
D[which.max(aux)]
##plot(D,aux,log="x",pch='.',cex=3)

Dmax <- min(D[which(aux>=max(aux)-2)])
Dmax <- D[which.max(aux)]
Dmax
cat("#Dmax:\t",Dmax,"\n");


## Find MLE for rho from the Beta-Binomial. 
Dmax2 <- Dmax
#aux <- optim(rep(0,sum(hetInd)),fn=logLikBetaBinomial2,gr=gLogLikBetaBinomial,D=Dmax2,R=ref[hetInd],A=alt[hetInd],method="L-BFGS-B")
auxLogis <- optim(rep(0,sum(hetInd)),fn=logLikBetaBinomial2,gr=gLogLikBetaBinomial,D=Dmax2,R=ref[hetInd],A=alt[hetInd],method="L-BFGS-B", hessian=TRUE)
rho3 <- plogis(auxLogis$par)
## Could I use rho2 instead?


## hard coded second derivative of the llk; identical to approximated result above
#rho.var.est <- g2LogLikBetaBinomial(auxLogis$par, D=Dmax2, R=ref[hetInd], A=alt[hetInd])


## Recaluclate Het LRT using the beta-bionomial 
lrt3 <- logLikBetaBinomialRhoEps(rho3,eps,Dmax2,ref[hetInd],alt[hetInd]) - logLikBetaBinomialRhoEps(0.5,eps,Dmax,ref[hetInd],alt[hetInd])
pval3 <- (1-pchisq(2*lrt3,df=1))

## simple rho calcuation using epsilon
rho2 <- rho
rho2[hetInd] <- rho3

aux <- pmax(logLikBetaBinomialRhoEps(0.0,eps,Dmax,ref,alt),
		logLikBetaBinomialRhoEps(1.0,eps,Dmax,ref,alt),
		logLikBetaBinomialRhoEps(0.5,eps,Dmax,ref,alt))
lrt2 <- logLikBetaBinomialRhoEps(rho2,eps,Dmax2,ref,alt) - aux;
pval2 <- (1-pchisq(2*lrt2,df=1))

qv2 <- qvalue(pval2[hetInd],pi0.method="bootstrap")
sum(qv2$qv<0.1)
100-qv2$pi0*100

## output data including beta estimates for the heterozygotes
out.dat.hets <- out.dat[hetInd, ]
out.dat.hets$betas <- round(auxLogis$par, digits=4)
out.dat.hets$se <- round(1/(diag(auxLogis$hessian)^.5), digits=4)
out.dat.hets$z <- round(out.dat.hets$betas/out.dat.hets$se, digits=4)
out.dat.hets$pvals <- round(2*pnorm(-abs(out.dat.hets$z)), digits=6)
out.dat.hets$qvals <- round(qv2$qv, digits=6)

write.table(out.dat.hets, file=gzfile(paste(outDir,'/het.data.txt.gz',sep='')), quote=F, row.names=F, col.names=T, sep='\t')

####################################

####################################

aux <- (pbeta(1-eps,1+alt,1+ref)-pbeta(eps,1+alt,1+ref))/(1-2*eps)
##bf <- aux/exp(asenull$loglik)
bf <- aux/exp(ase2$loglik)
##bf <- aux/exp(pmax(log.gt[,"g0"],log.gt[,"g1t0"],log.gt[,"g2"]))

##hist(log(bf[bf>4]),breaks=100)
##plot(-log10(pval2),log(bf[hetInd]))

aux <- cumsum(sort(bf))
aux <- aux/(1:length(aux))
p0 <- mean(aux<1.1)

########################################################################
########################################################################
########################################################################
########################################################################

cat("#rhoMeanCI:",round(mean(rho3)-sd(rho3)*2/sqrt(length(rho)),digits=4),
		round(mean(rho3)+sd(rho3)*2/sqrt(length(rho)),digits=4),"\n",sep="\t");


qval <- het*0+1
qval[hetInd] <- qv2$qval

pileup2 <- cbind(pileup[ind,c("chr","pos-1","pos","ref","alt","rsID","af","errors","ref.matches","alt.matches")],otac=otac[ind],het=het,ase2$gt,pval2=pval2,bf=bf,qval)
# pileup2[pileup2$qval<0.1,]
# write.table(file=gzfile(paste(fileName,'.all.bed.gz',sep='')),pileup2,row.names=F,col.names=F,quote=F,sep='\t')
# write.table(file=gzfile(paste(fileName,'.sig.bed.gz',sep='')),pileup2[pileup2$qval<0.2,],row.names=F,col.names=T,quote=F,sep='\t')


## Examine pval3
qv3 <- qvalue(pval3,pi0.method="bootstrap")
sum(qv3$qv<0.1)
100-qv3$pi0*100
qval3 <- het*0+1
qval3[hetInd] <- qv3$qval

# Adjust the out.dat.hets

se    = rep(NA, length(het))
betas = rep(NA, length(het))
z     = rep(NA, length(het))
se[hetInd]    <- out.dat.hets$se
betas[hetInd] <- out.dat.hets$betas
z[hetInd]     <- out.dat.hets$z

tmpPval <- rep(1,length(het))
tmpPval[hetInd] <- pval3
pileup3 <- cbind(pileup[ind,c("chr","pos-1","pos","ref","alt","rsID","af","errors","ref.matches","alt.matches")],otac=otac[ind],se=se,betas=betas,z=z,het=het,ase2$gt,bf=bf,pval2=pval2,qval.p2=qval,pval3=tmpPval,qval.p3=qval3)
write.table(file=gzfile(paste(fileName,'.all.bed.gz',sep='')),pileup3,row.names=F,col.names=T,quote=F,sep='\t')
write.table(file=gzfile(paste(fileName,'.hets.bed.gz',sep='')),pileup3[hetInd,],row.names=F,col.names=T,quote=F,sep='\t')


#######################
## Plotting the data ##
#######################
png(file=paste(fileName,'.qqP.z.png',sep=''))

rho.se <- 1/(diag(auxLogis$hessian)^.5)
rho.z <- (auxLogis$par)/rho.se
rho.pvals <- 2*pnorm(-abs(rho.z))

qq(rho.pvals)
title(main=paste(outDir, " p-values", sep=''))
my.points <- qqplot(-log10(ppoints(sum(hetInd))),-log10(pval3),plot.it=F)
points(my.points, cex=5, col='blue', pch='.')
abline(a=0, b=1, col='red')
legend("topleft", c(paste("Z-score"), paste("LLK ratio test")),
				fill=c("black","blue"))
dev.off()

########################################################################
png(file=paste(fileName,'.qqNorm.z.png',sep=''))
qqnorm(rho.z, main=paste("qqNorm: ", outDir, " Z-scores", sep=''))
qqline(rho.z)
dev.off()

########################################################################
## Plot quality values to see if they make any sense
png(file=paste(fileName,'.plots.quality.png',sep=''))

plot(as.numeric(names(qual.table)),as.numeric(qual.table),xlab="Base calling qualities",ylab="Freq.",pch='.',cex=6,log='y')
title(sName)
abline(h=0,lty=3)

dev.off()

########################################################################
png(file=paste(fileName,'.plots.ecdf.png',sep=''))
plot(ecdf(log10(otac[otac>0])),main=sName,xlab=expression(log[10](x)),ylab=expression(F[x>0](x)))
abline(v=log10(c(tr.l,tr.u)),lty=3)
dev.off()

########################################################################
##chrCol <- rainbow(length(chrList))
chrCol <- rep(c("orange","darkblue"),length(chrList))[1:length(chrList)]
names(chrCol) <- chrList

png(file=paste(fileName,'.plots.oraf.png',sep=''),width=800,height=400)
layout(t(c(1,1,1,2)))
par(cex=1.0)
oldmar <- par("mar")
mar <- oldmar
mar[4] <- 0.0
par(mar=mar)
ind2 <- ind & (abs(pileup$af-0.5)<0.4)
##x <- 1:sum(ind2)
plot(oraf[ind2],xlab="Chromosome order",ylab="Obs. reference allele Freq.",pch='.',cex=3,col=chrCol[pileup$chr[ind2]],axes=F)
axis(2)
x.at <- c(which(!duplicated(pileup$chr[ind2])),sum(ind2))
axis(1,at=x.at,labels=FALSE,cex=0.3)
abline(v=x.at,lty=3)
text((x.at[-1]+x.at[-length(x.at)])*0.5, -0.15, labels = chrList, srt = 45, pos =3, xpd = TRUE,cex=0.7)
title(sName)
abline(h=0.5,lty=3)
mar <- oldmar
mar[2] <- 0.0
par(mar=mar)
aux <- hist(jitter(oraf[ind2]),breaks=101,plot=F)
barplot(aux$counts+1,horiz=TRUE,log="x")
title(xlab="Freq.")
par(mar=oldmar)
dev.off()

########################################################################
png(file=paste(fileName,'.plots.rhoVsPhi.png',sep=''))

plot(phi[hetInd],rho3,pch='.',cex=2,xlab=expression(phi),ylab=expression(hat(rho)),main=sName,xlim=c(0,1),ylim=c(0,1))
lines(smooth.spline(phi[hetInd],rho3),col='blue',lwd=3)

## Get the mean value rho, and the upper and lower stderrs
rho.mean   <- mean(rho3)
rho.stderr <- sd(rho3)*2/sqrt(length(rho3))
rho.upper  <- rho.mean + rho.stderr
rho.lower  <- rho.mean - rho.stderr
cat('## rho:',rho.mean,rho.upper,rho.lower,'\n',sep='\t')

## Check the correlation between rho and phi
rho.cor <- cor.test(rho3, phi[hetInd])
cat('## phi vs rho:', rho.cor$estimate, rho.cor$p.value,'\n',sep='\t')

legend("topleft",
	paste("rho: ", round(rho.mean,digits=6)," +/- ",round(rho.stderr,digits=5),'\n',
		    "cor: ", round(rho.cor$estimate,digits=4), ", p-val: ", round(rho.cor$p.value,digits=4),
		    sep=''))
dev.off()


########################################################################
png(file=paste(fileName,'.plots.rho.png',sep=''))

plot(rho0,rho,pch='.',cex=2,xlab=expression(rho[0]),ylab=expression(hat(rho)),main=sName)
points(rho0[hetInd],rho3,pch='.',cex=2,col='blue')
abline(0,1,col='red',lty=2)

dev.off()


########################################################################
png(file=paste(fileName,'.plots.qqplot.png',sep=''))

qq(pval[hetInd])
qqp <- qqplot(-log10(ppoints(sum(hetInd))),-log10(pval3),plot.it=F)
points(qqp,pch='.',cex=4,col='blue')

qqp <- qqplot(-log10(ppoints(sum(hetInd))),-log10(pval2[hetInd]),plot.it=F)
points(qqp,pch='.',cex=7,col='darkgreen')
legend("topleft",c(paste("D = Inf, eps =",round(eps*100,digits=2),"%"),
				paste("Het 99%, D = ",round(Dmax,digits=2)),
				paste("Het unc., pi0 = ",round(qv2$pi0*100,digits=2),"%")),fill=c("black","blue","darkgreen"))

dev.off()


########################################################################
cat("End:",sName,"\n")
## THE-END

