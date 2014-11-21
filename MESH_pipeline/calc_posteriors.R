#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
parmfiles <- scan(pipe('ls *.hm_config.est'), character(0))
bffiles <- scan(pipe('ls *.hm_config.in'), character(0))
n.samples <- length(parmfiles)

for(ii in 1:n.samples){

#parmf <- "P1_18507_T15C1.hm_config.est"
#bffile <- "P1_18507_T15C1.hm_config.in"

parmf <- parmfiles[ii]
bffile <- bffiles[ii]   

pi0 <- as.numeric(scan(pipe(paste("less ", parmf, " | awk '$1==\"pi0:\" {print $2}'", sep='')),character(0)))
pi00 <- 1-pi0

config1 <- as.numeric(scan(pipe(paste("less ", parmf, " | awk '$1==\"config:\" {print $2}'", sep='')),character(0)))
config2 <- as.numeric(scan(pipe(paste("less ", parmf, " | awk '$1==\"config:\" {print $4}'", sep='')),character(0)))
config3 <- as.numeric(scan(pipe(paste("less ", parmf, " | awk '$1==\"config:\" {print $6}'", sep='')),character(0)))
config <- c(config1, config2, config3)

priors <- c(pi00 * config, pi0)
sum(priors)

grid1 <- as.numeric(scan(pipe(paste("less ", parmf, " | awk '$1==\"grid:\" {print $2}'", sep='')),character(0)))
grid2 <- as.numeric(scan(pipe(paste("less ", parmf, " | awk '$1==\"grid:\" {print $4}'", sep='')),character(0)))
grid3 <- as.numeric(scan(pipe(paste("less ", parmf, " | awk '$1==\"grid:\" {print $6}'", sep='')),character(0)))
grid <- c(grid1, grid2, grid3)

data <- read.table(bffile)
names(data) <- c('rsID', 'config', 'bf1', 'bf2', 'bf3')

snps <- as.character(unique(data[, 'rsID']))

bfs <- t(sapply(snps, function(snp){
  #snp <- snps[1]

  temp <- subset(data, rsID==snp, select=c('bf1', 'bf2', 'bf3'))

  bf1 <- temp[1, ]
  bf2 <- temp[2, ]
  bf3 <- temp[3, ]

  bf1 <- (10^bf1) %*% grid
  bf2 <- (10^bf2) %*% grid
  bf3 <- (10^bf3) %*% grid

  bfs <- c(bf1, bf2, bf3, 1)
}))

colnames(bfs) <- c('bf_c1', 'bf_c2', 'bf_c3', 'bf_c4')

n.snps <- dim(bfs)[1]

post <- t(sapply(1:n.snps, function(ii){
  #ii <- 1
  this.bf <- bfs[ii, c('bf_c1', 'bf_c2', 'bf_c3', 'bf_c4')]
  const <- 1/sum((this.bf*priors))
  post <- this.bf*priors*const
  post
}))

colnames(post) <- c('pos_c1', 'pos_c2', 'pos_c3', 'pos_c4')
post <- data.frame(snps, post, stringsAsFactors=FALSE, row.names=NULL)

bfs <- data.frame(snps, bfs, stringsAsFactors=FALSE, row.names=NULL)

system('mkdir -p posteriors_bfs')
root <- gsub('.hm_config.est', '', parmf)

write.table(post, file=paste('./posteriors_bfs/', root, '_posteriors.txt', sep=''), quote=FALSE, row.names=FALSE, sep='\t')
write.table(bfs, file=paste('./posteriors_bfs/', root, '_bayesFactors.txt', sep=''), quote=FALSE, row.names=FALSE, sep='\t')

}
