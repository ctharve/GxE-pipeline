###########################################################################
##
## CTH 111714
## helper functions for processing and plotting oputput from William Wen's MESH
##
###########################################################################
## capture full treatment names from treatment codes for system treatments
get_treat <- function(t.code){
  temp <- read.table('/wsu/home/groups/piquelab/charvey/GxE/derived_data/covariates/GxE_treatment_key.new.txt',
                     header=TRUE,
                     stringsAsFactors=FALSE)
  temp[temp$Treatment_ID==t.code, 'Short_Name'] 
}

## 
get_info <- function(ids){
  info.dat <- Reduce(rbind, lapply(ids, function(id){
    ## id <- betaTreat_dat$id[1]
    myLabels <- strsplit(gsub('.txt', '',  gsub('allPlates_', '', id)), ':')[[1]]
    direction <- strsplit(myLabels[1], "_")[[1]][1] 
    bin <- as.numeric(strsplit(myLabels[1], "_")[[1]][2])
    rsID <- strsplit(myLabels[2], "_")[[1]][1]
    myTreat <- get_treat(strsplit(myLabels[2], "_")[[1]][2])
    #treatment <- strsplit(myLabels[2], "_")[[1]][2]
    index <- as.numeric(strsplit(myLabels[2], "_")[[1]][3])
    out_dat <- data.frame(direction, bin, id=paste(rsID, myTreat, sep='_'),  stringsAsFactors=FALSE)
    out_dat
  }))
  info.dat
}

## Capture posterior samples from MESH models for violin plots
## fileList must have absolute path
meshf_capturePosteriorSamples <- function(fileList, filterLevel){  
  myFiles <- fileList
  if(filterLevel=='control'){
      filter_col <- 3
      filter_label <- 'control'
  } else if (filterLevel== 'treatment'){
      filter_col <- 4
      filter_label <- 'treatment'
  } else {
      stop('error: need to specify \'filterLevel\' control or treatment')
  }
  
  return_dat <- Reduce(rbind, lapply(myFiles, function(file){
    ## file <- myFiles[1]
    target <- paste0('less ', file, ' | tail -10000 | awk \'{print $', filter_col, '}\'')
    myDat <- as.numeric(scan(pipe(target), character(0)))
    n.samples <- length(myDat)
    fileName <- tail(strsplit(file, '/')[[1]], n=1) 
    myLabels <- strsplit(gsub('.hm_config.out', '',  gsub('allPlates_', '', fileName)), '_')[[1]]

    ## if original data, no-binning, use 'all_data' as subsampling label
    if(!(myLabels[3] %in% c('600', '1200', '2400'))){
      count_tag <- 'all_data'
    } else {
      count_tag <- myLabels[3]
    }
    
    out_dat <- data.frame(direction=rep(myLabels[1], n.samples), bin=rep(myLabels[2], n.samples), nSNPS=rep(count_tag, n.samples), samples=myDat, stringsAsFactors=FALSE)  
    out_dat$treat <- filter_label
    out_dat
  }))
  return_dat
}                             


meshf_median.quartile <- function(x){
    out <- quantile(x, probs = c(0.025,0.5,0.975))
      names(out) <- c("ymin","y","ymax")
      return(out)
  }
meshf_med <- function(x){
    out <- quantile(x, probs = c(0.5))
      names(out) <- c("y")
      return(out)
  }
