get_treat <- function(t.code){
  temp <- read.table('/wsu/home/groups/piquelab/charvey/GxE/derived_data/covariates/GxE_treatment_key.new.txt', header=TRUE, stringsAsFactors=FALSE)
  temp[temp$Treatment_ID==t.code, 'Short_Name'] 
}

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
