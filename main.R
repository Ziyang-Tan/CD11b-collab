library(dplyr)
library(flowCore)
library(mmR)
library(tibble)

read_and_filter_fcs = function(file){
  dat = file@exprs[,-1] 
  colnames(dat) = markernames(file)
  nonEmptyChannels = file@parameters@data[c("name", "desc")] %>% filter(name != desc)
  dat = dat[,colnames(dat) %in% nonEmptyChannels$desc]
  return(dat)
}

target_channel <- 'CD11b'
target_subpop <- c('Monocytes', 'Neutrophils')
exp_path <- '/Users/tan/cytof_data/EXP-19-CL8755'

namelist <- list.files('/Users/tan/cytof_data/EXP-19-CL8755/', recursive = T) %>% 
  basename() %>% 
  strsplit('.', fixed = T) %>% 
  sapply(function(x)x[1])
namelist <- namelist[namelist != 'abundance']

datlist <- lapply(namelist, function(x){
  file_path <- Sys.glob(paste0(exp_path, '/*/*/', x, '.fcs'))
  info_path <- Sys.glob(paste0(exp_path, '/*/*/', x, '.csv'))
  file <- read.FCS(file_path, truncate_max_range = F) %>%
    read_and_filter_fcs()
  info <- mm.fastread(info_path)
  dat <- info %>%
    add_column(target_channel = file[,target_channel]) %>%
    filter(level1 %in% target_subpop) %>%
    select(target_channel, level1, level2)
  dat <- dat %>%
    add_column(source = rep(strsplit(x, '_')[[1]][1], dim(dat)[1]))
  return(dat)
})
dat <- datlist %>%
  bind_rows()





