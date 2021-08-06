library(dplyr)
library(flowCore)
library(mmR)
library(tibble)
library(ggplot2)
source('query_from_bigquery.R')

read_and_filter_fcs = function(file){
  dat = file@exprs[,-1] 
  colnames(dat) = markernames(file)
  nonEmptyChannels = file@parameters@data[c("name", "desc")] %>% dplyr::filter(name != desc)
  dat = dat[,colnames(dat) %in% nonEmptyChannels$desc]
  return(dat)
}

# settings
target_channel <- 'CD11b'
target_subpop <- c('Monocytes', 'Neutrophils')
exp_id <- 'EXP-19-CL8757'
exp_path <- file.path('/Users/tan/cytof_data', exp_id)

# load data
namelist <- list.files(exp_path, recursive = T) %>% 
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
    dplyr::filter(level1 %in% target_subpop) %>%
    select(target_channel, level1, level2)
  dat <- dat %>%
    add_column(sample_id = rep(strsplit(x, '_')[[1]][1], dim(dat)[1]))
  return(dat)
})
dat <- datlist %>%
  bind_rows()

# fetch sample meta from bigquery
ids <- dat %>% select(sample_id) %>% distinct() %>% unlist(use.names = F)
sql <- paste0('select distinct sample_id, subject_type, study_id, family_id, baby_age 
              from cradle.all_samples where sample_id in (',
              paste(ids, collapse = ','), ')')
meta <- query_from_bigquery(sql)
meta <- meta[!base::duplicated(meta$sample_id, fromLast = T),]
sql2 <- paste0('select distinct `group`, study_id 
               from cradle.subjects_baby where study_id in (',
               paste(meta$study_id, collapse = ','), ')')  
group_info <- query_from_bigquery(sql2)
meta <- meta %>% 
  left_join(group_info, by = 'study_id') %>%
  mutate(sample_id = as.character(sample_id)) %>%
  mutate(group = case_when(
    subject_type %in% c('Father', 'Mother') ~ 'Adult',
    TRUE ~ group
  ))

# process 
dat.monocytes <- dat %>% 
  mutate(target_channel = asinh(target_channel)) %>% #transformation
  dplyr::filter(level1 == 'Monocytes') %>% 
  group_by(sample_id) %>% 
  summarise(exp = mean(target_channel)) %>%
  left_join(meta, by = 'sample_id')

dat.neutrophils <- dat %>% 
  mutate(target_channel = asinh(target_channel)) %>% #transformation
  dplyr::filter(level1 == 'Neutrophils') %>% 
  group_by(sample_id) %>% 
  summarise(exp = mean(target_channel)) %>%
  left_join(meta, by = 'sample_id')

mm.fastwrite(dat.monocytes, path = file.path('/Users/tan/CD11b-collab/processed_data/Monocytes', paste0(exp_id, '.csv')))
mm.fastwrite(dat.neutrophils, path = file.path('/Users/tan/CD11b-collab/processed_data/Neutrophils', paste0(exp_id, '.csv')))


