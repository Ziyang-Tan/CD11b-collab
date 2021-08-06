library(dplyr)
library(mmR)
library(tibble)
library(ggplot2)


dir_path <- '/Users/tan/CD11b-collab/processed_data/Neutrophils'
dat_list <- lapply(list.files(dir_path), function(x){
  batch <- strsplit(x, '.', fixed = T)[[1]][1]
  d <- mm.fastread(file.path(dir_path, x)) %>%
    add_column(batch = batch)
})

dat <- bind_rows(dat_list) %>%
  dplyr::filter(!is.na(study_id))

ggplot(dat, aes(x = group, y = exp)) +
  geom_boxplot()

t.test((dat %>% dplyr::filter(group == 'Preterm'))$exp, 
       (dat %>% dplyr::filter(group == 'Term'))$exp)

mm.fastwrite(dat, path = '/Users/tan/CD11b-collab/processed_data/CD11b_Neutrophils.csv')
