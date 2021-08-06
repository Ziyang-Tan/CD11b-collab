library(bigrquery)

query_from_bigquery <- function(sql){
  #sql <- "SELECT * FROM `cradle-259115.HLA.typing`"
  tb <- bq_project_query('cradle-259115', sql)
  return(bq_table_download(tb, bigint = 'integer64'))
}


#sql <- "select sample_id, subject_type from cradle.all_samples where sample_id in (2600100430, 2600100439)"
#tb <- bq_project_query('cradle-259115', sql)
#bq_table_download(tb, bigint = 'integer64')
