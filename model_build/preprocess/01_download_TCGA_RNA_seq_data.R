# This code will get all clinical indexed data from TCGA
library(data.table)
library(dplyr)
library(regexPipes)
library(TCGAbiolinks)
library(DT)

dir.create("./output")


projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl(projects,'^TCGA',perl=T)]
projects <- "TCGA-BRCA"
for(proj in projects){
  print(proj)

  query.exp.hg38 <- GDCquery(project = proj,
                             data.category = "Transcriptome Profiling",
                             data.type = "Gene Expression Quantification",
                             workflow.type = "HTSeq - FPKM"
  )
  GDCdownload(query.exp.hg38,directory = "./output")

  info_table <- query.exp.hg38$results[[1]]


  write.table(info_table,file = paste0("./output/",proj,"_sample_info.tsv"),sep = "\t",row.names = F)
}
 
