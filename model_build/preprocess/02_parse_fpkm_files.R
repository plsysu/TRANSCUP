data_dir <- "./output/"
result_dir <- "./output/01_raw_fpkm/"
dir.create(result_dir)

# load all tumor types: 33
tumor_types <- read.table("../../data/tumor_type.tsv",header = F,sep = "\t",stringsAsFactors = F)
tumor_types <- tumor_types$V1

all33_fpkm_df <- NULL
# load sample info
for(t in 1:length(tumor_types)){
  sample_info <- read.table(paste0(data_dir,tumor_types[t],"_sample_info.tsv"),header = T,sep = "\t",stringsAsFactors = F)
  all_fpkm_df <- NULL
  for(i in 1:dim(sample_info)[1]){
    fpkm_df <- read.table(gzfile(paste0(data _dir,tumor_types[t],
                                        "/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/",
                                        sample_info$file_id[i],
                                        "/",
                                        sample_info$file_name[i])),
                          header = F,
                          sep = "\t",
                          row.names = 1)
    colnames(fpkm_df) <- sample_info$cases[i]
    if(is.null(all_fpkm_df)){
      all_fpkm_df <- fpkm_df
    }else{
      all_fpkm_df <- cbind(all_fpkm_df,fpkm_df)
    }
    
  }
  
  write.csv(all_fpkm_df,file = paste0(result_dir,tumor_types[t],"_fpkm.csv"))
  if(is.null(all33_fpkm_df)){
    all33_fpkm_df <- all_fpkm_df
  }else{
    all33_fpkm_df <- cbind(all33_fpkm_df,all_fpkm_df)
  }
}

write.csv(all33_fpkm_df,file = paste0(result_dir,"all33_fpkm.csv"))