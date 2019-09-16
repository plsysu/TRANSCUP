data_dir <- "./output/"
result_dir <- "./output/01_raw_fpkm/"
dir.create(result_dir)
# load all tumor types: 33
tumor_types <- read.table("../../data/tumor_type.tsv",header = F,sep = "\t",stringsAsFactors = F)
tumor_types <- tumor_types$V1

all33_sample_info_df <- NULL
# load sample info
for(t in 1:length(tumor_types)){
  sample_info <- read.table(paste0(data_dir,tumor_types[t],"_sample_info.tsv"),header = T,sep = "\t",stringsAsFactors = F)
  sample_info$tumor_type <- unlist(strsplit(tumor_types[t],"-"))[2]
  if(is.null(all33_sample_info_df)){
    all33_sample_info_df <- sample_info
  }else{
    all33_sample_info_df <- rbind(all33_sample_info_df,sample_info)
  }
}
write.table(all33_sample_info_df,file = paste0(result_dir,"all33_sample_info.tsv"),sep = "\t",row.names = F)