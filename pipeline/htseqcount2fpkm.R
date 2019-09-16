library(dplyr)
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop(" input htseq-count, gene length table and output fpkm file must be supplied. ", call.=FALSE)
}
htseq.count <- read.table(args[1],header=F,sep="\t",stringsAsFactors=F)
gene_length <- read.table(args[2],header=T,sep="\t",stringsAsFactors=F)
colnames(htseq.count) <- c("gene_name","count")
gene_merge <- gene_length %>% left_join(htseq.count)
N <- sum(gene_merge$count)
gene_merge_fpkm <- transform(gene_merge, fpkm = count / length / N * 10^9)
write.table(gene_merge_fpkm,file=args[3],sep="\t",row.names=F,quote=F)
