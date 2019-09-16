# First, import the GTF-file that you have also used as input for htseq-count
library(GenomicFeatures)
args = commandArgs(trailingOnly=TRUE)
input_gtf <- args[1]
output_gene_length_file <- args[2]

txdb <- makeTxDbFromGFF(input_gtf,format="gtf")
# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
gene_length <- unlist(exonic.gene.sizes)
gene_length_df <- data.frame(gene_name=names(gene_length),length=as.numeric(gene_length))
write.table(gene_length_df,file=output_gene_length_file,sep="\t",row.names=F,quote=F)
