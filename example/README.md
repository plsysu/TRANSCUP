# TRANSCUP usage example

We use a public RNA-seq data set contains 20 samples to show the usage of TRANSCUP pipeline.

The accession number is [SRR2089755](https://www.ncbi.nlm.nih.gov/sra/?term=SRR2089755[accn]). 

The reference paper is:

> Lee JR et al., "Transcriptome analysis of paired primary colorectal carcinoma and liver metastases reveals fusion transcripts and similar gene expression profiles in primary carcinoma and liver metastases.", BMC Cancer, 2016 Jul 26;16:539
# Install software to download NCBI data
* [Aspera connect](https://downloads.asperasoft.com/connect2/)

* [sratoolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/)


## Step0: Download the SRA file
Use Aspera or other FTP method.
<pre>
mkdir ./sra

for i in `cat SRR2089755_SRR_Acc_List.txt`
do
	id1=`echo $i |cut -c1-3`
	id2=`echo $i |cut -c1-6`
	echo $id1
	echo $id2
	echo $i
	~/.aspera/connect/bin/ascp -i  ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -T -k 1 anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/${id1}/${id2}/${i}/${i}.sra ./sra
done
</pre>

## Step1: Convert sra  to fastq
<pre>
mkdir ./fq_raw
fastq-dump --split-3 --gzip  -O  ./fq_raw  ./sra/*.sra
</pre>
## Step3: Prepare sample_sheet.txt
Change every fastq's file path to the actual file's absolute pathname.
## Step4: Prepare data
<pre>
python prepare_data.py \
	--sample_sheet ./sample_sheet.txt
	--seq_type PE \
	--out_dir ./fq_in
</pre>
## Step5: Make a Snakefile
<pre>
python ./transcup.py \
	--fq_dir ./fq_in \
	--config ./config/database.software.resource.json \
	--cluster_config ./config/slurm.cluster.resource.json \
	--out_dir ./output
</pre>
## Step6: Run locally or  on  a  SLURM/SGE cluster
<pre>
cd ./output
nohup sh RUN_local.sh > log.local &
</pre>
<pre>
nohup sh RUN_slurm_cluster.sh > log.slurm &
</pre>
<pre>
nohup sh RUN_SGE_cluster.sh > log.sge &
</pre>