

OUTFOLDER = ROOTFOLDER + 'out/'
if not 'TRIMMOMATICIN' in globals():
    TRIMMOMATICIN = FASTQFOLDER
if not 'TRIMMOMATICOUT' in globals():
    TRIMMOMATICOUT = OUTFOLDER + 'cliptrim/'
ALIGNOUT = OUTFOLDER + 'align_star/'
TMPDIR = OUTFOLDER + 'tmp/'
STATSOUT = OUTFOLDER + 'stats/'
HTSEQCOUNTSOUT = OUTFOLDER + 'htseqcount/'
HTSEQCOUNT2FPKMSOUT = OUTFOLDER + 'htseqcount2fpkm/'
RNASEQCOUT = OUTFOLDER + 'rnaseqc/'
PREDICTOUT = OUTFOLDER + 'predict/'
'''
Collecting file names
'''
PAIREDFASTQFILES = [file.replace(FASTQFOLDER, '').replace('.fastq.gz','')for file in glob.glob(FASTQFOLDER + '*/PAIREDEND/*.fastq.gz')]
SINGLEFASTQFILES = [file.replace(FASTQFOLDER, '').replace('.fastq.gz','')for file in glob.glob(FASTQFOLDER + '*/SINGLEEND/*.fastq.gz')]
SAMPLENAMES = [samplename.replace(FASTQFOLDER,'').replace('/','')for samplename in glob.glob(FASTQFOLDER + '*/')]
PAIREDENDSAMPLENAMES = [token for token in set([samplename.split('/')[0].strip() for samplename in PAIREDFASTQFILES])]
SINGLEENDSAMPLENAMES = [token for token in set([samplename.split('/')[0].strip() for samplename in SINGLEFASTQFILES])]
#MAPPINGFILES = glob.glob(FASTQFOLDER + '/*/mapping.tsv')
#print(PAIREDFASTQFILES)
print(PAIREDENDSAMPLENAMES)
print(len(PAIREDENDSAMPLENAMES))

import os
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR)

'''
methods
'''

def getFastqsAfterTrimming(wildcards):
    '''
        get the fastqs after trimming
    '''
    paired = expand(TRIMMOMATICOUT + '{files}_PAIRED.fastq.gz', files=PAIREDFASTQFILES)
    single = expand(TRIMMOMATICOUT + '{files}.fastq.gz', files=SINGLEFASTQFILES)
    all = []
    for i in single:
        all.append(i)
    for i in paired:
        all.append(i)
    return all


def getFastQForSampleSingleEnd(wildcards):
    sample = wildcards.sample
    singles = expand(TRIMMOMATICOUT + '{files}.fastq.gz', files=SINGLEFASTQFILES)
    wc2single = []
    for single in singles:
        if '/' + sample+ '/' in single:
            wc2single.append(single)
    return wc2single
def getFastQForSampleSingleEndStarString(wildcards):
    fastqs = getFastQForSampleSingleEnd(wildcards)
    s = ''
    for fastq in fastqs:
        s = s + ',' + fastq
    return s[1:]

def getFastQForSamplePairedEnd(wildcards):
    sample = wildcards.sample
    paireds = expand(TRIMMOMATICOUT + '{files}_PAIRED.fastq.gz', files=PAIREDFASTQFILES)
    wc2paired = []
    for paired in paireds:
        if '/' + sample+ '/' in paired:
            wc2paired.append(paired)
    return wc2paired

def getFastQForSamplePairedEndStarString(wildcards):
    fastqs = getFastQForSamplePairedEnd(wildcards)
    left = []
    right = []
    for fastq in fastqs:
        if '_R1_PAIRED' in fastq:
            left.append(fastq)
            continue
        if '_R2_PAIRED' in fastq:
            right.append(fastq)
            continue
        raise ValueError('There re fastqs that do not fit the scheme for this sample: ' + fastq)
    if len(left) != len(right):
        raise ValueError('Not all paired end files are matched: ' + str(fastqs))
    left.sort()
    right.sort()

    return ','.join(left) + ' ' + ','.join(right)



def getAllFastqFiles(wildcards):
    infiles = expand(FASTQFOLDER + '{files}.fastq.gz', files=SINGLEFASTQFILES)
    single = expand(TRIMMOMATICOUT + '{files}.fastq.gz', files=SINGLEFASTQFILES)
    inpairedfiles = expand(FASTQFOLDER + '{files}.fastq.gz', files=PAIREDFASTQFILES)
    cliptrimpaired = expand(TRIMMOMATICOUT + '{files}_PAIRED.fastq.gz', files=PAIREDFASTQFILES)
    all = single + infiles + inpairedfiles + cliptrimpaired
    return all
def getAllFastqCountFiles(wildcards):
    fastqs = getAllFastqFiles(wildcards)
    counts = []
    for fastq in fastqs:
        counts.append(fastq.replace('fastq.gz', 'count'))
    return counts
def getAllFastqcFiles(wildcards):
    fastqs = getAllFastqFiles(wildcards)
    fastqcs = []
    for fastq in fastqs:
        fastqcs.append(fastq.replace('fastq.gz', 'fastqc.done'))
    return fastqcs
def getAllBamFiles(wildcards):
    alignpaired = expand(ALIGNOUT + 'PAIREDEND/{files}Aligned.out.sortedCor.bam', files=PAIREDENDSAMPLENAMES)
    alignsingle = expand(ALIGNOUT + 'SINGLEEND/{files}Aligned.out.sortedCor.bam', files=SINGLEENDSAMPLENAMES)
    align = alignpaired + alignsingle
    return align
def getBamHeaderFiles(wildcards):
    bams = getAllBamFiles(wildcards)
    headers = []
    for bam in bams:
        headers.append(bam.replace('.bam', '.header'))
    return headers
def getFlagstatFiles(wildcards):
    bams = getAllBamFiles(wildcards)
    flagstats = []
    for bam in bams:
        flagstats.append(bam.replace('.bam', '.flagstat'))
    return flagstats
def getAllStatsFilesForASample(wildcards):
    flagstats = getFlagstatFiles(None)
    counts = getAllFastqCountFiles(None)
    sample = wildcards.sample
    out = []
    for flagstat in flagstats:
        if sample in flagstat:
            out.append(flagstat)
    for count in counts:
        if sample in count:
            out.append(count)
    return out
def getAllAlignFiles(wildcards):
    return (expand(ALIGNOUT + 'PAIREDEND/{files}Aligned.out.sortedCor.bam', files=PAIREDENDSAMPLENAMES) +
            expand(ALIGNOUT + 'SINGLEEND/{files}Aligned.out.sortedCor.bam', files=SINGLEENDSAMPLENAMES))
def getAllHtseqFiles(wildcards):
    return expand(HTSEQCOUNTSOUT + 'PAIREDEND/{sample}.count.tsv', sample=PAIREDENDSAMPLENAMES) + expand(HTSEQCOUNTSOUT + 'SINGLEEND/{sample}.count.tsv', sample=SINGLEENDSAMPLENAMES)

def getAllFPKMFiles(wildcards):
    return expand(HTSEQCOUNT2FPKMSOUT + 'PAIREDEND/{sample}.fpkm.tsv', sample=PAIREDENDSAMPLENAMES) + expand(HTSEQCOUNT2FPKMSOUT + 'SINGLEEND/{sample}.fpkm.tsv', sample=SINGLEENDSAMPLENAMES)

def getAllrnaseqcFiles(wildcards):
    return expand(RNASEQCOUT + 'PAIREDEND/{sample}/{sample}Aligned.out.sortedCor.bam.metrics.tsv', sample=PAIREDENDSAMPLENAMES) + expand(RNASEQCOUT + 'SINGLEEND/{sample}/{sample}Aligned.out.sortedCor.bam.metrics.tsv', sample=SINGLEENDSAMPLENAMES)


def getAllRandomForestPredictFiles(wildcards):
    return expand(PREDICTOUT + 'PAIREDEND/{sample}_RF_pred_probs.csv', sample=PAIREDENDSAMPLENAMES) + expand(PREDICTOUT + 'SINGLEEND/{sample}_RF_pred_probs.csv', sample=SINGLEENDSAMPLENAMES)


'''
rules

'''
rule all:
    input:
        #getBamHeaderFiles,
        #getAllFastqcFiles,
        #getFlagstatFiles,
        #getAllFastqCountFiles,
        getAllHtseqFiles,
        getAllFPKMFiles,
        getAllAlignFiles,
        getAllrnaseqcFiles,
        getAllRandomForestPredictFiles
    output:
        touch(ROOTFOLDER + 'complete.txt')


rule trimmomatic_paired:
    input:
        forwardF = TRIMMOMATICIN + '{sample}/PAIREDEND/{fastq}_R1.fastq.gz',
        reverseF = TRIMMOMATICIN + '{sample}/PAIREDEND/{fastq}_R2.fastq.gz',
    output:
        forwardP = TRIMMOMATICOUT + '{sample}/PAIREDEND/{fastq}_R1_PAIRED.fastq.gz',
        forwardUP = TRIMMOMATICOUT + '{sample}/PAIREDEND/{fastq}_R1_UNPAIRED.fastq.gz',
        reverseP = TRIMMOMATICOUT + '{sample}/PAIREDEND/{fastq}_R2_PAIRED.fastq.gz',
        reverseUP = TRIMMOMATICOUT + '{sample}/PAIREDEND/{fastq}_R2_UNPAIRED.fastq.gz',
        #trimlog = TRIMMOMATICOUT + '{sample}/PAIREDEND/{fastq}_clipTrim.log.gz'
    params:
        #trimlog = temp(TRIMMOMATICOUT + '{sample}/PAIREDEND/{fastq}_clipTrim.log'),
        adapter = config['resources']['general']['sequencingAdapter'],
        slidingwindow = config['tools']['trimmomatic']['paired']['slidingwindow'],
        phred = config['tools']['trimmomatic']['paired']['phred'],
        mode = config['tools']['trimmomatic']['paired']['mode'],
        minQual = config['tools']['trimmomatic']['paired']['minQual'],
        leadingQual = config['tools']['trimmomatic']['paired']['leadingQual'],
        tailingQual = config['tools']['trimmomatic']['paired']['tailingQual'],
        seedmismatches = config['tools']['trimmomatic']['paired']['seedmismatches'],
        palindrom = config['tools']['trimmomatic']['paired']['palindrom'],
        score = config['tools']['trimmomatic']['paired']['score'],
        minlen = config['tools']['trimmomatic']['paired']['minlen'],
        min_adapt_len = config['tools']['trimmomatic']['paired']['min_adapt_len'],
        keep_both = config['tools']['trimmomatic']['paired']['keep_both'],
        lsfoutfile = TRIMMOMATICOUT + '{sample}/PAIREDEND/{fastq}.fastq.lsfout.log',
        lsferrfile = TRIMMOMATICOUT + '{sample}/PAIREDEND/{fastq}.fastq.lsferr.log',
        scratch = config['tools']['trimmomatic']['scratch'],
        mem = config['tools']['trimmomatic']['mem'],
        time = config['tools']['trimmomatic']['time']
    benchmark:
        TRIMMOMATICOUT + '{sample}/PAIREDEND/{fastq}.fastq.gz.benchmark'
    threads:
        config['tools']['trimmomatic']['paired']['threads']
    log:
        TRIMMOMATICOUT + '{sample}/PAIREDEND/{fastq}.fastq.gz.stdout.log'
    shell:
        ('{config[tools][trimmomatic][call]} ' +
        '{params.mode} ' +
        '{params.phred} ' +
        '-threads {threads} ' +
        '{input.forwardF} ' +
        '{input.reverseF} ' +
        '{output.forwardP} ' +
        '{output.forwardUP} ' +
        '{output.reverseP} ' +
        '{output.reverseUP} ' +
        'ILLUMINACLIP:{params.adapter}:{params.seedmismatches}:{params.palindrom}:{params.score}:{params.min_adapt_len}:{params.keep_both} ' +
        'SLIDINGWINDOW:{params.slidingwindow}:{params.minQual} ' +
        'LEADING:{params.leadingQual} ' +
        'TRAILING:{params.tailingQual} ' +
        'MINLEN:{params.minlen} ' +
        ' >{log}')

rule trimmomatic_single:
    input:
        fastq = TRIMMOMATICIN + '{sample}/SINGLEEND/{fastq}.fastq.gz',
        adapter = config['tools']['trimmomatic']['single']['adapterfile']
    output:
        TRIMMOMATICOUT + '{sample}/SINGLEEND/{fastq}.fastq.gz'
    params:
        slidingwindow = config['tools']['trimmomatic']['single']['slidingwindow'],
        phred = config['tools']['trimmomatic']['single']['phred'],
        mode = config['tools']['trimmomatic']['single']['mode'],
        minQual = config['tools']['trimmomatic']['single']['minQual'],
        leadingQual = config['tools']['trimmomatic']['single']['leadingQual'],
        tailingQual = config['tools']['trimmomatic']['single']['tailingQual'],
        seedmismatches = config['tools']['trimmomatic']['single']['seedmismatches'],
        palindrom = config['tools']['trimmomatic']['single']['palindrom'],
        score = config['tools']['trimmomatic']['single']['score'],
        minlen = config['tools']['trimmomatic']['single']['minlen'],
        lsfoutfile = TRIMMOMATICOUT + '{sample}/SINGLEEND/{fastq}_clipTrim.lsfout.log',
        lsferrfile = TRIMMOMATICOUT + '{sample}/SINGLEEND/{fastq}_clipTrim.lsferr.log',
        scratch = config['tools']['trimmomatic']['scratch'],
        mem = config['tools']['trimmomatic']['mem'],
        time = config['tools']['trimmomatic']['time']
    benchmark:
        TRIMMOMATICOUT + '{sample}/SINGLEEND/{fastq}.benchmark'
    threads:
        config['tools']['trimmomatic']['single']['threads']
    log:
        trimlog = TRIMMOMATICOUT + '{sample}/SINGLEEND/{fastq}_clipTrim.log',
        stdoutlog = TRIMMOMATICOUT + '{sample}/SINGLEEND/{fastq}_clipTrim.stdout.log'
    shell:
        ('{config[tools][trimmomatic][call]} ' +
        '{params.mode} ' +
        '{params.phred} ' +
        '-threads {threads} ' +
        '{input.fastq} ' +
        '{output} ' +
        'ILLUMINACLIP:{input.adapter}:{params.seedmismatches}:{params.palindrom}:{params.score} ' +
        'SLIDINGWINDOW:{params.slidingwindow}:{params.minQual} ' +
        'LEADING:{params.leadingQual} ' +
        'TRAILING:{params.tailingQual} ' +
        'MINLEN:{params.minlen} ' +
        ' > {log.stdoutlog}')


rule alignSingleFileStar:
    '''
        align all fastq files that are single
        '''
    input:
        fastq = getFastQForSampleSingleEnd
    output:
        ALIGNOUT + 'SINGLEEND/{sample}Aligned.out.bam',
        ALIGNOUT + 'SINGLEEND/{sample}Aligned.toTranscriptome.out.bam'
    params:
        fastq = getFastQForSampleSingleEndStarString,
        lsfoutfile = ALIGNOUT + 'SINGLEEND/{sample}_star.lsfout.log',
        lsferrfile = ALIGNOUT + 'SINGLEEND/{sample}_star.lsferr.log',
        scratch = config['tools']['star']['scratch'],
        mem = config['tools']['star']['mem'],
        time = config['tools']['star']['time'],
        params = config['tools']['star']['params'],
        index = config['tools']['star']['index'],
    benchmark:
        ALIGNOUT + 'SINGLEEND/{sample}_star.benchmark'
    threads:
        config['tools']['star']['threads']
    log:
        ALIGNOUT + 'SINGLEEND/{sample}_star.log'
    shell:
        '{config[tools][star][call]} --genomeDir {params.index} --readFilesIn {params.fastq} --runThreadN {threads} --outFileNamePrefix {params.out_prefix} {params.params} '


rule alignPairedFileStar:
    '''
        align all fastq files that are paired
        '''
    input:
        fastq = getFastQForSamplePairedEnd
    output:
        ALIGNOUT + 'PAIREDEND/{sample}Aligned.out.bam',
        ALIGNOUT + 'PAIREDEND/{sample}Aligned.toTranscriptome.out.bam'
    params:
        fastq = getFastQForSamplePairedEndStarString,
        lsfoutfile = ALIGNOUT + 'PAIREDEND/{sample}_star.lsfout.log',
        lsferrfile = ALIGNOUT + 'PAIREDEND/{sample}_star.lsferr.log',
        scratch = config['tools']['star']['scratch'],
        mem = config['tools']['star']['mem'],
        time = config['tools']['star']['time'],
        params = config['tools']['star']['params'],
        index = config['tools']['star']['index'],
        out_prefix = ALIGNOUT + 'PAIREDEND/{sample}'
    benchmark:
        ALIGNOUT + 'PAIREDEND/{sample}_star.benchmark'
    threads:
        config['tools']['star']['threads']
    log:
        ALIGNOUT + 'PAIREDEND/{sample}_star.log'
    shell:
        '{config[tools][star][call]} --genomeDir {params.index} --readFilesIn {params.fastq} --runThreadN {threads} --outFileNamePrefix {params.out_prefix} {params.params} '


rule samtools_sort_name_Single:
    input:
        bam = ALIGNOUT + 'SINGLEEND/{sample}Aligned.out.bam'
    output:
        bam = ALIGNOUT + 'SINGLEEND/{sample}Aligned.out.sortedName.bam'
    params:
        mem = config['tools']['samtools']['sortName']['mem'],
        scratch = config['tools']['samtools']['sortName']['scratch'],
        time = config['tools']['samtools']['sortName']['time'],
        params = config['tools']['samtools']['sortName']['params'],
        lsfoutfile = ALIGNOUT + 'SINGLEEND/{sample}.sortName.lsfout.log',
        lsferrfile = ALIGNOUT + 'SINGLEEND/{sample}.sortName.lsferr.log'
    threads:
        config['tools']['samtools']['sortName']['threads']
    shell:
        '{config[tools][samtools][call]} sort -n -O bam {params.params} --threads {threads} {input.bam} > {output.bam}'




rule samtools_sort_name_Paired:
    input:
        bam = ALIGNOUT + 'PAIREDEND/{sample}Aligned.out.bam'
    output:
        bam = ALIGNOUT + 'PAIREDEND/{sample}Aligned.out.sortedName.bam'
    params:
        mem = config['tools']['samtools']['sortName']['mem'],
        scratch = config['tools']['samtools']['sortName']['scratch'],
        time = config['tools']['samtools']['sortName']['time'],
        params = config['tools']['samtools']['sortName']['params'],
        lsfoutfile = ALIGNOUT + 'PAIREDEND/{sample}.sortName.lsfout.log',
        lsferrfile = ALIGNOUT + 'PAIREDEND/{sample}.sortName.lsferr.log'
    threads:
        config['tools']['samtools']['sortName']['threads']
    shell:
        '{config[tools][samtools][call]} sort -n -O bam {params.params} --threads {threads} {input.bam} > {output.bam}'


rule samtools_sort_pos_Single:
    input:
        bam = ALIGNOUT + 'SINGLEEND/{sample}Aligned.out.bam'
    output:
        bam = ALIGNOUT + 'SINGLEEND/{sample}Aligned.out.sortedCor.bam'
    params:
        mem = config['tools']['samtools']['sortCor']['mem'],
        scratch = config['tools']['samtools']['sortCor']['scratch'],
        time = config['tools']['samtools']['sortCor']['time'],
        params = config['tools']['samtools']['sortCor']['params'],
        lsfoutfile = ALIGNOUT + 'SINGLEEND/{sample}.sortCor.lsfout.log',
        lsferrfile = ALIGNOUT + 'SINGLEEND/{sample}.sortCor.lsferr.log'
    threads:
        config['tools']['samtools']['sortCor']['threads']
    shell:
        '{config[tools][samtools][call]} sort -O bam {params.params} --threads {threads} {input.bam} > {output.bam}'




rule samtools_sort_pos_Paired:
    input:
        bam = ALIGNOUT + 'PAIREDEND/{sample}Aligned.out.bam'
    output:
        bam = ALIGNOUT + 'PAIREDEND/{sample}Aligned.out.sortedCor.bam'
    params:
        mem = config['tools']['samtools']['sortCor']['mem'],
        scratch = config['tools']['samtools']['sortCor']['scratch'],
        time = config['tools']['samtools']['sortCor']['time'],
        params = config['tools']['samtools']['sortCor']['params'],
        lsfoutfile = ALIGNOUT + 'PAIREDEND/{sample}.sortCor.lsfout.log',
        lsferrfile = ALIGNOUT + 'PAIREDEND/{sample}.sortCor.lsferr.log'
    threads:
        config['tools']['samtools']['sortCor']['threads']
    shell:
        '{config[tools][samtools][call]} sort -O bam {params.params} --threads {threads} {input.bam} > {output.bam}'


rule picards_mark_PCR_duplicates_Single:
    input:
        bam=ALIGNOUT + 'SINGLEEND/{sample}Aligned.out.sortedCor.bam'
    output:
        bam=ALIGNOUT + 'SINGLEEND/{sample}Aligned.out.sorted.MarkDuplicates.bam'
    params:
        lsfoutfile = ALIGNOUT + 'SINGLEEND/{sample}.MarkDuplicates.lsfout.log',
        lsferrfile = ALIGNOUT + 'SINGLEEND/{sample}.MarkDuplicates.lsferr.log',
        scratch = config['tools']['picard']['markduplicates']['scratch'],
        mem = config['tools']['picard']['markduplicates']['mem'],
        time = config['tools']['picard']['markduplicates']['time'],
        params = config['tools']['picard']['markduplicates']['params']
    threads:
        config['tools']['picard']['markduplicates']['threads']
    benchmark:
        ALIGNOUT + 'SINGLEEND/{sample}Aligned.out.sorted.MarkDuplicates.benchmark'
    log:
        log = ALIGNOUT + 'SINGLEEND/{sample}Aligned.out.sorted.MarkDuplicates.log',
        metrics = ALIGNOUT + 'SINGLEEND/{sample}Aligned.out.sorted.MarkDuplicates.metrics'
    shell:
        ('{config[tools][picard][call]} ' +
        'MarkDuplicates ' +
        'INPUT={input.bam} ' +
        'OUTPUT={output.bam}  ' +
        '{params.params} ' +
        'TMP_DIR={TMPDIR} ' +
        'METRICS_FILE={log.metrics} ' +
        ' >{log.log} 2>&1')



rule picards_mark_PCR_duplicates_Paired:
    input:
        bam=ALIGNOUT + 'PAIREDEND/{sample}Aligned.out.sortedCor.bam'
    output:
        bam=ALIGNOUT + 'PAIREDEND/{sample}Aligned.out.sorted.MarkDuplicates.bam'
    params:
        lsfoutfile = ALIGNOUT + 'PAIREDEND/{sample}.MarkDuplicates.lsfout.log',
        lsferrfile = ALIGNOUT + 'PAIREDEND/{sample}.MarkDuplicates.lsferr.log',
        scratch = config['tools']['picard']['markduplicates']['scratch'],
        mem = config['tools']['picard']['markduplicates']['mem'],
        time = config['tools']['picard']['markduplicates']['time'],
        params = config['tools']['picard']['markduplicates']['params']
    threads:
        config['tools']['picard']['markduplicates']['threads']
    benchmark:
        ALIGNOUT + 'PAIREDEND/{sample}Aligned.out.sorted.MarkDuplicates.benchmark'
    log:
        log = ALIGNOUT + 'PAIREDEND/{sample}Aligned.out.sorted.MarkDuplicates.log',
        metrics = ALIGNOUT + 'PAIREDEND/{sample}Aligned.out.sorted.MarkDuplicates.metrics'
    shell:
        ('{config[tools][picard][call]} ' +
        'MarkDuplicates ' +
        'INPUT={input.bam} ' +
        'OUTPUT={output.bam}  ' +
        '{params.params} ' +
        'TMP_DIR={TMPDIR} ' +
        'METRICS_FILE={log.metrics} ' +
        ' >{log.log} 2>&1')

rule rnaseqc_Single:
    input:
        bam = ALIGNOUT + 'SINGLEEND/{sample}Aligned.out.sortedCor.bam',
        gtf = config['resources']['H_sapiens_hg38']['GTF_collapse']
    output:
        RNASEQCOUT + 'SINGLEEND/{sample}/{sample}Aligned.out.sortedCor.bam.metrics.tsv'
    params:
        lsfoutfile = RNASEQCOUT + '{sample}/{sample}.lsfout.log',
        lsferrfile = RNASEQCOUT + '{sample}/{sample}.lsferr.log',
        scratch = config['tools']['rnaseqc']['scratch'],
        mem = config['tools']['rnaseqc']['mem'],
        time = config['tools']['rnaseqc']['time'],
        params = config['tools']['rnaseqc']['params'],
        out_prefix = RNASEQCOUT + 'SINGLEEND/{sample}'
    benchmark:
        RNASEQCOUT + '{sample}/{sample}.benchmark'
    threads:
        config['tools']['rnaseqc']['threads']
    shell:
        '{config[tools][rnaseqc][call]} {input.gtf} {input.bam} {params.params} {params.out_prefix}'


rule rnaseqc_Paired:
    input:
        bam = ALIGNOUT + 'PAIREDEND/{sample}Aligned.out.sortedCor.bam',
        gtf = config['resources']['H_sapiens_hg38']['GTF_collapse']
    output:
        RNASEQCOUT + 'PAIREDEND/{sample}/{sample}Aligned.out.sortedCor.bam.metrics.tsv'
    params:
        lsfoutfile = RNASEQCOUT + 'PAIREDEND/{sample}/{sample}.lsfout.log',
        lsferrfile = RNASEQCOUT + 'PAIREDEND/{sample}/{sample}.lsferr.log',
        scratch = config['tools']['rnaseqc']['scratch'],
        mem = config['tools']['rnaseqc']['mem'],
        time = config['tools']['rnaseqc']['time'],
        params = config['tools']['rnaseqc']['params'],
        out_prefix = RNASEQCOUT + 'PAIREDEND/{sample}'
    benchmark:
        RNASEQCOUT + 'PAIREDEND/{sample}/{sample}.benchmark'
    threads:
        config['tools']['rnaseqc']['threads']
    shell:
        '{config[tools][rnaseqc][call]} {input.gtf} {input.bam} {params.params} {params.out_prefix}'




rule htseqcountSingle:
    input:
        bam = ALIGNOUT + 'SINGLEEND/{sample}Aligned.out.sortedName.bam',
        gtf = config['resources']['H_sapiens_hg38']['GTF']
    output:
        HTSEQCOUNTSOUT + 'SINGLEEND/{sample}.count.tsv'
    params:
        lsfoutfile = HTSEQCOUNTSOUT + 'SINGLEEND/{sample}.lsfout.log',
        lsferrfile = HTSEQCOUNTSOUT + 'SINGLEEND/{sample}.lsferr.log',
        scratch = config['tools']['htseqcount']['scratch'],
        mem = config['tools']['htseqcount']['mem'],
        time = config['tools']['htseqcount']['time'],
        params = config['tools']['htseqcount']['params']
    benchmark:
        HTSEQCOUNTSOUT + 'SINGLEEND/{sample}.benchmark'
    threads:
        config['tools']['htseqcount']['threads']
    shell:
        '{config[tools][htseqcount][call]} {params.params} {input.bam} {input.gtf} > {output}'




rule htseqcountPaired:
    input:
        bam = ALIGNOUT + 'PAIREDEND/{sample}Aligned.out.sortedName.bam',
        gtf = config['resources']['H_sapiens_hg38']['GTF']
    output:
        HTSEQCOUNTSOUT + 'PAIREDEND/{sample}.count.tsv'
    params:
        lsfoutfile = HTSEQCOUNTSOUT + 'PAIREDEND/{sample}.lsfout.log',
        lsferrfile = HTSEQCOUNTSOUT + 'PAIREDEND/{sample}.lsferr.log',
        scratch = config['tools']['htseqcount']['scratch'],
        mem = config['tools']['htseqcount']['mem'],
        time = config['tools']['htseqcount']['time'],
        params = config['tools']['htseqcount']['params']
    benchmark:
        HTSEQCOUNTSOUT + 'PAIREDEND/{sample}.benchmark'
    threads:
        config['tools']['htseqcount']['threads']
    shell:
        '{config[tools][htseqcount][call]} {params.params} {input.bam} {input.gtf} > {output}'


rule htseqcount2fpkmSingle:
    input:
        countF = HTSEQCOUNTSOUT + 'SINGLEEND/{sample}.count.tsv',
        gene_length = config['resources']['H_sapiens_hg38']['gene_length']
    output:
        HTSEQCOUNT2FPKMSOUT + 'SINGLEEND/{sample}.fpkm.tsv'
    params:
        lsfoutfile = HTSEQCOUNT2FPKMSOUT + 'SINGLEEND/{sample}.lsfout.log',
        lsferrfile = HTSEQCOUNT2FPKMSOUT + 'SINGLEEND/{sample}.lsferr.log',
        scratch = config['tools']['htseqcount2fpkm']['scratch'],
        mem = config['tools']['htseqcount2fpkm']['mem'],
        time = config['tools']['htseqcount2fpkm']['time']
    benchmark:
        HTSEQCOUNT2FPKMSOUT + 'SINGLEEND/{sample}.benchmark'
    threads:
        config['tools']['htseqcount2fpkm']['threads']
    shell:
        '{config[tools][htseqcount2fpkm][call]} {input.countF} {input.gene_length} {output}'


rule htseqcount2fpkmPaired:
    input:
        countF = HTSEQCOUNTSOUT + 'PAIREDEND/{sample}.count.tsv',
        gene_length = config['resources']['H_sapiens_hg38']['gene_length']
    output:
        HTSEQCOUNT2FPKMSOUT + 'PAIREDEND/{sample}.fpkm.tsv'
    params:
        lsfoutfile = HTSEQCOUNT2FPKMSOUT + 'PAIREDEND/{sample}.lsfout.log',
        lsferrfile = HTSEQCOUNT2FPKMSOUT + 'PAIREDEND/{sample}.lsferr.log',
        scratch = config['tools']['htseqcount2fpkm']['scratch'],
        mem = config['tools']['htseqcount2fpkm']['mem'],
        time = config['tools']['htseqcount2fpkm']['time']
    benchmark:
        HTSEQCOUNT2FPKMSOUT + 'PAIREDEND/{sample}.benchmark'
    threads:
        config['tools']['htseqcount2fpkm']['threads']
    shell:
        '{config[tools][htseqcount2fpkm][call]} {input.countF} {input.gene_length} {output}'


rule countReadsInFastq:
    input:
        fastq = '{sample}.fastq.gz'
    output:
        countF = '{sample}.count'
    params:
        scratch = config['tools']['countReadsInFastq']['scratch'],
        mem = config['tools']['countReadsInFastq']['mem'],
        time = config['tools']['countReadsInFastq']['time'],
        lsfoutfile = '{sample}.count.lsfout.log',
        lsferrfile = '{sample}.count.lsferr.log'
    threads:
        config['tools']['countReadsInFastq']['threads']
    shell:
        'zcat {input.fastq} | wc -l | awk \'{{print $1/4}}\' > {output.countF}'

rule runFastqcOnFastqs:
    input:
        fastq = '{sample}.fastq.gz'
    output:
        fastqc = touch('{sample}.fastqc.done')
    params:
        scratch = config['tools']['fastqc']['scratch'],
        mem = config['tools']['fastqc']['mem'],
        time = config['tools']['fastqc']['time'],
        params = config['tools']['fastqc']['params'],
        lsfoutfile = '{sample}.fastqc.lsfout.log',
        lsferrfile = '{sample}.fastqc.lsferr.log'
    threads:
        config['tools']['fastqc']['threads']
    shell:
        ('{config[tools][fastqc][call]} {params.params} -t {threads} {input.fastq}')

rule runFlagstat:
    input:
        bam = '{sample}.bam'
    output:
        flagstat = '{sample}.flagstat'
    params:
        mem = config['tools']['samtools']['flagstat']['mem'],
        scratch = config['tools']['samtools']['flagstat']['scratch'],
        time = config['tools']['samtools']['flagstat']['time'],
        lsfoutfile = '{sample}.flagstat.lsfout.log',
        lsferrfile = '{sample}.flagstat.lsferr.log'
    threads:
        config['tools']['samtools']['flagstat']['threads']
    shell:
        '{config[tools][samtools][call]} flagstat {input.bam} > {output.flagstat}'




rule runSamtoolsHeader:
    input:
        bam = '{sample}.bam'
    output:
        header = '{sample}.header'
    params:
        mem = config['tools']['samtools']['flagstat']['mem'],
        scratch = config['tools']['samtools']['flagstat']['scratch'],
        time = config['tools']['samtools']['flagstat']['time'],
        lsfoutfile = '{sample}.header.lsfout.log',
        lsferrfile = '{sample}.header.lsferr.log'
    threads:
        config['tools']['samtools']['flagstat']['threads']
    shell:
        '{config[tools][samtools][call]} view -H {input.bam} > {output.header}'




rule RandomForestPredictSingle:
    input:
        fpkm = HTSEQCOUNT2FPKMSOUT + 'SINGLEEND/{sample}.fpkm.tsv'
    output:
        type = PREDICTOUT + 'SINGLEEND/{sample}_RF_pred_tumor_type.csv',
        probs = PREDICTOUT + 'SINGLEEND/{sample}_RF_pred_probs.csv',
    params:
        lsfoutfile = PREDICTOUT + 'SINGLEEND/{sample}.RF.lsfout.log',
        lsferrfile = PREDICTOUT + 'SINGLEEND/{sample}.RF.lsferr.log',
        scratch = config['tools']['Random_forest_model_predict']['scratch'],
        model = config['tools']['Random_forest_model_predict']['model'],
        mem = config['tools']['Random_forest_model_predict']['mem'],
        time = config['tools']['Random_forest_model_predict']['time'],
        out_prefix = PREDICTOUT + 'SINGLEEND/{sample}'
    benchmark:
        PREDICTOUT + 'SINGLEEND/{sample}.RF.benchmark'
    threads:
        config['tools']['Random_forest_model_predict']['threads']
    shell:
        '{config[tools][Random_forest_model_predict][call]} -s {wildcards.sample} -d {input.fpkm} -m {params.model} -o {params.out_prefix}'


rule RandomForestPredictPaired:
    input:
        fpkm = HTSEQCOUNT2FPKMSOUT + 'PAIREDEND/{sample}.fpkm.tsv'
    output:
        type = PREDICTOUT + 'PAIREDEND/{sample}_RF_pred_tumor_type.csv',
        probs = PREDICTOUT + 'PAIREDEND/{sample}_RF_pred_probs.csv',
    params:
        lsfoutfile = PREDICTOUT + 'PAIREDEND/{sample}.RF.lsfout.log',
        lsferrfile = PREDICTOUT + 'PAIREDEND/{sample}.RF.lsferr.log',
        scratch = config['tools']['Random_forest_model_predict']['scratch'],
        model = config['tools']['Random_forest_model_predict']['model'],
        mem = config['tools']['Random_forest_model_predict']['mem'],
        time = config['tools']['Random_forest_model_predict']['time'],
        out_prefix = PREDICTOUT + 'PAIREDEND/{sample}'
    benchmark:
        PREDICTOUT + 'PAIREDEND/{sample}.RF.benchmark'
    threads:
        config['tools']['Random_forest_model_predict']['threads']
    shell:
        '{config[tools][Random_forest_model_predict][call]} -s {wildcards.sample} -d {input.fpkm} -m {params.model} -o {params.out_prefix}'


