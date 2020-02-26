#Run10
#snakemake --snakefile who_flu_noInd.py -j 48 --config ifq=/Data-RAID5/iseq/raw/illumina/iSeq10/ out=/Data-RAID5/iseq/output/iSeq10analysis/

import subprocess, sys, os, glob 
from os.path import join
from os.path import basename
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

#############################################################################
#                                                                           #
# Description                                                               #
#                                                                           #
#############################################################################

#Pipeline designed to obtain consensus sequences from raw fastq data using denovo assemblies.
#The steps are: QualityTrim -> Denovo assembly -> blastn (using a custom flu database or the entire database) to identify the contigs -> Use the references identified to map the raw reads (bwa) -> Generate consensus


# Input parameters  ------------------------------------------------------------------------
#

IFQ = config["ifq"]

workspace= config["out"]

threads=4

#NCBI databases
#blastDB="/Data-RAID5/iseq/databases/influenza_ncbi_2Dec19/influenza.fna"
#incomplete entries discarded DB
#blastDB="/Data-RAID5/iseq/databases/influenza_ncbi_2Dec19_fixNames/influenza2.fna"
#Database weekly updated
blastDB="/home/migrau/scratch/databases/flu_ncbi/influenza.fna"
#Full DB
#blastDB="/Data-RAID5/iseq/databases/ncbi/nt"

segmentsFlu=["HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2"]

## Functions -------------------------------------------------------------------

def multi2singleFasta(content):
    block=[]
    res=""
    dataList= content.split("\n")
    for line in dataList:
        if line.startswith('>'):
            if block:
                res+=''.join(block) + '\n'
                block = []
            res+=line+ '\n'
        else:
            block.append(line.strip())

    if block:
        res+=''.join(block) + '\n'
    return res

#Function which replaces the id-records (`>blastID|segment|`) for the sample's name and segment (`>NXXX_HA`). 
def fixNames(fafile):
    res=""
    segment=""
    listSeg=["(HA)","(MP)","(NA)","(NP)","(NS)","(PA)","(PB1)","(PB2)","(NEP)","(NS1)","(M1)","(M2)"," NEP ", " HA ", " MP ", " NA ", " NP ", " NS "," PA "," PB1 "," PB2 "," NS1 "," M1 "," M2 "]
    segNumber={"PB2":"1","PB1":"2","PA":"3","HA":"4","NP":"5","NA":"6","MP":"7","NS":"8",}
    sampleName=fafile.split("/")[-1].split("_")[0]
    for index,record in enumerate(SeqIO.parse(fafile, 'fasta')):
        for g in listSeg:
            tmpDescr=record.description.replace(":","").replace(","," ").replace("(NEP)","(NS)").replace("(NS1)","(NS)").replace("(M1)","(MP)").replace("(M2)","(MP)").replace(" NEP "," NS ").replace(" NS1 "," NS ").replace(" M1 "," MP ").replace(" M2 "," MP ").replace("SEGMENT 1","(PB2)").replace("SEGMENT 2","(PB1)").replace("SEGMENT 3","(PA)").replace("SEGMENT 4","(HA)").replace("SEGMENT 5","(NP)").replace("SEGMENT 6","(NA)").replace("SEGMENT 7","(MP)").replace("SEGMENT 8","(NS)").replace("SEG 1","(PB2)").replace("SEG 2","(PB1)").replace("SEG 3","(PA)").replace("SEG 4","(HA)").replace("SEG 5","(NP)").replace("SEG 6","(NA)").replace("SEG 7","(MP)").replace("SEG 8","(NS)")
            if g in tmpDescr:
                segment=g.replace("(","").replace(")","").replace(" ","")
        #res+=">"+sampleName+"."+segNumber[str(segment)]+"\n"+str(record.seq)+"\n"
        res+=">"+sampleName+"_"+str(segment)+"\n"+str(record.seq)+"\n"
    return res

#Prepare final output separated by subtype
def groupBy (filesin, oFolder):
    listSeg=["HA","MP","NA","NP","NS","PA","PB1","PB2"]
    segments={}
    for file in files.split(" "):
        for index,record in enumerate(SeqIO.parse(file, 'fasta')):
            for seg in listSeg:
                if seg in record.id:
                    segments[seg]+=record.id+record.seq
    for k in segments.keys():
        with open (oFolder+k+".fa",'w') as fo:
            fo.write()
 
        
# Rules ------------------------------------------------------------------------
# 

#PAIRED READS
SAMPLES, PAIR= glob_wildcards(IFQ+"/{sample}_L001_{pair}_001.fastq.gz")

rule all:
    input:
        expand(workspace+'assemblies/{sample}/contigs.fasta', sample=SAMPLES),
        expand(workspace+'alignments/{sample}.noInd.bam', sample=SAMPLES),
        expand(workspace+'blast/{sample}/results.fasta', sample=SAMPLES),
        expand(workspace+'blast/{sample}/results.out', sample=SAMPLES),
        expand(workspace+'consensus/pre/{sample}.fa', sample=SAMPLES),
        expand(workspace+"qualtrim/{sample}.R1.paired.fastq", sample=SAMPLES),
        expand(workspace+"qualtrim/{sample}.R2.paired.fastq", sample=SAMPLES),
        expand(workspace+'alignments/post/{sample}.noInd.bam', sample=SAMPLES),
        expand(workspace+'consensus/post/{sample}.fa',  sample=SAMPLES),
        expand(workspace+'alignments/post/{sample}.ind.vcf', sample=SAMPLES),
        expand(workspace+'consensus/post/rename/{sample}.fa', sample=SAMPLES),
        expand(workspace+'circos/{sample}-covX.txt', sample=SAMPLES),
        expand(workspace+"circos/{sample}-coverage.png", sample=SAMPLES),
        expand(workspace+"circos/png/{sample}-coverage.png", sample=SAMPLES),
        expand(workspace+"circos/tmp/{sample}_depth.coverage",sample=SAMPLES)

#QUALITY FILTER
rule filter:
   input:
        faR1=expand(IFQ+"{{sample}}_L001_{pair}_001.fastq.gz", pair=["R1"]),
        faR2=expand(IFQ+"{{sample}}_L001_{pair}_001.fastq.gz", pair=["R2"])
   output:
        R1out=workspace+"qualtrim/{sample}.R1.paired.fastq",
        R2out=workspace+"qualtrim/{sample}.R2.paired.fastq",
        R1out_unpaired=workspace+"qualtrim/{sample}.R1.unpaired.fastq",
        R2out_unpaired=workspace+"qualtrim/{sample}.R2.unpaired.fastq"
   shell:"""
      
      java -jar /home/migrau/apps/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 {input.faR1} {input.faR2} {output.R1out} {output.R1out_unpaired} {output.R2out} {output.R2out_unpaired} ILLUMINACLIP:/home/migrau/apps/Trimmomatic-0.39/adapters/TruSeq2-SE.fa:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 MINLEN:100 HEADCROP:10 TRAILING:5
  """ 

#Denovo assembly using SPADES paired
# rule spades:
#     input:
#         R1out=workspace+"qualtrim/{sample}.R1.paired.fastq",
#         R2out=workspace+"qualtrim/{sample}.R2.paired.fastq"
#     output:
#         contigs=workspace+'assemblies/{sample}/contigs.fasta'
#     params:
#         folder=workspace+'assemblies/{sample}/',
#         spadesFolder=spadesFolder
#     shell:"""
#         python {params.spadesFolder}/spades.py --careful --pe1-1 {input.R1out} --pe1-2 {input.R2out} -t 8 -k 21,33,55,77 -o {params.folder}
#     """


#Denovo assembly using IVA
rule iva:
    input:
        R1out=workspace+"qualtrim/{sample}.R1.paired.fastq",
        R2out=workspace+"qualtrim/{sample}.R2.paired.fastq"
    output:
        contigs=workspace+'assemblies/{sample}/contigs.fasta'
    params:
        folder=workspace+'assemblies/{sample}/'
    shell:"""
        rm -r {params.folder}
        iva --seed_ext_min_cov 3 --seed_min_kmer_cov 10 -f {input.R1out} -r {input.R2out} {params.folder}
    """

#Blast contigs from IVA
rule blastn:
    input:
        contigs=workspace+'assemblies/{sample}/contigs.fasta'
    output:
        longContigs=workspace+'assemblies/{sample}/longContigs.fasta',
        blastout=workspace+'blast/{sample}/results.out'
    params:
        blast_database=blastDB
    shell:"""
        bioawk -c fastx '{{ if(length($seq) > 200) {{ print ">"$name; print $seq }}}}' {input.contigs} > {output.longContigs}
        blastn -db {params.blast_database} -query {output.longContigs} -out {output.blastout} -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle length mismatch sstart send' -evalue 0.01 -num_threads 2
    
    """

#Extract best case references from DB
rule extract_references:
    input:
        blastout=workspace+'blast/{sample}/results.out'
    output:
        blastFasta=workspace+'blast/{sample}/results.fasta'
    params:
        blast_database=blastDB,
        segmentsFlu=segmentsFlu
    shell:"""
        #awk '!seen[$1]++' {input.blastout} > {input.blastout}_best
        #Remove partial results first
        sed '/partial/d' {input.blastout} | awk '!seen[$1]++' > {input.blastout}_best
        python tools/bestPerGene_segment.py {input.blastout}_best > {input.blastout}_bestGenes
        cut -d$'\t' -f6 {input.blastout}_bestGenes | sort | uniq > {input.blastout}_uniq
        blastdbcmd -db {params.blast_database} -dbtype nucl -entry_batch {input.blastout}_uniq -outfmt "%f" -out {output.blastFasta}
    """

#Mapping reads to the new refs
rule mapping:
    input:
        blastFasta=workspace+'blast/{sample}/results.fasta',
        R1out=workspace+"qualtrim/{sample}.R1.paired.fastq",
        R2out=workspace+"qualtrim/{sample}.R2.paired.fastq"
    output:
        sam=workspace+'alignments/{sample}.sam',
        bam=workspace+'alignments/{sample}.bam',
        bamNoInd=workspace+'alignments/{sample}.noInd.bam'
    params:
        threads=threads
    shell:"""
        awk '{{ print toupper($0) }}' {input.blastFasta} > {input.blastFasta}_UP
        rgid=$(echo {input.R1out}  | md5sum | cut -d " " -f1)
        rgpu=${{rgid}}.PU
        #In case FLU not present
        touch {input.blastFasta}_UP
        bwa index {input.blastFasta}_UP
        if [ -s {input.blastFasta}_UP ]
        then
            bwa mem -R "@RG\\tID:$rgid\\tPL:illumina\\tPU:$rgpu\\tSM:{input.R1out}" {input.blastFasta}_UP -t {params.threads} {input.R1out} {input.R2out} > {output.sam} 
        else
            #In case FLU not present
            touch {output.sam}
        fi
        samtools view -bT {input.blastFasta}_UP {output.sam} | samtools sort > {output.bam}
        #Remove indels from bam file
        samtools view -h {output.bam} | awk '$1 ~ "^@" || $6 !~ "I|D"' | samtools view -b - > {output.bamNoInd}
    """

#PlanB
#kindel consensus --min-depth 1 {output.bam} > {output.preConsensus}

#Obtain consensus
rule consensus:
    input:
        bamNoInd=workspace+'alignments/{sample}.noInd.bam',
        blastFasta=workspace+'blast/{sample}/results.fasta',
    output:
        vcf=workspace+'alignments/{sample}.vcf',
        preConsensus=workspace+'consensus/pre/{sample}.fa',
    shell:"""
        if [ -s {input.blastFasta}_UP ]
        then
            samtools mpileup -A -uf {input.blastFasta}_UP {input.bamNoInd} | bcftools call -mv -Oz -o {output.vcf}
            tabix {output.vcf}
            cat {input.blastFasta}_UP | bcftools consensus {output.vcf} > {output.preConsensus}
        else
            touch {output.vcf}
            touch {output.preConsensus}
        fi
    """
#Repeat two last rules, polishing the final consensus
rule re_mapping:
    input:
        preConsensus=workspace+'consensus/pre/{sample}.fa',
        R1out=workspace+"qualtrim/{sample}.R1.paired.fastq",
        R2out=workspace+"qualtrim/{sample}.R2.paired.fastq"
    output:
        sam=workspace+'alignments/post/{sample}.sam',
        bam=workspace+'alignments/post/{sample}.bam',
        bamNoInd=workspace+'alignments/post/{sample}.noInd.bam',
        depths=workspace+'alignments/post/{sample}.noInd.depth'
    params:
        threads=threads
    shell:"""
        rgid=$(echo {input.R1out} | md5sum | cut -d " " -f1)
        rgpu=${{rgid}}.PU
        bwa index {input.preConsensus}
        if [ -s {input.preConsensus} ]
        then
            bwa mem -R "@RG\\tID:$rgid\\tPL:illumina\\tPU:$rgpu\\tSM:{input.R1out}" {input.preConsensus} -t {params.threads} {input.R1out} {input.R2out} > {output.sam}
        else
            touch {output.sam} 
        fi
        samtools view -bT {input.preConsensus} {output.sam} | samtools sort > {output.bam}
        samtools view -h {output.bam} | awk '$1 ~ "^@" || $6 !~ "I|D"' | samtools view -b - > {output.bamNoInd}
        samtools depth {output.bamNoInd} -a | awk '$3 == 0' > {output.depths}
    """

#Final consensus discarding ALL the indels from the bam file (normally the case of FLUA)
rule FINALconsensus:
    input:
        bamNoInd=workspace+'alignments/post/{sample}.noInd.bam',
        bam=workspace+'alignments/post/{sample}.bam',
        preConsensus=workspace+'consensus/pre/{sample}.fa',
        depths=workspace+'alignments/post/{sample}.noInd.depth'
    output:
        vcfIndel=workspace+'alignments/post/{sample}.ind.vcf',
        vcf=workspace+'alignments/post/{sample}.vcf',
        vcfilt=workspace+'alignments/post/{sample}.filt.vcf',
        postConsensus=workspace+'consensus/post/{sample}.fa',
        bed=workspace+'alignments/post/{sample}.bed',
        postConsensusMasked=workspace+'consensus/post/{sample}.masked.fasta',
        postConsensusMaskedRenamed=workspace+'consensus/post/{sample}.masked.renamed.fasta'
    shell:"""
        awk -v OFS='\t' '{{print $1,$2,$2+1}}' {input.depths} > {output.bed}
        if [ -s {input.preConsensus} ]
        then        
            bcftools mpileup -A -Ov -f {input.preConsensus} {input.bam} | bcftools call -mv -Ov -o {output.vcfIndel}
            bcftools mpileup -A -Ov -f {input.preConsensus} {input.bamNoInd} | bcftools call -mv -Ov -o {output.vcf}
            #Filter
            bgzip -c {output.vcf} > {output.vcf}.gz
            #bcftools view -e '(DP4[2]+DP4[3])<(DP4[0]+DP4[1])' {output.vcf}.gz -o {output.vcfilt}
            python tools/filterVCF_noInd.py {output.vcf} 20 > {output.vcfilt}
            sed -i '/^$/d' {output.vcfilt}
            bgzip -c {output.vcfilt} > {output.vcfilt}.gz
            tabix -p vcf {output.vcfilt}.gz
            bedtools maskfasta -fi {input.preConsensus} -bed {output.bed} -fo {output.postConsensusMasked}
            python tools/fix_bedtool_id.py {input.preConsensus} {output.postConsensusMasked} > {output.postConsensusMaskedRenamed}
            #Added -I in bcftools call to include IUPAC code
            cat {output.postConsensusMaskedRenamed} | bcftools consensus -I {output.vcfilt}.gz > {output.postConsensus}
        else
            touch {output.postConsensus}
            touch {output.vcfIndel}
            touch {output.vcf}
            touch {output.vcfilt}
            touch {output.postConsensusMasked}
            touch {output.postConsensusMaskedRenamed}
        fi
    """

#Rename output consensus
rule rename:
    input:
        fasta=workspace+'consensus/post/{sample}.fa'
    output:
        fasta=workspace+'consensus/post/rename/{sample}.fa'
    run:
        newFasta = fixNames(input.fasta)
        with open(output.fasta,'w') as fo:
            fo.write(newFasta)

#TODO group consensus by gene
#rule group:
#    input:
#        files=expand(workspace+'consensus/post/rename/{sample}.fa', sample= SAMPLES)
#    output:
#        segment=expand(workspace+'consensus/post/segment/{segment}.fa', segment= ["HA","MP","NA","NS","NP","PA","PB1","PB2"])
#    params:
#        segmentFolder=workspace+'consensus/post/segment/'
#    run:   
#        groupBy(input.files, params.segmentFolder)


#Group by subtype and create depth plots
#it requires some parts from Uma pipe and circos.
rule subTypeAndDepthPlots:
    input:
        fasta=workspace+'consensus/post/rename/{sample}.fa',
        R1=workspace+"qualtrim/{sample}.R1.paired.fastq",
        R2=workspace+"qualtrim/{sample}.R2.paired.fastq"
    output:
        covX=workspace+"circos/{sample}-covX.txt",
        png=workspace+"circos/{sample}-coverage.png",
        pngend=workspace+"circos/png/{sample}-coverage.png",
        tmpend=workspace+"circos/tmp/{sample}_depth.coverage"
    params:
        tmp=workspace+'consensus/post/rename/{sample}-covX.txt',
        alltxt=workspace+'consensus/post/rename/{sample}*.txt',
        circos=workspace+"circos/",
        circosTMP=workspace+"circos/tmp/",
        consensus=workspace+'consensus/post/rename/{sample}',
        tmpdepth=workspace+"circos/{sample}_depth.coverage",
        tmppng=workspace+"circos/{sample}-coverage.png",
        fasta=workspace+'consensus/post/'
    shell:"""
        #identify subtypes
        #python tools/groupSubtype.py {params.fasta}
        python tools/groupSubtype_1file.py {input.fasta}
        #depthPlots
        python /home/yi-mo/pipelines/miguel/paired/tools/createGraphfiles_Full_17Nov.py {input.fasta} {input.R1} {input.R2}
        mv {params.tmp} {output.covX}
        mv {params.alltxt} {params.circos}
        cp /home/migrau/pipelines/miguel/paired/tools/*conf {params.circos}
        if [ -s {params.tmpdepth} ] 
        then
            python /home/yi-mo/pipelines/miguel/paired/tools/generate_covplot.py {output.covX}
        else
            touch {params.tmpdepth}
            touch {params.tmppng}
            python /home/yi-mo/pipelines/miguel/paired/tools/generate_covplot.py {output.covX}
        fi
        mv {params.consensus}_depth.* {params.circosTMP}
        mv {params.consensus}_reftemp* {params.circosTMP}
        cp {output.png} {output.pngend}
    """
