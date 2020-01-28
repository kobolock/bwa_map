
##### load config and sample sheets #####
configfile: "config/config.yaml"

def getFastq(wildcards):
  return config["samples"][wildcards.sample]


## A partir d'un fichier .tsv
##===========================
#import pandas as pd
#
#samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
#
#units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
#units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])
#
#def getFastq(wildcards):
#  return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
##===========================

rule all:
  input:
     "calls/all.vcf"

rule bwa_map:
  input:
     ref = "/data2/fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa",
#     R1=lambda wildcards: config["samples"][wildcards.sample]["R1"],
#     R2=lambda wildcards: config["samples"][wildcards.sample]["R2"]
     fastq = getFastq 
  output:
     temp("mapped_reads/{sample}.sam")
  params:
     name="alignement", 
     nthread= 8,
     rg=r"@RG\tID:{sample}\tSM:{sample}"
  log:
     "logs/bwa_mem/{sample}.log"
  threads: 8
  conda:
     "env/bwa_mem.yaml"
  shell:
     "(bwa mem -R '{params.rg}' -t {threads} {input.ref} {input.fastq} > {output}) 2> {log}"
 
rule samtobam:
  input:
     "mapped_reads/{sample}.sam"
  output:
     "mapped_reads/{sample}.bam"
  params:
     name="SamToBam",
     nthread=4
  conda:
     "env/samtools.yaml"
  shell:
     "samtools view -Sb {input} > {output}"   

rule samtools_sort:
  input:
     "mapped_reads/{sample}.bam"
  output:
     "sorted_reads/{sample}.bam"
  params:
     name="sorted",
     nthread=4
  conda:
     "env/samtools.yaml"
  shell:
     "samtools sort -T sorted_reads/{wildcards.sample} -O bam {input} >  {output}"   

rule samtools_index:
  input:
     "sorted_reads/{sample}.bam"
  output:
     "sorted_reads/{sample}.bam.bai"
  params:
     name="samtools_index",
     nthread=4
  conda:
     "env/samtools.yaml"
  shell:
     "samtools index {input}"

rule bcftools_call:
  input:
     fa="/data2/fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa",
     bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
     bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
  output:
     "calls/all.vcf"
  params:
     name="bcftools",
     nthread=4
  log:
     "logs/bcftools_call.log"
  conda:
     "env/samtools.yaml"
  shell:
     "(samtools mpileup -g -f {input.fa} {input.bam} | "
     "bcftools call -mv - > {output}) 2> {log}"

#rule stats:
#  input:
#     "calls/all.vcf"
#  output:
#     "plots/quals.svg"
#  params:
#     name="stats",
#     nthread=2
#  conda:
#     "env/stats.yaml"
#  script:
#     "scripts/plot-quals.py"
 
