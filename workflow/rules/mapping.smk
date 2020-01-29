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

rule bwa_map:
    input:
        ref = "/data2/fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa",
        fastq = getFastq 
    output:
        temp("data/mapped_reads/{sample}.bam")
    params:
        name="alignement", 
        nthread= 8,
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    log:
        "data/logs/bwa_mem/{sample}.log"
    threads: 8
    conda:
        "workflow/env/mapping.yaml"
    shell:
        "bwa mem -R '{params.rg}' -t {threads} {input.ref} {input.fastq} | "
        "samtools view -Sb - > {output}"
 
rule samtools_sort:
    input:
        "data/mapped_reads/{sample}.bam"
    output:
        "data/sorted_reads/{sample}.bam"
    params:
        name="sorted",
        nthread=4
    conda:
        "workflow/env/mapping.yaml"
    shell:
        "samtools sort -o {output} {input}"   

#rule samtools_index:
#  input:
#     "sorted_reads/{sample}.bam"
#  output:
#     "sorted_reads/{sample}.bam.bai"
#  params:
#     name="samtools_index",
#     nthread=4
#  conda:
#     "env/samtools.yaml"
#  shell:
#     "samtools index {input}"

