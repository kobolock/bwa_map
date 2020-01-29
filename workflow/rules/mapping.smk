configfile: "config/config.yaml"

def getFastq(wildcards):
    return config["samples"][wildcards.sample]

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
        "../env/mapping.yaml"
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
        "../env/mapping.yaml"
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
#     "workflow/env/samtools.yaml"
#  shell:
#     "samtools index {input}"

