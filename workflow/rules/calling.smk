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
 
