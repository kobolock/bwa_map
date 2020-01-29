
##### load config and sample sheets #####
configfile: "config/config.yaml"

#SAMPLES = config["samples"][wildcards.sample]


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

include: "workflow/rules/mapping.smk"


def getTargets(wildcards):
    targets = []
    targets.extend(expand("data/sorted_reads/{sample}.bam", sample=config["samples"]))
    return targets
   
rule all:
    input: getTargets 

 
