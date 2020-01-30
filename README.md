# Test du tutorial de l'outil snakemake
=======================================

Lors de l'installation initiale sur le cluster, la version de snakemake était la
v.3.13.0. Or, pour l'utilisation de la librairie pandas et du module validate, il
est nécessaire d'avoir une version >= à la v.5.1.

Pour créer un environnement avec une version spécifique d'un outil :  

$ conda create -n smk_env snakemake=5.10.0
