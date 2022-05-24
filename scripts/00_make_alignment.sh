#!/bin/bash -l
#$ -S /bin/bash
#$ -o /SAN/breuerlab/pathseq1/oc/22-7_ebv_popstruc
#$ -e /SAN/breuerlab/pathseq1/oc/22-7_ebv_popstruc
#$ -l h_rt=48:00:00
#$ -l tmem=200G,h_vmem=200G
#$ -N mafft
#$ -wd /SAN/breuerlab/pathseq1/oc/22-7_ebv_popstruc
#$ -V
#$ -R y
#$ -pe mpi 16



core_msa=data/core_ebv_msa.fa
new_to_align=data/ebv_raw.fa
msa_out=data/ebv_msa
mafft=/share/apps/genomics/mafft-7.453/bin/mafft



#---------------------- 1 get sequence data and align it
# download formatted EBV sequences from genbank, formayt, same to multi-fasta
$python3 scripts/02_accession_list_rename_download.py  > $new_to_align

# align sequences to core msa
# // todo this may be something we have to do before, high memory requirement
$mafft --keeplength --thread 16 --addfragments $new_to_align $core_msa > ${msa_out}.fa
