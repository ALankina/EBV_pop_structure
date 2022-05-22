#!/bin/bash
# qrsh -l h_vmem = 10g, tmem=10g -pe smp 16
# script to run all the EBV analyses






### run vars
generate_msa_from_core=T
core_msa=data/core_ebv_msa.fa
new_to_align=data/ebv_raw.fa
msa_out=data/ebv_msa
threads=4



### set the environment for the cscluster
module load R # sets R and Rscript
mafft=/share/apps/genomics/mafft-7.453/bin/mafft
snpsites=/share/apps/genomics/snp-sites-2.4.1/bin/snp-sites
admixture=/share/apps/genomics/admixture_linux-1.3.0/admixture
python3=/share/apps/python-3.9.0-shared/bin/python3
pyenv=/SAN/breuerlab/pathseq1/oc/envs/pybio/bin/activate
plink=/share/apps/genomics/plink-2.0/bin/plink2
iqtree=/share/apps/genomics/iqtree-1.6.12/bin/iqtree

# activate envs
source $pyenv




#---------------------- 0 get sequence data and align it
# download formatted EBV sequences from genbank, formayt, same to multi-fasta
$python3 scripts/02_accession_list_rename_download.py  > $new_to_align

# align sequences to core msa
echo $mafft  --keeplength --thread $threads --add $new_to_align $core_msa > ${msa_out}.fa

# dirs
mkdir analysis
mkdir analysis/1-msa
mkdir analysis/2-phylo
mkdir analysis/3-clustering
mkdir analysis/4-admix


#---------------------- 1 find regions of structure




#---------------------- 2 phylogeny
source scripts/20_phylo.sh



#---------------------- 3 clustering




#---------------------- 4 admixture
# msa -> vcf
${snpsites} -v ${msa_out}.fa -o analysis/4-admix/msa.vcf

# vcf -> bed format
${plink} --vcf analysis/4-admix/msa.vcf --make-bed --double-id --max-alleles 2 --rm-dup --allow-extra-chr --maf 0.01 --out analysis/4-admix/msa

cd analysis/4-admix

for k in $(seq 1 3); 
do 
    echo $k
    # run admixture
    ${admixture} msa.bed ${k} --cv=20 \
    -s time `#random seed` | tee ${k}.log
    #-j ${threads} `#threads` \

    # get error
    if [[ "${k}" = "1" ]]
    then
        echo a
        echo k err > all.err
        echo ${k} `grep -h CV ${k}.log | grep -Po '0.[0-9]..'` >> all.err
    else
        echo b
        echo ${k} `grep -h CV ${k}.log | grep -Po '0.[0-9]..'` >> all.err
    fi

done

# generate plots


#---------------------- 5 LD