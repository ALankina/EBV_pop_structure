#!/bin/bash
# qrsh -l h_vmem = 10g, tmem=10g -pe smp 16
# script to run all the EBV analyses






### run vars
msa_out=data/ebv_msa
threads=4



### set the environment for the cscluster
module load R # sets R and Rscript
snpsites=/share/apps/genomics/snp-sites-2.4.1/bin/snp-sites
admixture=/share/apps/genomics/admixture_linux-1.3.0/admixture
python3=/share/apps/python-3.9.0-shared/bin/python3
pyenv=/SAN/breuerlab/pathseq1/oc/envs/pybio/bin/activate
plink=/share/apps/genomics/plink-2.0/bin/plink2
iqtree=/share/apps/genomics/iqtree-1.6.12/bin/iqtree

# activate envs
source $pyenv





# dirs
mkdir analysis
mkdir analysis/1-msa
mkdir analysis/2-phylo
mkdir analysis/3-clustering
mkdir analysis/4-admix


#---------------------- 0 generate alignment
# this process can be very memory and time intensive and should be checked manually in any case
# please run scripts/00 beforehand


#---------------------- 1 find regions of structure
Rscript scripts/10_run_hmmcluster.R data/EBV_MSA_April2022.fa data/NC_007605.1.fa NC007605


#---------------------- 2 phylogeny
source scripts/20_phylo.sh



#---------------------- 3 clustering




#---------------------- 4 admixture
# msa -> vcf
${snpsites} -v ${msa_out}.fa -o analysis/4-admix/msa.vcf

# vcf -> bed format
${plink} --vcf analysis/4-admix/msa.vcf --make-bed --double-id --max-alleles 2 --rm-dup --allow-extra-chr --maf 0.01 --indep-pairwise 50 10 0.1 --out analysis/4-admix/msa

cd analysis/4-admix

for k in $(seq 1 10); 
do 
    echo $k
    # run admixture
    ${admixture} msa.bed ${k} -j${threads} --cv=20 \
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

cd ../..


### generate plots
# plot errors
Rscript scripts/40_plot_err_boxplot.R

# find the optimal K and plot the admix fractions


#---------------------- 5 LD