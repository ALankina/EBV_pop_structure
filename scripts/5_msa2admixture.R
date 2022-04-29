# this converts the main msa into plink format


# inputs
msa = "data/EBV_MSA_April2022.fasta"




#------------------- 1 - MSA 2 PLINK
#fasta -> vcf
outfile = paste0("analysis/5_admix","/",basename(msa),".vcf")
system(paste0("snp-sites -v ", msa, " -o ", outfile))

#vcf -> plink
# removes local LD within 50bp as recommended in manual for admixture
# EBV seems to have multi-allelic variants, so forcing biallelic only
system(paste0("plink2 --vcf ", outfile," --make-bed --double-id --max-alleles 2 --allow-extra-chr --indep-pairwise 50 10 0.1 --out ",outfile))











#------------------- 2 - PLINK 2 ADMIXTURE
