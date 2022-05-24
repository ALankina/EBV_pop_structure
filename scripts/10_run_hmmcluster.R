# this is a script that runs hmmcluster
# OC 2022
# Rscript scripts/10_run_hmmcluster.R data/EBV_MSA_April2022.fa data/NC_007605.1.fa NC007605 4 analysis/1-hmm/

args = commandArgs(trailingOnly=TRUE)
in_msa = args[1] # data/EBV_MSA_April2022.fa
ref_fa = args[2] # data/NC_007605.1.fa
ref_pattern = args[3] # NC_007605.1
threads = args[4]
#out_dir = args[5]

#setwd(out_dir)

# remotes::install_github("https://github.com/ucl-pathgenomics/hmmcluster")
library(hmmcluster)



get_regions_parallel(in_msa, ref_fa, ref_pattern, run_width = 200, run_java = "java", run_model = "AIC", run_cores = threads)

