# this converts the main msa into plink format


# req:
# admixture https://dalexander.github.io/admixture/download.html
# snp-sites

# inputs
msa = "data/EBV_MSA_April2022.fasta"
max_k = 10




#------------------- 1 - MSA 2 PLINK
vcffile = paste0("analysis/5_admix","/",basename(msa),".vcf")
plinkfile = paste0("analysis/5_admix","/",tools::file_path_sans_ext(basename(msa)))
#fasta -> vcf

system(paste0("snp-sites -v ", msa, " -o ", vcffile))

#vcf -> plink
# removes local LD within 50bp as recommended in manual for admixture
# EBV seems to have multi-allelic variants, so forcing biallelic only
#system(paste0("plink2 --vcf ", vcffile," --make-bed --double-id --max-alleles 2 --rm-dup --allow-extra-chr --maf 0.01 --indep-pairwise 50 10 0.1 --out ",plinkfile))
system(paste0("plink2 --vcf ", vcffile," --make-bed --double-id --max-alleles 2 --rm-dup --allow-extra-chr --maf 0.01 --out ",plinkfile))











#------------------- 2 - PLINK 2 ADMIXTURE
#runs admixture for random seed for each k of interest ad.dir times


for(k in c(1:max_k)){

    # ----------- run admix and store cross validation error in table
    command = paste0("admixture --cv=20 -s time -j 4  ",plinkfile,".bed ",k," | tee ", k,".error")
    command
    system(command,wait = T,ignore.stdout = T, ignore.stderr = T)
    t = system(paste0("grep -h CV 92_admix/even_labs/",sample,"/",sample,".",K,".out"), intern = T)
    t = stringr::str_split(t," ",simplify = T)[,4]
    t = readr::parse_number(t)
    even.df = rbind(even.df, data.frame(sample = sample, k = K, error = t)) # keep track of all runs


    # -------- plot this admix
    dat=read.table(paste0("92_admix/even_labs/",sample,"/",sample,".vcf.",K,".Q"))
    real_labs = c(paste0("africa_", 1:10), paste0("eu_", 11:20), paste0("Asia_", 1:2))
    dat$individual = real_labs
    dat = reshape2::melt(dat, id.vars = "individual")
    dat$group = str_split(dat$individual, "_",simplify = T)[,1]
    g = ggplot(dat, aes(fill=variable, y=value, x=individual)) +
        geom_bar(position="stack", stat="identity") +
        facet_grid(~fct_inorder(group), switch = "x", scales = "free", space = "free") +
        theme_classic() +
        labs(title = paste("clusters:", K)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_x_discrete(label = function(x) stringr::str_trunc(x, 30)) +
        scale_fill_gdocs(guide = FALSE) +
        theme(
        panel.spacing.x = unit(0.1, "lines"),
        #axis.text.x = element_blank(),
        axis.text=element_text(size=6),
        panel.grid = element_blank()
        ) +
        scale_fill_manual(values = pal)
    ggsave(filename = paste0("92_admix/even_labs/", sample,"/", sample, ".", K, ".png"),plot = g,device = "png")

}
  
  
  # find optimal
  write.csv(even.df, paste0("92_admix/even_labs/",sample,"/error.csv"))
  #optimK = which.min(even.df$error)
  #optimK = even.df$k[optimK]
  
  #file.copy(paste0("92_admix/even_labs/", sample,"/", sample, ".", optimK, ".png"),
  #          paste0("92_admix/even_labs/", sample,"/", sample, ".best.png"))
  
}





colnames(even.df) = c("k", "ad.dir", "cross_val_error")
for(ad.dir in 1:5){
  dat = even.df[even.df$ad.dir == ad.dir,]
  ggplot(data = dat, aes(x = k, y = cross_val_error)) +
    geom_line() +
    theme_minimal() +
    labs(title = paste(fasta, "- random seed run:", ad.dir))
  ggsave(paste0("92_admix/even_labs/random-seed", ad.dir, ".png"), device = "png")
}
  

