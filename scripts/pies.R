rm(list=ls())

library("ape")
library ("gtools")
library("tidyverse")
library("tidyr")
library("reshape2")
library("ggrepel")
library("plotly")
library("ggpubr")
library("vegan")


setwd("~/Documents/LIDo_Y2/EBV/HMM")
assignment<- read.csv("gt_ass_T1T2.csv", header = T)
info <- read_csv("info.csv")
nreg <-length(unique(gt_ass_T1T2$region))

pie_maker_basic<- function(data,n){
  
  region <- data[which(data$region==n),]
  region <- left_join(region,info, by="sample")
  
  ### basic pie
  nallele1 <-length(which(region$genotype==1))
  nallele2 <-length(which(region$genotype==2))
  nallele3 <-length(which(region$genotype==3))
  
  num_seqs<-c(nallele1, nallele2, nallele3)
  allele<-c("nallele1", "nallele2", "nallele3")
  region_piechart <- data_frame(allele, num_seqs)
  
  pie_default <- ggplot(region_piechart, aes(x="", y=num_seqs, fill=allele)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0)+
    theme_minimal()+ theme(legend.position='none',
                           axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.line.y = element_blank(),
                           panel.grid = element_blank(),
                           axis.text.x = element_blank()) +
  xlab(" ") +ylab(paste0("Allele distribution for region ", n))
  
  return (pie_default)
}
pie_maker_type<- function(data,n){
  
  region <- data[which(data$region==n),]
  region <- left_join(region,info, by="sample")
  region_T1 <- data[which(region$type=="T1"),]
  region_T2 <- data[which(region$type=="T2"),]
  
  ### pie for type 1
  nallele1 <-length(which(region$genotype==1))
  nallele2 <-length(which(region$genotype==2))
  nallele3 <-length(which(region$genotype==3))
  
  num_seqs<-c(nallele1, nallele2, nallele3)
  allele<-c("nallele1", "nallele2", "nallele3")
  region_piechart <- data_frame(allele, num_seqs)
  pie_number<-paste0(n)
  
  pie_T1 <- ggplot(region_piechart, aes(x="", y=num_seqs, fill=allele)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0)+
    theme_minimal()+ theme(legend.position='none',
                           axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.line.y = element_blank(),
                           panel.grid = element_blank(),
                           axis.text.x = element_blank()) +
    xlab(" ") +ylab(paste0("Allele distribution for T1 seqs in region ", n))
  
  ### pie for type 2
  nallele1 <-length(which(region_T2$genotype==1))
  nallele2 <-length(which(region_T2$genotype==2))
  nallele3 <-length(which(region_T2$genotype==3))
  
  num_seqs<-c(nallele1, nallele2, nallele3)
  allele<-c("nallele1", "nallele2", "nallele3")
  region_piechart <- data_frame(allele, num_seqs)
  pie_number<-paste0(n)
  
  pie_T2 <- ggplot(region_piechart, aes(x="", y=num_seqs, fill=allele)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0)+
    theme_minimal()+ theme(legend.position='none',
                           axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.line.y = element_blank(),
                           panel.grid = element_blank(),
                           axis.text.x = element_blank()) +
    xlab(" ") +ylab(paste0("Allele distribution for T2 seqs in region ", n))
  
pies_types<-pie_T1+pie_T2

}
pie_maker_geo<- function(data,n){
  
  region <- data[which(data$region==n),]
  region <- left_join(region,info, by="sample")
  region_GE <- data[which(region$geo=="GE"),]
  region_GA<- data[which(region$geo=="GA"),]
  region_GF<- data[which(region$geo=="GF"),]
  
  ### pie for GE
  nallele1 <-length(which(region_GE$genotype==1))
  nallele2 <-length(which(region_GE$genotype==2))
  nallele3 <-length(which(region_GE$genotype==3))
  
  num_seqs<-c(nallele1, nallele2, nallele3)
  allele<-c("nallele1", "nallele2", "nallele3")
  region_piechart <- data_frame(allele, num_seqs)
  pie_number<-paste0(n)
  
  pie_GE <- ggplot(region_piechart, aes(x="", y=num_seqs, fill=allele)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0)+
    theme_minimal()+ theme(legend.position='none',
                           axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.line.y = element_blank(),
                           panel.grid = element_blank(),
                           axis.text.x = element_blank()) +
    xlab(" ") +ylab(paste0("Allele distribution for European seqs in region ", n))
  
  ### pie for GA
  nallele1 <-length(which(region_GA$genotype==1))
  nallele2 <-length(which(region_GA$genotype==2))
  nallele3 <-length(which(region_GA$genotype==3))
  
  num_seqs<-c(nallele1, nallele2, nallele3)
  allele<-c("nallele1", "nallele2", "nallele3")
  region_piechart <- data_frame(allele, num_seqs)
  pie_number<-paste0(n)
  
  pie_GA <- ggplot(region_piechart, aes(x="", y=num_seqs, fill=allele)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0)+
    theme_minimal()+ theme(legend.position='none',
                           axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.line.y = element_blank(),
                           panel.grid = element_blank(),
                           axis.text.x = element_blank()) +
    xlab(" ") +ylab(paste0("Allele distribution for Asian seqs in region ", n))
  
  ### pie for GF
  nallele1 <-length(which(region_GF$genotype==1))
  nallele2 <-length(which(region_GF$genotype==2))
  nallele3 <-length(which(region_GF$genotype==3))
  
  num_seqs<-c(nallele1, nallele2, nallele3)
  allele<-c("nallele1", "nallele2", "nallele3")
  region_piechart <- data_frame(allele, num_seqs)
  pie_number<-paste0(n)
  
  pie_GF <- ggplot(region_piechart, aes(x="", y=num_seqs, fill=allele)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0)+
    theme_minimal()+ theme(legend.position='none',
                           axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.line.y = element_blank(),
                           panel.grid = element_blank(),
                           axis.text.x = element_blank()) +
    xlab(" ") +ylab(paste0("Allele distribution for African seqs in region ", n))
  
  pies_geos<-pie_GE+pie_GA+pie_GF
  
}


for (i in 1:nreg){
  p_b <- pie_maker_basic(assignment, i)
  p_t <- pie_maker_type(assignment, i)
  p_g <- pie_maker_geo(assignment, i)
  ggsave(p_b, file=paste0("/plots/pies_T1T2/allele_distr_", i, ".png") )
  ggsave(p_t, file=paste0("/plots/pies_T1T2/allele_distr_TYPE_", i, ".png") )
  ggsave(p_g, file=paste0("/plots/pies_T1T2/allele_distr_GEO_", i, ".png") )
}