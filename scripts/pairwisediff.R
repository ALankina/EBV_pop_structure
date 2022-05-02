rm(list = ls())

library(ape)
install.packages("tidyverse")
library("tidyverse")
install.packages("reshape2")
library("reshape2")
install.packages("ggrepel")
library("ggrepel")
#these packages are not necessary, but I like interactive clustering where I can click on points and see the sample names! 
install.packages("ggpubr")
library("ggpubr")

##read file
data <- read.dna("input_fastas/EBV_MSA_April2022_nojunk.fasta",format="fasta") 
check <- dist.dna(data, model="TN93",as.matrix = TRUE,pairwise.deletion=TRUE) #TN93 model bc literature
check1 <- melt(check) #melt df with pairwise differences 
check1_c <- na.omit(check1) #get rid of NAs if there
check1_c_n <- check1_c %>% filter(Var1!=Var2) #get rid of comparisons between themselves i.e. Sequences1 vs Sequences - which gives 0 as difference 

check1_c_n$difference <-0
for (i in 1:nrow(check1_c_n)){
  print(i)
  if ((str_detect(check1_c_n$Var1[i], "_T2"))&(str_detect(check1_c_n$Var2[i], "_T2")))
  {check1_c_n$difference[i]<-"Within Type 2"}
  else{
  if ((str_detect(check1_c_n$Var1[i], "_T1"))&(str_detect(check1_c_n$Var2[i], "_T1")))
  {check1_c_n$difference[i]<-"Within Type 1"}
  else
  {check1_c_n$difference[i]<-"Between types"}}
}  

col_vector1 <- c("mediumorchid3", "green", "red")
  
#Plot 
p2<- ggplot(check1_c_n,aes(x=difference, y=value, fill=difference)) +
geom_violin(alpha=0.5) + geom_boxplot(width=0.2, fill="Ivory")+
scale_fill_manual(values=col_vector1) +
theme_pubr() + #from package ggpubr
ylab("Pairwise proportion differences") +
xlab(" ")

ggsave(filename = "DistDifferences_T1T2.pdf", width = 7, height = 5, p2,useDingbats=FALSE)