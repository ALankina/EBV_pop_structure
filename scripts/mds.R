rm(list = ls())

library(ape)
install.packages("tidyverse")
library("tidyverse")
install.packages("reshape2")
library("reshape2")
install.packages("ggrepel")
library("ggrepel")
#these packages are not necessary, but I like interactive clustering where I can click on points and see the sample names! 
install.packages("plotly")
library("plotly")
install.packages("htmlwidgets")
library("htmlwidgets")
#this is for plotting, just because I like the style
install.packages("ggpubr")
library("ggpubr")

##read file
data <- read.dna("input_fastas/EBV_MSA_April2022_nojunk.fasta",format="fasta") 
check <- dist.dna(data, model="TN93",as.matrix = TRUE,pairwise.deletion=TRUE) #TN93 model bc literature
fit <- cmdscale(check,eig=TRUE, k=4)
mydata <- as.data.frame(fit$points)
mydata$sample <- row.names(mydata) 

#read metadata - this is important if you want to colour your points by group (geography, disease etc etc)
info <- read_csv("info.csv")
#Join mydata and metadata together 
mydatainfo <- left_join(mydata,info, by="sample")

#palletes
col_vector2 <- c('royalblue2', 'firebrick1','gold','cyan','darkorange','seagreen', 'mediumorchid3', 'gray75')
shape_vector <-c(17,19,15)

#Plotting MDS 
p <- ggplot(mydatainfo,aes(x=V1,y=V2,text=sample, colour=geo, shape=type)) +
  geom_point(size=4, alpha=0.4) +
  theme_pubr() +
  scale_shape_manual(values=shape_vector, name="EBV type") +
  scale_colour_manual(values = col_vector2, name="Geo tag") +
  #theme(legend.position = "None") +
  xlab('MDS1') + ylab("MDS2")
p

#for interactive plot
p1 <- ggplotly(p,tooltip = "text")
p1

#save plot as pdf 
ggsave(filename = "wholegenome_MDS.pdf", width = 7, height = 5, p,useDingbats=FALSE)