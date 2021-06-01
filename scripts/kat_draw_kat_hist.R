
args<-commandArgs(TRUE)
(file=args[1])
(step=args[2])
(xmax=as.numeric(args[3]))
(ymax=as.numeric(args[4]))

#setwd("/home/scum/N/main/kat/test_guosong_strain")

#file <- "CBS5682.kat.hist"

df <- read.table(file,sep = " ",header = F,skip = 6)
title_list = unlist(strsplit(file,split = "/"))
my_length = length(title_list)
title = title_list[[my_length]]
pdf(paste(file,"_",step,".pdf",sep=""))
plot(df$V1,df$V2,xlim=c(0,xmax),ylim=c(0,ymax),type="l",xlab="27-mer frequency",ylab="#. of distinct kmer",main = title)
dev.off()

