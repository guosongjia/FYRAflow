
args<-commandArgs(TRUE)
(file=args[1])
(xmax=as.numeric(args[2]))
(ymax=as.numeric(args[3]))

#setwd("/home/scum/N/main/kat/test_guosong_strain")

#file <- "CBS5682.kat.hist"

df <- read.table(file,sep = " ",header = F,skip = 6)
pdf(paste(file,".pdf",sep=""))
plot(df$V1,df$V2,xlim=c(0,xmax),ylim=c(0,ymax),type="l",xlab="27-mer frequency",ylab="#. of distinct kmer",main = file)
dev.off()

