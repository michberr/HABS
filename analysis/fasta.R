r<-read.table("~/habs/genomics/raw_fastq/MiSeq_final_fastq/OTU_Building_97/OTUmap.uc",sep="\t")
mc1.index<-r$V10=="OTU_8"
mc1<-r[mc1.index,]
mc2.index<-r$V10=="OTU_8291"
mc2<-r[mc2.index,]

v4<-mc2$V4
t<-v4!=97.2
table(t) #110996 = 97.2, 1214= not 97.2
not97.2<-mc2[t,]


derep<-read.table("~/habs/genomics/raw_fastq/MiSeq_final_fastq/OTU_Building_97/FinalOTUs.fasta",sep="\t")
odds<-derep[ seq(1, nrow(derep), by = 2),] 
evens<-derep[seq(2,nrow(derep),by=2),]
fasta<-data.frame(odds,evens)

foo <- data.frame(do.call('rbind', strsplit(as.character(fasta$odds),';;',fixed=TRUE)))

newfasta<-data.frame(foo$X1,evens)




