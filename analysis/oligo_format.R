

setwd("~/habs/genomics/data/oligotype/")

# Read in oligotype counts file
oligo <- read.csv("matrix_counts.txt",sep="\t")

# remove rows with less than sum of 10 and control samples
#oligo.pruned<-oligo[rowSums(oligo[,-1])>9,]
#oligo.pruned<-oligo.pruned[1:150,]
oligo.pruned <-oligo[1:303,]
row.names(oligo.pruned)<-oligo.pruned$samples


# Scale oligotype data by dividing by sample read sums
sums <- sample_sums(moth_good)
m<- merge(oligo.pruned,sums,by="row.names")
scale.oligo <-m[,c(3:5)]/m[,6]
rownames(scale.oligo) <- m$Row.names

#  Merge the sample metadata with oligotype data
d<-data.frame(map)
dd<-subset(d, Type=="sample" & Intensive=="n")
dd <- dd[,c("SampleID","Station","Date","Fraction")]
oligo.merge<-merge(dd,scale.oligo,by.x="SampleID",by.y="row.names",all=TRUE)

write.table(oligo.merge, "oligotype/oligo_formatted.txt")

# Fix NAs
missingsamp <- c("")


wcna <- which(o$Date=="9/29" & o$Fraction=="CNA")
w100 <- which(o$Date=="9/29" & o$Fraction=="100LTR")
w53 <- which(is.na(o$TG) & o$Fraction=="53LTR")
w3 <- which(is.na(o$TG) & o$Fraction=="3NA")
w22 <- which(is.na(o$TG) & o$Fraction=="22NA")

