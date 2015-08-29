#Load libraries
library(phyloseq)
packageVersion("phyloseq")
library(ggplot2)
library(ape)
library(vegan)
library(plyr)
library(scales)
library(grid)
library(reshape2)
theme_set(theme_bw())

# Set working directory
setwd("~/habs/genomics/r/data/")

# Uparse import of good run
otu<-read.csv("uparse.otutable",sep="\t")
tax<-read.csv("uparse60.tax",sep="\t")
map<-read.csv("metadata.txt",sep="\t")
rownames(map)<-map$SampleID
phyOTU<-otu_table(otu,taxa_are_rows=TRUE)
tax<-as.matrix(tax)
phyTax<-(tax_table(tax))
phySam<-sample_data(map)
physeq<-phyloseq(phyOTU,phyTax,phySam)

good<-subset_taxa(physeq,Class!="Chloroplast")
good<-subset_taxa(physeq,Kingdom=="Bacteria")
good<-subset_taxa(good,Family!="mitochondria")

# Uparse import of bad run 
botu<-read.csv("uparse.otutable",sep="\t")
btax<-read.csv("besttax.txt",sep="\t")
rownames(btax)<- rownames(botu)
bmap<-read.csv("newmap.txt",sep="\t")
rownames(bmap)<b-map$SampleID
bphyOTU<-otu_table(otu,taxa_are_rows=TRUE)
btax<-as.matrix(btax)
bphyTax<-(tax_table(btax))
bphySam<-sample_data(bmap)
bphyseq<-phyloseq(bphyOTU,bphyTax,bphySam)
