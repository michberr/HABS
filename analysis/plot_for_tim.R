# Data for Tim

timdata<-subset_samples(moth_good,Date %in% c("8/4/14","9/15/14"))

moth.long<-transform_and_melt(physeq = timdata,taxrank = "Phylum",prune = .01)

make_tax_barplot(df = moth.long,x = "Date",y="Abundance",tax = "Phylum",facet = "Fraction~Station",
                 title = "",
                 colors=phylum.colors,
                 xlab="",
                 ylab="Relative Abundance (taxa>1%)")


cyano.tim<-subset_samples(cyano,Date %in% c("8/4/14","9/15/14"))

# We will 
moth.long<-transform_and_melt(physeq = cyano.tim,taxrank = "Phylum",prune = .01)

names(moth.long)[30]<-"Genus"

brewcols<-c(brewer.pal(5,"Set3"),"grey")

mothcols<-c(Anabaena=brewcols[4],"Microcystis"=brewcols[1],Pseudanabaena=brewcols[2],Synechococcus="#fdb462",unclassified=brewcols[5],"Non-Cyano"="#4D4D4D")


moth.long$Genus <- factor(moth.long$Genus, levels = c("Anabaena","Microcystis","Pseudanabaena","Synechococcus","unclassified","Non-Cyano"))

moth.long<-arrange(moth.long,Genus)

make_tax_barplot(df = moth.long,x = "Date",y="Abundance",tax = "Genus",facet = "Fraction~Station",
                 title = "",
                 colors=mothcols,
                 xlab="",
                 ylab="Relative Abundance (taxa>1%)")

