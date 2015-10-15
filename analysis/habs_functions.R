# Functions specific to HABs formatting

##### FUNCTIONS ########

order_dates <- function(df){
  df$Date <- factor(df$Date, 
    levels = c("6/16","6/30","7/8","7/14","7/21",
      "7/29","8/4","8/11","8/18","8/25","9/2","9/8","9/15",
      "9/23","9/29","10/6","10/15","10/20","10/27"))
  return(df)
}



habs_format <- function(phy.long){

  # Reorder factor levels for dates
  phy.long <- order_dates(phy.long)

  # Set all samples from 8/11 to 0 because of bad sampling on that date
  w <- which(phy.long$Date == "8/11")
  phy.long$Abundance[w] <- 0
  
  # Reorder factor levers for stations
  phy.long$Station <- factor(phy.long$Station,
    levels = c("Toledo", "nearshore", "offshore"))
  
  # Reorder factor levels for Fraction. Change Fraction factor names
  phy.long$Fraction <- factor(phy.long$Fraction, 
    levels = c("CNA", "100LTR", "53LTR", "3NA", "22NA"))
  levels(phy.long$Fraction) <- c("Full", "100um", "53um","3um","0.22um")
  
  return (phy.long)

}



habs_ord_format <- function(physeq){

  # Refactor stations
  sample_data(physeq)$Station <- factor(sample_data(physeq)$Station,
    levels = c("WE2", "WE4", "WE12"))
  
  # Refactor months
  sample_data(physeq)$Month <- factor(sample_data(physeq)$Month,
    levels = c("May", "June", "July", "August", 
      "September", "October", "November"))

  return(physeq)

}



oligo_prep <- function(oligo.data, sample.data ){
  
  # merge oligotype data and sample data 
  oligo.merge <- merge(oligo.data, sample.data, by.x = "samples", by.y = "SampleID")
  
  # Scale oligotype data
  just.oligos <- subset(oligo.merge, 
                        select = -c(samples, Date, Station, Fraction, Sums))
  not.oligos <- subset(oligo.merge, 
                       select = c(samples, Date, Station, Fraction))
  oligos.scale <- just.oligos/oligo.merge$Sums
  oligo.remerge <- data.frame(oligos.scale, not.oligos)
  
  # Melt the oligotypes 
  oligo.melt <- melt(
    value.name = "count", 
    variable.name = "Oligotype", 
    id.vars = c("samples", "Station", "Date", "Fraction"),
    oligo.remerge
  )
  
  
  # Order dates
  oligo.melt <- order_dates(oligo.melt)
  
  # Reorder factor levers for stations
  oligo.melt$Station <- factor(oligo.melt$Station, 
    levels = c("WE2", "WE12", "WE4"))
  
  # Reorder factor levels for Fraction. Change Fraction factor names
  oligo.melt$Fraction <- factor(oligo.melt$Fraction, 
    levels = c("CNA", "100LTR", "53LTR", "3NA", "22NA"))
  
  levels(oligo.melt$Fraction) <- c("Full", "100um", "53um", "3um", ".22um")
  
  return(oligo.melt)
  
}

oligoplot <- function(oligo.melt, title){
  ggplot(oligo.melt,aes_string(x="Date",y="count",fill="Oligotype")) +
    facet_grid(Fraction~Station,scales = "free_y")+
    geom_bar(stat="identity")+
    scale_fill_manual(values = c("grey26","chartreuse3", "darkorange","royalblue", "red", 
                                 "cyan2", "darkgreen",phylum.colors,"white",phylum.colors,"white",phylum.colors)) + 
    scale_x_discrete(breaks=c("6/10","7/8",
                              "8/4","9/2", "10/6", "11/3"),
                     labels=c("Jun", "Jul",
                              "Aug", "Sep", "Oct", "Nov"),
                     drop = FALSE
    )+
    theme(axis.title.x = element_text(size=16, face="bold"),
          axis.text.x = element_text(angle=50, colour = "black", 
                                     vjust=1, hjust = 1, size=12,face="bold"),
          axis.text.y = element_text(colour = "black", size=12),
          axis.title.y = element_text(face="bold", size=16),
          plot.title = element_text(size = 18),
          legend.title = element_text(size=18),
          legend.text = element_text(size = 16),
          legend.position="right",
          strip.text.x = element_text(size=18, face="bold"),
          strip.text.y = element_text(size=16, face="bold"),
          strip.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill=NA, color="black"),
          axis.line = element_line(colour = "black"),
          panel.margin = unit(1, "lines")
    ) + 
    guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1))+   
    xlab("")+
    ylab("Relative Abundance \n")+
    ggtitle(title)
}

cyano_melt <- function (physeq) {
  physeq.long <- psmelt(physeq)
  physeq.long <- arrange(physeq.long, Samplenum, Species)
  
  # Reorder factor levels for dates
  physeq.long <- order_dates(physeq.long)
  
  # Reorder factor levers for stations
  physeq.long$Station <- factor(physeq.long$Station,
                                levels = c("WE2","WE12","WE4"))
  
  # Reorder factor levels for Fraction. Change Fraction factor names
  physeq.long$Fraction <- factor(physeq.long$Fraction, 
                                 levels = c("CNA", "100LTR", "53LTR", "3NA", "22NA"))
  levels(physeq.long$Fraction) <- c("Full", "100um", "53um","3um","0.22um")
  return(physeq.long)
}


# This function takes  a data frame in long format
# (such as the output from transform_and_melt)
# and produces a stacked barplot of community composition.
make_tax_barplot <- function(df, x, tax, facet, title, colors, xlab, ylab, relative, outline, guide) {
    ggplot(df, aes_string(x = x, y = "Abundance", fill = tax)) +
    facet_grid(reformulate(facet), scales="free_y") +
    geom_bar(stat = "identity", show_guide = guide) +
    scale_fill_manual(values = colors) +
    scale_x_discrete(
    breaks = c("6/10", "7/8", "8/4", "9/2", "10/6", "11/3"),
    labels = c("Jun", "Jul", "Aug", "Sep", "Oct", "Nov"), drop = FALSE
    ) +
    theme(
    axis.title.x = element_text(size = 16,face = "bold"),
    axis.text.x = element_text(angle = 50, colour = "black", vjust = 1, hjust = 1, size = 13),
    axis.text.y = element_text(colour = "black", size = 10),
    axis.title.y = element_text(face = "bold", size = 16),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    legend.position = "right",
    strip.text.x = element_text(size = 16, face = "bold"),
    strip.text.y = element_text(size = 16, face = "bold"),
    strip.background = element_rect(color = "white",size = 2, fill = NA),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
    panel.margin = unit(1, "lines")
    ) +
    guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(title) +
    if (relative) {
        geom_bar(
        stat = "identity",position = "fill", colour = outline,show_guide = FALSE
        )
    } else {
        geom_bar(stat = "identity",colour = outline,show_guide = FALSE)
    }
    
}

## run correlation analysis
test_htroph_corr <- function(htroph.df, cyano.df, i, physeq.scale){
  corr <- corr.test(cyano.df[,i], htroph.df, use = "everything", method = "spearman",adjust = "bonferroni")
  
  pvalue <- corr$p
  r <- corr$r
  
  sig <- which(pvalue < 0.05)
  otu.names <- colnames(pvalue)[sig]
  
  sig.tax <- data.frame(tax_table(physeq.scale))[otu.names,]
  sig.p <- pvalue[sig]
  sig.r <- r[sig]
  sig.dat<-data.frame(sig.p,sig.r,sig.tax)
  names(sig.dat)[1:2] <- c("pvalue","r")
  sig.dat <- arrange(sig.dat, r)
  sig.dat$Genus <- droplevels(sig.dat$Genus)
  return(sig.dat)
}

stackbar_theme <- theme(
  axis.title.x = element_text(size = 16,face = "bold"),
  axis.text.x = element_text(angle = 50, 
                             face = "bold", 
                             vjust = 1, 
                             hjust = 1, 
                             size = 15),
  axis.text.y = element_text(colour = "black", size = 10),
  axis.title.y = element_text(face = "bold", size = 16),
  plot.title = element_text(face = "bold", size = 22),
  legend.title = element_text(face = "bold", size = 16),
  legend.text = element_text(size = 14),
  strip.text.x = element_text(size = 16, face = "bold"),
  strip.text.y = element_text(size = 16, face = "bold"),
  strip.background = element_rect(color = "white",size = 2, fill = NA),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
  panel.margin = unit(1, "lines")
) 

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}


# Scales OTU relative abundance and changes the names of the OTUs to include
# Inputs are a phyloseq object and the taxonomic rank to combine with OTU # for the label 
scale_otu <- function(data, taxrank, datenames) {
  
  data.otu <- t(otu_table(data))
  
  # scale the realtive abundanec of each OTU
  data.otu <- scale(data.otu)
  
  # Change the taxa names to be taxrank + OTU number
  colnames(data.otu) <- paste(tax_table(data)[ ,taxrank], 
                              tax_table(data)[ ,"Species"])

  # Fix rownames and column names
  rownames(data.otu) <- datenames
  colnames(data.otu) <- gsub("[.]", " ", colnames(data.otu))
  
  # Return scaled otu table with taxa as rows
  return(t(data.otu))
  
}


# Computes euclidean distance between OTUs and performs 
# hierarchical clustering with complete linkage.
# Cuts tree into number of classes given by treecut
cluster_otu <- function(data, treecut) {
  
  d <- vegdist(t(data), method = "euclidean")
  
  data.clust <- hclust(d = d, method = "complete")
  
  cols <- cutree(data.clust, cut)
  
  data <- data.frame(data.otu)
  data <- data[ ,data.clust$order]
  
  newcols <- cols[data.clust$order]
  
  return(list(data = data, cols = newcols))

} 


