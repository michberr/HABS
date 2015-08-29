
##### FUNCTIONS ########

# Scale reads using mcmurdie holmes method to a given library size of n
scale_reads <- function(physeq, n) {
  physeq.scale <-
    transform_sample_counts(physeq, function(x) {
      (n * x/sum(x))
    })
  otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}

# Order dates correctly in a data frame
order_dates <- function(df) {
  df$Date <-
    factor(
      df$Date, levels = c(
        "5/27","6/10","6/16","6/30","7/8","7/14","7/21",
        "7/29","8/4","8/11","8/18","8/25","9/2","9/8","9/15",
        "9/23","9/29","10/6","10/15","10/20","10/27","11/3"
      )
    )
  return(df)
}

# Function to run adonis test on a physeq object and a variable from metadata 
doadonis <- function(physeq, category) {
  physeq.scale <- scale_reads(physeq, min(sample_sums(physeq)))
  bdist <- phyloseq::distance(physeq.scale, "bray")
  col <- as(sample_data(physeq), "data.frame")[ ,category]
  
  # Adonis test
  adonis.bdist <- adonis(bdist ~ col)
  print(adonis.bdist)
  
  # Homogeneity of dispersion test
  betatax = betadisper(bdist,col)
  p = permutest(betatax)
  print("The homogeneity of dispersion:")
  print(p$tab)
}


# Function to make a pcoa plot with bray-cutis distance where color variable is discrete
make_pcoa <- function(physeq) {
  physeq.scale <- scale_reads(physeq, min(sample_sums(physeq)))
  sample_data(physeq.scale) <- order_dates(sample_data(physeq.scale))
  physeq.bray.nmds <- ordinate(physeq.scale, method = "PCoA","bray")
}


# This function takes a phyloseq object, agglomerates OTUs to the desired taxonomic rank, 
# prunes out OTUs below a certain relative proportion in a sample (ie 1% ) 
# and melts the phyloseq object into long format.
transform_and_melt <- function(physeq, taxrank, prune) {
  # Agglomerate all otu's by phylum level
  physeq_taxrank <- tax_glom(physeq,taxrank = taxrank)
  
  # Create a new phyloseq object which removes taxa from each sample below the prune parameter
  physeq_taxrank.prune <- transform_sample_counts(physeq_taxrank,function(x) {x/sum(x)})
  otu_table(physeq_taxrank.prune)[otu_table(physeq_taxrank.prune) < prune] <- 0
  physeq_taxrank.prune <- prune_taxa(taxa_sums(physeq_taxrank.prune) > 0, physeq_taxrank.prune)
  
  # Set all samples from 8/11 to 0 because of bad sampling on that date
  w <- which(sample_data(physeq_taxrank.prune)$Date == "8/11")
  otu_table(physeq_taxrank.prune)[,w] <- 0
  
  # Melt into long format and sort by samplenum and phylum
  physeq_taxrank.long <- psmelt(physeq_taxrank.prune)
  physeq_taxrank.long <- arrange(physeq_taxrank.long, Samplenum, Phylum)
  
  # Reorder factor levels for dates
  physeq_taxrank.long <- order_dates(physeq_taxrank.long)
  
  # Reorder factor levers for stations
  physeq_taxrank.long$Station <- factor(physeq_taxrank.long$Station,
                                        levels = c("WE2","WE12","WE4"))
  
  # Reorder factor levels for Fraction. Change Fraction factor names
  physeq_taxrank.long$Fraction <- factor(physeq_taxrank.long$Fraction, 
                                         levels = c("CNA", "100LTR", "53LTR", "3NA", "22NA"))
  levels(physeq_taxrank.long$Fraction) <- c("Full", "100um", "53um","3um","0.22um")
  
  # Return long data frame
  return(physeq_taxrank.long)
  
}

# This function takes  a data frame in long format 
# (such as the output from transform_and_melt) 
# and produces a stacked barplot of community composition.
make_tax_barplot <- function(df, x, y, tax, facet, title, colors, xlab, ylab, relative, outline, guide) {
  ggplot(df, aes_string(x = x, y = y, fill = tax)) +
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
      panel.margin = unit(1, "lines"),
      panel.background = element_rect(fill = "#a8a8a8")
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

myord <- function(physeq, fraction, method, distance, title){
  physeq.scale <- scale_sub(physeq, fraction)
  physeq.ordinate <- ordinate(physeq.scale, method = method, distance = distance)
  # Plot the ordination
  myordplot(physeq.scale, physeq.ordinate, title)
}

scale_sub <- function (physeq, fraction) {
  
  # Subset samples to a given filter fraction
  # Note: this dumb code is because you can't pass parameters to 
  # subset_samples inside functions
  subidx <- as.vector(sample_data(physeq)[,"Fraction"] == fraction)
  phy.prune <- prune_samples(subidx,physeq)
  phy.prune <- prune_taxa(taxa_sums(phy.prune)>0, phy.prune)
  
  # Scale reads to minimum sample size
  phy.scale <- scale_reads(phy.prune, min(sample_sums(phy.prune)))
  
  # Refactor stations
  sample_data(phy.scale)$Station <- factor(sample_data(phy.scale)$Station,
                                           levels = c("WE2", "WE4", "WE12"))
  
  # Refactor months
  sample_data(phy.scale)$Month <- factor(sample_data(phy.scale)$Month,
                                         levels = c("May", "June", "July", "August", "September", "October", 
                                                    "November"))
  
  return(phy.scale)
  
}

myordplot <- function (physeq, ordination, title){
  plot_ordination(
    physeq = physeq, 
    ordination = ordination, 
    color = "Month", 
    shape = "Station", 
    axes = c(1,2)) + 
    geom_point(
      aes(colour = Month), 
      alpha = 0.7, 
      size = 4) +
    geom_point(
      colour="grey90", 
      size = 1.5) + 
    scale_color_manual(
      values = c("#a65628", "red", "#ffae19", "#4daf4a", 
                 "#1919ff", "darkorchid3", "magenta")) +
    theme(plot.margin = unit(c(2, 2, 2, 2), "lines")) +
    ggtitle(title)
  
}

oligo_prep <- function(oligo.data, sample.data ){
  
  # merge oligotype data and sample data
  oligo.merge <- merge(oligo.data, sample.data, by.x="samples", by.y="SampleID")
  
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
                                 "cyan2", "darkgreen",phylum.colors)) + 
    scale_x_discrete(breaks=c("6/10","7/8",
                              "8/4","9/2","10/6","11/3"),
                     labels=c("Jun", "Jul",
                              "Aug", "Sep",  "Oct","Nov"),
                     drop=FALSE
    )+
    theme(axis.title.x = element_text(size=16,face="bold"),
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
          axis.line = element_line(colour = "black")
    )+
    guides(fill = guide_legend(reverse= TRUE,keywidth=1,keyheight=1))+   
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
