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



refactor_months_and_stations  <- function(physeq) {
  
  # Refactor stations
  sample_data(physeq)$Station <- factor(sample_data(physeq)$Station,
    levels = c("nearshore", "Toledo", "offshore"))

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
  ggplot(oligo.melt,aes_string(x = "Date", y = "count", fill = "Oligotype")) +
    facet_grid(Fraction~Station, scales = "free_y") +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("grey26","chartreuse3", "darkorange","royalblue", "red", 
                                 "cyan2", "darkgreen", phylum.colors,"white", phylum.colors, "white", phylum.colors)) + 
    scale_x_discrete(breaks = c("6/10","7/8", "8/4","9/2", "10/6", "11/3"),
                     labels = c("Jun", "Jul","Aug", "Sep", "Oct", "Nov"),
                     drop = FALSE
    ) +
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
  corr <- corr.test(cyano.df[ ,i], htroph.df, use = "everything", method = "spearman",adjust = "bonferroni")
  
  pvalue <- corr$p
  r <- corr$r
  
  sig <- which(pvalue < 0.05)
  otu.names <- colnames(pvalue)[sig]
  
  sig.tax <- data.frame(tax_table(physeq.scale))[otu.names,]
  sig.p <- pvalue[sig]
  sig.r <- r[sig]
  sig.dat <- data.frame(sig.p, sig.r, sig.tax)
  names(sig.dat)[1:2] <- c("pvalue", "r")
  sig.dat <- arrange(sig.dat, r)
  sig.dat$Genus <- droplevels(sig.dat$Genus)
  return(sig.dat)
}

stackbar_theme <- theme(
  axis.title.x = element_text(size = 16, face = "bold"),
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
  g <- ggplotGrob(plots[[1]] + theme(legend.position = "bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position = "none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}


# Scales OTU relative abundance and changes the names of the OTUs to include
# Inputs are a phyloseq object and the taxonomic rank to combine with OTU # for the label 
scale_otu <- function(physeq, taxrank) {
  
  # Scale relative abundance of each OTU
  # Assumes taxa are rows in original physeq object
  scaled.otu <- scale(t(otu_table(physeq)))
  
  # Change the taxa names to be taxrank + OTU number
  colnames(scaled.otu) <- paste(tax_table(physeq)[ ,taxrank], 
                                tax_table(physeq)[ ,"Species"])

  # Fix rownames and column names
  rownames(scaled.otu) <- sample_data(physeq)$Date
  colnames(scaled.otu) <- gsub("[.]", " ", colnames(scaled.otu))

  return(scaled.otu)
} 


# Computes euclidean distance between OTUs and performs 
# hierarchical clustering with complete linkage.
cluster_otu <- function(data, dist.method) {
  if (dist.method == "spearman") {
    data.sim <- cor(data, method = "spearman") 
    data.clust <- hclust(d = as.dist(1 - abs(data.sim)), method = "complete")
  } else {
    d <- dist(t(data), method = dist.method)
    data.clust <- hclust(d = d, method = "complete")
  }
  
}

heat_wrapper <- function(physeq, taxrank, dist.method, k) {
  
  # Scale OTUs
  scaled.otu <- scale_otu(
    physeq = physeq, 
    taxrank = taxrank
  )
  
  # hierarchical clustering of OTUs using complete linkage
  clustered.otu <- cluster_otu(scaled.otu, dist.method)
  
  # cut the dendrogram into groups
  groups <- cutree(clustered.otu, k = k)[clustered.otu$order]
  
  # Melt to long format
  melted.otu <- melt(scaled.otu, varnames = c("Date", "OTU"))
  
  melted.otu$OTU <- revalue(
    melted.otu$OTU, 
    c("Microcystis Otu00005" = "Microcystis", 
      "Synechococcus Otu00007" = "Syne. 1",
    "Pseudanabaena Otu00037" = "Pseuda.",
    "unclassified Otu00049" = "unclass",
    "Synechococcus Otu00044" = "Syne. 2")
    
  )
  
  
  # Order dates
  melted.otu <- order_dates(melted.otu)
  
  # Reorder levels of OTUs to match hierarchical clustering
  #otu.order <- levels(melted.otu$OTU)[clustered.otu$order]
  otu.order = c("Microcystis", "unclass", "Pseuda.", "Syne. 2", "Syne. 1")
  melted.otu$OTU <- factor(melted.otu$OTU, levels = otu.order)
  
  
  return(list(clusters = groups, data = melted.otu))
  

}


# heat_wrapper <- function(physeq, taxrank, dist.method, k) {
#   
#   # Scale OTUs
#   scaled.otu <- scale_otu(
#     physeq = physeq, 
#     taxrank = taxrank
#   )
#   
#   # hierarchical clustering of OTUs using complete linkage
#   clustered.otu <- cluster_otu(scaled.otu, dist.method)
#   
#   # cut the dendrogram into groups
#   groups <- cutree(clustered.otu, k = k)[clustered.otu$order]
#   
#   # Melt to long format
#   melted.otu <- melt(scaled.otu, varnames = c("Date", "OTU"))
#   
#   melted.otu$OTU <- revalue(
#     melted.otu$OTU, 
#     c("Microcystis Otu00005" = "Microcystis", 
#       "Synechococcus Otu00007" = "Syne. 1",
#       "Pseudanabaena Otu00037" = "Pseuda.",
#       "unclassified Otu00049" = "unclassified",
#       "Synechococcus Otu00044" = "Syne. 2",
#       "Synechococcus Otu00147" = "Syne. 3",
#       "Synechococcus Otu00177" = "Syne. 4")
#     
#   )
#   
#   
#   # Order dates
#   melted.otu <- order_dates(melted.otu)
#   
#   # Reorder levels of OTUs to match hierarchical clustering
#   #otu.order <- levels(melted.otu$OTU)[clustered.otu$order]
#   otu.order = c("Microcystis", "unclassified", "Pseuda.", "Syne. 4", "Syne. 2", "Syne. 1", "Syne. 3")
#   melted.otu$OTU <- factor(melted.otu$OTU, levels = otu.order)
#   
#   
#   return(list(clusters = groups, data = melted.otu))
#   
#   
# }


heat_log_wrapper <- function(physeq, taxrank, dist.method, k) {
  
  # Scale OTUs
  scaled.otu <- scale_otu_log(
    physeq = physeq, 
    taxrank = taxrank
  )
  
  # Only keep columns without NAs
  scaled.otu <- scaled.otu[ ,colSums(is.na(scaled.otu)) == 0]
  
  # hierarchical clustering of OTUs using complete linkage
  clustered.otu <- cluster_otu(scaled.otu, dist.method)
  
  # cut the dendrogram into groups
  groups <- cutree(clustered.otu, k = k)[clustered.otu$order]
  
  # Melt to long format
  melted.otu <- melt(scaled.otu, varnames = c("Date", "OTU"))
  
  # Order dates
  melted.otu <- order_dates(melted.otu)
  
  # Reorder levels of OTUs to match hierarchical clustering
  otu.order <- levels(melted.otu$OTU)[clustered.otu$order]
  melted.otu$OTU <- factor(melted.otu$OTU, levels = otu.order)
  
  return(list(clusters = groups, data = melted.otu))
  
  
}

# Inputs are a phyloseq object and the taxonomic rank to combine with OTU # for the label 
scale_otu_log <- function(physeq, taxrank) {
  
  # Scale relative abundance of each OTU by adding median, loging and then z-score
  # Assumes taxa are rows in original physeq object
  otu.t <- t(otu_table(physeq))
  otu.log <- apply(otu.t, 2, function(x) {log(x + median(x))})
  # hist(otu.log[ ,1]) yay this looks more normal - but still zero inflated
  scaled.otu <- scale(otu.log) 
  
  # Change the taxa names to be taxrank + OTU number
  colnames(scaled.otu) <- paste(tax_table(physeq)[ ,taxrank], 
                                tax_table(physeq)[ ,"Species"])
  
  # Fix rownames and column names
  rownames(scaled.otu) <- sample_data(physeq)$Date
  colnames(scaled.otu) <- gsub("[.]", " ", colnames(scaled.otu))
  
  return(scaled.otu)
}  



make_station_bray_heatmap <- function(physeq, self = FALSE, range, xlab, ylab, title) {
  
  
  # Calculate bray
  braydist <- as.matrix(phyloseq::distance(physeq, method = "bray"))
  
  # rows will be evens columns will be odds
  if (!self) {
    rows <- seq(from = 2, to = ncol(braydist), by = 2)
    cols <- seq(from = 1, to = ncol(braydist), by = 2)
    braydist <- braydist[rows, cols]
  }
  
  
  
  # rename the columns by date instead of sample number
  colnames(braydist) <- unique(sample_data(physeq)$Date)
  rownames(braydist) <- unique(sample_data(physeq)$Date)
  
  # Melt
  station.bray.melt <- melt(braydist, varnames = c("st1", "st2"))
  
  
  # ggplot
  g <- ggplot(station.bray.melt, aes(x = st1, y = st2) ) +
    geom_tile(aes(fill = value)) + 
    theme(
      axis.text.x = element_text(
        angle = 45, 
        vjust = 1, 
        hjust = 1, 
        size = 11, 
        face = "bold"
      ),
      axis.text.y = element_text(face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(size = 16, face = "bold")
    ) + 
    xlab(xlab) +
    ylab(ylab) + 
    ggtitle(title)
  
  if (!self) {
    g = g + scale_fill_gradient(low = "red", high = "white", limits = range)
  } else {
    g = g + scale_fill_gradient(low = "blue", high = "white", limits = range)
  }
  
  
}




nutplot_unpooled <- function(df, nutrient, title) {
  ggplot(df, aes_string(
    x = "Date",
    y = nutrient,
    group = "Station", 
    shape = "Station", 
    color = "Station")
  ) +
    geom_line(size = 0.8) +
    geom_point(size = 1.8) +
    ggtitle(title) + 
    xlab("") + ylab("ug/L") +       
    scale_color_manual(values = c("#E96446","#52361b","#87CEFA")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9, face = "bold"),     
          axis.title.y = element_text(size = 14),
          plot.title = element_text(size = 16, face = "bold"))
}


scale_otu_lineplots <- function(physeq) {
  
  # Scale relative abundance of each OTU
  # Assumes taxa are rows in original physeq object
  scaled.otu <- scale(t(otu_table(physeq)))
  
  # Make new data frame
  scaled.otu <- data.frame(
    Days = sample_data(physeq)$Days, 
    Date = sample_data(physeq)$Date, 
    Station = sample_data(physeq)$Station, 
    Abundance = scaled.otu, 
    row.names = NULL
  )
  colnames(scaled.otu)[4] <- "Abundance"
  scaled.otu <- order_dates(scaled.otu)
  
  scaled.otu$Station <- factor(scaled.otu$Station, levels = c("nearshore", "Toledo", "offshore"))
  
  return(scaled.otu)
} 

plot_otu_lineplots <- function(df, title) {
  ggplot(df, aes(x = Date, y = Abundance, group = Station, color = Station )) +
    geom_line(size = 1.2) + 
    scale_color_manual(values = c("#E96446", "#52361b", "#87CEFA")) +
    xlab("") + 
    ylab(title) +
    theme(
      axis.text.x = element_blank(),
      axis.title.y = element_text(size = 14),
      plot.margin = unit(c(1, 1.5, -0.2, 1.5), "cm")
    ) 
  
}