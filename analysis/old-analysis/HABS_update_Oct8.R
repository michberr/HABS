---
  title: 'Lake Erie HABs Time Series'
output: html_document
---
  


```{r load libraries, warning = FALSE, message = FALSE, echo = FALSE}
#Load libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(scales)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(magrittr)
source("habs_functions.R")
source("miseqR.R")
library(mvtsplot)
theme_set(theme_bw())
```

```{r mothur import, echo = FALSE}

# Import mothur files and sample metadata
sharedfile = "mothur/allhabs.shared"
taxfile = "mothur/allhabs.taxonomy"
mapfile = "other/habs_metadata.csv"

mothurdata = import_mothur(mothur_shared_file = sharedfile, 
                           mothur_constaxonomy_file = taxfile)

# Add the OTU number as a column in the taxonomy file
tax_table(mothurdata) <- cbind(tax_table(mothurdata), 
                               row.names(tax_table(mothurdata)))

# Rename the taxonomy columns
colnames(tax_table(mothurdata)) <- c("Kingdom", "Phylum", "Class", "Order", 
                                     "Family", "Genus", "Species")

# Import sample metadata
map <- read.csv(mapfile)
map <- sample_data(map)
rownames(map) <- map$SampleID

# Merge mothurdata object with sample metadata
moth_merge <- merge_phyloseq(mothurdata, map)

# Filter out non-samples (i.e. water,mock) and samples from intensive cruises.
# Also prune out taxa which were only present in removed samples
moth_chloro <- subset_samples(moth_merge,Type == "sample")
moth_chloro <- prune_taxa(taxa_sums(moth_chloro) > 0, moth_chloro)

# Filter out non-bacteria, chloroplasts and mitochondria
moth_chloro %>%
  subset_taxa(Kingdom == "Bacteria" & Family != "mitochondria" & Class != "Chloroplast") -> moth_good

phylum.colors <- c("#CBD588", "#5F7FC7", "orange", "#DA5724", "#89C5DA",
                   "#508578", "#CD9BCD", "#74D944", "#D7C1B1", "#AD6F3B", "#673770","#D14285",
                   "#689030", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F",
                   "#D1A33D", "#8A7C64", "#599861"
)

```

```{r echo=FALSE}
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
```