################# Load Libraries ################
setwd("~/chabs/miseq_may2015/analysis")
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(car)
library(RColorBrewer)
library(gtable)
source("habs_functions.R")
theme_set(theme_bw())


############################ Load Data #############################

## Nutrient data
# Import metadata file with nutrients, pigments and toxin
nutrient <- read.csv("other/nutrients.csv")


# Format nutrient data
nutrient <-
  nutrient %>%
  filter(!(Date %in% c("5/27", "6/10", "11/3"))) %>%
  order_dates() %>%
  mutate(H2O2 = H2O2/1000) %>%
  select(Date, Station, Phycocyanin, ParMC, Nitrate, SRP, H2O2) %>%
  melt(id.vars = c("Date", "Station"), variable.name = "Nutrient")


#########
## Metagenome  and qPCR data
  
  
## qPCR data
# Import qPCR data
qPCR <- read.csv("~/chabs/qPCR/qPCR_data_nored.csv", stringsAsFactors = FALSE)

# Format qPCR data
qPCR <-
  qPCR %>%
    mutate(Percent = Percent * 0.01) %>%
    rename("qpcr" = Percent) 

# Import metagenome data
mcy <- read.csv("other/mcy.csv", stringsAsFactors = FALSE)
names(mcy)[20] <- "mcyD"
names(mcy)[21] <- "mcyE"
mcy$Station <- gsub("L", "", mcy$Station)


# Format fraction
Fractionlabel <-  
  mcy$Description %>%
  as.character() %>%
  strsplit(split = "-") %>%
  unlist() %>%
  recode(
    "'100LTR' = '100um'; 
    '53LTR' = '53um';
    '3NA' = '3um'"
  ) %>%
  matrix(ncol = 2, byrow = TRUE)

mcy$Fraction <- Fractionlabel[ ,2]


mcy.joined <-
  mcy %>%
  select(Station, Date, Fraction, mcyD, mcyE) %>%
  full_join(qPCR, by = c("Date", "Station", "Fraction")) %>%
  melt(
    variable.name = "Gene", 
    id.vars = c("Station", "Date", "Fraction")
  ) %>%
  order_dates()

mcy.joined$Gene <- factor(mcy.joined$Gene, levels = c("mcyD", "qpcr", "mcyE"))



############################ Plotting functions ################################

plot_one_nutrient <- function(nutrient.df, nutrient, ylabel, i) {

  # Filter by nutrient
  data <- 
    nutrient.df %>%
    filter(Nutrient %in% nutrient)
  
  g <- 
    ggplot(data, aes(x = Date, y = value, group = Nutrient)) +
    geom_point(
      aes(colour = Nutrient), 
      alpha = 0.6, 
      size = 2.25
    ) + 
    geom_point(
      colour = "grey90", 
      size = 1.2,
      alpha = 0.6
    ) + 
    geom_line(
      aes(color = Nutrient),
      size = 0.7
    ) +
    xlab("") + 
    ylab(ylabel) +
    scale_color_manual(values = colors[i]) + 
    theme(
      axis.title.y = element_text(size = 10, face = "bold", vjust = 0.7),
      legend.key.size = unit(0.4, "cm"),
      axis.text.x = element_blank(),
      axis.ticks = element_blank() 
    )
  
  #If its the first plot, make the margins smaller on bottom
  if(i == 1) {
    g <- g + 
      theme(
        plot.margin = unit(c(1, 1.5, -0.3, 1.5), "cm")
      ) 
  # Else make margins on top and bottom smaller
  } else {
    g <- g +
      theme(
        plot.margin = unit(c(-0.3, 1.5, -0.3, 1.5), "cm")
      )
  }
  
  # If its the 5th plot, show the legend, otherwise don't
  if (i == 5) {
    g <- g +
      guides(
        color = guide_legend(title = NULL)
      ) 
  } else {
    g <- g +
      guides(color = FALSE)
  }
  
}



make_mcyplot <- function(data, fraction, ylab, i, sample.dates) {
  
  data.fraction <-
    data %>%
    filter(Fraction == fraction)
  
  g <- 
    ggplot(data.fraction, 
           aes(x = Date, y = value, group = Gene, fill = Gene)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_x_discrete(
      breaks = levels(data[ , "Date"]),
      labels = levels(data[ , "Date"]),
      drop = FALSE
    ) + 
    scale_y_continuous(lim = c(0, 1)) +
    scale_fill_manual(NULL, labels = c("metaG mcyD", "qPCR mcyD", "metaG mcyE"), values = c("#1f78b4", "#e31a1c", "#a6cee3")) + 
    ylab(ylab) +
    xlab("") + 
    theme(
      axis.text.x = element_text(
        angle = 45, 
        vjust = 1, 
        hjust = 1,
        size = 8, 
        face = "bold",
        color = sample.dates
      ),
      axis.title.y = element_text(size = 10, face = "bold", vjust = 0.7),
      legend.key.size = unit(0.4, "cm"),
      plot.margin = unit(c(-0.3, 1.5, -0.3, 1.5), "cm")
    ) +
    guides(
      fill = guide_legend(title = NULL)
    ) 
  
}


################## Set Variables #################

# Set station
S = "WE12"
title = "Station 2: nearshore"



#1f78b4 light blue
#a6cee3 dark blue 
#e31a1c   red 


sample.dates.we2 = 
  c(
    rep("black", 4), 
      "#e31a1c", 
      "#e31a1c",
      "#1f78b4",
      "black",
      "#e31a1c",
      "black",
      "#e31a1c",
      "black",
      "#e31a1c",
      "#e31a1c",
      "#1f78b4",
      "#e31a1c",
      "#e31a1c",
      "#1f78b4",
      "black"
  )

     
sample.dates.we12 = 
  c(
    rep("black", 2),
        "#1f78b4", # 7/8 
        "black",   # 7/14
        "#e31a1c", # 7/21
        "#e31a1c",
        "#1f78b4",
        "black",
    rep("#e31a1c", 6), # 8/18 - 9/23
        "#1f78b4",
        "#e31a1c",
        "#e31a1c",
        "#1f78b4",
        "black"
  )


          

#################### Nutrient plots ##################


data.station <- nutrient %>% 
  filter(Station == S)

colors <- brewer.pal(7, "Set1")
colors <- colors[c(1, 3, 4, 7, 2, 5)]


# Nitrogen
nitr <- plot_one_nutrient(
  nutrient.df = data.station, 
  nutrient = "Nitrate",
  ylabel = "Nitrate",
  i = 1
) + ggtitle("Station 12")

# Phosphorus
srp <- plot_one_nutrient(
  nutrient.df = data.station, 
  nutrient = "SRP",
  ylabel = "SRP",
  i = 2
)

h2o2 <- plot_one_nutrient(
  nutrient.df = data.station, 
  nutrient = "H2O2",
  ylabel = "H2O2",
  i = 3
)

phyco <- plot_one_nutrient(
  nutrient.df = data.station, 
  nutrient = "Phycocyanin",
  ylabel = "PC",
  i = 4
)

parmc <- plot_one_nutrient(
  nutrient.df = data.station, 
  nutrient = c("ParMC", "DissMC"),
  ylabel = "MC",
  i = c(5,6)
)

########## MCY plots ####################

mcy.station <- 
  mcy.joined %>%
    filter(Station == S)


Full <- make_mcyplot(
  data = mcy.station, 
  fraction = "CNA", 
  ylab = "Full",
  i = 1,
  sample.dates = sample.dates.we2
)

F100 <- make_mcyplot(
  data = mcy.station, 
  fraction = "100um", 
  ylab = "ratio toxin:16s",
  i = 2,
  sample.dates = sample.dates.we2
)

F53 <- make_mcyplot(
  data = mcy.station, 
  fraction = "53LTR", 
  ylab = "53um",
  i = 3,
  red.dates = R
)

F3 <- make_mcyplot(
  data = mcy.station, 
  fraction = "3NA", 
  ylab = "3um",
  i = 4,
  red.dates = R
)


########### Grobs ##################


g1 = ggplotGrob(nitr)
g2 = ggplotGrob(srp)
g3 = ggplotGrob(h2o2)
g4 = ggplotGrob(phyco)
g5 = ggplotGrob(parmc)
g6 = ggplotGrob(Full)
g7 = ggplotGrob(F100)
g8 = ggplotGrob(F53)
g9 = ggplotGrob(F3)


g1 <- gtable_add_cols(g1, g7$widths[6])
g2 <- gtable_add_cols(g2, g7$widths[6])
g3 <- gtable_add_cols(g3, g7$widths[6])
g4 <- gtable_add_cols(g4, g7$widths[6])
#g5 <- gtable_add_cols(g5, g7$widths[6])
#g7 <- gtable_add_cols(g7, g7$widths[6])
#g8 <- gtable_add_cols(g8, g7$widths[6])
#g9 <- gtable_add_cols(g9, g7$widths[6])

# Draw everything
#grid.draw(rbind(g1, g2, g3, g4, g5, g6, g7, g8, g9, size = "last"))

gg <- grid.draw(rbind(g1, g2, g3, g4, g5, g7, size = "last"))

ggsave(gg,  "")



##################

g10 <- gtable_add_cols(g10, g7$widths[6])
grid.draw(rbind(g7,g10, size = "first"))
