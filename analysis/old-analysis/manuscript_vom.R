```{r cyanolineplots, fig.height = 12, fig.width = 9}

cyano.otus <- c("Otu00005", "Otu00007", "Otu00037", "Otu00044", "Otu00049")
titles <- c("Microcystis", "Synechococcus", "Pseudanabaena", "Synechococcus", "Unclassified")
plot.list <- list()


for (i in 1:length(cyano.otus)) {
  # Subset taxa to OTU of interest
  cyano.otu <- 
    cyanos %>%
    subset_taxa(Species == cyano.otus[i])
  
  # Z-score scaling
  cyano.scale <- scale_otu_lineplots(cyano.otu)
  
  # Plot all three stations on one plot
  plot <- plot_otu_lineplots(df = cyano.scale, title = titles[i])
  
  # If its not the first plot remove the legend
  if (i != 1) {plot <- plot + theme(legend.position = "none")}
  
  # If its not the last plot, remove x-axis
  if (i != length(cyano.otus)) {plot <- plot + theme(axis.text.x = element_blank())}
  
  plot.list[[i]] <- plot
}

plot.list <- sapply(plot.list, function(x) {ggplotGrob(x)} )

for (i in 2:5) {
  plot.list[[i]] <- gtable_add_cols(plot.list[[i]], plot.list[[1]]$widths[6])
}

plots <- rbind(plot.list[[1]], plot.list[[2]], plot.list[[3]], plot.list[[4]], plot.list[[5]], size = "first")

#dev.off()
grid.draw(plots)

```
**Figure Legend** These are all the cyanobacterial OTUs with a mean relative abundance > 0.005 (75 reads per sample). I applied a z-score transformation to each OTU and plotted its abundance over time.


```{r eval = FALSE}

pH <- 
  ggplot(braypc.htrophs.sampledata, aes(x = PC1, y = pH, group = Station, color = Station)) +
  geom_point(size = 4) + 
  scale_color_brewer(palette = "Set2") +
  xlab("") + 
  theme(
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 14), 
    legend.position = "none",
    plot.margin = unit(c(-0.2, 1.5, -0.2, 1.5), "cm")
  )

turb <- 
  ggplot(braypc.htrophs.sampledata, aes(x = PC1, y = Turbidity, group = Station, color = Station)) +
  geom_point(size = 4) + 
  scale_color_brewer(palette = "Set2") +
  theme(
    axis.text.x = element_text(
      angle = 45, 
      vjust = 1, 
      hjust = 1, 
      size = 12, 
      face = "bold"
    ), 
    axis.title.y = element_text(size = 14),
    legend.position = "none",
    plot.margin = unit(c(-0.2, 1.5, 1, 1.5), "cm")
  )




ggplot(braypc.sampledata, aes(x = PC1, y = Turbidity, group = Station, color = Station)) +
  geom_point(size = 4) + 
  scale_color_brewer(palette = "Set2") +
  theme(
    axis.text.x = element_text(
      angle = 45, 
      vjust = 1, 
      hjust = 1, 
      size = 12, 
      face = "bold"
    ), 
    axis.title.y = element_text(size = 14),
    legend.position = "none",
    plot.margin = unit(c(-0.2, 1.5, 1, 1.5), "cm")
  )

ggplot(braypc.sampledata, aes(x = PC1, y = pH, group = Station, color = Station)) +
  geom_point(size = 4) + 
  scale_color_brewer(palette = "Set2") +
  theme(
    axis.text.x = element_text(
      angle = 45, 
      vjust = 1, 
      hjust = 1, 
      size = 12, 
      face = "bold"
    ), 
    axis.title.y = element_text(size = 14),
    legend.position = "none",
    plot.margin = unit(c(-0.2, 1.5, 1, 1.5), "cm")
  )

```

```{r include = FALSE}
# library(longitudinalData)
# 
# get_frech_dist <- function(otu.df) {
#   
#   otu.we2 <- filter(otu.df, Station == "nearshore")
#   otu.we12 <- filter(otu.df, Station == "Toledo")
#   otu.we4 <- filter(otu.df, Station == "offshore")
# 
#   otu.we2.we12 <- distFrechet(Px = otu.we2$Days, Py = otu.we2$Abundance, Qx = otu.we12$Days, Qy = otu.we12$Abundance, FrechetSumOrMax = "sum")
#   otu.we2.we4 <- distFrechet(Px = otu.we2$Days, Py = otu.we2$Abundance, Qx = otu.we4$Days, Qy = otu.we4$Abundance, FrechetSumOrMax = "sum")
#   otu.we4.we12 <- distFrechet(Px = otu.we4$Days, Py = otu.we4$Abundance, Qx = otu.we12$Days, Qy = otu.we12$Abundance, FrechetSumOrMax = "sum")
#   
#   d <- data.frame(we2.we12 = otu.we2.we12, we2.we4 = otu.we2.we4, we4.we12 = otu.we4.we12)
#   
#   return(d)
# }
# 
# cyano.otus <- c("Otu00005", "Otu00007", "Otu00037", "Otu00044", "Otu00049")
# titles <- c("Microcystis", "Synechococcus", "Pseudanabaena", "Synechococcus", "Unclassified")
# frech.dists <- data.frame()
# 
# 
# for (i in 1:length(cyano.otus)) {
#   
#   # Subset taxa to OTU of interest
#   cyano.otu <- 
#     cyanos %>%
#       subset_taxa(Species == cyano.otus[1])
#   
#   # Scale OTU across all sites
#   cyano.scale <- scale_otu_lineplots(cyano.otu)
#   
#   # Z-score scaling
#   cyano.frech <- get_frech_dist(cyano.scale)
#   
#   frech.dists[1,] <- cyano.frech
# }
#  
# 
# frech_dists <- rbind(mc, syne44, syne7, pseu, un)
# frech_dists$OTU <- c("Microcystis", "Syne44", "Syne7", "Pseuda.", "Unclass.")
# 
# frech_dists.melt <- melt(frech_dists, value.name = "Frechet")
# ```
# 
# ```{r}
# 
# ggplot(frech_dists.melt, aes(x = variable, y = OTU)) +
#   geom_tile(aes(fill = Frechet)) +
#   scale_fill_gradient(low = "#FF3300", high = "#000033") +
#   xlab("") +
#   ylab("")
#   

#** Figure Legend: ** I computed the Frechet distance for each OTU pairwise between the sites. This is a heatmap summarizing these distances

```


```{r}
# autocorrelation of error terms
lm.mode <- lm(formula = PC1 ~ pH, data = braypc.htrophs.sampledata)

plot(lm.mode)


df <- data.frame(days = braypc.htrophs.sampledata$Days[1:50], resid = lm.mode$residuals, station = braypc.htrophs.sampledata$Station[1:50])
ggplot(df, aes(x = days, y = resid, group = station, color = station)) + geom_line() + ggtitle("PC1~pH model")

acf(df$resid)
pacf(df$resid)
fit <- arima(df$resid, order = c(4,0,0))

library(astsa)
sarima(df$resid, 4, 0, 0)

#gls(PC1 ~ pH, data = braypc.htrophs.sampledata, correlation = corARMA(p=4, q=0))

```

```{r env variables PC2}




```
