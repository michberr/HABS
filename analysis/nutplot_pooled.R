## Plot pooled

# Make a summary
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

```{r plot pooled, fig.height=10, fig.width=10, warning=FALSE}

nutplot_pooled <- function(df,nutrient,title){
  nut_sum <- summarySE(df,measurevar = nutrient,groupvars = "Date",na.rm=TRUE)
  nut_sum$ymax <- nut_sum[,nutrient]+nut_sum[,"se"]
  nut_sum$ymin <- nut_sum[,nutrient]-nut_sum[,"se"]
  ggplot(nut_sum,aes_string(x="Date",y=nutrient,group=1))+
    geom_line(size=.5)+
    geom_errorbar(aes(ymax=ymax,ymin=ymin), width=.3)+
    geom_point()+
    ggtitle(title)+ scale_shape(solid = FALSE) + xlab("")+
    ylab("ug/L")+
    theme(axis.text.x = element_text(angle=45,
                                     vjust=1, hjust = 1, size=10, face="bold"),     
          axis.title.y = element_text(size=14),
          plot.title = element_text(size=16,face="bold")
    )
} 

p_NO3 <- nutplot_pooled(df = nutrient,nutrient = "Nitrate",title = "Nitrate")         

p_NH4 <- nutplot_pooled(df = nutrient,nutrient = "Ammonia",title = "Ammonia") 

p_SRP <- nutplot_pooled(df = nutrient,nutrient = "SRP",title = "Soluble Reactive Phosphorus")

p_parmc <- nutplot_pooled(df = nutrient,nutrient = "ParMC",title = "Particulate Microcystin")

p_phyco <- nutplot_pooled(df = nutrient,nutrient = "Phycocyanin",title = "Phycocyanin")

p_chla <- nutplot_pooled(df = nutrient,nutrient = "Chla",title = "Chlorophyll A")

theme_set(theme_bw())
multiplot(p_NO3,p_NH4,p_SRP,p_parmc,p_phyco,p_chla,cols = 2)
```