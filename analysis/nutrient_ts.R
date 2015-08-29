library(astsa)
# Subset stations
WE2.dat <- subset(nutrient,Station=="WE2")[-3,c(1,3,15,18)]
WE2.dat$NitrateDiff <- c(NA,WE2.dat$Nitrate[-1]-WE2.dat$Nitrate[-20])
WE2.dat$ParMCDiff <- c(NA,WE2.dat$ParMC[-1]-WE2.dat$ParMC[-20])
cWE2<-ccf(WE2.dat$NitrateDiff[-1],WE2.dat$ParMCDiff[-1])
cWE2<-ccf(WE2.dat$Nitrate,WE2.dat$ParMC)


WE4.dat <- subset(nutrient,Station=="WE4")[-3,c(1,3,15,18)]
WE4.dat$NitrateDiff <- c(NA,WE4.dat$Nitrate[-1]-WE4.dat$Nitrate[-19])
WE4.dat$ParMCDiff <- c(NA,WE4.dat$ParMC[-1]-WE4.dat$ParMC[-19])
cWE4<-ccf(WE4.dat$NitrateDiff[-1],WE4.dat$ParMCDiff[-1])
cWE4<-ccf(WE4.dat$Nitrate,WE4.dat$ParMC)

WE12.dat <- subset(nutrient,Station=="WE12")[-2,c(1,3,15,18)]
WE12.dat$NitrateDiff <- c(NA,
                          WE12.dat$Nitrate[-1]-WE12.dat$Nitrate[-19])
WE12.dat$ParMCDiff <- c(NA,WE12.dat$ParMC[-1]-WE12.dat$ParMC[-19])
cWE12<-ccf(WE12.dat$NitrateDiff[-1],WE12.dat$ParMCDiff[-1])
cWE12<-ccf(WE12.dat$Nitrate,WE12.dat$ParMC)

# WE2 ParMC ~ H2O2
WE2.H2O2.ts<-ts(WE2.dat$H2O2)
WE2ParMC<-ts(WE2.dat$ParMC)


WE2.H2O2.lag1<-lag(WE2.H2O2.ts,-1)
WE2.H2O2.lag2<-lag(WE2.H2O2.ts,-2)
WE2.H2O2.lag3<-lag(WE2.H2O2.ts,-3)

WE2.ts<-ts.intersect(WE2ParMC,WE2.H2O2.lag1,WE2.H2O2.lag2,WE2.H2O2.lag3)

WE2.lag1<-lm(WE2ParMC~WE2.H2O2.lag1,data=WE2.ts) #ns
WE2.lag2<-lm(WE2ParMC~WE2.H2O2.lag2,data=WE2.ts) #sig ***
WE2.lag3<-lm(WE2ParMC~WE2.H2O2.lag3,data=WE2.ts) #ns 

# WE4 H202
WE4.H2O2.ts<-ts(WE4.dat$H2O2)
WE4ParMC<-ts(WE4.dat$ParMC)


WE4.H2O2.lag1<-lag(WE4.H2O2.ts,-1)
WE4.H2O2.lag2<-lag(WE4.H2O2.ts,-2)
WE4.H2O2.lag3<-lag(WE4.H2O2.ts,-3)

WE4.ts<-ts.intersect(WE4ParMC,WE4.H2O2.lag1,WE4.H2O2.lag2,WE4.H2O2.lag3)

WE4.lag1<-lm(WE4ParMC~WE4.H2O2.lag1,data=WE4.ts) #ns
WE4.lag2<-lm(WE4ParMC~WE4.H2O2.lag2,data=WE4.ts) #sig ***
WE4.lag3<-lm(WE4ParMC~WE4.H2O2.lag3,data=WE4.ts) #ns 


# WE12 H202
WE12.ts<-ts(WE12.dat)
WE12.H2O2.ts<-ts(WE12.dat$H2O2)
WE12ParMC<-ts(WE12.dat$ParMC)


WE12.H2O2.lag1<-lag(WE12.H2O2.ts,-1)
WE12.H2O2.lag2<-lag(WE12.H2O2.ts,-2)
WE12.H2O2.lag3<-lag(WE12.H2O2.ts,-3)

WE12.ts<-ts.intersect(WE12ParMC,WE12.H2O2.lag1,WE12.H2O2.lag2,WE12.H2O2.lag3)

WE12.lag1<-lm(WE12ParMC~WE12.H2O2.lag1,data=WE12.ts) #sig ***
WE12.lag2<-lm(WE12ParMC~WE12.H2O2.lag2,data=WE12.ts) #ns
WE12.lag3<-lm(WE12ParMC~WE12.H2O2.lag3,data=WE12.ts) #ns 




# WE2 ParMC~Nitrate

WE2.N.ts<-ts(WE2.dat$Nitrate)
WE2.ParMC.ts<-ts(WE2.dat$ParMC)


WE2.N.lag4<-lag(WE2.H2O2.ts,-4)
WE2.N.lag7<-lag(WE2.H2O2.ts,-7)

WE2.ts<-ts.intersect(WE2.ParMC.ts,WE2.N.lag4,WE2.N.lag7,WE2.N.ts)

WE2.lag4<-lm(WE2.ParMC.ts~WE2.N.lag4,data=WE2.ts) 
WE2.lag7<-lm(WE2.ParMC.ts~WE2.N.lag7,data=WE2.ts)
WE2.lag0<-lm(WE2.ParMC.ts~WE2.N.ts,data=WE2.ts)

WE2.NO3diff.ts<-ts(WE2.dat$NitrateDiff)


WE2.ParMCdiff.ts<-ts(WE2.dat$ParMCDiff)

WE2.NO3diff.lagg1<-lag(WE2.NO3diff.ts,1)
WE2.NO3diff.lag1<-lag(WE2.NO3diff.ts,-1)
WE2.NO3diff.lag2<-lag(WE2.NO3diff.ts,-2)
WE2.NO3diff.lag3<-lag(WE2.NO3diff.ts,-3)
WE2.NO3diff.lag4<-lag(WE2.NO3diff.ts,-4)
WE2.NO3diff.lag5<-lag(WE2.NO3diff.ts,-5)
WE2.NO3diff.lag6<-lag(WE2.NO3diff.ts,-6)


WE2.ts<-ts.intersect(WE2.ParMCdiff.ts,WE2.NO3diff.lagg1,WE2.NO3diff.lag1,WE2.NO3diff.lag2,WE2.NO3diff.ts,WE2.NO3diff.lag3,WE2.NO3diff.lag4,WE2.NO3diff.lag5,WE2.NO3diff.lag6)

WE2.lag1<-lm(WE2.ParMCdiff.ts~WE2.NO3diff.lag1,data=WE2.ts) #ns
WE2.lag2<-lm(WE2.ParMCdiff.ts~WE2.NO3diff.lag2,data=WE2.ts) #ns
WE2.lagg1<-lm(WE2.ParMCdiff.ts~WE2.NO3diff.lagg1,data=WE2.ts) #ns
WE2.lag0<-lm(WE2.ParMCdiff.ts~WE2.NO3diff.ts,data=WE2.ts) #ns
WE2.lag3<-lm(WE2.ParMCdiff.ts~WE2.NO3diff.lag3,data=WE2.ts)
WE2.lag4<-lm(WE2.ParMCdiff.ts~WE2.NO3diff.lag4,data=WE2.ts)
WE2.lag5<-lm(WE2.ParMCdiff.ts~WE2.NO3diff.lag5,data=WE2.ts)
WE2.lag6<-lm(WE2.ParMCdiff.ts~WE2.NO3diff.lag6,data=WE2.ts)


# WE4 ParMC~Nitrate

WE4.mc.no3 <- lm(ParMCDiff~NitrateDiff,data=WE4.dat)
summary(WE4.mc.no3) # sig * 

# WE12 ParMC~Nitrate
WE12.mc.no3 <- lm(ParMCDiff~NitrateDiff,data=WE12.dat)
summary(WE12.mc.no3) # sig **