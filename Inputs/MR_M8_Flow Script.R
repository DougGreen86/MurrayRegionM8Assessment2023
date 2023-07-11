# Script for the assessment of the flow data for the Murray Region Matter 8 assessment



#### Set environment ####

#clear working environment
rm(list=ls())

setwd("C:/Ranalysis/MR_M8_Reporting/Inputs")

#Load Libraries needed

library(vegan)
library(mgcv)
library(lattice)
library(xts)
library(zoo)
library(ggplot2)
library(lme4)
library(MASS)
library(tidyverse)
library(fs)
library(lubridate)
library("rstan")
library("rstanarm")
library("doParallel")
library(naniar)
library(EGRET)
library(scales)
library(ggpubr)



#Months that comprise the low flow season
lfsMonths <- c(12,1:4)
#Months that comprise the T1 seasons
T1Months <- 5:6
# look up table for trend assessment
luTrend <- read_csv("C:/Ranalysis/MR_M8_Reporting/Inputs/luTrend.csv")

# A function to get the trend from p and slope
trend_set <- function(x) {
  
  if(x >= 0.99) luTrend$Trend[9] else            # Virtually certain decrease
    if(between(x, 0.95,0.9899)) luTrend$Trend[8] else      # Extremely likely decrease
      if(between(x,0.9,0.9499)) luTrend$Trend[7] else       #Very likely decrease
        if(between(x,0.66,0.8999)) luTrend$Trend[6] else    #Likely decrease
          if(between(x,0.33,0.6599)) luTrend$Trend[5] else    #About as likely as not increase
            if(between(x,0.10,0.3299)) luTrend$Trend[4] else   #Likely increase
              if(between(x,0.05,0.0999)) luTrend$Trend[3] else  #Very likely increase
                if(between(x,0.0101,0.0499)) luTrend$Trend[2] else  #Extremely likely increase
                  if(x<=0.01) luTrend$Trend[1] else                 #Virtually certain increase
                    if(is.na(x)) luTrend$Trend[10] else   # not applicable
                      luTrend$Trend[11]                        # unknown
  
}


dat <- read_csv("WorldsEnd_flow.csv") %>% #gathers the data together into a single tibble
  dplyr::mutate(Flow = A4261148
                , Date = dmy(Date)
                , month = month(Date)
                , year = year(Date)
                , flowYear = if_else(month<12,year,year+1)
                , lfs = if_else(month %in% lfsMonths,"Yes","No")
                , T1 = if_else(month %in% T1Months,"Yes","No")
                , Site = "A4261148"
  ) 

datInter <- dat %>% #creates a seperate tibble to play with
  drop_na() %>%
  dplyr::group_by(Site,flowYear) %>% #groups data by type (current/no dams), site and flow year
  dplyr::summarise(ZeroFlowDays = sum(Flow < 0.05)
                   ,inter = sum(Flow < 0.05)/365) %>% # counts the number of days with flow over 0.05ML/day
  dplyr::ungroup() 

## Low flow fresh --------
#A fresh is defined in the WAPs as 2 time the median non-zero flow
#In this script we use the no dams data to calculate it for a site

lfsFreshThreshold <- dat %>% #creates a seperate tibble to play with
  dplyr::filter(lfs == "Yes"
                , Flow > 0.05
  ) %>% # filters the tibble to only include the no dams data
  dplyr::group_by(Site) %>%
  dplyr::summarise(lfsFreshThresh = 2*median(Flow)) #calculates the LFS fresh threshold value


datLFS <- dat %>% #creates a seperate tibble to play with
  drop_na() %>%
  dplyr::left_join(lfsFreshThreshold) %>%
  dplyr::filter(lfs == "Yes"
                , Flow > lfsFreshThresh
  ) %>% #pulls out the data for each year that is in the LFS
  dplyr::group_by(flowYear,Site) %>% #groups data into year, type and site
  dplyr::summarise(lfsFreshDays = n()
                   ,lfsFresh = n()/sum(days_in_month(lfsMonths))) %>% #sums the number of days over the threshold
  dplyr::ungroup()


## T1 Fresh ---------

T1FreshThreshold <- dat %>% #creates a seperate tibble to play with
  dplyr::filter(T1 == "Yes"
                , Flow > 0.05
  ) %>% #filters by no dam scenario and if datais part of T1
  dplyr::group_by(Site) %>%
  dplyr::summarise(T1FreshThresh = 2*median(Flow)) %>% #calculates the T1 fresh flow threshold
  dplyr::ungroup()

datT1 <- dat %>% #creates a seperate tibble to play with
  drop_na() %>%
  dplyr::left_join(T1FreshThreshold) %>%
  dplyr::filter(T1 == "Yes"
                , Flow > T1FreshThresh
  ) %>% #pulls out the T1 season data
  dplyr::group_by(flowYear,Site) %>% #Groups data 
  dplyr::summarise(T1FreshDays = n()
                   ,T1Fresh = n()/sum(days_in_month(T1Months))) %>% #sums the number of days over the threshold
  dplyr::ungroup()

datMedian <- dat %>% #creates a seperate tibble to play with
  drop_na() %>%
  dplyr::group_by(flowYear,Site) %>% #Groups data 
  dplyr::summarise(medianFlow = median(Flow),
                   meanFlow = mean(Flow)) %>% 
  dplyr::ungroup()

datMissing <- dat %>%
  dplyr::group_by(flowYear,Site) %>%
  dplyr::summarise(missingDays = sum(is.na(Flow)))



# Summary data tibble

datSum <- datInter %>% #creates a seperate tibble to play with
  dplyr::left_join(datLFS) %>%
  dplyr::left_join(datT1) %>%
  dplyr::mutate(lfsFreshDays = replace_na(lfsFreshDays,0)
                ,lfsFresh = replace_na(lfsFresh,0)
                ,T1FreshDays = replace_na(T1FreshDays,0)
                ,T1Fresh = replace_na(T1Fresh,0)) %>%
  dplyr::left_join(datMedian) %>%
  dplyr::left_join(datMissing)

datSum %>%
  ggplot() +
  geom_boxplot(aes (Site,inter)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Site") +
  ylab("Intermittency")

ggsave("../plots/Inter_boxplot.png")

datSum %>%
  ggplot() +
  geom_boxplot(aes (Site,lfsFresh)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Site") +
  ylab("LFS Fresh")

ggsave("../plots/LFSFresh_boxplot.png")

datSum %>%
  ggplot() +
  geom_boxplot(aes (Site,T1Fresh)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Site") +
  ylab("T1 Fresh")

ggsave("../plots/T1Fresh_boxplot.png")

datSum %>%
  ggplot() +
  geom_boxplot(aes (Site,medianFlow)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Site") +
  ylab("Median Flow")

ggsave("../plots/Median_boxplot.png")

datSum %>%
  ggplot() +
  geom_boxplot(aes (Site,meanFlow)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Site") +
  ylab("Mean Flow")

ggsave("../plots/Mean_boxplot.png")

datSum %>%
  ggplot() +
  geom_point(aes (flowYear,medianFlow)) +
  geom_smooth(aes (flowYear,medianFlow),method = "lm")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Site") +
  ylab("Median Flow")

ggsave("../plots/Median_PointPlot.png")

datSum %>%
  ggplot() +
  geom_point(aes (flowYear,meanFlow)) +
  geom_smooth(aes (flowYear,meanFlow),method = "lm")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Site") +
  ylab("Mean Flow")

ggsave("../plots/Mean_pointPlot.png")



### Runoff Analysis ###

# read in Burra rainfall data and convert to monthly totals 
datRain <- read_csv("Burra_rainfallSUB.csv") %>%
  dplyr::mutate(Date = dmy(paste(Day,"-",Month,"-",Year)),
                DataType = "Rainfall",
                value = Rainfall)%>%
  dplyr::select(., -c(Station, Year, Month, Day, Rainfall)) %>%
  drop_na()

datRainASS <- datRain %>%
  dplyr::filter(Year >= 2014)


datRainLT <- read_csv("Burra_rainfall.csv") %>%
  dplyr::mutate(Date = dmy(paste(Day,"-",Month,"-",Year)),
                DataType = "Long term rainfall",
                value = Rainfall)%>%
  dplyr::select(., -c(Station, Year, Month, Day, Rainfall)) %>%
  drop_na()

datFlow <- dat %>%
  dplyr::mutate(DataType = "A4261148",
                value = Flow) %>%
  dplyr::select(., -c(A4261148, month, year, flowYear, lfs, T1, Site,Flow))%>%
  drop_na()

datFlowLT <- read_csv("WorldsEnd_LongTerm.csv") %>% #gathers the data together into a single tibble
  dplyr::mutate(FlowLT = A4260536
                , Date = dmy(Date)
                , month = month(Date)
                , year = year(Date)
                , flowYear = if_else(month<12,year,year+1)
                , lfs = if_else(month %in% lfsMonths,"Yes","No")
                , T1 = if_else(month %in% T1Months,"Yes","No")
                , Site = "A4260536") %>%
  dplyr::mutate(DataType = "A4260536",
                value = FlowLT) %>%
  dplyr::select(., -c(A4260536, month, year, flowYear, lfs, T1, Site,FlowLT))%>%
  drop_na()

datRO <- datFlow %>%
  dplyr::left_join(datRain, by = "Date") %>%
  dplyr::mutate(Flow = value.x,
                Rainfall = value.y) %>%
  dplyr::select(., -c(DataType.x, DataType.y, value.x, value.y)) %>%
  dplyr::mutate_if(is.numeric, list(~na_if(., Inf))) %>%
  drop_na() %>%
  dplyr::mutate(month = month(Date)
                , year = year(Date)) %>%
  dplyr::group_by(month, year) %>%
  dplyr::summarise(monthyFlow = sum(Flow)
                ,monthlyRain = sum(Rainfall)) %>%
  dplyr::mutate(Date = dmy(paste(1,"-",month,"-",year))) %>%
  dplyr::mutate(RO = monthyFlow/monthlyRain) %>%  
  dplyr::mutate_if(is.numeric, list(~na_if(., Inf))) %>%
  drop_na()


datTall <- rbind(datFlow,datRain)
datTallFlow <- rbind(datFlow,datFlowLT)


ggplot(datFlow, aes(x=Date, y = value)) +
  geom_line() +
  geom_smooth(color = "red") +
  geom_smooth(method = "lm") +
  ylab("Daily Flow (ML/day)")+
  scale_y_continuous(trans = "log2", labels = comma) 
ggsave("../plots/log daily flow.png")


ggplot(datRain, aes(x=Date, y = log(value)))+
  geom_line() +
  geom_smooth(method = "lm")+
  geom_smooth(color = "red")+
  ylab("Log Daily Rainfall (mm)")
ggsave("../plots/log daily rain.png")


ggplot(datFlowLT, aes(x=Date, y = value)) +
  geom_line() +
  geom_smooth(method = "lm")+
  geom_smooth(color = "red")+
  ylab("Daily Flow (ML/day)") +
  scale_y_continuous(trans = "log2", labels = comma) 
ggsave("../plots/log daily flow - long term.png")

ggplot(datRainLT, aes(x=Date, y = log(value)))+
  geom_line() +
  geom_smooth(method = "lm")+
  geom_smooth(color = "red")+
  ylab("Log Daily Rainfall (mm)")
ggsave("../plots/log daily rain.png")

ggplot(datRO, aes(x=Date, y = RO)) +
  geom_line() +
  geom_smooth(color = "red") +
  geom_smooth(method = "lm") +
  ylab("Monthy Runoff (ML/day/Km2)")+
  scale_y_continuous(trans = "log2", labels = comma)
ggsave("../plots/log monthy runoff.png")



ggplot(datTall, aes(x=Date, y = value, color = (DataType), alpha(0.2))) +
  geom_line() +
  geom_smooth(method = "lm") +
  scale_color_manual(values=c("red", "blue"))+
  scale_y_continuous(trans = "log2", labels = comma)  +
  ylab("Daily flow (ML/d)")
  
ggplot(datTallFlow, aes(x=Date, y = value, color = (DataType), alpha(0.2))) +
    geom_line() +
    geom_smooth(method = "lm") +
  scale_color_manual(values=c("orange", "purple"))+
  labs(color = "Flow Data") +
  ylab("Daily flow (ML/d)")+
  scale_y_continuous(trans = "log2", labels = comma) 

ggsave("../plots/FlowBoth.png")

modRainASS <- stan_glm(value ~ Date 
                    , data = datRain
                    , family = gaussian
                    , cores = 12
                    , iter = 4000)

modFlow <- stan_glm(value ~ Date 
                    , data = datFlow
                    , family = gaussian
                    , cores = 12
                    , iter = 4000)

modFlowLT <- stan_glm(value ~ Date 
                    , data = datFlowLT
                    , family = gaussian
                    , cores = 12
                    , iter = 4000)

modRainLT <- stan_glm(value ~ Date 
                      , data = datRainLT
                      , family = gaussian
                      , cores = 12
                      , iter = 4000)

modRO <- stan_glm(RO ~ Date 
                      , data = datRO
                      , family = gaussian
                      , cores = 12
                      , iter = 4000)

# Model check plots to be looked at to ensure suitablility
ggsave("../plots/modRain.png", stan_trace(modRain))
ggsave("../plots/modRain.png", plot(modRain, "rhat_hist"))
ggsave("../plots/modFlow.png", stan_trace(modFlow))
ggsave("../plots/modFlow.png", plot(modFlow, "rhat_hist"))
ggsave("../plots/modFlowLT.png", stan_trace(modFlowLT))
ggsave("../plots/modFlowLT.png", plot(modFlowLT, "rhat_hist"))
ggsave("../plots/modRainLT.png", stan_trace(modRainLT))
ggsave("../plots/modRainLT.png", plot(modRainLT, "rhat_hist"))
ggsave("../plots/modRO.png", stan_trace(modRO))
ggsave("../plots/modRO.png", plot(modRO, "rhat_hist"))

## Rainfall data processing
divideBy <- nrow(as_tibble(modRain))

# provide estimates of slope
modRainRes <- as_tibble(modRain) %>%
  dplyr::select(Date) %>%
  write_csv("../tbl/modRainRes.csv")

trendRain <- modRainRes %>%
  dplyr::summarise(negSlope = sum(Date<0)/divideBy #count slopes less than 0
                   , posSlope = sum(Date>0)/divideBy
                   , ci90loSlope = quantile(Date, 0.05)
                   , ci90upSlope = quantile(Date, 0.95)
                   , meanSlope = mean(Date)
  ) %>%
  dplyr::mutate(Trend = map(negSlope,trend_set)
                , intervalSlope = paste0(round(ci90loSlope,3)," to ",round(ci90upSlope,3))
  ) %>%
  tidyr::unnest() %>% #brings out as character but retains data
  dplyr::mutate(Trend = fct_explicit_na(Trend, na_level = "Not applicable"),
                Metric = "Rainfall") %>%
   write_csv("../tbl/trendRain.csv")

## rainASSfall data processing
divideBy <- nrow(as_tibble(modRainASS))

# provide estimates of slope
modrainASSRes <- as_tibble(modRainASS) %>%
  dplyr::select(Date) %>%
  write_csv("../tbl/modrainASSRes.csv")

trendrainASS <- modrainASSRes %>%
  dplyr::summarise(negSlope = sum(Date<0)/divideBy #count slopes less than 0
                   , posSlope = sum(Date>0)/divideBy
                   , ci90loSlope = quantile(Date, 0.05)
                   , ci90upSlope = quantile(Date, 0.95)
                   , meanSlope = mean(Date)
  ) %>%
  dplyr::mutate(Trend = map(negSlope,trend_set)
                , intervalSlope = paste0(round(ci90loSlope,3)," to ",round(ci90upSlope,3))
  ) %>%
  tidyr::unnest() %>% #brings out as character but retains data
  dplyr::mutate(Trend = fct_explicit_na(Trend, na_level = "Not applicable"),
                Metric = "rainASSfall") %>%
  write_csv("../tbl/trendrainASS.csv")

## Flow data processing
divideBy <- nrow(as_tibble(modFlow))

# provide estimates of slope
modFlowRes <- as_tibble(modFlow) %>%
  dplyr::select(Date) %>%
  write_csv("../tbl/modFlowRes.csv")

trendFlow <- modFlowRes %>%
  dplyr::summarise(negSlope = sum(Date<0)/divideBy #count slopes less than 0
                   , posSlope = sum(Date>0)/divideBy
                   , ci90loSlope = quantile(Date, 0.05)
                   , ci90upSlope = quantile(Date, 0.95)
                   , meanSlope = mean(Date)
  ) %>%
  dplyr::mutate(Trend = map(negSlope,trend_set)
                , intervalSlope = paste0(round(ci90loSlope,3)," to ",round(ci90upSlope,3))
  ) %>%
  tidyr::unnest() %>% #brings out as character but retains data
  dplyr::mutate(Trend = fct_explicit_na(Trend, na_level = "Not applicable"),
                Metric = "Flow") %>%
   write_csv("../tbl/trendFlow.csv")


## Flow data processing
divideBy <- nrow(as_tibble(modFlowLT))

# provide estimates of slope
modFlowResLT <- as_tibble(modFlowLT) %>%
  dplyr::select(Date) %>%
  write_csv("../tbl/modFlowResLT.csv")

trendFlowLT <- modFlowResLT %>%
  dplyr::summarise(negSlope = sum(Date<0)/divideBy #count slopes less than 0
                   , posSlope = sum(Date>0)/divideBy
                   , ci90loSlope = quantile(Date, 0.05)
                   , ci90upSlope = quantile(Date, 0.95)
                   , meanSlope = mean(Date)
  ) %>%
  dplyr::mutate(Trend = map(negSlope,trend_set)
                , intervalSlope = paste0(round(ci90loSlope,3)," to ",round(ci90upSlope,3))
  ) %>%
  tidyr::unnest() %>% #brings out as character but retains data
  dplyr::mutate(Trend = fct_explicit_na(Trend, na_level = "Not applicable"),
                Metric = "Long Term Flow") %>%
   write_csv("../tbl/trendFlowLT.csv")

## Rainfall data processing
divideBy <- nrow(as_tibble(modRainLT))

# provide estimates of slope
modRainResLT <- as_tibble(modRainLT) %>%
  dplyr::select(Date) %>%
  write_csv("../tbl/modRainResLT.csv")

trendRainLT <- modRainResLT %>%
  dplyr::summarise(negSlope = sum(Date<0)/divideBy #count slopes less than 0
                   , posSlope = sum(Date>0)/divideBy
                   , ci90loSlope = quantile(Date, 0.05)
                   , ci90upSlope = quantile(Date, 0.95)
                   , meanSlope = mean(Date)
  ) %>%
  dplyr::mutate(Trend = map(negSlope,trend_set)
                , intervalSlope = paste0(round(ci90loSlope,3)," to ",round(ci90upSlope,3))
  ) %>%
  tidyr::unnest() %>% #brings out as character but retains data
  dplyr::mutate(Trend = fct_explicit_na(Trend, na_level = "Not applicable"),
                Metric = "Long Term Rainfall") %>%
  write_csv("../tbl/trendRainLT.csv")

## RO data processing
divideBy <- nrow(as_tibble(modRO))

# provide estimates of slope
modRORes <- as_tibble(modRO) %>%
  dplyr::select(Date) %>%
  write_csv("../tbl/modRO.csv")

trendRO <- modRORes %>%
  dplyr::summarise(negSlope = sum(Date<0)/divideBy #count slopes less than 0
                   , posSlope = sum(Date>0)/divideBy
                   , ci90loSlope = quantile(Date, 0.05)
                   , ci90upSlope = quantile(Date, 0.95)
                   , meanSlope = mean(Date)
  ) %>%
  dplyr::mutate(Trend = map(negSlope,trend_set)
                , intervalSlope = paste0(round(ci90loSlope,3)," to ",round(ci90upSlope,3))
  ) %>%
  tidyr::unnest() %>% #brings out as character but retains data
  dplyr::mutate(Trend = fct_explicit_na(Trend, na_level = "Not applicable"),
                Metric = "Monthly Runoff") %>%
  write_csv("../tbl/trendRO.csv")



datBig <- dat %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(BigFlow = sum(Flow >20)) %>%
  dplyr::mutate(BigFlow = replace_na(BigFlow,0),
                Site = "Worlds End") %>% 
  dplyr::filter(year != "2008")  %>%
  dplyr::filter(year != "2009")

datBig %>%
  ggplot(aes (year,BigFlow)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Site") +
  ylab("Flows > 20 ML/d")

datBigLT <- datFlowLT %>%
  dplyr::mutate(year = year(Date)) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(BigFlow = sum(value >20)) %>%
  dplyr::mutate(Site = "Worlds End")

datBigCombo <- rbind(datBigLT,datBig)

datBigCombo %>%
  ggplot() +
  geom_boxplot(aes (Site,BigFlow)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Site") +
  ylab("Flows > 20 ML/d")

datBigCombo %>%
  ggplot(aes (year,BigFlow)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Year") +
  ylab("Number of days > 20 ML/d")

ggsave("../plots/BigFlow.PNG")

modBigFlow <- stan_glm(BigFlow ~ year 
                  , data = datBigCombo
                  , family = neg_binomial_2()
                  , cores = 12
                  , iter = 4000)

# Model check plots to be looked at to ensure suitablility
ggsave("../plots/modBigFlow.png", stan_trace(modBigFlow))
ggsave("../plots/modBigFlow.png", plot(modBigFlow, "rhat_hist"))

## Rainfall data processing
divideBy <- nrow(as_tibble(modBigFlow))

# provide estimates of slope
modBigFlowRes <- as_tibble(modBigFlow) %>%
  dplyr::select(year) %>%
  write_csv("../tbl/modBigFlowRes.csv")

trendBigFlow <- modBigFlowRes %>%
  dplyr::summarise(negSlope = sum(year<0)/divideBy #count slopes less than 0
                   , posSlope = sum(year>0)/divideBy
                   , ci90loSlope = quantile(year, 0.05)
                   , ci90upSlope = quantile(year, 0.95)
                   , meanSlope = mean(year)
  ) %>%
  dplyr::mutate(Trend = map(negSlope,trend_set)
                , intervalSlope = paste0(round(ci90loSlope,3)," to ",round(ci90upSlope,3))
  ) %>%
  tidyr::unnest() %>% #brings out as character but retains data
  dplyr::mutate(Trend = fct_explicit_na(Trend, na_level = "Not applicable"),
                Metric = "Big Flows") %>%
  write_csv("../tbl/trendBigFlows.csv")

OverallTrendLT <- bind_rows(trendFlowLT,
                            trendFlow,
                            trendRain,
                            trendRainLT,
                            trendRO,
                            trendBigFlow) %>%
  write_csv("../tbl/Overall trend results.csv")



### Rainfall Flow Runoff Figure  ###


rainGG <- ggplot(datRain, aes(x=Date, y = value))+
  geom_line() +
  geom_smooth(color = "red")+
  geom_smooth(method = "lm")+
  ylab("Daily Rainfall (mm)")+
  scale_y_continuous(trans = "log2", labels = comma)

flowGG <- ggplot(datFlow, aes(x=Date, y = value)) +
  geom_line() +
  geom_smooth(color = "red") +
  geom_smooth(method = "lm") +
  ylab("Daily Flow (ML/day)")+
  scale_y_continuous(trans = "log2", labels = comma) 


ROGG <- ggplot(datRO, aes(x=Date, y = RO)) +
  geom_line() +
  geom_smooth(color = "red") +
  geom_smooth(method = "lm") +
  ylab("Monthy Runoff (ML/mm)")+
  scale_y_continuous(trans = "log2", labels = comma)

ggarrange(rainGG, flowGG, ROGG,
          ncol = 1, nrow = 3, align = "v",
          widths = 0.5)
ggsave("../plots/Recent Combo graph.PNG")




datRainannual <- datRainLT %>%
  dplyr::mutate(Year = year(Date)) %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(AnnualRain = sum(value)) %>%
dplyr::filter(Year != "2020")

ggplot(datRainannual, aes(x=Year, y = AnnualRain)) +
  geom_point() +
  geom_smooth(method = "lm")

annRailLM <- lm(AnnualRain ~ Year, datRainannual)
