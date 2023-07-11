# Script for the Evaluation of environmental outcomes for the Murray Region Matter 8 assessment
# Assessment uses flow data from the worlds end gauge station (both data sets)
# Assessment uses rainfall data from Burra Community School
# Impact metric calculations at end of script

# Author Douglas Green douglas.green1@sa.gov.au
#### Set environment ####

#clear working environment
rm(list=ls())

setwd("C:/Ranalysis/MR_M8_Reporting/Inputs")

#Load Libraries needed

library(mgcv)
library(lattice)
library(xts)
library(ggplot2)
library(MASS)
library(tidyverse)
library(fs)
library(lubridate)
library("rstan")
library("rstanarm")
library("doParallel")
library(scales)
library(ggpubr)


#Trend assessment based on classes from IPCC ratings for trends
# look up table for trend assessment
luTrend <- read_csv("C:/Ranalysis/MR_M8_Reporting/Inputs/luTrend.csv")

# A function to get the trend from p and slope
# function requires the proportion of negative slopes from the Bayesian modelling. 
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

# Read in Flow data (short term flow data)
dat <- read_csv("WorldsEnd_flow.csv") %>% #gathers the data together into a single tibble
  dplyr::mutate(Flow = A4261148
                , Date = dmy(Date)
                , month = month(Date)
                , year = year(Date)
                , flowYear = if_else(month<12,year,year+1)
                , Site = "A4261148"
  ) 

# Read in flow data (from longer term but discontuinued gauge)
datFlowLT <- read_csv("WorldsEnd_LongTerm.csv") %>% #gathers the data together into a single tibble
  dplyr::mutate(FlowLT = A4260536
                , Date = dmy(Date)
                , month = month(Date)
                , year = year(Date)
                , flowYear = if_else(month<12,year,year+1)
                , Site = "A4260536") %>%
  dplyr::mutate(DataType = "A4260536",
                value = FlowLT) %>%
  dplyr::select(., -c(A4260536, month, year, flowYear, Site,FlowLT))%>%
  drop_na()

# Read in long term rain data
datRainLT <- read_csv("Burra_rainfall.csv") %>%
  dplyr::mutate(Date = dmy(paste(Day,"-",Month,"-",Year)),
                DataType = "Long term rainfall",
                value = Rainfall)%>%
  dplyr::select(., -c(Station, Year, Month, Day, Rainfall)) %>%
  drop_na()

# read in  rainfall data subset to same window as short term flow data
datRain <- read_csv("Burra_rainfallSUB.csv") %>%
  dplyr::mutate(Date = dmy(paste(Day,"-",Month,"-",Year)),
                DataType = "Rainfall",
                value = Rainfall)%>%
  dplyr::select(., -c(Station, Year, Month, Day, Rainfall)) %>%
  drop_na()


# Calculate zero flow days
datInter <- dat %>% #creates a seperate tibble to play with
  drop_na() %>%
  dplyr::group_by(Site,flowYear) %>% #groups data by type (current/no dams), site and flow year
  dplyr::summarise(ZeroFlowDays = sum(Flow < 0.05)
                   ,inter = sum(Flow < 0.05)/365) %>% # counts the number of days with flow over 0.05ML/day
  dplyr::ungroup() 


# Create  Figure 3 (Flow trends)

datFlow <- dat %>%
  dplyr::mutate(DataType = "A4261148",
                value = Flow) %>%
  dplyr::select(., -c(A4261148, month, year, flowYear, Site,Flow))%>%
  drop_na()

datTallFlow <- rbind(datFlow,datFlowLT)

ggplot(datTallFlow, aes(x=Date, y = value, color = (DataType), alpha(0.2))) +
  geom_line() +
  geom_smooth(method = "lm") +
  scale_color_manual(values=c("orange", "purple"))+
  labs(color = "Flow Data") +
  ylab("Daily flow (ML/d)")+
  scale_y_continuous(trans = "log2", labels = comma) 

#save plot to plots folder
ggsave("../plots/FlowBoth.png")


### Rainfall Flow Runoff Figure (Figure 4)

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


# Plotting Higher Flows (Figure 5)

datBig <- dat %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(BigFlow = sum(Flow >20)) %>%
  dplyr::mutate(BigFlow = replace_na(BigFlow,0),
                Site = "Worlds End") %>% 
  dplyr::filter(year != "2008")  %>%
  dplyr::filter(year != "2009")
datBigLT <- datFlowLT %>%
  dplyr::mutate(year = year(Date)) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(BigFlow = sum(value >20)) %>%
  dplyr::mutate(Site = "Worlds End")
datBigCombo <- rbind(datBigLT,datBig)

datBigCombo %>%
  ggplot(aes (year,BigFlow)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Year") +
  ylab("Number of days > 20 ML/d")

ggsave("../plots/BigFlow.PNG")


## Trend assessments 

# Bayesian models for trend assessment 
# Standard glm, no random factors
# Gaussian as the data in continuous
# Cores = 12 for faster processing - not overly important in this case
# iter = 4000 means 4000 model runs

modRain <- stan_glm(value ~ Date 
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

modBigFlow <- stan_glm(BigFlow ~ year 
                       , data = datBigCombo
                       , family = neg_binomial_2()
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
ggsave("../plots/modBigFlow.png", stan_trace(modBigFlow))
ggsave("../plots/modBigFlow.png", plot(modBigFlow, "rhat_hist"))

# Processing the model results 
# Three steps
# 1) set length of results
# 2) extract slopes from each run
# 3) apply function to assess the slopes and produce summary table

## Rainfall (2008-2020) data processing
divideBy <- nrow(as_tibble(modRain))
modRainRes <- as_tibble(modRain) %>%
  dplyr::select(Date) %>%
  write_csv("../tbl/modRainRes.csv")
trendRain <- modRainRes %>%
  dplyr::summarise(negSlope = sum(Date<0)/divideBy #count slopes less than 0
                   , posSlope = sum(Date>0)/divideBy
                   , ci90loSlope = quantile(Date, 0.05)
                   , ci90upSlope = quantile(Date, 0.95)
                   , meanSlope = mean(Date)) %>%
  dplyr::mutate(Trend = map(negSlope,trend_set)
                , intervalSlope = paste0(round(ci90loSlope,3)," to ",round(ci90upSlope,3))) %>%
  tidyr::unnest(cols = c(Trend)) %>% #brings out as character but retains data
  dplyr::mutate(Trend = fct_explicit_na(Trend, na_level = "Not applicable"),
                Metric = "Rainfall") %>%
  write_csv("../tbl/trendRain.csv")


## Flow data (2008-2020) processing
divideBy <- nrow(as_tibble(modFlow))
modFlowRes <- as_tibble(modFlow) %>%
  dplyr::select(Date) %>%
  write_csv("../tbl/modFlowRes.csv")
trendFlow <- modFlowRes %>%
  dplyr::summarise(negSlope = sum(Date<0)/divideBy #count slopes less than 0
                   , posSlope = sum(Date>0)/divideBy
                   , ci90loSlope = quantile(Date, 0.05)
                   , ci90upSlope = quantile(Date, 0.95)
                   , meanSlope = mean(Date)) %>%
  dplyr::mutate(Trend = map(negSlope,trend_set)
                , intervalSlope = paste0(round(ci90loSlope,3)," to ",round(ci90upSlope,3))) %>%
  tidyr::unnest(cols = c(Trend)) %>% #brings out as character but retains data
  dplyr::mutate(Trend = fct_explicit_na(Trend, na_level = "Not applicable"),
                Metric = "Flow") %>%
  write_csv("../tbl/trendFlow.csv")


## Flow data (long term) processing
divideBy <- nrow(as_tibble(modFlowLT))
modFlowResLT <- as_tibble(modFlowLT) %>%
  dplyr::select(Date) %>%
  write_csv("../tbl/modFlowResLT.csv")
trendFlowLT <- modFlowResLT %>%
  dplyr::summarise(negSlope = sum(Date<0)/divideBy #count slopes less than 0
                   , posSlope = sum(Date>0)/divideBy
                   , ci90loSlope = quantile(Date, 0.05)
                   , ci90upSlope = quantile(Date, 0.95)
                   , meanSlope = mean(Date)) %>%
  dplyr::mutate(Trend = map(negSlope,trend_set)
                , intervalSlope = paste0(round(ci90loSlope,3)," to ",round(ci90upSlope,3))) %>%
  tidyr::unnest(cols = c(Trend))%>% #brings out as character but retains data
  dplyr::mutate(Trend = fct_explicit_na(Trend, na_level = "Not applicable"),
                Metric = "Long Term Flow") %>%
  write_csv("../tbl/trendFlowLT.csv")

## Rainfall (long term) data processing
divideBy <- nrow(as_tibble(modRainLT))
modRainResLT <- as_tibble(modRainLT) %>%
  dplyr::select(Date) %>%
  write_csv("../tbl/modRainResLT.csv")
trendRainLT <- modRainResLT %>%
  dplyr::summarise(negSlope = sum(Date<0)/divideBy #count slopes less than 0
                   , posSlope = sum(Date>0)/divideBy
                   , ci90loSlope = quantile(Date, 0.05)
                   , ci90upSlope = quantile(Date, 0.95)
                   , meanSlope = mean(Date)) %>%
  dplyr::mutate(Trend = map(negSlope,trend_set)
                , intervalSlope = paste0(round(ci90loSlope,3)," to ",round(ci90upSlope,3))) %>%
  tidyr::unnest(cols = c(Trend)) %>% #brings out as character but retains data
  dplyr::mutate(Trend = fct_explicit_na(Trend, na_level = "Not applicable"),
                Metric = "Long Term Rainfall") %>%
  write_csv("../tbl/trendRainLT.csv")

## RO data processing
divideBy <- nrow(as_tibble(modRO))
modRORes <- as_tibble(modRO) %>%
  dplyr::select(Date) %>%
  write_csv("../tbl/modRO.csv")
trendRO <- modRORes %>%
  dplyr::summarise(negSlope = sum(Date<0)/divideBy #count slopes less than 0
                   , posSlope = sum(Date>0)/divideBy
                   , ci90loSlope = quantile(Date, 0.05)
                   , ci90upSlope = quantile(Date, 0.95)
                   , meanSlope = mean(Date)) %>%
  dplyr::mutate(Trend = map(negSlope,trend_set)
                , intervalSlope = paste0(round(ci90loSlope,3)," to ",round(ci90upSlope,3))) %>%
  tidyr::unnest(cols = c(Trend)) %>% #brings out as character but retains data
  dplyr::mutate(Trend = fct_explicit_na(Trend, na_level = "Not applicable"),
                Metric = "Monthly Runoff") %>%
  write_csv("../tbl/trendRO.csv")

## RHigh flow days data processing
divideBy <- nrow(as_tibble(modBigFlow))
modBigFlowRes <- as_tibble(modBigFlow) %>%
  dplyr::select(year) %>%
  write_csv("../tbl/modBigFlowRes.csv")
trendBigFlow <- modBigFlowRes %>%
  dplyr::summarise(negSlope = sum(year<0)/divideBy #count slopes less than 0
                   , posSlope = sum(year>0)/divideBy
                   , ci90loSlope = quantile(year, 0.05)
                   , ci90upSlope = quantile(year, 0.95)
                   , meanSlope = mean(year)) %>%
  dplyr::mutate(Trend = map(negSlope,trend_set)
                , intervalSlope = paste0(round(ci90loSlope,3)," to ",round(ci90upSlope,3))) %>%
  tidyr::unnest(cols = c(Trend)) %>% #brings out as character but retains data
  dplyr::mutate(Trend = fct_explicit_na(Trend, na_level = "Not applicable"),
                Metric = "Big Flows") %>%
  write_csv("../tbl/trendBigFlows.csv")

# Combine all trend assessments into single table
OverallTrendLT <- bind_rows(trendFlowLT,
                            trendFlow,
                            trendRain,
                            trendRainLT,
                            trendRO,
                            trendBigFlow) %>%
  write_csv("../tbl/Overall trend results.csv")


### Impact Metric calculations

### Impact metric model development are undertaken in a seperate script
### This script simply uses the model
### predictions done against whole dataset for comparison

# Set working directory
setwd("C:/Ranalysis/BRT")

#Load libraries 
library(gbm)
library(ggplot2)
library(vegan)
library(mgcv)
library(lattice)
library(xts)
library(zoo)
library(lme4)
library(MASS)
library(tidyverse)
library(fs)
library(lubridate)
library("rstan")
library("rstanarm")
library("doParallel")
library(naniar)


# get functions for BRT modelling
source("brt.functions.R")

#read in spatial data to be predicted against
Pred.data <- read.csv("C:/Ranalysis/BRT/Pred_Data_incBurra.csv", as.is=T)
# linking file for river names 
RiverLink <- read_csv("RiverLink.csv")

#####Pred Data transform
#use the same transform as for the input data
#create dataframe for transformed data
Pred.data.transform <- Pred.data

Pred.data.transform$PV0305_SUM_Far <- log(as.numeric(Pred.data$PV0305_SUM_Far))
Pred.data.transform$PV0608_MEAN_Far <- log(Pred.data$PV0608_MEAN_Far)
Pred.data.transform$PV0608_SUM_Far <- log(Pred.data$PV0608_SUM_Far)
Pred.data.transform$PV0911_MEAN_Far <- Pred.data$PV0911_MEAN_Far^3
Pred.data.transform$PV0911_SUM_Far <- log(Pred.data$PV0911_SUM_Far)
Pred.data.transform$PV1202_SUM_Far <- log(Pred.data$PV1202_SUM_Far)
Pred.data.transform$PV2008.2017_MEAN_Far <- Pred.data$PV2008.2017_MEAN_Far^2
Pred.data.transform$PV2008.2017_SUM_Far <- log(Pred.data$PV2008.2017_SUM_Far)
Pred.data.transform$SA_Evap_SUM_Far <- log(Pred.data$SA_Evap_SUM_Far)
Pred.data.transform$Soil_pH_MEAN_Far <- log(Pred.data$Soil_pH_MEAN_Far)
Pred.data.transform$Soil_pH_SUM_Far <- log(Pred.data$Soil_pH_SUM_Far)
Pred.data.transform$Rain_SUM_Far <- log(Pred.data$Rain_SUM_Far)
Pred.data.transform$Soil_txt_SUM_Far <- log(Pred.data$Soil_txt_SUM_Far)
Pred.data.transform$Temp_SUM_Far <- log(Pred.data$Temp_SUM_Far)
Pred.data.transform$Water_erros_0305_MEAN_Far <- sqrt(Pred.data$Water_erros_0305_MEAN_Far)
Pred.data.transform$Water_erros_0305_SUM_Far <- log(Pred.data$Water_erros_0305_SUM_Far)
Pred.data.transform$Water_erros_0608_MEAN_Far <- sqrt(Pred.data$Water_erros_0608_MEAN_Far)
Pred.data.transform$Water_erros_0608_SUM_Far <- log(Pred.data$Water_erros_0608_SUM_Far)
Pred.data.transform$Water_erros_0911_MEAN_Far <- sqrt(Pred.data$Water_erros_0911_MEAN_Far)
Pred.data.transform$Water_erros_0911_SUM_Far <- log(Pred.data$Water_erros_0911_SUM_Far)
Pred.data.transform$Water_erros_1202_SUM_Far <- sqrt(Pred.data$Water_erros_1202_SUM_Far)
Pred.data.transform$PV0305_SUM_Near <- log(Pred.data$PV0305_SUM_Near)
Pred.data.transform$PV0608_MEAN_Near <- log(Pred.data$PV0608_MEAN_Near)
Pred.data.transform$PV0608_SUM_Near <- log(Pred.data$PV0608_SUM_Near)
Pred.data.transform$PV0911_MEAN_Near <- Pred.data$PV0911_MEAN_Near^3
Pred.data.transform$PV0911_SUM_Near <- log(Pred.data$PV0911_SUM_Near)
Pred.data.transform$PV1202_SUM_Near <- log(Pred.data$PV1202_SUM_Near)
Pred.data.transform$PV2008.2017_MEAN_Near <- Pred.data$PV2008.2017_MEAN_Near^2
Pred.data.transform$PV2008.2017_SUM_Near <- log(Pred.data$PV2008.2017_SUM_Near)
Pred.data.transform$SA_Evap_SUM_Near <- log(Pred.data$SA_Evap_SUM_Near)
Pred.data.transform$Soil_pH_MEAN_Near <- log(Pred.data$Soil_pH_MEAN_Near)
Pred.data.transform$Soil_pH_SUM_Near <- log(Pred.data$Soil_pH_SUM_Near)
Pred.data.transform$Rain_SUM_Near <- log(Pred.data$Rain_SUM_Near)
Pred.data.transform$Soil_txt_SUM_Near <- log(Pred.data$Soil_txt_SUM_Near)
Pred.data.transform$Temp_SUM_Near <- log(Pred.data$Temp_SUM_Near)
Pred.data.transform$Water_erros_0305_MEAN_Near <- sqrt(Pred.data$Water_erros_0305_MEAN_Near)
Pred.data.transform$Water_erros_0305_SUM_Near <- log(Pred.data$Water_erros_0305_SUM_Near)
Pred.data.transform$Water_erros_0608_MEAN_Near <- sqrt(Pred.data$Water_erros_0608_MEAN_Near)
Pred.data.transform$Water_erros_0608_SUM_Near <- log(Pred.data$Water_erros_0608_SUM_Near)
Pred.data.transform$Water_erros_0911_MEAN_Near <- sqrt(Pred.data$Water_erros_0911_MEAN_Near)
Pred.data.transform$Water_erros_0911_SUM_Near <- log(Pred.data$Water_erros_0911_SUM_Near)
Pred.data.transform$Water_erros_1202_SUM_Near <- sqrt(Pred.data$Water_erros_1202_SUM_Near)
Pred.data.transform$Water_erros_2008_2017_MEAN_Near <- sqrt(Pred.data$Water_erros_2008_2017_MEAN_Near)
Pred.data.transform$DamVol_ML <- log(Pred.data$DamVol_ML)
Pred.data.transform$DamDens <- sqrt(Pred.data$DamDens)


#load final model for calculations

load(file = "winning.model.rda") #will be called KIaddon.model.simp

# Check model imported ok
KIaddon.model.simp$contributions
KIaddon.model.simp$n.trees
KIaddon.model.simp$cv.statistics
KIaddon.model.simp$self.statistics

#### Predict #####

#use the model to predict the results agaisnt the rest of the dataset
preds <- predict.gbm(KIaddon.model.simp, Pred.data.transform, 
                     n.trees=KIaddon.model.simp$gbm.call$best.trees, 
                     type="response")
#NA in data error is OK

#correct tranform of data and rescale to 0 to 1
preds <- (preds^2)/(sqrt(3))
Scaled.IMP.MET <- Pred.data.transform$IMP_Met/(sqrt(3))

pred.export <- cbind(Pred.data.transform,Scaled.IMP.MET,preds)

#linear regression for correction
Lin_Reg <- lm(preds~Scaled.IMP.MET, pred.export)

#save slope and intercept as values for correction
slope<- Lin_Reg$coefficients[2]
inter<- Lin_Reg$coefficients[1]

#apply correction
preds.corr <- (preds/slope) - inter

#remove negitive values
preds.corr[preds.corr<0]<-0
pred.export <- cbind(Pred.data.transform,Scaled.IMP.MET,preds,preds.corr)

#Export final values 
write.csv(pred.export, file = "C:/Ranalysis/BRT/Imp_Met_predictions_incBurra.csv")


#### figure for reporting for M8 #####

datSum <- as_tibble(pred.export) %>%
  dplyr::left_join(RiverLink)
datSum %>%
  ggplot() +
  geom_boxplot(aes(River,preds.corr, color = Region)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #  facet_grid(~Region) +
  xlab("River") +
  ylab("Impact Metric") +
  labs(color = "Flow Data")

ggsave("Impact Metric boxplot for MR M8.PNG")




