# library ----

library(tidyverse)
library(lubridate)
library(survival)
library(magrttr)
# Import data

DATA = read.csv("C:/Users/conducta/Dropbox/PHD/DATA/Csv/RawData4.5.csv",sep =",",header = T)
source("C:/Users/conducta/Dropbox/PHD/Custom functions/Censoring.R")

# Cox Matrix
Cox_Matrix <- DATA %>%
  rename(Clutch = PUESTA, Brood = NIDADA) %>% 
  mutate_at(vars(ECLOSION1, ECLOSION2, FECHFINAL1, FECHFINAL2, ESTECLOS2), ymd) %>% 
  mutate(DateInitial = ESTECLOS2) %>% # Ref date
  filter(SEMANIPULO == "f") %>%
  filter(!FECHFINAL1 < ESTECLOS2) %>% # filter only nest with cohabitation
  filter(Brood == 2) %>% filter(Clutch == 2) %>%
  filter(!(MURIO2 %in% c("t") & LastAge2 %in% c(0:5) & SENCONTRO2 %in% c("f"))) %>% # Remove depredation
  filter(!(MURIO1 %in% c("t") & LastAge1 %in% c(0:5) & SENCONTRO1 %in% c("f"))) %>%
  Censoring(., 60) %>%
  rowwise() %>% 
  mutate(LastDateBroodReduc = case_when(Death601 == 1 & Death602 == 0 ~ (DateInitial + Time602),
                                        Death601 == 0 & Death602 == 1 ~ (DateInitial + Time601),
                                        Death601 == 1 & Death602 == 1 ~ (max(DateInitial + Time601, DateInitial + Time602)),
                                        Death601 == 0 & Death602 == 0 ~ (DateInitial + 60))) %>% 
  mutate(FirstDateBroodReduc = case_when(Death601 == 1 & Death602 == 0 ~ FECHFINAL1,
                                         Death601 == 0 & Death602 == 1 ~ FECHFINAL2,
                                         Death601 == 1 & Death602 == 1 ~ (min(FECHFINAL1, FECHFINAL2)),
                                         Death601 == 0 & Death602 == 0 ~ (DateInitial + 60))) %>% 
  ungroup() %>% 
  dplyr::mutate(LastTime = as.numeric(abs(difftime(DateInitial, LastDateBroodReduc , units = "days")))) %>% 
  dplyr::mutate(FirstTime = as.numeric(abs(difftime(DateInitial, FirstDateBroodReduc , units = "days")))) %>% 
  mutate(FirstTime = if_else(FirstTime == 0, 1, FirstTime)) %>% # Fix issue that seems like chick 2 born death 
  mutate(BroodReduction = case_when(Death601 == 1 | Death602 == 1 ~ 1,
                                    Death601 == 0 & Death602 == 0 ~ 0)) %>% 
  drop_na(BroodReduction, LastTime, FirstTime) %>% glimpse()

Cox_Matrix_add <- Cox_Matrix %>% mutate(HatchAsync = ECLOSION2 - ECLOSION1) %>% # add some others important variables for the analysis
  mutate_at(vars(HatchAsync), as.numeric) %>%
  filter(!(HatchAsync < 0 | HatchAsync > 25)) %>%
  filter(!(LastAge1 < 0 | LastAge2 < 0)) %>%
  filter(!FECHFINAL1 < ESTECLOS2) %>%
  mutate(LayingRoundRank = round(LayRank, digits = 1)) %>%
  mutate(LayingCatRank = case_when(
    LayingRoundRank %in% c(0.0,0.1) ~ 1,
    LayingRoundRank %in% c(0.2) ~ 2,
    LayingRoundRank %in% c(0.3) ~ 3,
    LayingRoundRank %in% c(0.4,0.5) ~ 4,
    LayingRoundRank %in% c(0.6:1) ~ 5)) %>%
  mutate(HatchingRoundRank = round(HatchRank, digits = 1)) %>%
  mutate(HatchingCatRank = case_when(
    HatchingRoundRank %in% c(0.0,0.3) ~ 1,
    HatchingRoundRank %in% c(0.3) ~ 2,
    HatchingRoundRank %in% c(0.4) ~ 3,
    HatchingRoundRank %in% c(0.5,0.6) ~ 4,
    HatchingRoundRank %in% c(0.7:1) ~ 5)) %>%
  mutate(HatchCat = case_when(
    HatchAsync %in% c(0:2) ~ 1,
    HatchAsync %in% c(3:5) ~ 2,
    HatchAsync %in% c(6:8) ~ 3,
    HatchAsync %in% c(9:11) ~ 4,
    HatchAsync %in% c(12:14) ~ 5)) %>%
  select(WORKYEAR, BroodReduction, LastDateBroodReduc, FirstDateBroodReduc, LastTime, FirstTime,
         NIDO, DateInitial, HatchAsync,
         AGEHEMB, ANILLOHEMB, AGEMACH, ANILLOMACH, CoupleExpLocal, ESTECLOS1,
         HatchRank, HatchCat,HatchingRoundRank, HatchingCatRank,
         PROPORTIONALRANK, LayingCatRank, LayingRoundRank, LayRank, CONFIRMHEMB) %>%
  mutate(across(where(is.Date), as.character))

# Analysis of multiple failure-time data (Unordered failure events of the same type) @Prentice, Williams, and Peterson (1981)

LastCountingProcess = Cox_Matrix_add %>% survSplit(formula = Surv(LastTime, BroodReduction) ~ .,
                                              cut = seq(from = 0, to = 60, by = 7), event = "BroodReduction",
                                              start = "Start_time", end = "Stop_time") %>% select(-FirstTime)

FirstCountingProcess = Cox_Matrix_add %>% survSplit(formula = Surv(FirstTime, BroodReduction) ~ .,
                                               cut = seq(from = 0, to = 60, by = 7), event = "BroodReduction",
                                               start = "Start_time", end = "Stop_time") %>% select(NIDO, Start_time, WORKYEAR, BroodReduction)



DataCountingProcess <- left_join(LastCountingProcess, FirstCountingProcess, by = c("NIDO", "Start_time", "WORKYEAR")) %>%  # merge both counting data
  rowwise() %>%
  mutate(BroodReduction = max(c_across(c(BroodReduction.x, BroodReduction.y)), na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(NIDO) %>% mutate(Time = 1:n()) %>%
  mutate_at(vars(DateInitial), ymd) %>%
  mutate(DATE = min(DateInitial) + Stop_time) %>%
  group_by(NIDO) %>% 
  mutate(Event = lag(cumsum(BroodReduction))) %>% 
  mutate(Event = case_when(Start_time == 0 ~ 0,
                           is.na(Event) ~ lead(Event),
                           TRUE ~ Event)) %>%
  ungroup() %>% 
  group_by(NIDO, Event) %>% 
  mutate(Gap = 1:n()) %>% 
  mutate(DateForPropRank = DATE - 7) %>%
  mutate(Time0 = Time) %>% mutate(Time = Time + 1) %>%
  group_by(WORKYEAR) %>%
  mutate(ESTECLOS1 = ymd(ESTECLOS1)) %>%
  mutate(Ndays = as.numeric(abs(difftime(min(ESTECLOS1, na.rm = T), DateForPropRank, units = "days")))) %>%
  mutate(CountingPropRank = Ndays/max(Ndays, na.rm = T))




write.csv(DataCountingProcess, file = "ClimMultiFailCox.csv", row.names = F)




# Diagnostic data in cox model
  

  
# Test data in model
DataCountingProcess %>% names()
DataCountingProcess %>% mutate(across(NIDO, as.factor)) 
baseline <- coxph(Surv(Time0, Time, BroodReduction) ~ HatchAsync + PROPORTIONALRANK +
                    Event + strata(Event), data = DataCountingProcess)

summary(baseline)  
# Diagnostic


# check_collinearity(Coxme)
rms::vif(baseline)


library(survminer)

zph <- cox.zph(baseline) ; zph
ggcoxzph(zph)

aa_fit <-aareg(Surv(Time0, Time, BroodReduction) ~
                 Gap + Event,
               data = DataCountingProcess)



# autoplot(aa_fit)
# Final Checking influential observations

# type = "dfbeta" , type = "deviance"

survminer::ggcoxdiagnostics(baseline, type = , linear.predictions = TRUE)

ggcoxdiagnostics(baseline, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())

ggcoxdiagnostics(survminer, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw()) # We see more negative values correspond to individual that "lived too long". Only problem with diagnostic

# Identify outliers
library(coxrobust)
res_mart = resid(Cox,type="martingale") %>% stack()# martingale or deviance
res_mart %>% #dplyr::slice_max(values, n = 100) %>%
  ggplot(aes(ind, values)) + geom_point()

# Final Checking non linearity

ggcoxfunctional(Surv(Time0, Time, BroodReduction) ~
                  PROPORTIONALRANK +
                  # HatchAsync +
                  RefittedChl +
                  I(RefittedChl^2) +
                  RefittedSST +
                  frailty(WORKYEAR),
                data = CoxData)
  
  
  


# old
# 
# LastDataCountingProcess = Data %>%
#   timeSplitter(by = 7, # Instead of 3 I need to choose 7 because i'm using climwin as weeks and not as days
#                event_var = "BroodReduction",
#                time_var = "LastTime",
#                event_start_status = "0") %>% select(-FirstTime)
# 
# FirstDataCountingProcess = Data %>%
#   timeSplitter(by = 7, # Instead of 3 I need to choose 7 because i'm using climwin as weeks and not as days
#                event_var = "BroodReduction",
#                time_var = "FirstTime",
#                event_start_status = "0") %>% select(-LastTime, -DateInitial, -ESTECLOS1)
