# ClimWin Cox Create data

# library ----

library(tidyverse)
library(lubridate)
library(Greg) # timeSplitter

# Import data
path = "conducta"

DATA = read.csv(str_glue("C:/Users/{path}/Dropbox/PHD/DATA/Csv/RawData4.5.csv"),sep =",",header = T)
source(str_glue("C:/Users/{path}/Dropbox/PHD/Custom functions/Censoring.R"))
# Check created proprank
DATA %>% mutate(ESTECLOS1 = ymd(ESTECLOS1)) %>% filter(WORKYEAR == 2019) %>% ggplot(aes(HatchRank, ESTECLOS1)) + geom_point()

# Create new data ----

# Cox Matrix
Cox_Matrix <- DATA %>%
  rename(Clutch = PUESTA, Brood = NIDADA) %>% 
  mutate_at(vars(ECLOSION1, ECLOSION2, FECHFINAL1, FECHFINAL2, ESTECLOS2), ymd) %>% 
  mutate(DateInitial = ESTECLOS2) %>% # Ref date
  filter(SEMANIPULO == "f") %>%
  filter(!FECHFINAL1 < ESTECLOS2) %>% # filter only nest with cohabitation
  filter(Brood == 2) %>% filter(Clutch == 2) %>%
  filter(!(MURIO2 %in% c("t") & LastAge2 %in% c(0:5) & SENCONTRO2 %in% c("f"))) %>%
  filter(!(MURIO1 %in% c("t") & LastAge1 %in% c(0:5) & SENCONTRO1 %in% c("f"))) %>%
  mutate(Death601 = case_when(MURIO1 == "f" & LastAge1 >= 60 ~ 0, # Right censoring
                              MURIO1 == "t" & LastAge1 < 60 ~ 1,
                              MURIO1 == "t" & LastAge1 >= 60 ~ 0,
                              MURIO1 == "f" & LastAge1 < 60 ~ NA_real_,
                              TRUE ~ NA_real_)) %>%
  mutate(Death602 = case_when(MURIO2 == "f" & LastAge2 >= 60 ~ 0,
                              MURIO2 == "t" & LastAge2 < 60 ~ 1,
                              MURIO2 == "t" & LastAge2 >= 60 ~ 0,
                              MURIO2 == "f" & LastAge2 < 60 ~ NA_real_,
                              TRUE ~ NA_real_)) %>%
  filter(!(Death601 == "t" & Death602 == "t")) %>% # CHECK THIS!!! If Adding maybe a variable first
  # mutate(N1DateBroodReduction = ifelse(FECHFINAL1 < FECHFINAL2, as.character(FECHFINAL1), 60)) %>%
  # mutate_at(vars(N1DateBroodReduction), ymd) %>%
  mutate(BroodReduction = ifelse(xor(Death601 == 1, Death602 == 1), 1, 0)) %>%
  mutate(DateBroodReduc1 = case_when(BroodReduction == 1 & Death601 == 1 ~ as.character(FECHFINAL1),
                                     BroodReduction == 1 & Death602 == 1 ~ as.character(FECHFINAL2),
                                     BroodReduction == 0 ~ as.character(DateInitial + 60),
                                     TRUE ~ NA_character_)) %>%
  # mutate(N2DateBroodReduction = ifelse(FECHFINAL1 < FECHFINAL2, as.character(FECHFINAL2), as.character(FECHFINAL1))) %>%
  # mutate(N1BroodReduction = ifelse(FECHFINAL1 < FECHFINAL2 & , 0, 1)) %>%
  # mutate(N2BroodReduction = ifelse(Brood == FLEDGED, 0, 1)) %>%
  mutate(HatchAsync = ECLOSION2 - ECLOSION1) %>%
  mutate_at(vars(HatchAsync), as.numeric) %>%
  filter(!(HatchAsync < 0 | HatchAsync > 25)) %>%
  filter(!(LastAge1 < 0 | LastAge2 < 0)) %>%
  filter(!FECHFINAL1 < ESTECLOS2) %>%
  mutate(DateBroodReduc1 = ymd(DateBroodReduc1)) %>%
  mutate(Time = as.numeric(abs(difftime(DateInitial, DateBroodReduc1 , units = "days")))) %>% # Check this!!!
  mutate_at(vars(DateInitial, DateBroodReduc1, BroodReduction), as.character) %>%
  filter_at(vars(BroodReduction, Time), all_vars(!is.na(.))) %>%
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
  select(WORKYEAR, BroodReduction, DateBroodReduc1, Time, NIDO, DateInitial, HatchAsync,
         AGEHEMB, ANILLOHEMB, AGEMACH, ANILLOMACH, CoupleExpLocal, ESTECLOS1,
         HatchRank, HatchCat,HatchingRoundRank, HatchingCatRank,
         PROPORTIONALRANK, LayingCatRank, LayingRoundRank, LayRank, CONFIRMHEMB) %>%
  mutate(ESTECLOS1= as.character(ESTECLOS1)) %>% glimpse()

# Plot to know who to create categories
# Data %>% ggplot(aes(HatchingRoundRank)) + geom_bar() + scale_x_continuous(n.breaks = 10)
# Data %>% filter(WORKYEAR == 2018) %>% ggplot(aes(HatchingRank, ESTECLOS1)) + geom_point()
# Transform data as counting process
DataCountingProcess = Cox_Matrix %>% mutate(Death = factor(BroodReduction,
                                                     levels = 0:1,
                                                     labels = c("Alive", "Death"))) %>%
  timeSplitter(by = 7, # Instead of 3 I need to choose 7 because i'm using climwin as weeks and not as days
               event_var = "BroodReduction",
               time_var = "Time",
               event_start_status = "0") %>% # timeSplitter doesn't accept dates so need to transform into characters
  group_by(NIDO) %>% mutate(Time = 1:n()) %>%
  mutate_at(vars(DateInitial), ymd) %>%
  mutate(DATE = min(DateInitial) + Stop_time) %>%
  mutate(DateForPropRank = DATE - 7) %>%
  mutate(Time0 = Time) %>% mutate(Time = Time + 1) %>%
  group_by(WORKYEAR) %>%
  mutate(ESTECLOS1 = ymd(ESTECLOS1)) %>%
  mutate(Ndays = as.numeric(abs(difftime(min(ESTECLOS1, na.rm = T), DateForPropRank, units = "days")))) %>%
  mutate(CountingPropRank = Ndays/max(Ndays, na.rm = T)) %>%
  glimpse()

write.csv(DataCountingProcess, file = "C:/Users/ivan/Dropbox/PHD/DATA/Csv/ClimCox.csv", row.names = F)



c <- read.csv("C:/Users/conducta/Dropbox/PHD/DATA/Csv/ClimCox.csv")
waldo::compare(c$ECLOSION2, as.character(Data$ECLOSION2))
waldo::compare(c$Time, as.integer(DataCountingProcess$Time))
