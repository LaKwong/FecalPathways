### Figure out which distribution to use
# normal probably fits pretty well, but must be bounded at 0 because neither SA nor conc or load can be neg


######### Simple soil model ##############

### Always calculate for ONE hand

### Since the mouthing frequencies do not follow distinct age group categories, group by month of age instead. This will create small n in many categories ###
## How to account for correlation within the same child? For now, ignore

library(readxl)
library(tidyr)
library(mc2d)
library(stats)
library(foreign)
library(triangle)
library(ggplot2)
library(data.table)
library(stringr)
library(reshape)
library(fitdistrplus)
library(logspline)
library(plyr)
library(dplyr)
#library(compositions)

# Load soil mass data

soilMass.base <- read_excel("C:/Users/Laura Kwong/Box Sync/VO Soil/Mass of soil on children's hands/Soil mass_RO1_analysis.xlsx", sheet = "final")
soilMass.base <- soilMass.base[1:272,c("hh", "RO1.round", "survey.day", "survey.month", "rep.dup", "motherOrChild", "volTotal.ml", "volSample.ml", "filterMass.g", "filterSoilMass.g")]
soilMass.base$survey.yr <- 2016
soilMass.base$surveyDate.soil <- dmy(paste(as.numeric(soilMass.base$survey.day), as.numeric(soilMass.base$survey.month), soilMass.base$survey.yr, sep = "-"))
soilMass.base$soil.g <- soilMass.base$filterSoilMass.g - soilMass.base$filterMass.g
soilMass.base$soil.mg <- (soilMass.base$filterSoilMass.g - soilMass.base$filterMass.g)*1000

#### Think about the detection limit! The precision of the scale is 5 mg = 0.005 g, so anything less than 0.005 g should be treated as...
## THIS ASSUMPTION WILL MAKE A BIG DIFFERENCE
## How many samples are below the detection limit of 0.005 g?
qplot(soilMass.base$soil.mg, geom="histogram", binwidth = 1)

sum(soilMass.base$soil.mg <= 5)/length(soilMass.base$soil.mg) # 57.0% of observations (including the reps) are less than the LOD

## Can try 1/2 the detection limit, though I remember hearing that there are better ways to handle this. 

##################################################################
###################################################################
## For now I will use this simple method and dig more into it later. 
soilMass.base$soil.mg.halfLOD <- ifelse(soilMass.base$soil.mg <= 5, 2.5, soilMass.base$soil.mg)
qplot(soilMass.base$soil.mg.halfLOD, geom="histogram", binwidth = 1)

#####################################################################
####################################################################


# AFTER replacing the values below the detection limit with 1/2 the detection limit, average the reps
soilMass.base.reps <- soilMass.base %>%
  filter(!is.na(rep.dup))

soilMass.rep.means <- soilMass.base.reps %>%
  group_by(hh, motherOrChild) %>%
  #group_by_(setdiff(names(soilMass.base), "rep.dup")) %>% # uses "group_by_" so I can use quoted columns; groups by all col except rep.dup
  summarise(soil.mg.halfLOD = mean(soil.mg.halfLOD))

# in soilMass.base, replace the reps with the rep mean
# keep only REP 1, set the rep.dup col to null, remove the soil.mg col and merge with the soilMass.rep.means dataset by hh and motherOrChild to add the soil.mg col based on the means of the reps
soilMass.reps <- soilMass.base.reps %>%
  filter(rep.dup == "REP1") %>%
  select(-soil.mg.halfLOD) %>%
  left_join(soilMass.rep.means[,c("hh", "motherOrChild", "soil.mg.halfLOD")], by = c("hh", "motherOrChild"))
soilMass.reps$rep.dup = "NA"

soilMass.NOreps <- soilMass.base %>%
  filter(is.na(rep.dup))

# Add the rep means back to the dataset
soilMass <- rbind(soilMass.NOreps, soilMass.reps) %>%
  arrange(hh)

# In case the rep mean was <0.005g, replace again with 1/2 LOD
##################################################################
###################################################################
## For now I will use this simple method and dig more into it later. 
soilMass$soil.mg.halfLOD <- ifelse(soilMass$soil.mg < 5, 2.5, soilMass$soil.mg)
qplot(soilMass$soil.mg.halfLOD, geom="histogram", binwidth = 1)

#####################################################################
####################################################################

# Calculate the mass of soil on one hand by finding the mass / ml of sample tested, multiplying by the total sample volume (250 mL for children and 350 mL for mothers) and dividing by 2 to get the mass on one hand
soilMass$soilOneHand.mg <- ((soilMass$soil.mg.halfLOD/soilMass$volSample.ml) * soilMass$volTotal.ml)/2
qplot(soilMass$soilOneHand.mg, geom="histogram", binwidth = 1)

soilMass.narrow <- soilMass[,c("hh", "RO1.round", "surveyDate.soil", "motherOrChild", "soilOneHand.mg")]


# Separate into mother and child mass
soilMass.child.Round7 <- soilMass.narrow %>%
  filter(motherOrChild == "CH", RO1.round == 7) %>%
  select(hh, RO1.round, surveyDate.soil, soilOneHand.mg)

soilMass.child.Round8 <- soilMass.narrow %>%
  filter(motherOrChild == "CH", RO1.round == 8) %>%
  select(hh, RO1.round, surveyDate.soil, soilOneHand.mg)

soilMass.mom.Round7 <- soilMass.narrow %>%
  filter(motherOrChild == "MH", RO1.round == 7) %>%
  select(hh, RO1.round, surveyDate.soil, soilOneHand.mg)

soilMass.mom.Round8 <- soilMass.narrow %>%
  filter(motherOrChild == "MH", RO1.round == 8) %>%
  select(hh, RO1.round, surveyDate.soil, soilOneHand.mg)

soilMass.wide <- soilMass.child.Round7 %>%
  left_join(soilMass.child.Round8, by = "hh") %>%
  left_join(soilMass.mom.Round7, by = "hh") %>%
  left_join(soilMass.mom.Round8, by = "hh")

names(soilMass.wide) <- c("hh", "RO1.c7", "surveyDate.soil.child.7", "soilOneHand.mg.child.7",
                          "RO1.c8", "surveyDate.soil.child.8", "soilOneHand.mg.child.8",
                          "RO1.m7", "surveyDate.soil.mom.7", "soilOneHand.mg.mom.7",
                          "RO1.m8", "surveyDate.soil.mom.8", "soilOneHand.mg.mom.8")
soilMass.wide <- soilMass.wide[, c("hh", "surveyDate.soil.child.7", "soilOneHand.mg.child.7",
                                   "surveyDate.soil.child.8", "soilOneHand.mg.child.8",
                                   "surveyDate.soil.mom.7", "soilOneHand.mg.mom.7",
                                   "surveyDate.soil.mom.8", "soilOneHand.mg.mom.8")]

###################################################################
###################################################################

####### In the file "WASHB anthro" #############
# Load anthropometric data
# Calculate a hand surface area based on body mass
# Calc age difference between when anthropometic data was collected and when soil mass was collected
# If the difference is small - assume that growth between the two time periods was negligible and use the hand surface area caclulated from the same child's measure height and weight 
# This is likely not a good measurement because young children grow rapidly
# If the difference is large - select the surface area from the distribution of children the same number of months old

anth.soil.wide <- read.csv("C:/Users/Laura Kwong/Box Sync/VO Fecal Intake Model/anth.soil.wide.csv")

anth.soil.wide <- anth.soil.wide %>%
  filter(arm == "Control")

anth.soilMass.wide <- left_join(anth.soil.wide, soilMass.wide, by = "hh")
# The dates had the correct format when saved but the format was not imported so need to change from factor to date
anth.soilMass.wide$unique_dob <- ymd(anth.soilMass.wide$unique_dob)
anth.soilMass.wide$surveyDate.mid <- ymd(anth.soilMass.wide$surveyDate.mid)
anth.soilMass.wide$surveyDate.end <- ymd(anth.soilMass.wide$surveyDate.end)

anth.soilMass.wide$age.soil.child.7 <- anth.soilMass.wide$surveyDate.soil.child.7 - anth.soilMass.wide$unique_dob
anth.soilMass.wide$age.soil.child.8 <- anth.soilMass.wide$surveyDate.soil.child.8 - anth.soilMass.wide$unique_dob
anth.soilMass.wide$round7DaysFromEnd <- anth.soilMass.wide$surveyDate.soil.child.7 - anth.soilMass.wide$surveyDate.end
anth.soilMass.wide$round8DaysFromEnd <- anth.soilMass.wide$surveyDate.soil.child.8 - anth.soilMass.wide$surveyDate.end


anth.soilMass.wide %>%
  select(age.soil.child.7) %>%
  filter(!is.na(age.soil.child.7)) / 365*12 # ~ most kids about 32-38 months old when sampled

anth.soilMass.wide %>%
  select(round8DaysFromEnd) %>%
  filter(!is.na(round8DaysFromEnd))
######################################################
######################################################
########## Bummer, the kids whos hand we sampled are WAY older than the kids in WASHB, so can't use the WASHB anthro data to est g/cm2
### How else can we get the SA of the hands - I think Ayse did not collect anthropometry data

### Maybe DHS data? - nope, birthweights only --> Ask Steve, kishor, mahbub?

######### For now, let's use the surface area of the oldest children in the WASHB dataset
## There are only 7 Control arm kids > 30 months at WASHB Endline
## Since there are only 7 kids, we'll just use the mean for now
anth.C.over30.hand.SA.WASHB <- anth.soilMass.wide %>%
  filter(survey.age.mo.end > 30) %>%
  select(C.hand.SA.end)

C.over30.hand.SA.WASHB <- anth.C.over30.hand.SA.WASHB %>%
  summarise(over30.hand.SA.mean = mean(C.hand.SA.end))

# Would be better if we can get month-specific est for hand surface area, but right now we can't so just average the mass of soil on child's hand over rounds 7 and 8
soilMass.wide$soilOneHand.mg.child.mean <- rowMeans(soilMass.wide[, c("soilOneHand.mg.child.7", "soilOneHand.mg.child.8")], na.rm = TRUE)

# Calc grams soil/cm2 using SA from WASHB data
soilMass.wide$C.hand.SA.est <- C.over30.hand.SA.WASHB[[1]]
soilMass.wide$C.soil.conc.mg.cm2 <- soilMass.wide$soilOneHand.mg.child.mean/soilMass.wide$C.hand.SA.est


######################################################################################
## We CAN use the WASHB anthropometry measurements for mothers. 
## Use endline values ONLY FOR THE MOTHERS WHOSE HANDS WERE SAMPLED, which are closer in date to the when soil on hands were measured

# # For the hand SA of all mothers in the WASHB dataset
# boxplot(anth.soilMass.wide[,"M.hand.SA.end"])  # clearly there are outliers (all values above 1000)
# 
# anth.M.hand.SA <- anth.soilMass.wide %>%
#   select(M.hand.SA.end) %>%
#   filter(M.hand.SA.end < 1000)
# 
# # Fitting of the distribution ' norm ' by maximum likelihood 
# # Parameters:
# #   estimate Std. Error
# # mean 322.94354  0.8323223
# # sd    27.04727  0.5885408

# For only the mothers in the WASHB dataset whose hands were sampled for soil
boxplot(anth.soilMass.wide[anth.soilMass.wide$hh %in% soilMass.wide$hh,"M.hand.SA.end"])  # clearly there are outliers (all values above 1000)

anth.M.hand.SA <- anth.soilMass.wide %>%
  filter(hh %in% soilMass.wide$hh) %>% 
  select(M.hand.SA.end)

M.hand.SA.end.dist <- fitdist(anth.M.hand.SA$M.hand.SA.end, "norm", method = "mle")
# mean = M.hand.SA.end.dist$estimate[1] 
# sd = M.hand.SA.end.dist$estimate[2] 

# Fitting of the distribution ' norm ' by maximum likelihood 
# Parameters:
#   estimate Std. Error
# mean 314.5447   4.412041
# sd    25.7264   3.119784

# Average the mass of soil on mother's hand over rounds 7 and 8
soilMass.wide$soilOneHand.mg.mom.mean <- rowMeans(soilMass.wide[, c("soilOneHand.mg.mom.7", "soilOneHand.mg.mom.8")], na.rm = TRUE)

# Calc grams soil/cm2 using SA from WASHB data
soilMass.wide$M.hand.SA.random <- rnorm(nrow(soilMass.wide), mean = M.hand.SA.end.dist$estimate[1], sd = M.hand.SA.end.dist$estimate[2] )
soilMass.wide$M.soil.conc.mg.cm2 <- soilMass.wide$soilOneHand.mg.mom.mean/soilMass.wide$M.hand.SA.random

################## Calc distribution of loading (mg) ###########
# Probablistic multiply by conc (g/cm2) by SA (cm2) of observed children to get est load g soil 
# -> If this method then need to # Determine a hand size distribution for children 3-6 mo old, 6-12, 12-24, 24-36 months old

C.hand.conc.mg.cm2.dist <- fitdist(soilMass.wide$C.soil.conc.mg.cm2, "norm", method = "mle")
M.hand.conc.mg.cm2.dist <- fitdist(soilMass.wide$M.soil.conc.mg.cm2, "norm", method = "mle")

# Use same conc mg/cm2 on hands for midline and endline load calc because assuming conc does not change with age
anth.soilMass.wide$C.hand.conc.random <- rnorm(nrow(anth.soilMass.wide), mean = C.hand.conc.mg.cm2.dist$estimate[1], sd = C.hand.conc.mg.cm2.dist$estimate[2] )
anth.soilMass.wide$C.hand.load.mid <- anth.soilMass.wide$C.hand.SA.mid * anth.soilMass.wide$C.hand.conc.random
anth.soilMass.wide$C.hand.load.end <- anth.soilMass.wide$C.hand.SA.end * anth.soilMass.wide$C.hand.conc.random

qplot(anth.soilMass.wide$C.hand.conc.random, geom = "histogram", binwidth = 0.005)
qplot(anth.soilMass.wide$M.hand.conc.random, geom = "histogram", binwidth = 0.005)

anth.soilMass.wide$M.hand.conc.random <- rnorm(nrow(anth.soilMass.wide), mean = M.hand.conc.mg.cm2.dist$estimate[1], sd = M.hand.conc.mg.cm2.dist$estimate[2] )
anth.soilMass.wide$M.hand.load.mid <- anth.soilMass.wide$M.hand.SA.mid * anth.soilMass.wide$M.hand.conc.random
anth.soilMass.wide$M.hand.load.end <- anth.soilMass.wide$M.hand.SA.end * anth.soilMass.wide$M.hand.conc.random


qplot(anth.soilMass.wide$C.hand.load.end, geom="histogram", binwidth = 1, xlim = c(-25, 80))
qplot(anth.soilMass.wide$M.hand.load.end, geom="histogram", binwidth = 1, xlim = c(-25, 80))

# scatter plot of child age (months) vs load on one hand (mg)
C.load.mid <- anth.soilMass.wide[, c("survey.age.mo.mid", "C.hand.load.mid")]
names(C.load.mid) <- c("age", "load")
C.load.end <- anth.soilMass.wide[, c("survey.age.mo.end", "C.hand.load.end")]
names(C.load.end) <- c("age", "load")

C.load.age <- rbind(C.load.mid, C.load.end)
scatter.smooth(C.load.age)

### Doesn't make sense to track mother age with loading - she is a month older but probably not grown


# Load HM frequencies and soil ingestion frequencies for each individual, not hh
## Structured Observation 
### Using data from all arms, justify by not sig diff (Kwong, 2016)

############################################
### Go back and figure out how in SO I calculated Soil_h_m_tot_freq --> what can I use this for?
###########################################3
SO.allObs <- read.csv("C:/Users/Laura Kwong/Dropbox/Fecal pathways/Structured obs/Analysis/SO_160310_submitted to IJERPH 160313/dem_obs_full_freq.csv")
names(SO.allObs)
SO.HM.SM <- SO.allObs %>%
  filter(location == 3) %>%
  mutate(round = "so.1", Mouth_hands_d = NA, Mouth_hands_nd = NA) %>%
  select(round, participant_id,  age_SO_mo, age_SO_group, hand_m_tot_freq, Mouth_hands_d, Mouth_hands_nd, soil_m_tot_freq, soil_h_m_tot_freq)
  
names(SO.HM.SM) <- c("round", "hh", "age", "age.group", "HM", "HM_d", "HM_nd", "SM", "SHM")

# ############ SO Data if I want to use the age group ################3
# ## From SO of 149 kids
# SO.HM.base <- read.csv("C:/Users/Laura Kwong/Dropbox/Fecal pathways/Structured obs/Analysis/SO_160310_submitted to IJERPH 160313/shape_scale_mean_med_SAVE.csv")
# # SO scale and shape for f6
# HM.SO.f6.scale <- SO.HM.base[SO.HM.base$X == "h.f6", "scale.all"]
# HM.SO.f6.shape <- SO.HM.base[SO.HM.base$X == "h.f6", "shape.all"]
# 
# # SO scale and shape for u6_12
# HM.SO.u6_12.scale <- SO.HM.base[SO.HM.base$X == "h.u6_12", "scale.all"]
# HM.SO.u6_12.shape <- SO.HM.base[SO.HM.base$X == "h.u6_12", "shape.all"]
# 
# # SO scale and shape for u24
# HM.SO.u24.scale <- SO.HM.base[SO.HM.base$X == "h.u24", "scale.all"]
# HM.SO.u24.shape <- SO.HM.base[SO.HM.base$X == "h.u24", "shape.all"]




# Load HM frequencies and soil ingestion frequencies for each individual, not hh
## Video Observations 
vo123.objclass.base <- read.csv("C:/Users/Laura Kwong/Box Sync/VO R123/vo.11.objclass.csv")
VO.HM.SM <- vo123.objclass.base %>%
  filter(actobj.class %in% c("Mouth_hands", "Mouth_hands_d", "Mouth_hands_nd", "Mouth_soil")) %>%
  select(vo.num, hh, age.vo, age.vo.group, actobj.class, freq) %>%
  spread(actobj.class, freq) %>%
  mutate(soil_h_m_tot_freq = NA)

names(VO.HM.SM) <- c("round", "hh", "age", "age.group", "HM", "HM_d", "HM_nd", "SM", "SHM")


HM.SM <- rbind(SO.HM.SM, VO.HM.SM)

scatter.smooth(HM.SM[,c("age", "HM")])
scatter.smooth(HM.SM[,c("age", "SM")])

################### Will it change anything if I do fractional months?

########### Assign a hand mouthing frequency to each child in the WASHB data set based ONLY on age in months ###############


###### Double check that in SO I didn't find any child or hh characteristics except age to be imp #####

#### Alternative methods to decide how the hand mouthing frequency should be set for children of diff ages (mo)

# Create a weibull distribution for each age.mo --> I'll try this method first
# But make wieull dist for children <6 mo because not many data point and this is the cutoff for exclusive breatfeeding - should have very little hands_d here. 
# And make a weibull dist for children >= 24 mo because lack data and expect more similar after this age

 ################### Don't have the data to show more similar from 24-36 mo old

# Alternative methods:
# Create a weibull distriution for each age.group and apply to each age group 
# Create a weibull distriution for each age.group and apply to each age month
# Assign the median or median freq by age.mo 

HM.SM.f6.dist <- fitdist(HM.SM[HM.SM$age < 6, "HM"], "weibull", method = "mle")
HM.SM.o24.dist <- fitdist(HM.SM[HM.SM$age >= 24, "HM"], "weibull", method = "mle")

# For the month-specific distributions, can I make weibull for each month?
HM.SM.dist <- function(x) {fitdist(HM.SM[HM.SM$age == x, "HM"], "weibull", method = "mle")}

for(i in 6:23){ #max(HM.SM$age)
  assign(paste("HM.SM.", i, ".dist", sep = ""), HM.SM.dist(i))
}

anth.soilMass.wide <- anth.soilMass.wide %>%
  select(-c(X.1, X)) %>%
  filter(!is.na(hh))

str(anth.soilMass.wide)
anth.soilMass.wide$HM.freq.mid <- NA
anth.soilMass.wide$HM.freq.end <- NA

# For some reason there is an NA row that is produced in the anth.soilMass.wide file that I cannot figure out whenceforth it comes 
# or how to get rid of it except by using & !is.na(survey.age.mo.mid)
attach(anth.soilMass.wide)

for(i in 6:23){
  anth.soilMass.wide[survey.age.mo.mid > i & survey.age.mo.mid < (i + 1) & !is.na(survey.age.mo.mid), "HM.freq.mid"] <- rweibull(nrow(anth.soilMass.wide[survey.age.mo.mid > i & survey.age.mo.mid < (i + 1) & !is.na(survey.age.mo.mid),]), shape = get(paste("HM.SM.", i, ".dist", sep = ""))$estimate[[1]], scale = get(paste("HM.SM.", i, ".dist", sep = ""))$estimate[[2]] )
  anth.soilMass.wide[survey.age.mo.end > i & survey.age.mo.end < (i + 1) & !is.na(survey.age.mo.mid), "HM.freq.end"] <- rweibull(nrow(anth.soilMass.wide[survey.age.mo.end > i & survey.age.mo.end < (i + 1) & !is.na(survey.age.mo.mid),]), shape = get(paste("HM.SM.", i, ".dist", sep = ""))$estimate[[1]], scale = get(paste("HM.SM.", i, ".dist", sep = ""))$estimate[[2]] )
}

# for(i in 6:23){
#   anth.soilMass.wide[anth.soilMass.wide$survey.age.mo.mid == i, "HM.freq.mid"] <- rweibull(nrow(anth.soilMass.wide[anth.soilMass.wide$survey.age.mo.mid == i,]), shape = get(paste("HM.SM.", i, ".dist", sep = ""))$estimate[[1]], scale = get(paste("HM.SM.", i, ".dist", sep = ""))$estimate[[2]] )
#   anth.soilMass.wide[anth.soilMass.wide$survey.age.mo.end == i, "HM.freq.end"] <- rweibull(nrow(anth.soilMass.wide[anth.soilMass.wide$survey.age.mo.end == i,]), shape = get(paste("HM.SM.", i, ".dist", sep = ""))$estimate[[1]], scale = get(paste("HM.SM.", i, ".dist", sep = ""))$estimate[[2]] )
# }

# 45 children < 6 mo at midline
anth.soilMass.wide[anth.soilMass.wide$survey.age.mo.mid < 6 & !is.na(survey.age.mo.mid), "HM.freq.mid"] <- rweibull(nrow(anth.soilMass.wide[anth.soilMass.wide$survey.age.mo.mid < 6 & !is.na(survey.age.mo.mid), ]), shape = HM.SM.f6.dist$estimate[[1]], scale = HM.SM.f6.dist$estimate[[2]] )
# 0 children < 6 mo at endline
anth.soilMass.wide[anth.soilMass.wide$survey.age.mo.end < 6 & !is.na(survey.age.mo.end), "HM.freq.end"] <- rweibull(nrow(anth.soilMass.wide[anth.soilMass.wide$survey.age.mo.end < 6 & !is.na(survey.age.mo.end), ]), shape = HM.SM.f6.dist$estimate[[1]], scale = HM.SM.f6.dist$estimate[[2]] )

# 1 child >= 24 mo at midline
anth.soilMass.wide[anth.soilMass.wide$survey.age.mo.mid >= 24 & !is.na(survey.age.mo.mid), "HM.freq.mid"] <- rweibull(nrow(anth.soilMass.wide[anth.soilMass.wide$survey.age.mo.mid >= 24 & !is.na(survey.age.mo.mid), ]), shape = HM.SM.o24.dist$estimate[[1]], scale = HM.SM.o24.dist$estimate[[2]] )
#210 children >= 24 mo at endline
anth.soilMass.wide[anth.soilMass.wide$survey.age.mo.end >= 24 & !is.na(survey.age.mo.end), "HM.freq.end"] <- rweibull(nrow(anth.soilMass.wide[anth.soilMass.wide$survey.age.mo.end >= 24 & !is.na(survey.age.mo.end), ]), shape = HM.SM.o24.dist$estimate[[1]], scale = HM.SM.o24.dist$estimate[[2]] )

detach(anth.soilMass.wide)

##################### Multiple HM freq by amt of soil on hands ######################
REally you should be using mc2d for this rather than assigning in an excel sheet. /......




###### Translate from ingestion per hour to ingestion per day based on number of hours awake per day ###############
######## This makes the assumption that the frequency of HM and SM soil mouthing averaged over the duration of observation represents the average over the entire day
##### Check this assumption by.......
  #### Rationalizing that we captured at least one (of how many?) feeding periods.



# Load hours awake per day (VO.R3 Status_161210)
obs.stat.v3.raw <- read_excel("C:/Users/Laura Kwong/Box Sync/VO R3/VO.R3 Status_170221.xlsx")
# Calc hours awake per day at diff age groups
#### Ignores within-child correlation of awake hours at different ages

hr.awake.base <- obs.stat.v3.raw %>%
  select(calc.age.mo, awake.day.hr, minus12.months, minus12.awake.time, minus24.months, minus24.awake.time)

soil.survey.v3 <- obs.stat.v3.raw %>%
  select(soil.today, soil.yest, soil.sevenday, soil.typicalday, soil.reaction, soil.healthy, soil.parent.eat, soil.parent.eat.day, soil.parent.eat.wk, soil.parent.eat.preg)

# The average number of hours awake was collected by recall on during the third observation
# Recall error could cause 