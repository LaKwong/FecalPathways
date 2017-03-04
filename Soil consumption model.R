### Figure out which distribution to use
# normal probably fits pretty well, but must be bounded at 0 because neither SA nor conc or load can be neg


### Exclude obj mouthing?

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

ndvar <- 101 #1001, 10001
ndunc <- 10
seed <- 123

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


### Measured soil load on children's hands
soil.C.load <- data.frame(combine(soilMass.wide$soilOneHand.mg.child.7, soilMass.wide$soilOneHand.mg.child.8))
names(soil.C.load) <- "soilOneHand.mg.child"
soil.C.load <- soil.C.load[!is.na(soil.C.load$soilOneHand.mg.child),]
qplot(soil.C.load, geom = "histogram", binwidth = 1)
qplot(log(soil.C.load), geom = "histogram", binwidth = 0.1)

# ### No distribution fits so instead I'll draw from the empirical data
# soil.C.load.dist <- fitdist(soil.C.load, "lnorm", method = "mle")
# soil.C.load.mcstoc <- mcstoc(rlnorm, type="V", meanlog = soil.C.load.dist$estimate[[1]], sdlog = soil.C.load.dist$estimate[[2]], rtrunc=TRUE, linf=0)

soil.C.load.mcdata <- mcdata(sample(soil.C.load, size=ndvar(), replace=TRUE),type="V")

# p0 = qplot(soil.C.load, geom = 'blank') +   
#   geom_line(aes(y = ..density.., colour = 'Empirical'), stat = 'density') +  
#   stat_function(fun = dnorm, aes(colour = 'Normal')) +                       
#   geom_histogram(aes(y = ..density..), alpha = 0.4) +                        
#   scale_colour_manual(name = 'Density', values = c('red', 'blue')) + 
#   theme(legend.position = c(0.85, 0.85))
# 
# print(p0)

### Measured soil load on mom's hands
soil.M.load <- data.frame(combine(soilMass.wide$soilOneHand.mg.mom.7, soilMass.wide$soilOneHand.mg.mom.8))
names(soil.M.load) <- "soilOneHand.mg.mom"
soil.M.load <- soil.M.load[!is.na(soil.M.load$soilOneHand.mg.mom),]
qplot(soil.M.load, geom = "histogram", binwidth = 1)

# ### No distribution fits so instead I'll draw from the empirical data
# soil.M.load.dist <- fitdist(soil.M.load, "lnorm", method = "mle")

soil.M.load.mcdata <- mcdata(sample(soil.M.load, size=ndvar(), replace=TRUE),type="V")


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
# Hand SA distribution for WASHB children OVER 30 months at endline
## There are only 7 Control arm kids > 30 months at WASHB Endline
anth.C.over30.hand.SA.WASHB <- anth.soilMass.wide %>%
  filter(survey.age.mo.end > 30) %>%
  select(C.hand.SA.end)
qplot(anth.C.over30.hand.SA.WASHB$C.hand.SA.end, geom = "histogram", binwidth = 1)

# # Dist doesn't fit so use empirical
# anth.C.over30.hand.SA.WASHB.dist <- fitdist(anth.C.over30.hand.SA.WASHB$C.hand.SA.end, "norm", method = "mle")
# C.over30.hand.SA.WASHB.mcstoc <- mcstoc(rweibull, type="V", shape = anth.C.over30.hand.SA.WASHB.dist$estimate[[1]], scale = anth.C.over30.hand.SA.WASHB.dist$estimate[[2]], rtrunc=TRUE, linf=0)

anth.C.over30.hand.SA.WASHB.mcdata <- mcdata(sample(anth.C.over30.hand.SA.WASHB$C.hand.SA.end, size=ndvar(), replace=TRUE),type="V")


# Hand SA distribution for WASHB mothers at endline
anth.M.hand.SA.WASHB.end <- anth.soilMass.wide %>%
  select(M.hand.SA.end) %>%
  filter(M.hand.SA.end < 1000) ## Filter out the unreasonable hand SA values - there are not outliers but rather mistakes with height or weight measurement
qplot(anth.M.hand.SA.WASHB.end$M.hand.SA.end, geom = "histogram", binwidth = 1)

# Dist DOES fit so use empirical
anth.M.hand.SA.WASHB.end.dist <- fitdist(anth.M.hand.SA.WASHB.end$M.hand.SA.end, "norm", method = "mle")
M.hand.SA.WASHB.end.mcstoc <- mcstoc(rnorm, type="V", mean = anth.M.hand.SA.WASHB.end.dist$estimate[[1]], sd = anth.M.hand.SA.WASHB.end.dist$estimate[[2]], rtrunc=TRUE, linf=0)

############## Est soil concentration on hands based on measured soil load for RO1 children and mothers and WASHB data for hand surface area ################
soil.C.conc <- soil.C.load.mcdata / anth.C.over30.hand.SA.WASHB.mcdata
soil.M.conc <- soil.M.load.mcdata / M.hand.SA.WASHB.end.mcstoc



##### To estimate the age month-specific load on the hands of WASHB kids #####
## Determine the month-specific hand surface area distribution for kids
C.hand.SA.WASHB <- anth.soilMass.wide %>%
  #select(survey.age.mo.mid, C.hand.SA.mid, survey.age.mo.end, C.hand.SA.end) %>%
  unite(mid, survey.age.mo.mid, C.hand.SA.mid) %>%
  unite(end, survey.age.mo.end, C.hand.SA.end) %>%
  select(hh, mid, end) %>%
  gather(key = round, value = age_hand.SA, mid, end) %>%
  separate(age_hand.SA, c("age", "C.hand.SA"), sep = "_", convert = TRUE) %>%
  mutate(age.mo = floor(age)) %>%
  select(hh, age.mo, C.hand.SA) %>%
  filter(!is.na(C.hand.SA))

scatter.smooth(C.hand.SA.WASHB$age, C.hand.SA.WASHB$C.hand.SA)

# ############# How about instead of making a distribution, you just use the empricial data
# # For the month-specific distributions, can I make normal/weibull for each month?
# C.hand.SA.WASHB.f6.dist <- fitdist(C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo < 6 & !is.na(C.hand.SA.WASHB$age.mo), "C.hand.SA"], "norm", method = "mle")
# C.hand.SA.WASHB.o24.dist <- fitdist(C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo >= 24 & !is.na(C.hand.SA.WASHB$age.mo), "C.hand.SA"], "norm", method = "mle")
# 
# hand.SA.find.dist <- function(data, x, dist) {fitdist(data[data$age.mo == x, "C.hand.SA"], dist, method = "mle")}
# 
# for(i in c(6:13, 18:23)){ #max(HM.SM$age)
#   assign(paste("C.hand.SA.WASHB.", i, sep = ""), hand.SA.find.dist(C.hand.SA.WASHB, i, "norm"))
# }


for(i in c(4:40)){ #max(HM.SM$age)
  assign(paste("C.hand.SA.WASHB.", i, sep = ""), C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo == i, "C.hand.SA"])
  assign(paste("C.hand.SA.WASHB.", i, ".mcdata", sep = ""), mcdata(sample(C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo == i, "C.hand.SA"], size = ndvar(), replace = TRUE), type = "V"))
}


## Determine the all-ages hand surface area distribution for mothers
M.hand.SA.WASHB <- anth.soilMass.wide %>%
  select(hh, M.hand.SA.mid, M.hand.SA.end) %>%
  gather(key = round, value = M.hand.SA, M.hand.SA.mid, M.hand.SA.end) %>%
  select(hh, M.hand.SA) %>%
  filter(M.hand.SA < 1000 & !is.na(M.hand.SA))

qplot(M.hand.SA.WASHB$M.hand.SA, geom = "histogram", binwidth = 10)

# Dist DOES fit so use empirical
M.hand.SA.WASHB.dist <- fitdist(M.hand.SA.WASHB$M.hand.SA, "norm", method = "mle")
M.hand.SA.WASHB.mcstoc <- mcstoc(rnorm, type="V", mean = M.hand.SA.WASHB.dist$estimate[[1]], sd = M.hand.SA.WASHB.dist$estimate[[2]], rtrunc=TRUE, linf=0)
#M.hand.SA.WASHB.mcdata <- mcdata(sample(M.hand.SA.WASHB$M.hand.SA, size=ndvar(), replace=TRUE),type="V")

################## Calc distribution of loading (mg) ###########
##### Est load on hands of all WASHB children using the est soil concentration on RO1 hands and the dist of hand SA of the WASHB hands
# Probablistic multiply by conc (g/cm2) by SA (cm2) of observed children to get est load g soil 
# Age month-specific for children
for(i in c(4:40)){
  assign(paste("C.soil.load.", i, sep = ""), soil.C.conc * get(paste("C.hand.SA.WASHB.", i, ".mcdata", sep = "")))
}

M.soil.load <- soil.M.conc * M.hand.SA.WASHB.mcstoc



######################## Add Hand-to-mouth contact data ################################
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
  select(participant_id,  round, age_SO_mo, age_SO_group, hand_m_tot_freq, Mouth_hands_d, Mouth_hands_nd, soil_m_tot_freq, soil_h_m_tot_freq)

names(SO.HM.SM) <- c("hh", "round", "age", "age.group", "HM", "HM_d", "HM_nd", "SM", "SHM")

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
  select(hh, vo.num, age.vo, age.vo.group, actobj.class, freq) %>%
  spread(actobj.class, freq) %>%
  mutate(soil_h_m_tot_freq = NA)

names(VO.HM.SM) <- c("hh", "round", "age", "age.group", "HM", "HM_d", "HM_nd", "SM", "SHM")


HM.SM <- rbind(SO.HM.SM, VO.HM.SM) %>%
  arrange(age, round, hh)

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

# ## For all hand_mouthing (use for round SO)
# HM.SM.f6.dist <- fitdist(HM.SM[HM.SM$age < 6, "HM"], "weibull", method = "mle")
# HM.SM.o24.dist <- fitdist(HM.SM[HM.SM$age >= 24, "HM"], "weibull", method = "mle")
# 
# # # For the month-specific distributions, can I make weibull for each month?
# find.dist <- function(data, x, dist) {fitdist(data[data$age == x, "HM"], dist, method = "mle")}
# 
# for(i in 6:23){ #max(HM.SM$age)
#   assign(paste("HM.SM.", i, ".dist", sep = ""), find.dist(HM.SM, i, "weibull"))
#   assign(paste("HM.SM.", i, ".mcdata", sep = ""), mcdata(sample(HM.SM[HM.SM$age == i, "HM"], size = ndvar(), replace = TRUE), type = "V"))
# }


#### Do NOT calc HM_d_child, HM_d_mother, HM_nd_child, HM_nd_mother on only the children in the video obs 
## Instead use HM * MONTH-SPECIFIC % HM_d_child, HM * MONTH-SPECIFIC % HM_d_mother, etc to estimate these values for all Structured observation kids. 
# Problem will be that ages don't align


## For now, use HM ONLY and do 25% mom, 75% child just to get an idea of the magnitude
## Probably best I can do is split by age group and take mean % d/nd, child/mom

Frac_HM_dueTo_Child <- 0.75
###########################
### Must change to reflect that most mom HM contacts do not have recontam so really only want number of indep feeding events
Frac_HM_dueTo_Mom <- 0.25 
####################


HM.SM$HM_child <- HM.SM$HM * Frac_HM_dueTo_Child
HM.SM$HM_mom <- HM.SM$HM * Frac_HM_dueTo_Mom

## For all HM_child
HM.SM.f6.child.dist <- fitdist(HM.SM[HM.SM$age < 6, "HM_child"], "weibull", method = "mle")
HM.SM.f6.child.mcstoc <- mcstoc(rweibull, type="V", shape = HM.SM.f6.child.dist$estimate[[1]], scale = HM.SM.f6.child.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.SM.o24.child.dist <- fitdist(HM.SM[HM.SM$age >= 24, "HM_child"], "weibull", method = "mle")
HM.SM.o24.child.mcstoc <- mcstoc(rweibull, type="V", shape = HM.SM.o24.child.dist$estimate[[1]], scale = HM.SM.o24.child.dist$estimate[[2]], rtrunc=TRUE, linf=0)

# # For the month-specific distributions, can I make weibull for each month?
find.dist <- function(data, x, dist) {fitdist(data[data$age == x, "HM_child"], dist, method = "mle")}

for(i in 6:23){ #max(HM.SM$age)
  assign(paste("HM.SM.", i, ".child.dist", sep = ""), find.dist(HM.SM, i, "weibull"))
  assign(paste("HM.SM.", i, ".child.mcdata", sep = ""), mcdata(sample(HM.SM[HM.SM$age == i, "HM_child"], size = ndvar(), replace = TRUE), type = "V"))
}

## For all HM_mom
HM.SM.f6.mom.dist <- fitdist(HM.SM[HM.SM$age < 6, "HM_mom"], "weibull", method = "mle")
HM.SM.f6.mom.mcstoc <- mcstoc(rweibull, type="V", shape = HM.SM.f6.mom.dist$estimate[[1]], scale = HM.SM.f6.mom.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.SM.o24.mom.dist <- fitdist(HM.SM[HM.SM$age >= 24, "HM_mom"], "weibull", method = "mle")
HM.SM.o24.mom.mcstoc <- mcstoc(rweibull, type="V", shape = HM.SM.o24.mom.dist$estimate[[1]], scale = HM.SM.o24.mom.dist$estimate[[2]], rtrunc=TRUE, linf=0)

# # For the month-specific distributions, can I make weibull for each month?
find.dist <- function(data, x, dist) {fitdist(data[data$age == x, "HM_mom"], dist, method = "mle")}

for(i in 6:23){ #max(HM.SM$age)
  assign(paste("HM.SM.", i, ".mom.dist", sep = ""), find.dist(HM.SM, i, "weibull"))
  assign(paste("HM.SM.", i, ".mom.mcdata", sep = ""), mcdata(sample(HM.SM[HM.SM$age == i, "HM_mom"], size = ndvar(), replace = TRUE), type = "V"))
}

#################### Combine load_child/load_mom with HM_child/HM_mom with mouth SA dist, removal efficiency dist

#HF	Child own hand fracion mouthed/event	-	day	beta	3.7	25				Zartarian 2005 as presented in Ozkaynak 2011	Zartarian 2005
HF.child.mcstoc <- mcstoc(rbeta, type="V", shape1 = 3.7, shape2 = 25, rtrunc=TRUE, linf=0)

# Children will mouth a smaller portion of mom's hand than of own hand 
############### I have no idea how much smaller, lets just say 50% less so if they mouth 12% of own hand, only mouth 6% of mom's hand
Frac_HF_for_mom <- 0.5
HF.mom.mcstoc <- HF.mcstoc * Frac_HF_for_mom

# HMRE	Hand mouthing removal(transfer) eff.	-	day	beta	2	8				Cohen Hubal 2008 as presented in Ozkaynak 2011
HMRE.mcstoc <- mcstoc(rbeta, type="V", shape1 = 2, shape2 = 8, rtrunc=TRUE, linf=0)

### For now, for children under 6 mo, use the soil loading est for children 4 mo old
C.soil.simple.f6 <- ((C.soil.load.4 * HM.SM.f6.child.mcstoc * HF.child.mcstoc) + (M.soil.load * HM.SM.f6.mom.mcstoc * HF.mom.mcstoc))  * HMRE.mcstoc
C.soil.simple.o24 <- ((C.soil.load.4 * HM.SM.o24.child.mcstoc * HF.child.mcstoc) + (M.soil.load * HM.SM.o24.mom.mcstoc * HF.mom.mcstoc))  * HMRE.mcstoc

for(i in 6:23){
  assign(paste("C.soil.simple.", i, sep = ""), ((get(paste("C.soil.load.", i, sep = "")) * get(paste("HM.SM.", i, ".child.mcdata", sep = ""))) + 
                                                  (M.soil.load * get(paste("HM.SM.", i, ".mom.mcdata", sep = "")))) * HF.mcstoc * HMRE.mcstoc)
  #C.soil.simple.10 <- ((C.soil.load.10 * HM.SM.10.child.mcdata) + (M.soil.load * HM.SM.10.mom.mcdata)) * HF.mcstoc * HMRE.mcstoc
}

########## Eval the mc model ###########
soil.test1 <- mcmodel({ 
  C.soil.simple.f6 <- ((C.soil.load.4 * HM.SM.f6.child.mcstoc * HF.child.mcstoc) + (M.soil.load * HM.SM.f6.mom.mcstoc * HF.mom.mcstoc))  * HMRE.mcstoc
  # C.soil.simple.o24 <- ((C.soil.load.4 * HM.SM.o24.child.mcstoc) + (M.soil.load * HM.SM.o24.mom.mcstoc)) * HF.mcstoc * HMRE.mcstoc
  # 
  # for(i in 6:23){
  #   assign(paste("C.soil.simple.", i, sep = ""), ((get(paste("C.soil.load.", i, sep = "")) * get(paste("HM.SM.", i, ".child.mcdata", sep = ""))) + 
  #                                                   (M.soil.load * get(paste("HM.SM.", i, ".mom.mcdata", sep = "")))) * HF.mcstoc * HMRE.mcstoc)
  #   #C.soil.simple.10 <- ((C.soil.load.10 * HM.SM.10.child.mcdata) + (M.soil.load * HM.SM.10.mom.mcdata)) * HF.mcstoc * HMRE.mcstoc
  # }
  res <- mc(C.soil.load.4, HM.SM.f6.child.mcstoc, HF.child.mcstoc, M.soil.load, HM.SM.f6.mom.mcstoc, HF.mom.mcstoc, HMRE.mcstoc, C.soil.simple.f6)
})

expr <- soil.test1
res <- evalmcmod(expr, nsv=ndvar(), nsu=ndunc(), seed=seed)
plot(res, prec = 0.001, stat = c("median", "mean"), lim = c(0.025, 0.25, 0.75, 0.975), na.rm =
       TRUE, griddim = NULL, xlab = NULL, ylab = "Fn(x)", main = "", draw = TRUE, paint = TRUE)
tor.res <- tornado(res)
plot(tor.res)
summary(res)
hist(res)


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





