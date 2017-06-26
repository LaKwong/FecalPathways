### Figure out which distribution to use
# normal probably fits pretty well, but must be bounded at 0 because neither SA nor conc or load can be neg

# What about HM and OM distributions that fit terribly (either nothing fits, non-weibull fits better but I feel that I should be consistent and use the same [weibull] for all)
# ex. line 980 ### all of the distributions fit terribly, normal sort of fits the best
# ex. HM.child.24_36 and HM.child 36_48

# How do I know the ratio of variable to uncertainty and use mcratio?
# 
# I think  I have  used mcdata wrong, instead of mcdata(data) I need to generate an 
# empirical distribution out of the data-- then it selects not just from the data, but from values within the distribution. How does it generate the distribution??

# in mcstoc(rempiricalD...) when do I need to specific nsv = and when not? why don't I need to specify nsu?


########## Ask Jenna or Steve
# ## Show the WASHB survey results for soil.today are comparable to the observed results to give confidence that we can trust the seven day results
# ## Observed consumption vs surveyed
# <6 mo:     2/24 =  8.3% vs today: 0.0%,  yest: 0.0%,  twoday: 0.0%,  sevenday: 0.0% ##### Review the <6 contacts and ensure they are actually mouthing - why didn't mother's report?
# 6-12 mo: 33/123 = 26.8% vs today: 38.5%, yest: 52.2%, twoday: 52.9%, sevenday: 55.2%
# 12-24 mo: 14/54 = 25.9% vs today: 35.0%, yest: 47.5%, twoday: 52.5%, sevenday: 60.0%
# 24-36 mo:   0/5 =  0.0% vs today: 10.6%, yest: 26.6%, twoday: 28.2%, sevenday: 42.3%
# 
# # Could use sevenday parental report as long-term average fraction of children consuming soil, BUT I am only modeling for ONE day so I should use only the one-day results. 
# # Otherwise I will overestimate the daily (though be closer to the long-term average?) 



#### To replace na using dplyr! ######333
# df %>% replace(is.na(.), 0)


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
library(lubridate)
library(irr)
library(dunn.test)
library(goftest)
library(plyr)
library(dplyr)
#library(compositions)
'%notin%' <- function(x, y) { !(x %in% y)}

ndvar <- 1001 #1001, 10001 = number of simulations in the variability dimension
ndunc <- 10 # = number of simulations in the uncertainty dimension
seed <- 888



#####################################################################################
#####################################################################################

################## Soil ingestion from contacts with hands #########################

#####################################################################################
######################################################################################


# Load soil mass data

soilMass.base <- read_excel("C:/Users/Tareq/Box Sync/VO Soil/Mass of soil on children's hands/Soil mass_RO1_analysis.xlsx", sheet = "final")
soilMass.base <- soilMass.base[1:272,c("hh", "RO1.round", "survey.day", "survey.month", "rep.dup", "motherOrChild", "volTotal.ml", "volSample.ml", "filterMass.g", "filterSoilMass.g")]
soilMass.base$survey.yr <- 2016
soilMass.base$surveyDate.soil <- dmy(paste(as.numeric(soilMass.base$survey.day), as.numeric(soilMass.base$survey.month), soilMass.base$survey.yr, sep = "-"))
soilMass.base$soil.g <- soilMass.base$filterSoilMass.g - soilMass.base$filterMass.g
soilMass.base$soil.mg <- (soilMass.base$filterSoilMass.g - soilMass.base$filterMass.g)*1000
soilMass.base$soil.mg.adjforLOD <- NA

#### Think about the detection limit! The precision of the scale is 5 mg = 0.005 g, so anything less than 0.005 g should be treated as...
## THIS ASSUMPTION WILL MAKE A BIG DIFFERENCE
## How many samples are below the detection limit of 0.005 g?
qplot(soilMass.base$soil.mg, geom="histogram", binwidth = 1)

sum(soilMass.base$soil.mg <= 5)/length(soilMass.base$soil.mg) # 57.0% of observations (including the reps) are less than the LOD

##################################################################
###################################################################
## Conduct sensitivity analysis to det effect of choosing what happens to values < LOD

# Option 1: Replace values < LOD with 1/2 the LOD
# belowLOD_replacementValue <- 2.5
# soilMass.base$soil.mg.halfLOD <- ifelse(soilMass.base$soil.mg <= 5, belowLOD_replacementValue, soilMass.base$soil.mg) # could also do < 5
# qplot(soilMass.base$soil.mg.halfLOD, geom="histogram", binwidth = 1)

# Option 2: Replace values < LOD with uniform selection 0 to LOD 
test <- soilMass.base %>%
  mutate(soil.mg.adjforLOD = ifelse(soil.mg <= 5 & motherOrChild == "CH" & volSample.ml == 50, round(runif(1, 0, 5), digits = 0),1000))
         
soilMass.base[soilMass.base$soil.mg <= 5, c("soil.mg.adjforLOD")] <- runif(nrow(soilMass.base[soilMass.base$soil.mg <= 5, c("soil.mg.halfLOD")]), 0, 5)
# code for making the LOD replacement - specific to Child vs Mother and aliquot tested - don't need to do this because this is captured in the soil on one hand calc where I divide by the size of the aliquot and multiply by the total sample vol - if the adjLOD for the sample is 5 then I can get from 
# soilMass.base[soilMass.base$soil.mg <= 5 & soilMass.base$motherOrChild == "CH" & soilMass.base$volSample.ml == 50, c("soil.mg.adjforLOD")] <- 
#   round(runif(nrow(soilMass.base[soilMass.base$soil.mg <= 5 & soilMass.base$motherOrChild == "CH" & soilMass.base$volSample.ml == 50, c("soil.mg.halfLOD")]), 0, 5), digits = 0)
# soilMass.base[soilMass.base$soil.mg <= 5 & soilMass.base$motherOrChild == "CH" & soilMass.base$volSample.ml == 200, c("soil.mg.adjforLOD")] <- 
#   round(runif(nrow(soilMass.base[soilMass.base$soil.mg <= 5 & soilMass.base$motherOrChild == "CH" & soilMass.base$volSample.ml == 200, c("soil.mg.halfLOD")]), 0, 5), digits = 0)
# soilMass.base[soilMass.base$soil.mg <= 5 & soilMass.base$motherOrChild == "MH" & soilMass.base$volSample.ml == 50, c("soil.mg.adjforLOD")] <- 
#   round(runif(nrow(soilMass.base[soilMass.base$soil.mg <= 5 & soilMass.base$motherOrChild == "MH" & soilMass.base$volSample.ml == 50, c("soil.mg.halfLOD")]), 0, 5), digits = 0)
# soilMass.base[soilMass.base$soil.mg <= 5 & soilMass.base$motherOrChild == "MH" & soilMass.base$volSample.ml == 200, c("soil.mg.adjforLOD")] <- 
#   round(runif(nrow(soilMass.base[soilMass.base$soil.mg <= 5 & soilMass.base$motherOrChild == "MH" & soilMass.base$volSample.ml == 200, c("soil.mg.halfLOD")]), 0, 5), digits = 0)
soilMass.base[soilMass.base$soil.mg > 5, c("soil.mg.adjforLOD")] <- soilMass.base[soilMass.base$soil.mg > 5, c("soil.mg")]
soilMass.base$soil.mg.adjforLOD
qplot(soilMass.base$soil.mg.adjforLOD, geom="histogram", binwidth = 1)

#####################################################################
####################################################################
# BEFORE replacing the values below the detection limit, average the reps - doesn't make sense to average a value that don't at all trust was accurate. Intead just take the (more) accurate value (the one > LOD)
soilMass.base.reps <- soilMass.base %>%
  filter(!is.na(rep.dup))

soilMass.rep.means <- soilMass.base.reps %>%
  dplyr::group_by(hh, motherOrChild) %>%
  #group_by_(setdiff(names(soilMass.base), "rep.dup")) %>% # uses "group_by_" so I can use quoted columns; groups by all col except rep.dup
  dplyr::summarise(soil.mg.adjforLOD = mean(soil.mg.adjforLOD, na.rm = TRUE))

# in soilMass.base, replace the reps with the rep mean
# keep only REP 1, set the rep.dup col to null, remove the soil.mg col and merge with the soilMass.rep.means dataset by hh and motherOrChild to add the soil.mg col based on the means of the reps
soilMass.reps <- soilMass.base.reps %>%
  filter(rep.dup == "REP1") %>%
  select(-soil.mg.adjforLOD) %>%
  left_join(soilMass.rep.means[,c("hh", "motherOrChild", "soil.mg.adjforLOD")], by = c("hh", "motherOrChild")) %>%
  mutate(rep.dup = "NA")

soilMass.NOreps <- soilMass.base %>%
  filter(is.na(rep.dup))

# Add the rep means back to the dataset
soilMass <- rbind(soilMass.NOreps, soilMass.reps) %>%
  arrange(hh)

# Calculate the mass of soil on one hand by finding the mass / ml of sample tested, multiplying by the total sample volume (250 mL for children and 350 mL for mothers) and dividing by 2 to get the mass on one hand
soilMass$soilOneHand.mg <- ((soilMass$soil.mg.adjforLOD/soilMass$volSample.ml) * soilMass$volTotal.ml)/2

## Now replace the samples that were < LOD with replacements that are specific to Child vs Mother and aliquot of 50 vs 200
soilMass.base[soilMass.base$soil.mg <= 5 & soilMass.base$motherOrChild == "CH" & soilMass.base$volSample.ml == 50, c("soil.mg.adjforLOD")] <-
  round(runif(nrow(soilMass.base[soilMass.base$soil.mg <= 5 & soilMass.base$motherOrChild == "CH" & soilMass.base$volSample.ml == 50, c("soil.mg.halfLOD")]), 0, 5), digits = 0)
soilMass.base[soilMass.base$soil.mg <= 5 & soilMass.base$motherOrChild == "CH" & soilMass.base$volSample.ml == 200, c("soil.mg.adjforLOD")] <-
  round(runif(nrow(soilMass.base[soilMass.base$soil.mg <= 5 & soilMass.base$motherOrChild == "CH" & soilMass.base$volSample.ml == 200, c("soil.mg.halfLOD")]), 0, 5), digits = 0)
soilMass.base[soilMass.base$soil.mg <= 5 & soilMass.base$motherOrChild == "MH" & soilMass.base$volSample.ml == 50, c("soil.mg.adjforLOD")] <-
  round(runif(nrow(soilMass.base[soilMass.base$soil.mg <= 5 & soilMass.base$motherOrChild == "MH" & soilMass.base$volSample.ml == 50, c("soil.mg.halfLOD")]), 0, 5), digits = 0)
soilMass.base[soilMass.base$soil.mg <= 5 & soilMass.base$motherOrChild == "MH" & soilMass.base$volSample.ml == 200, c("soil.mg.adjforLOD")] <-
  round(runif(nrow(soilMass.base[soilMass.base$soil.mg <= 5 & soilMass.base$motherOrChild == "MH" & soilMass.base$volSample.ml == 200, c("soil.mg.halfLOD")]), 0, 5), digits = 0)

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


############################################################################################################ 
## Instead of finding a dist for load on hands and then a dist for surface area of hands and dividing these dist to get concentration, 
## use empirical load divided by modeled hand SA based on age = concentration and then find the dist of that concentration   ################
#soil.C.load <- soilMass.wide[, c("hh", "child.hand.SA.age.soil.7", "soilOneHand.mg.child.7", "child.hand.SA.age.soil.7", "soilMass.wide$soilOneHand.mg.child.8")]
soil.C.load <- soilMass.wide[, c("hh", "soilOneHand.mg.child.7", "soilOneHand.mg.child.8")]

# ### Measured soil load on children's hands
# soil.C.load <- data.frame(combine(soilMass.wide$soilOneHand.mg.child.7, soilMass.wide$soilOneHand.mg.child.8))
# names(soil.C.load) <- "soilOneHand.mg.child"
# soil.C.load <- soil.C.load[!is.na(soil.C.load$soilOneHand.mg.child),]
# qplot(soil.C.load, geom = "histogram", binwidth = 1)
# qplot(log(soil.C.load), geom = "histogram", binwidth = 0.1)
# 
# # ### No distribution fits so instead I'll draw from the empirical data
# # soil.C.load.dist <- fitdist(soil.C.load, "lnorm", method = "mle")
# # soil.C.load.mcstoc <- mcstoc(rlnorm, type="V", meanlog = soil.C.load.dist$estimate[[1]], sdlog = soil.C.load.dist$estimate[[2]], rtrunc=TRUE, linf=0)
# 
# soil.C.load.mcdata <- mcdata(sample(soil.C.load, size=ndvar(), replace=TRUE),type="V")

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

length(soil.M.load) #130
summary(soil.M.load)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   3.500   7.000   8.461  11.380  49.000

qplot(soil.M.load, geom = "histogram", binwidth = 1)

# ### No distribution fits so instead I'll draw from the empirical data
# soil.M.load.dist <- fitdist(soil.M.load, "lnorm", method = "mle")

soil.M.load.mcstoc <- mcstoc(rempiricalD, values = soil.M.load, type="V", nsv = ndvar()) # I think this does the same thing as sample(soil.M.load, size=ndvar(), replace=TRUE)
# I think I can't make a continuous function out of this because it doesn't fit into a distribution
# soil.M.load.mcdata <- mcstoc(rempiricalC(n = ndvar(), min = min(soil.M.load), max = max(soil.M.load), values = soil.M.load), type="V")
#soil.M.load.mcdata <- mcdata(sample(soil.M.load, size=ndvar(), replace=TRUE),type="V")


###################################################################
###################################################################

######## Det the Hand SA of the children and mother's at the time they gave the hand rinse sample for soil ###########

##### First determine how old the children were when they gave the hand rinse sample for soil 

####### In the file "WASHB anthro" #############
# Load anthropometric data
# Calculate a hand surface area based on body mass
# Calc age difference between when anthropometic data was collected and when asoil mass was collected
# If the difference is small - assume that growth between the two time periods was negligible and use the hand surface area caclulated from the same child's measure height and weight 
# This is likely not a good measurement because young children grow rapidly
# If the difference is large - select the surface area from the distribution of children the same number of months old

#anth.soil.wide <- read.csv("C:/Users/Tareq/Box Sync/VO Fecal Intake Model/anth.soil.wide.csv")
anth.soil.wide <- read.csv("C:/Users/Tareq/Box Sync/VO Fecal Intake Model/anth.soil.wide_const_HandSAratio.csv")



anth.soil.wide <- anth.soil.wide %>%
  #filter(arm %in% c ("Control"))
  filter(arm %in% c ("Control", "Sanitation", "Handwashing", "Water", "WSH")) 
# I did't want to use all the sanitation arm children for fear that they have different growth than the control arm children (though I guess Steve's study shows that they don't have sig dif growth...)
# So I guess I will use all the data possible, just not the nutrition kids

missingWASHBdata <- setdiff(soilMass.wide$hh, anth.soil.wide$hh)
# There are 31 hh that do not have any WASHB midline or endline data so I don't know how old the kids were at the time of soil sampling. 
# [1] 36302 40502 41003 41004 41005 41007 42302 42305 42307 42906 43401 43404 43405 43407 45202 45903 45905 45907 45908 46701 46707 46708 47901 47903 47905 47907 48701 48704 48708 51005 51007
write.csv(missingWASHBdata, "C:/Users/Tareq/Box Sync/VO Soil/Mass of soil on children's hands/soilsample hh that are missing WASHB data.csv")

allsoilsamplehh <- unique(soilMass.wide$hh)
write.csv(allsoilsamplehh, "C:/Users/Tareq/Box Sync/VO Soil/Mass of soil on children's hands/all soil sample hh.csv")





anth.soilMass.wide <- left_join(anth.soil.wide, soilMass.wide, by = "hh")
# The dates had the correct format when saved but the format was not imported so need to change from factor to date
anth.soilMass.wide$unique_dob <- ymd(anth.soilMass.wide$unique_dob)
anth.soilMass.wide$surveyDate.mid <- ymd(anth.soilMass.wide$surveyDate.mid)
anth.soilMass.wide$surveyDate.end <- ymd(anth.soilMass.wide$surveyDate.end)

length(unique(anth.soil.wide$hh)) #3314 hh in the ("Control", "Sanitation", "Handwashing", "Water", "WSH") arms
sd(anth.soil.wide$survey.age.mo.mid, na.rm = TRUE)
sd(anth.soil.wide$survey.age.mo.end, na.rm = TRUE)

# # The code below shows that the kids who's hand we sampled are WAY older than the kids in WASHB, 
# # so can't use the WASHB anthro data to est g/cm2
# 
anth.soilMass.wide$age.soil.child.7 <- (anth.soilMass.wide$surveyDate.soil.child.7 - anth.soilMass.wide$unique_dob) / 365 * 12
anth.soilMass.wide$age.soil.child.8 <- (anth.soilMass.wide$surveyDate.soil.child.8 - anth.soilMass.wide$unique_dob) / 365 * 12
# anth.soilMass.wide$round7DaysFromEnd <- anth.soilMass.wide$surveyDate.soil.child.7 - anth.soilMass.wide$surveyDate.end # at the time the soil round 7 sample was collected, children were 360 year older than the were at WASHB endline 
# anth.soilMass.wide$round8DaysFromEnd <- anth.soilMass.wide$surveyDate.soil.child.8 - anth.soilMass.wide$surveyDate.end# at the time the soil round 8 sample was collected, children were 460 year older than the were at WASHB endline 
# 
# anth.soilMass.wide %>%
#   select(age.soil.child.7) %>%
#   filter(!is.na(age.soil.child.7)) / 365*12 # ~ most kids about 32-38 months old when sampled
# 
# anth.soilMass.wide %>%
#   select(round8DaysFromEnd) %>%
#   filter(!is.na(round8DaysFromEnd)) / 365*12


############### To get the child hand SA of the appropriate-age children,
## Could use the surface area of the oldest children in the WASHB dataset, but there are only 7 kids in the control arm
# # Hand SA distribution for WASHB 7 control children OVER 30 months at endline
# anth.C.over30.hand.SA.WASHB <- anth.soilMass.wide %>%
#   filter(survey.age.mo.end > 30) %>%
#   select(C.hand.SA.end)
# qplot(anth.C.over30.hand.SA.WASHB$C.hand.SA.end, geom = "histogram", binwidth = 1)
# # # Dist doesn't fit so use empirical
# # anth.C.over30.hand.SA.WASHB.dist <- fitdist(anth.C.over30.hand.SA.WASHB$C.hand.SA.end, "norm", method = "mle")
# # C.over30.hand.SA.WASHB.mcstoc <- mcstoc(rweibull, type="V", shape = anth.C.over30.hand.SA.WASHB.dist$estimate[[1]], scale = anth.C.over30.hand.SA.WASHB.dist$estimate[[2]], rtrunc=TRUE, linf=0)
# 
# anth.C.over30.hand.SA.WASHB.mcdata <- mcdata(sample(anth.C.over30.hand.SA.WASHB$C.hand.SA.end, size=ndvar(), replace=TRUE),type="V")

# ## Instead est the hand SA for the kids at the age when they gave the hand soil sample by using the trendline of age ~ hand SA from all the control arm WASHB kids
# #### We can analyze all the children in WASHB and see that there's a pretty linear trend between hand size and age so try modeling hand surface area based on age
# C.hand.SA.WASHB <- anth.soilMass.wide %>%
#   #select(survey.age.mo.mid, C.hand.SA.mid, survey.age.mo.end, C.hand.SA.end) %>%
#   unite(mid, survey.age.mo.mid, C.hand.SA.mid) %>%
#   unite(end, survey.age.mo.end, C.hand.SA.end) %>%
#   select(hh, mid, end) %>%
#   gather(key = round, value = age_hand.SA, mid, end) %>%
#   separate(age_hand.SA, c("age", "C.hand.SA"), sep = "_", convert = TRUE) %>%
#   mutate(age.mo = age) %>%#mutate(age.mo = floor(age)) %>%
#   select(hh, age.mo, C.hand.SA) %>%
#   filter(!is.na(C.hand.SA))
# 
# scatter.smooth(C.hand.SA.WASHB$age, C.hand.SA.WASHB$C.hand.SA)
# C.hand.SA.WASHB.lm <- lm(C.hand.SA ~ age.mo, C.hand.SA.WASHB)
# summary(C.hand.SA.WASHB.lm)$r.squared # 0.613 # yes, same as adj r squared value
# summary(C.hand.SA.WASHB.lm)$adj.r.squared # 0.613
# 

##### The r-squared of the linear trend is pretty low (0.613), so better to use the slope for the particular child. 
# True, that they could have changed their growth trajectory between endline and the RO1 samples, but I think this will be more accurate than averaging all the kids
C.hand.SA.WASHB <- anth.soilMass.wide %>%
  #select(survey.age.mo.mid, C.hand.SA.mid, survey.age.mo.end, C.hand.SA.end) %>%
  unite(mid, survey.age.mo.mid, C.hand.SA.mid) %>%
  unite(end, survey.age.mo.end, C.hand.SA.end) %>%
  select(hh, mid, end) %>%
  gather(key = round, value = age_hand.SA, mid, end) %>%
  separate(age_hand.SA, c("age", "C.hand.SA"), sep = "_", convert = TRUE) %>%
  mutate(age.mo = age) %>%#mutate(age.mo = floor(age)) %>%
  select(hh, age.mo, C.hand.SA) %>%
  filter(!is.na(C.hand.SA)) %>%
  # filter on the kids who gave hand rinse samples
  filter(hh %in% soilMass.wide$hh)

C.hand.SA.WASHB$hh <- as.factor(C.hand.SA.WASHB$hh)
ggplot(data = C.hand.SA.WASHB, 
       aes(x = age.mo, y = C.hand.SA)) +
  geom_point(aes(colour = hh)) + 
  geom_line(aes(colour = hh)) +
  ggtitle("Surface Area by Age in Months for Children Who Provided a Hand Rinse Sample for Soil") +
  theme(legend.position="none") +
  labs(x="Child age (months)", y="Surface Area (cm2)") 

anth.soilMass.wide$age.soil.child.7 <- as.numeric(anth.soilMass.wide$age.soil.child.7)
anth.soilMass.wide$age.soil.child.8 <- as.numeric(anth.soilMass.wide$age.soil.child.8)

anth.soilMass.wide$child.hand.SA.age.soil.7 <- predict(lm(C.hand.SA~age.mo, data = C.hand.SA.WASHB), newdata = data.frame(age.mo = anth.soilMass.wide$age.soil.child.7))
anth.soilMass.wide$child.hand.SA.age.soil.8 <- predict(lm(C.hand.SA~age.mo, data = C.hand.SA.WASHB), newdata = data.frame(age.mo = anth.soilMass.wide$age.soil.child.8))

C.hand.SA <- data.frame(c(anth.soilMass.wide[!is.na(anth.soilMass.wide$child.hand.SA.age.soil.7), "child.hand.SA.age.soil.7"], 
                          anth.soilMass.wide[!is.na(anth.soilMass.wide$child.hand.SA.age.soil.8), "child.hand.SA.age.soil.8"]))
names(C.hand.SA) <- c("C.hand.SA")
qplot(C.hand.SA$C.hand.SA, geom = "histogram", binwidth = 1)
# These hand surface areas are pretty close to what I'd expect! ~150
# # Exposure Factors Handbook, 2011 for 6-12 mo old
# # 6 to < 12 mo	0.51	5.3	0.02703	270.3	0.024	0.027	120	135 - ONE HAND mean 120, 95th pc 135 - my children are WAY smaller 


summary(c(anth.soilMass.wide$age.soil.child.7, anth.soilMass.wide$age.soil.child.8), na.rm = TRUE)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   32.02   35.31   36.82   36.82   38.30   41.72    6737 

############### To get the mother hand SA of the appropriate-age children, 
# # could use the hand SA distribution for WASHB mothers at endline
# anth.M.hand.SA.WASHB.end <- anth.soilMass.wide %>%
#   select(M.hand.SA.end) %>%
#   filter(M.hand.SA.end < 1000) ## Filter out the unreasonable hand SA values - there are not outliers but rather mistakes with height or weight measurement
# qplot(anth.M.hand.SA.WASHB.end$M.hand.SA.end, geom = "histogram", binwidth = 1)
# # In Yu (2009), XS slim women had a Measured Body SA of 15,220.56 cm2, assume ONE HAND is 2.26% = 343.9847 cm2 - my values are MUCH higher. Why?
# 
# # Dist DOES fit so use empirical
# anth.M.hand.SA.WASHB.end.dist <- fitdist(anth.M.hand.SA.WASHB.end$M.hand.SA.end, "norm", method = "mle")
# M.hand.SA.WASHB.end.mcstoc <- mcstoc(rnorm, type="V", mean = anth.M.hand.SA.WASHB.end.dist$estimate[[1]], sd = anth.M.hand.SA.WASHB.end.dist$estimate[[2]], rtrunc=TRUE, linf=0)
# 
# #### However, it is more parsimonious to use values from midline and endline combined # M.hand.SA.WASHB.mcstoc
# # This is not necessarily less accurate because we don't know if mothers gained or lost weight from endline to soil sample time)

## Determine the hand surface area distribution for WASHB mothers, pooling midline and endline data
M.hand.SA.WASHB <- anth.soilMass.wide %>%
  select(hh, M.hand.SA.mid, M.hand.SA.end) %>%
  gather(key = round, value = M.hand.SA, M.hand.SA.mid, M.hand.SA.end) %>%
  select(hh, M.hand.SA) %>%
  filter(M.hand.SA < 1000 & !is.na(M.hand.SA))

qplot(M.hand.SA.WASHB$M.hand.SA, geom = "histogram", binwidth = 1)

# Dist DOES fit so use dist
M.hand.SA.WASHB.dist <- fitdist(M.hand.SA.WASHB$M.hand.SA, "norm", method = "mle")

# Fitting of the distribution ' norm ' by maximum likelihood 
# Parameters:
#       estimate Std. Error
# mean 320.9020  0.5920519
# sd    26.2847  0.4186439
M.hand.SA.WASHB.mcstoc <- mcstoc(rnorm, type="V", mean = M.hand.SA.WASHB.dist$estimate[[1]], sd = M.hand.SA.WASHB.dist$estimate[[2]], rtrunc=TRUE, linf=0)
#M.hand.SA.WASHB.mcdata <- mcdata(sample(M.hand.SA.WASHB$M.hand.SA, size=ndvar(), replace=TRUE),type="V")

# In any case the values for mothers at endline and pooled midline and endline are very close
# # > M.hand.SA.WASHB.end.mcstoc
# # node    mode  nsv nsu nva variate min mean median max Nas type outm
# # 1    x numeric 1001   1   1       1 234  322    322 404   0    V each
# # > M.hand.SA.WASHB.mcstoc
# # node    mode  nsv nsu nva variate min mean median max Nas type outm
# # 1    x numeric 1001   1   1       1 232  321    321 388   0    V each


############## Est soil concentration on hands based on measured soil load for RO1 children and mothers and WASHB data for hand surface area ################
# Use the actual load divided by the est hand.SA based on # 

 ####### To est the soil adherence ON THE PORTION OF SKIN THAT IS CARRYING SOIL AND NOT ON THE OTHER PORTION, divide not by the total hand surface area but only on the "contaminated" surface area on the front of the hand 
# get the contaminated surface area on the front of the hand from Leckie_2000 page 489 - floor 23%, ceiling 35%

anth.soilMass.wide <- anth.soilMass.wide %>%
  mutate(frac.C.hand.contaminated = runif(nrow(anth.soilMass.wide), 0.23, 0.35)) %>% #, frac.M.hand.contaminated = runif(nrow(anth.soilMass.wide), 0.23, 0.35)) %>% 
  mutate(SA.C.hand.contaminated.7 = child.hand.SA.age.soil.7 * frac.C.hand.contaminated, SA.C.hand.contaminated.8 = child.hand.SA.age.soil.8 * frac.C.hand.contaminated)
  
anth.soilMass.wide[anth.soilMass.wide$hh %in% soilMass.wide$hh, c("child.hand.SA.age.soil.7", "soilOneHand.mg.child.7", "child.hand.SA.age.soil.8", "soilOneHand.mg.child.8")]

anth.soilMass.wide <- anth.soilMass.wide %>%
  mutate(soil.C.conc.7 = soilOneHand.mg.child.7 / SA.C.hand.contaminated.7, soil.C.conc.8 = soilOneHand.mg.child.8 / SA.C.hand.contaminated.8) #%>%
soil.C.conc.long <- data.frame(combine(anth.soilMass.wide$soil.C.conc.7, anth.soilMass.wide$soil.C.conc.8))
names(soil.C.conc.long) <- "soil.C.conc"


length(soil.C.conc.long) #68 observations
summary(soil.C.conc.long$soil.C.conc) 
# MassDEP uses 1.5 mg/cm2 as conservative; lowest mentioned in Finley is Duggan 1985 with 0.12 mg/cm2 - mine are even lower than this!
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.00000 0.08089 0.20240 0.29570 0.41790 3.53500       3

soil.C.empirical.36mo.conc <- mcstoc(rempiricalD, values = soil.C.conc.long$soil.C.conc, type="V", nsv = ndvar())
#soil.C.conc <- mcdata(sample(soil.C.conc.long$soil.C.conc, size=ndvar(), replace=TRUE),type="V")

##################################################################################################
 ### Figure out the soil adherence for younger children, who contact the soil more frequently ######

# If we want to base it on the weibull distribution:
# load frequency of soil contacts
voso.11.objclass.full <- read.csv("C:/Users/Tareq/Box Sync/VO R123/voso.11.objclass.full.csv")

# Touching FREQUENCY
vo.11.objclass.summary <- voso.11.objclass.full %>%
  filter(round != "so1", actobj.class == "Touch_soil") %>%
  group_by(age.group) %>%
  summarise(freq.med = median(freq, na.rm = TRUE), freq.mean = mean(freq, na.rm = TRUE), freq.sd = sd(freq, na.rm = TRUE))

# 
# # What distribution describes the touching freq for each age group? (This code is from VO_Longitudinal Analysis_170614)
# ## The code below shows that Weibull fits very well for all the data combined
# ## I will use non-parametric tests for analysis. However, I think it would be appropriate to use parametric tests since the fit to normal is so close. 
# Kwong.tbl <- voso.11.objclass.full %>% ungroup() %>% filter(actobj.class == "Touch_soil") %>% select(freq) # age.group %in% c("6-12 months", "12-24 months", "24-36 months", "36-48 months"))
# Kwong <- as.vector(unlist(Kwong.tbl))
# qplot(Kwong, geom = "histogram", binwidth = 10)
# plot(ecdf(Kwong),xlab="", ylab="",main="")
# 
# #Kolmogorov-Smirnov, Cramer-von-Mises and Anderson Darling, Chi-square
# norm  <-  fitdistr(na.omit(Kwong), densfun="normal")
# curve(pnorm(x, mean=norm$estimate[1], sd=norm$estimate[2]),from=0, to=200, add=TRUE,col="blue",lwd=2)
# ks.test(Kwong,rnorm(68, mean=norm$estimate[1], sd=norm$estimate[2]),alternative="two.sided") #D=0. ,p=0.
# cvm.test(Kwong,"pnorm",mean=norm$estimate[1], sd=norm$estimate[2]) #omega2=0. ,p=0.
# ad.test(Kwong,"pnorm",mean=norm$estimate[1], sd=norm$estimate[2]) #An=0. ,p=0. 
# chisq.test(Kwong,rnorm(68, mean=norm$estimate[1], sd=norm$estimate[2]))#x= 
# 
# lognorm  <-  fitdistr(na.omit(Kwong[Kwong != 0]), densfun="lognormal")
# curve(plnorm(x, meanlog=lognorm$estimate[1], sdlog=lognorm$estimate[2]),from=0, to=200, add=TRUE,col="red",lwd=2)
# ks.test(Kwong,rlnorm(68, meanlog=lognorm$estimate[1], sdlog=lognorm$estimate[2]),alternative="two.sided") #D=0. ,p=0.
# cvm.test(Kwong,"plnorm",meanlog=lognorm$estimate[1], sdlog=lognorm$estimate[2]) #omega2=0. ,p=0.
# ad.test(Kwong,"plnorm",meanlog=lognorm$estimate[1], sdlog=lognorm$estimate[2]) #An=0. ,p=0.
# chisq.test(Kwong,rlnorm(68, meanlog=lognorm$estimate[1], sdlog=lognorm$estimate [2])) #x= ,p=0.
# 
# weib  <-  fitdistr(na.omit(Kwong[Kwong != 0]), densfun=dweibull,start=list(scale=10,shape=1), method = "Nelder-Mead")
# curve(pweibull(x, scale=weib$estimate[1], shape=weib$estimate[2]),from=0, to=200, add=TRUE,col="green",lwd=2)
# ks.test(Kwong,rweibull(68, scale=weib$estimate[1], shape=weib$estimate[2]),alternative="two.sided") #D=0. ,p=0.
# cvm.test(Kwong,"pweibull",scale=weib$estimate[1], shape=weib$estimate[2]) #omega2=0. ,p=0. 
# ad.test(Kwong,"pweibull",scale=weib$estimate[1], shape=weib$estimate[2]) #An=0. , p=0.
# chisq.test(Kwong,rweibull(68, scale=weib$estimate[1], shape=weib$estimate[2]))#x= ,p=0.
# 
# legend("topleft",legend=c("normal","lognormal","Weibull"),col=c("blue","red","green"),lty=c(1,1,1),lwd=c(2,2,2))
# title(main="Soil touching by Bangladeshi children", xlab="Frequency (contacts/hr)",ylab="Percentile")
# 
# fitdistr(na.omit(Kwong.tbl[Kwong.tbl$freq != 0 & Kwong.tbl$age.group %in% "6-12 months", "freq"]), densfun=dweibull,start=list(scale=10,shape=1), method = "Nelder-Mead")
# 
# vo.11.loc.objclass <- read.csv("C:/Users/Tareq/Box Sync/VO R123/vo.11.loc.objclass.csv")
# names(vo.11.loc.objclass)
# vo.11.loc.objclass <- vo.11.loc.objclass %>%
#   select(-X)
# 
# touch.soil.weib <- matrix(NA, nrow = 2, ncol = 10)
# m = 1
# n = 1
# # I can't get the loop below to save var to .GlobalEnv
# for (actobj.type in c( "Touch_soil")){ #"Mouth_hands_d",
#   for (location in c("Inside", "Outside")){
#     for (ages in 4:7){
#       assign(paste("vo.11", ages, location, actobj.type, sep = "."), vo.11.loc.objclass %>%
#                filter(age.group == ages, loc == location, actobj.class == actobj.type)) #, envir = .GlobalEnv)
#       assign(paste("weib", ages, location, actobj.type, sep = "."), fitdistr(na.omit(get(paste("vo.11", ages, location, actobj.type, sep = "."))$freq), densfun=dweibull, start=list(scale=10, shape=1), method="Nelder-Mead"), envir = .GlobalEnv) # if get error of non-finite finite differences, try changing the start values or using method="Nelder-Mead"
#       weib <- fitdistr(na.omit(get(paste("vo.11", ages, location, actobj.type, sep = "."))$freq), densfun=dweibull, start=list(scale=10, shape=1), method="Nelder-Mead") # if get error of non-finite finite differences, try changing the start values or using method="Nelder-Mead"
#       #assign(paste("shape", age.group, location, actobj.type,sep = "."), as.numeric(weib$estimate[2])) #, envir = .GlobalEnv)
#       shape <- as.numeric(weib$estimate[2]) #, envir = .GlobalEnv)
#       #assign(paste("scale", age.group, location, actobj.type,sep = "."), as.numeric(weib$estimate[1])) #, envir = .GlobalEnv)
#       scale <- as.numeric(weib$estimate[1]) #, envir = .GlobalEnv)
#       touch.soil.weib[m, n] <- shape
#       touch.soil.weib[m, n+1] <- scale
#       n = n+2
#     }
#     n = 1
#     m = m+1
#   }
# }
# touch.soil.weib
# colnames(touch.soil.weib) <- c("4.shape", "4.scale", "5.shape", "5.scale", "6.shape", "6.scale", "7.shape", "7.scale", "extra", "extra")
# row.names(touch.soil.weib) <- c("touch.soil.in", "touch.soil.out") #"hands_d.in", "hands_d.out", 
# write.csv(touch.soil.weib, "C:/Users/Tareq/Box Sync/VO Soil/Touching.soil.weib.shape.scale.csv")

# If we just want to base it on the summaries: Frequency of touching soil: 
#   2 of 2 children <6 mo, med 35.6, mean 35.9 (inside contacts med 0.4 by child 71105, outside contacts mean 80.7 by child 66506);
# 19 of 19 children 6-12 mo; med 128.2, mean 130.7 contacts/hr; 
# 28 of 28 children 12-24 mo, med 24.55, mean 38.27 contacts/hr, 
# 22 of 23 24-36 mo, med 11.4, mean 20.6 contacts/hr, 6 of 6 children 36-48 mo, med 4, mean 3.63 contacts/hr

# But the number of CONSECUTIVE contacts is much lower and very similar for children. 
# Non-consecutive contacts will likely reduce the load of soil on children's hands and allow them to have high reloading when they do touch soil again. 
vo.9.loc <- read.csv("C:/Users/Tareq/Box Sync/VO R123/vo.9.analyze.csv") #1202536
vo.9.loc.base <- vo.9.loc %>%
  mutate(hh = as.factor(hh), round = as.factor(vo.num), age.mo = age.vo, age.group = age.vo.group)
vo.9.long.loc<- vo.9.loc.base %>% select(-X, -vo.num, -age.vo, -age.vo.group)

# https://stackoverflow.com/questions/19998836/create-a-sequential-number-counter-within-each-contiguous-run-of-equal-values
# sequence(rle(as.character(dataset$input))$lengths)
vo.9.touch.soil <- vo.9.long.loc %>%
  arrange(round, hh, clip, coder, hand, realtime) %>%
  filter(act == "Touch", obj %notin% c("Nothing", "NotInView"), actobj.class != "NA_NA")
vo.9.touch.soil$lengths <- sequence(rle(as.character(vo.9.touch.soil$actobj.class))$lengths)
vo.9.touch.soil %>%
  filter(actobj.class == "Touch_soil") %>%
  group_by(age.group) %>% #hand
  summarise(median = median(lengths), mean = mean(lengths), sd = sd(lengths), max = max(lengths))
# The median (max) number of consecutive soil contacts among children <6 mo was 2 (14), among children 6-12 mo 3 (37), 
# among children 12-24 mo 2 (59), among children 24-36 mo 1 (17), among children 36-48 mo 1 (3)

### BUT DO NOT ASSUME THAT THE ADHERENCE FACTOR SCALES LINEARALY - RODES 2001 finds that consecutive presses decreased the transfer factor by 3, requiring ~100 presses to reach equilibrium
# Assume that we took the hand rinse sample when the hand was at pseudo-equilibrium, which Rodes says is reached after ~100 contacts
# Tim Julian suggests there is a similar equilibrium with rotavirus after 10 min of contacts (and there is approx 33 seconds between each contact so 10 min = 20 contacts)
c = 0.65
total = 0
for(i in 0:36){
  total = total + sum(c/3^i)
}
total
# Since the drop-off is so quick, the load after 4 contacts is approximately equal to the load after 100 contacts

# For the median = 0.132, initial contact transfers 0.09, after 4 contacts load = 0.1311, 100 contacts = 0.135
# For the mean = 0.293, initial contact transfers 0.20 mg/cm2, after 4 contacts load = 0.291, 100 contacts = 0.300
# For the min = 0.0141, initial contact transfers 0.01 mg/cm2, after 4 contacts load = 0.0145, 100 contacts = 0.015
# For the max = 0.951, initial contact transfers 0.65 mg/cm2, after 4 contacts load = 0.946, 11 contacts = 0.975, 25 contacts = 0.975, 36 contacts = 0.975, 100 contacts = 0.975, 128 contacts = 0.975

# Based on data in Rodes_2001_Figure 5
dermal.frac.transferred <- data.frame(x = seq(1, 56, 5), y = c(64, 46, 72, 56, 40, 48, 28, 42, 18, 42, 9, 37)) #, rep(NA, 9)
model <- lm(log(dermal.frac.transferred$y) ~ dermal.frac.transferred$x)
plot(dermal.frac.transferred$x, dermal.frac.transferred$y, type = "p", lwd = 3)
lines(dermal.frac.transferred$x, exp(model$fit), col = "red")
x = 3
exp(4.2185 + (-0.02136)*x)
# After 1 contacts = 66.5
# After 2 contacts = 65.1
# After 3 contacts = 63.7
# After 100 contacts = 8.04
# After 130 contacts = 4.23
# After 150 contacts = 2.7
# After 200 contacts = 0.949

# with three additional datapoints that Rodes made up in order to get things to level off:
# dermal.frac.transferred <- data.frame(x = seq(1, 71, 5), y = c(64, 46, 72, 56, 40, 48, 28, 42, 18, 42, 9, 37, 24, 22, 20)) #, rep(NA, 9)
# result in the equation: exp(4.14 + (-0.01779)*x)
# After 100 contacts = 8.04
# After 130 contacts = 4.23
# After 150 contacts = 2.7
# After 200 contacts = 0.949


# The difference in concentration between the surface and the hand drives the transfer, so at some point an equilibrium will be reached. 
# In in rural Bangladesh (at least outdoors) the dislodgable residue is very, very high.

# So on the hand we have 0.08 mg/cm2, which we assume is the long-term equilibrium



soil.C.conc.f6 <- soil.C.empirical.36mo.conc * as.numeric(vo.11.objclass.summary[vo.11.objclass.summary$age.group == "3-6 months", "freq.med"]/vo.11.objclass.summary[vo.11.objclass.summary$age.group == "36-48 months", "freq.med"])
soil.C.conc.6_12 <- soil.C.empirical.36mo.conc * as.numeric(vo.11.objclass.summary[vo.11.objclass.summary$age.group == "6-12 months", "freq.med"]/vo.11.objclass.summary[vo.11.objclass.summary$age.group == "36-48 months", "freq.med"])
soil.C.conc.12_24 <- soil.C.empirical.36mo.conc * as.numeric(vo.11.objclass.summary[vo.11.objclass.summary$age.group == "12-24 months", "freq.med"]/vo.11.objclass.summary[vo.11.objclass.summary$age.group == "36-48 months", "freq.med"])
soil.C.conc.24_36 <- soil.C.empirical.36mo.conc * as.numeric(vo.11.objclass.summary[vo.11.objclass.summary$age.group == "24-36 months", "freq.med"]/vo.11.objclass.summary[vo.11.objclass.summary$age.group == "36-48 months", "freq.med"])
soil.C.conc.36_48 <- soil.C.empirical.36mo.conc

soil.M.conc <- soil.M.load.mcstoc / (M.hand.SA.WASHB.mcstoc * runif(ndvar(), 0.23, 0.35)) #M.hand.SA.WASHB.end.mcstoc 
#soil.M.conc <- soil.M.load.mcdata / M.hand.SA.WASHB.mcstoc #M.hand.SA.WASHB.end.mcstoc 

# Ozkaynak Soil adherence to hands mg/cm2 : mcstoc(rlnorm, type = "V", meanlog = log(0.11), sdlog = log(2.0)) * Oz reports logrithmic dist with geometric mean = 0.11 and geometric sd = 2.0
# Produces values ~ double my estimates...not sure what to do about this....
# My est
# > soil.C.conc_24_36
# node    mode  nsv nsu nva variate     min   mean median   max  Nas type outm
# 1    x numeric 1001 101   1       1 0.00987 0.0811 0.0421 0.867 4453   VU each

# Ozkaynak est
# > mcstoc(rlnorm, type = "V", meanlog = log(0.11), sdlog = log(2.0))
# node    mode  nsv nsu nva variate    min  mean median  max Nas type outm
# 1    x numeric 1001   1   1       1 0.0111 0.139  0.109 1.01   0    V each



############# Get the empirical hand SA data #########################
## About we found the hand SA of all children in the WASHB dataset: C.hand.SA.WASHB 

# By age group
C.hand.SA.WASHB.f6 <- C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo < 6, "C.hand.SA"] # This variable isn't actually used
C.hand.SA.WASHB.f6.mcstoc <- mcstoc(rempiricalD, values = C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo < 6, "C.hand.SA"], type = "V", nsv = ndvar())
#C.hand.SA.WASHB.f6.mcdata <- mcdata(sample(C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo < 6, "C.hand.SA"], size = ndvar(), replace = TRUE), type = "V")

C.hand.SA.WASHB.6_12 <- C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo >= 6 & C.hand.SA.WASHB$age.mo < 12, "C.hand.SA"] # This variable isn't actually used
C.hand.SA.WASHB.6_12.mcstoc <- mcstoc(rempiricalD, values = C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo >= 6 & C.hand.SA.WASHB$age.mo < 12, "C.hand.SA"], type = "V", nsv = ndvar())
#C.hand.SA.WASHB.6_12.mcdata <- mcdata(sample(C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo >= 6 & C.hand.SA.WASHB$age.mo < 12, "C.hand.SA"], size = ndvar(), replace = TRUE), type = "V")

C.hand.SA.WASHB.12_24 <- C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo >= 12 & C.hand.SA.WASHB$age.mo < 24, "C.hand.SA"] # This variable isn't actually used
C.hand.SA.WASHB.12_24.mcstoc <- mcstoc(rempiricalD, values = C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo >= 12 & C.hand.SA.WASHB$age.mo < 24, "C.hand.SA"], type = "V", nsv = ndvar())
#C.hand.SA.WASHB.12_24.mcdata <- mcdata(sample(C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo >= 12 & C.hand.SA.WASHB$age.mo < 24, "C.hand.SA"], size = ndvar(), replace = TRUE), type = "V")

C.hand.SA.WASHB.24_36 <- C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo >= 24, "C.hand.SA"] # This variable isn't actually used
C.hand.SA.WASHB.24_36.mcstoc <- mcstoc(rempiricalD, values = C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo >= 24, "C.hand.SA"], type = "V", nsv = ndvar())
#C.hand.SA.WASHB.24_36.mcdata <- mcdata(sample(C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo >= 24, "C.hand.SA"], size = ndvar(), replace = TRUE), type = "V")

qplot(C.hand.SA.WASHB.24_36, geom = "histogram", binwidth = 1)


# # By month-specific age
# # For the month-specific distributions, can I make normal/weibull for each month?
# C.hand.SA.WASHB.f6.dist <- fitdist(C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo < 6 & !is.na(C.hand.SA.WASHB$age.mo), "C.hand.SA"], "norm", method = "mle")
# C.hand.SA.WASHB.24_36.dist <- fitdist(C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo >= 24 & !is.na(C.hand.SA.WASHB$age.mo), "C.hand.SA"], "norm", method = "mle")
# 
# hand.SA.find.dist <- function(data, x, dist) {fitdist(data[data$age.mo == x, "C.hand.SA"], dist, method = "mle")}
# 
# for(i in c(6:23)){ # may need to avoid 23 months bc there are not children this age: for(i in c(6:13, 18:23))#max(HM.SM$age)
#   assign(paste("C.hand.SA.WASHB.", i, sep = ""), C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo == i, "C.hand.SA"])
#   assign(paste("C.hand.SA.WASHB.", i, ".mcstoc", sep = ""), mcstoc(rempiricalD, values = C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo == i, "C.hand.SA"], type = "V", nsv = ndvar()))
#   # assign(paste("C.hand.SA.WASHB.", i, ".mcdata", sep = ""), mcdata(sample(C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo == i, "C.hand.SA"], size = ndvar(), replace = TRUE), type = "V"))
# }

## We determined the dist for WASHB mother's hand SA in line 300 above: M.hand.SA.WASHB.mcstoc


############# Define some distributional tests ##################
## Code to determine the best-fitting distribution --> Weibull is the best for all 
bestFitDist <- function(data, title, xaxis){
  Kwong  <- data #"hand_m_nd_tot_freq"]
  max(Kwong)
  breaks <- seq(0,120,by=5)
  HMfreq <- table(cut(Kwong,breaks,right=FALSE))
  cum_HMfreq0 <- c(0,cumsum(HMfreq)/sum(HMfreq))
  plot(breaks, cum_HMfreq0, main = title, xlab = xaxis, ylab = "Percentile") #Percentile = cumulative density
  plot(ecdf(Kwong),add=TRUE)
  
  #Kolmogorov-Smirnov, Cramer-von-Mises and Anderson Darling, Chi-square
  lognorm  <-  fitdistr(na.omit(Kwong),densfun="log-normal")
  curve(plnorm(x, meanlog=lognorm$estimate[1], sdlog=lognorm$estimate[2]),from=0, to=120, add=TRUE,col="red")
  ks.test(Kwong, rlnorm(ndvar(), meanlog=lognorm$estimate[1], sdlog=lognorm$estimate[2]),alternative="two.sided") #D=0.11722, p=0.61203
  cvm.test(Kwong, "plnorm",meanlog=lognorm$estimate[1], sdlog=lognorm$estimate[2]) #omega2=0.043878,0.9168
  ad.test(Kwong, "plnorm",meanlog=lognorm$estimate[1], sdlog=lognorm$estimate[2]) #An=0.33139, p=9119
  chisq.test(Kwong, rlnorm(ndvar(), meanlog=lognorm$estimate[1], sdlog=lognorm$estimate[2]))
  
  weib  <-  fitdistr(na.omit(Kwong),densfun=dweibull,start=list(scale=20,shape=1))
  curve(pweibull(x, scale=weib$estimate[1], shape=weib$estimate[2]),from=0, to=120, add=TRUE,col="green")
  ks.test(Kwong,rweibull(ndvar(), scale=weib$estimate[1], shape=weib$estimate[2]),alternative="two.sided") #D=0.11722, p=0.61203
  cvm.test(Kwong,"pweibull",scale=weib$estimate[1], shape=weib$estimate[2]) #omega2=0.043878,0.9168
  ad.test(Kwong,"pweibull",scale=weib$estimate[1], shape=weib$estimate[2]) #An=0.33139, p=9119
  chisq.test(Kwong,rweibull(ndvar(), scale=weib$estimate[1], shape=weib$estimate[2]))
  
  norm  <-  fitdistr(na.omit(Kwong),densfun="normal")
  curve(pnorm(x, mean=norm$estimate[1], sd=norm$estimate[2]),from=0, to=120, add=TRUE,col="blue")
  ks.test(Kwong,rnorm(ndvar(), mean=norm$estimate[1], sd=norm$estimate[2]),alternative="two.sided") #D=0.11722, p=0.61203
  cvm.test(Kwong,"pnorm",mean=norm$estimate[1], sd=norm$estimate[2]) #omega2=0.043878,0.9168
  ad.test(Kwong,"pnorm",mean=norm$estimate[1], sd=norm$estimate[2]) #An=0.33139, p=9119
  chisq.test(Kwong,rnorm(ndvar(), mean=norm$estimate[1], sd=norm$estimate[2]))
  
  # beta  <-  fitdistr(na.omit(Kwong), densfun="beta", start = list(shape1 = 2, shape2 = 4))
  # curve(pbeta(x, mean=beta$estimate[1], sd=beta$estimate[2]),from=0, to=95, add=TRUE,col="blue")
  # ks.test(Kwong,rbeta(ndvar(), mean=beta$estimate[1], sd=beta$estimate[2]),alternative="two.sided") #D=0.11722, p=0.6953
  # cvm.test(Kwong,"pbeta",mean=beta$estimate[1], sd=beta$estimate[2]) #omega2=0.043878,0.9168
  # ad.test(Kwong,"pbeta",mean=beta$estimate[1], sd=beta$estimate[2]) #An=0.33139, p=9119
  # chisq.test(Kwong,rbeta(ndvar(), mean=beta$estimate[1], sd=beta$estimate[2]))
  # legend("topleft",legend=c("log-normal","Weibull","normal", "beta"), col=c("red","green","blue", "purple"),lty=c(1,1,1,1))
  
  legend("topleft",legend=c("log-normal","Weibull","normal"), col=c("red","green","blue"),lty=c(1,1,1))
}

find.lognorm.dist <- function(mcnode){
  sample <- sample(mcnode, ndvar())
  sample <- sample[!is.na(sample)]
  mcnode.dist <- fitdist(sample, "lnorm", method = "mle")
  print(mcnode.dist)
  plot(mcnode.dist)
}

find.norm.dist <- function(mcnode){
  sample <- sample(mcnode, ndvar())
  mcnode.dist <- fitdist(sample, "norm", method = "mle") # lower = c(0, 0),
  print(mcnode.dist)
  plot(mcnode.dist, xlim = c(50, 150))
}
#find.norm.dist(C.hand.SA.WASHB.6_12.mcstoc)

find.weib.dist <- function(mcnode){
  sample <- sample(mcnode, ndvar())
  sample <- sample[!is.na(sample)]
  mcnode.dist <- fitdist(sample, "weibull", method = "mle")
  print(mcnode.dist)
  plot(mcnode.dist)
}

find.unif.dist <- function(mcnode){
  sample <- sample(mcnode, ndvar())
  sample <- sample[!is.na(sample)]
  mcnode.dist <- fitdist(sample, "unif", method = "mle")
  print(mcnode.dist)
  plot(mcnode.dist)
}

# fitW <- fitdist(serving, "weibull")
# fitg <- fitdist(serving, "gamma")
# fitln <- fitdist(serving, "lnorm")
# summary(fitW)
# summary(fitg)
# cdfcomp(list(fitW, fitg, fitln), legendtext=c("Weibull", "gamma", "lognormal"))
# denscomp(list(fitW, fitg, fitln), legendtext=c("Weibull", "gamma", "lognormal"))
# qqcomp(list(fitW, fitg, fitln), legendtext=c("Weibull", "gamma", "lognormal"))
# ppcomp(list(fitW, fitg, fitln), legendtext=c("Weibull", "gamma", "lognormal"))
# gofstat(list(fitW, fitg, fitln), fitnames=c("Weibull", "gamma", "lognormal"))

find.lognorm.dist(amt.soil.ingested.mg.day.u6)
find.lognorm.dist(amt.soil.ingested.mg.day.6_12)
find.lognorm.dist(amt.soil.ingested.mg.day.12_24)
find.lognorm.dist(amt.soil.ingested.mg.day.24_36)

# These should likely al lbe normal not weibull but I get error 100 when trying to use the normal dist, perhaps bcause the sd is neg?
# Lognormal actually fits really well, so I will use that 

hist(C.hand.SA.WASHB.f6.mcstoc)
find.lognorm.dist(C.hand.SA.WASHB.f6.mcstoc)
find.lognorm.dist(C.hand.SA.WASHB.6_12.mcstoc)
find.lognorm.dist(C.hand.SA.WASHB.12_24.mcstoc)
find.lognorm.dist(C.hand.SA.WASHB.24_36.mcstoc)

######################## Add Hand-to-mouth contact data ################################
# Load HM frequencies and soil ingestion frequencies for each individual, not hh
## Structured Observation 
### Using data from all arms, justify by not sig diff (Kwong, 2016)

############################################
### Go back and figure out how in SO I calculated Soil_h_m_tot_freq --> what can I use this for?
###########################################3

# This file is pre-analysis in the R code for SO analysis that results in SO.HM.SM
#dem_obs_full <- read.csv("C:/Users/Tareq/Box Sync/VO Longitudinal/count_mat_food_dem_fromSO_submitted to IJERPH 160313.csv")
SO.allObs <- read.csv("C:/Users/Tareq/Dropbox/Fecal pathways/Structured obs/Analysis/SO_160310_submitted to IJERPH 160313/dem_obs_full_freq.csv")
names(SO.allObs)
# Create a table with the FREQ of hand-to-mouth contacts, starting with hand_m_tot_freq
SO.HM.SM <- SO.allObs %>%
  filter(location == 3) %>%
  mutate(round = "so1", Mouth_hands_child = 0, Mouth_hands_child_nd = 0, Mouth_hands_child_d = 0, Mouth_hands_mom_nd = 0, Mouth_hands_mom_d = 0) %>%
  select(participant_id,  round, age_SO_mo, age_SO_group, hand_m_tot, hand_m_tot_freq, Mouth_hands_child, Mouth_hands_child_nd, Mouth_hands_child_d, Mouth_hands_mom_nd, Mouth_hands_mom_d,  soil_m_tot_freq, soil_h_m_tot_freq)

names(SO.HM.SM) <- c("hh", "round", "age", "age.group", "HM_count", "HM_freq", "HM_c_freq", "HM_c_nd_freq", "HM_c_d_freq", "HM_m_nd_freq", "HM_m_d_freq",  "SM_freq", "SHM_freq")

# ############ SO Data if I want to use the age group ################3
# ## From SO of 149 kids
# SO.HM.base <- read.csv("C:/Users/Tareq/Dropbox/Fecal pathways/Structured obs/Analysis/SO_160310_submitted to IJERPH 160313/shape_scale_mean_med_SAVE.csv")

# Load HM frequencies and soil ingestion frequencies for each individual, not hh
## Video Observations 
hands_summary <- read.csv("C:/Users/Tareq/Box Sync/VO R123/hands_summary.csv")
hands_summary <- hands_summary %>%
  mutate(hh = as.factor(hh), round = vo.num, age = age.vo, age.group = age.vo.group) %>%
  replace(is.na(.), 0)

VO.HM.SM.handscontacts <- hands_summary %>%
  mutate(m_feeding_freq = 0, soil_m_tot_freq = 0, soil_h_m_tot_freq = 0) %>% # soil contacts are not recorded in hands_summary so soil_m_tot_freq and soil_h_m_tot_freq are NA
  select(hh, round, age, age.group, count_allhands, freq_allhands, freq_child_hands, freq_child_hands_nd, freq_child_hands_d, freq_other_hands_nd, freq_other_hands_d, soil_m_tot_freq, soil_h_m_tot_freq)
names(VO.HM.SM.handscontacts) <- c("hh", "round", "age", "age.group", "HM_count", "HM_freq", "HM_c_freq", "HM_c_nd_freq", "HM_c_d_freq", "HM_m_nd_freq", "HM_m_d_freq", "SM_freq", "SHM_freq")

# to get the soil contacts, use vo123.objclass.base
vo123.objclass.base <- read.csv("C:/Users/Tareq/Box Sync/VO R123/vo.11.objclass.csv")
VO.HM.SM.soilcontacts <- vo123.objclass.base %>%
  filter(actobj.class %in% c("Mouth_soil")) %>%
  select(hh, round, age.mo, age.group, actobj.class, freq) %>%
  spread(actobj.class, freq) %>%
  mutate(soil_h_m_tot_freq = 0) %>%
  replace(is.na(.), 0) %>%
  select(hh, round, Mouth_soil, soil_h_m_tot_freq)
names(VO.HM.SM.soilcontacts) <- c("hh", "round", "SM_freq", "SHM_freq")
VO.HM.SM.soilcontacts$hh <- as.factor(VO.HM.SM.soilcontacts$hh)

VO.HM.SM <- left_join(VO.HM.SM.handscontacts[,c("hh", "round", "age", "age.group", "HM_count", "HM_freq", "HM_c_freq", "HM_c_nd_freq", "HM_c_d_freq", "HM_m_nd_freq", "HM_m_d_freq")], VO.HM.SM.soilcontacts[,c("hh", "round", "SM_freq", "SHM_freq")], by = c("hh", "round"))

HM.SM.prefeedingevents <- rbind(SO.HM.SM, VO.HM.SM)

scatter.smooth(HM.SM[,c("age", "HM_freq")])
scatter.smooth(HM.SM[,c("age", "SM_freq")])


# For all the kids that only have total HM values from the SO, need to est child, mom_d, mom_nd
# Est using the VO results for children 6-12 mo
# For the other children, use the actual child, mom_d, and mom_nd recorded in the VO


### Must change to reflect that most mom HM contacts do not have recontam so really only want number of indep feeding events
feeding.mom.freq <- read.csv("C:/Users/Tareq/Box Sync/VO R123/feeding.mom.freq.csv")
feeding.mom.freq <- feeding.mom.freq %>%
  mutate(hh = as.factor(hh)) %>%
  select(-X, -X.1) 
names(feeding.mom.freq) <- c("round", "hh", "age", "hh_feeding_events", "awake.h", "m_feeding_freq")

HM.SM <- left_join(HM.SM.prefeedingevents, feeding.mom.freq[,c("hh", "round", "m_feeding_freq")], by = c("hh", "round")) %>%
  replace(is.na(.), 0) %>%
  arrange(age, round, hh)

HM.proportions <- left_join(hands_summary, feeding.mom.freq, by = c("round", "hh", "age")) %>%
  select(round, hh, age, age.group, pf_child_hands_nd, pf_child_hands_d, pf_other_hands_nd, pf_other_hands_d)

scatter.smooth(HM.SM$age, HM.SM$m_feeding_freq)
# There is no trend by age, so use feeding freq from any age for a child of any age

########### FOR THE STRUCTURED OBSERVATION ONLY, ESTIMATE THE BREAKDOWN OF CONTACTS BETWEEN MOTHER AND CHILD, DIETARY and NON #########
## For the video observations, I can use the raw data

# Make an array of vectors, where each vector is the child_nd, child_d, mom_nd, feeding.freq for each child in the age group
# Then create a random sample of the indicies
# Select the random index of the entire row and get the child_nd, child_d, mom_nd, or feeding.freq col and save as mcdata

# age.group = 4 is 12-24 mo, age.group = 5 is 12-24 mo, age.group = 6 is 24-36 mo
## Before 26 May 2017 I was calculatin ghtis WRONG - I had the proportion but did not multiply by the actual FREQUENCY of mouthing....

# There are so few kids under 6 months that the proportion shouldn't be generalized, use the kids <6 and 6_12 mo to determine the proportions
HM.proportions.f6 <- HM.proportions %>% filter(age.group %in% c(2, 3,4)) %>% select(pf_child_hands_nd, pf_child_hands_d, pf_other_hands_nd, pf_other_hands_d) /100
HM.proportions.index.f6 <- sample(c(1:nrow(HM.proportions.f6)), size = nrow(HM.SM[HM.SM$round == "so1" & HM.SM$age.group %in% c(2, 3),]), replace = TRUE)
HM.SM[HM.SM$round == "so1" & HM.SM$age.group %in% c(2, 3), "HM_c_nd_freq"] <- HM.proportions.f6[HM.proportions.index.f6, "pf_child_hands_nd"] * HM.SM[HM.SM$round == "so1" & HM.SM$age.group %in% c(2, 3),"HM_freq"]
HM.SM[HM.SM$round == "so1" & HM.SM$age.group %in% c(2, 3), "HM_c_d_freq"] <- HM.proportions.f6[HM.proportions.index.f6, "pf_child_hands_d"] * HM.SM[HM.SM$round == "so1" & HM.SM$age.group %in% c(2, 3),"HM_freq"]
HM.SM[HM.SM$round == "so1" & HM.SM$age.group %in% c(2, 3), "HM_m_nd_freq"] <- HM.proportions.f6[HM.proportions.index.f6, "pf_other_hands_nd"] * HM.SM[HM.SM$round == "so1" & HM.SM$age.group %in% c(2, 3),"HM_freq"]
HM.SM[HM.SM$round == "so1" & HM.SM$age.group %in% c(2, 3), "HM_m_d_freq"] <- HM.proportions.f6[HM.proportions.index.f6, "pf_other_hands_d"] * HM.SM[HM.SM$round == "so1" & HM.SM$age.group %in% c(2, 3),"HM_freq"]
HM.SM[HM.SM$round == "so1" & HM.SM$age.group %in% c(2, 3), "HM_c_freq"] <- rowSums(HM.SM[HM.SM$round == "so1" & HM.SM$age.group %in% c(2, 3), c("HM_c_nd_freq", "HM_c_d_freq")]) 

HM.proportions.6_12 <- HM.proportions %>% filter(age.group == 4) %>% select(pf_child_hands_nd, pf_child_hands_d, pf_other_hands_nd, pf_other_hands_d) /100
HM.proportions.index.6_12 <- sample(c(1:nrow(HM.proportions.6_12)), size = nrow(HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 4,]), replace = TRUE)
HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 4, "HM_c_nd_freq"] <- HM.proportions.6_12[HM.proportions.index.6_12, "pf_child_hands_nd"] * HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 4,"HM_freq"]
HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 4, "HM_c_d_freq"] <- HM.proportions.6_12[HM.proportions.index.6_12, "pf_child_hands_d"] * HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 4,"HM_freq"]
HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 4, "HM_m_nd_freq"] <- HM.proportions.6_12[HM.proportions.index.6_12, "pf_other_hands_nd"] * HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 4,"HM_freq"]
HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 4, "HM_m_d_freq"] <- HM.proportions.6_12[HM.proportions.index.6_12, "pf_other_hands_d"] * HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 4,"HM_freq"]
HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 4, "HM_c_freq"] <- rowSums(HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 4, c("HM_c_nd_freq", "HM_c_d_freq")]) 

HM.proportions.12_24 <- HM.proportions %>% filter(age.group == 5) %>% select(pf_child_hands_nd, pf_child_hands_d, pf_other_hands_nd, pf_other_hands_d) /100
HM.proportions.index.12_24 <- sample(c(1:nrow(HM.proportions.12_24)), size = nrow(HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 5,]), replace = TRUE)
HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 5, "HM_c_nd_freq"] <- HM.proportions.12_24[HM.proportions.index.12_24, "pf_child_hands_nd"] * HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 5,"HM_freq"]
HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 5, "HM_c_d_freq"] <- HM.proportions.12_24[HM.proportions.index.12_24, "pf_child_hands_d"] * HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 5,"HM_freq"]
HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 5, "HM_m_nd_freq"] <- HM.proportions.12_24[HM.proportions.index.12_24, "pf_other_hands_nd"] * HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 5,"HM_freq"]
HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 5, "HM_m_d_freq"] <- HM.proportions.12_24[HM.proportions.index.12_24, "pf_other_hands_d"] * HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 5,"HM_freq"]
HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 5, "HM_c_freq"] <- rowSums(HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 5, c("HM_c_nd_freq", "HM_c_d_freq")]) 

HM.proportions.24_36 <- HM.proportions %>% filter(age.group == 6) %>% select(pf_child_hands_nd, pf_child_hands_d, pf_other_hands_nd, pf_other_hands_d) /100
HM.proportions.index.24_36 <- sample(c(1:nrow(HM.proportions.24_36)), size = nrow(HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 6,]), replace = TRUE)
HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 6, "HM_c_nd_freq"] <- HM.proportions.24_36[HM.proportions.index.24_36, "pf_child_hands_nd"] * HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 6,"HM_freq"]
HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 6, "HM_c_d_freq"] <- HM.proportions.24_36[HM.proportions.index.24_36, "pf_child_hands_d"] * HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 6,"HM_freq"]
HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 6, "HM_m_nd_freq"] <- HM.proportions.24_36[HM.proportions.index.24_36, "pf_other_hands_nd"] * HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 6,"HM_freq"]
HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 6, "HM_m_d_freq"] <- HM.proportions.24_36[HM.proportions.index.24_36, "pf_other_hands_d"] * HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 6,"HM_freq"]
HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 6, "HM_c_freq"] <- rowSums(HM.SM[HM.SM$round == "so1" & HM.SM$age.group == 6, c("HM_c_nd_freq", "HM_c_d_freq")]) 

##### Create mcnodes made using mcstoc(rempirical) from all the hand contact data that is summarized in HM.SM ######
## In the code that follows this block, Will be replaced by DISTRIBUTIONS for this data, if distributions exist
HM.child.f6.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(2,3), "HM_c_freq"], type = "V", nsv = ndvar())
HM.child.nd.f6.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(2,3), "HM_c_nd_freq"], type = "V", nsv = ndvar())
HM.child.d.f6.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(2,3), "HM_c_d_freq"], type = "V", nsv = ndvar())
HM.mom.nd.f6.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(2,3), "HM_m_nd_freq"], type = "V", nsv = ndvar())
HM.mom.d.events.f6.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(2,3), "m_feeding_freq"], type = "V", nsv = ndvar())

HM.child.6_12.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(4), "HM_c_freq"], type = "V", nsv = ndvar())
HM.child.nd.6_12.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(4), "HM_c_nd_freq"], type = "V", nsv = ndvar())
HM.child.d.6_12.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(4), "HM_c_d_freq"], type = "V", nsv = ndvar())
HM.mom.nd.6_12.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(4), "HM_m_nd_freq"], type = "V", nsv = ndvar())
HM.mom.d.events.6_12.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(4), "m_feeding_freq"], type = "V", nsv = ndvar())

HM.child.12_24.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(5), "HM_c_freq"], type = "V", nsv = ndvar())
HM.child.nd.12_24.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(5), "HM_c_nd_freq"], type = "V", nsv = ndvar())
HM.child.d.12_24.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(5), "HM_c_d_freq"], type = "V", nsv = ndvar())
HM.mom.nd.12_24.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(5), "HM_m_nd_freq"], type = "V", nsv = ndvar())
HM.mom.d.events.12_24.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(5), "m_feeding_freq"], type = "V", nsv = ndvar())

HM.child.24_36.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(6), "HM_c_freq"], type = "V", nsv = ndvar())
HM.child.nd.24_36.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(6), "HM_c_nd_freq"], type = "V", nsv = ndvar())
HM.child.d.24_36.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(6), "HM_c_d_freq"], type = "V", nsv = ndvar())
HM.mom.nd.24_36.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(6), "HM_m_nd_freq"], type = "V", nsv = ndvar())
HM.mom.d.events.24_36.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(6), "m_feeding_freq"], type = "V", nsv = ndvar())

HM.child.36_48.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(7), "HM_c_freq"], type = "V", nsv = ndvar())
HM.child.nd.36_48.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(7), "HM_c_nd_freq"], type = "V", nsv = ndvar())
HM.child.d.36_48.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(7), "HM_c_d_freq"], type = "V", nsv = ndvar())
HM.mom.nd.36_48.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(7), "HM_m_nd_freq"], type = "V", nsv = ndvar())
HM.mom.d.events.36_48.mcstoc <- mcstoc(rempiricalD, values = HM.SM[HM.SM$age.group %in% c(7), "m_feeding_freq"], type = "V", nsv = ndvar())

# Does HM.child follow a weibull dist? 
# for age.grouup = 4, the weibull fits all right, but for the others, there are too few datapoint for a good fit
# --> for now use the empirical distributions, use the weibull if I get more data later
descdist(unmc(HM.child.f6.mcstoc), discrete = FALSE)
bestFitDist(data = unmc(HM.child.f6.mcstoc)+ 10e-4, "Hand-to-mouth freq/hr", "Hand-to-mouth freq/hr")
descdist(unmc(HM.child.6_12.mcstoc), discrete = FALSE)
bestFitDist(data = unmc(HM.child.6_12.mcstoc)+ 10e-4, "Hand-to-mouth freq/hr", "Hand-to-mouth freq/hr")
descdist(unmc(HM.child.12_24.mcstoc), discrete = FALSE)
bestFitDist(data = unmc(HM.child.12_24.mcstoc)+ 10e-4, "Hand-to-mouth freq/hr", "Hand-to-mouth freq/hr")
descdist(unmc(HM.child.24_36.mcstoc), discrete = FALSE)
bestFitDist(data = unmc(HM.child.24_36.mcstoc)+ 10e-4, "Hand-to-mouth freq/hr", "Hand-to-mouth freq/hr")
# The distributional fit for HM.child.24_36 is pretty poor 
descdist(unmc(HM.child.36_48.mcstoc), discrete = FALSE)
bestFitDist(data = unmc(HM.child.36_48.mcstoc)+ 10e-4, "Hand-to-mouth freq/hr", "Hand-to-mouth freq/hr")
# The distributional fit for HM.child.36_48 is pretty poor 

########### Create DISTRIBUTIONS to represent mouthing frequency PER HOUR - will later conver to per day ############## 
### Use this distributions to REPLACE the mcnodes made above using mcstoc(empirical D) IF THE DISTRIBUTOIN EXISTS (will simply error and not overwrite if the distribution does not exist)
## For all HM_child
HM.child.f6.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(2,3), "HM_c_freq"], "weibull", method = "mle")
HM.child.f6.mcstoc <- mcstoc(rweibull, type="V", shape = HM.child.f6.dist$estimate[[1]], scale = HM.child.f6.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.child.nd.f6.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(2,3), "HM_c_nd_freq"], "weibull", method = "mle")
HM.child.nd.f6.mcstoc <- mcstoc(rweibull, type="V", shape = HM.child.nd.f6.dist$estimate[[1]], scale = HM.child.nd.f6.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.child.d.f6.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(2,3), "HM_c_d_freq"], "weibull", method = "mle")
HM.child.d.f6.mcstoc <- mcstoc(rweibull, type="V", shape = HM.child.d.f6.dist$estimate[[1]], scale = HM.child.d.f6.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.mom.nd.f6.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(2,3), "HM_m_nd_freq"], "weibull", method = "mle")
HM.mom.nd.f6.mcstoc <- mcstoc(rweibull, type="V", shape = HM.mom.nd.f6.dist$estimate[[1]], scale = HM.mom.nd.f6.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.mom.d.events.f6.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(2,3), "m_feeding_freq"], "weibull", method = "mle")
HM.mom.d.events.f6.mcstoc <- mcstoc(rweibull, type="V", shape = HM.mom.d.events.f6.dist$estimate[[1]], scale = HM.mom.d.events.f6.dist$estimate[[2]], rtrunc=TRUE, linf=0)

HM.child.6_12.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(4), "HM_c_freq"], "weibull", method = "mle")
HM.child.6_12.mcstoc <- mcstoc(rweibull, type="V", shape = HM.child.6_12.dist$estimate[[1]], scale = HM.child.6_12.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.child.nd.6_12.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(4), "HM_c_nd_freq"], "weibull", method = "mle")
HM.child.nd.6_12.mcstoc <- mcstoc(rweibull, type="V", shape = HM.child.nd.6_12.dist$estimate[[1]], scale = HM.child.nd.6_12.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.child.d.6_12.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(4), "HM_c_d_freq"], "weibull", method = "mle")
HM.child.d.6_12.mcstoc <- mcstoc(rweibull, type="V", shape = HM.child.d.6_12.dist$estimate[[1]], scale = HM.child.d.6_12.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.mom.nd.6_12.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(4), "HM_m_nd_freq"], "weibull", method = "mle")
HM.mom.nd.6_12.mcstoc <- mcstoc(rweibull, type="V", shape = HM.mom.nd.6_12.dist$estimate[[1]], scale = HM.mom.nd.6_12.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.mom.d.events.6_12.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(4), "m_feeding_freq"], "weibull", method = "mle")
HM.mom.d.events.6_12.mcstoc <- mcstoc(rweibull, type="V", shape = HM.mom.d.events.6_12.dist$estimate[[1]], scale = HM.mom.d.events.6_12.dist$estimate[[2]], rtrunc=TRUE, linf=0)

HM.child.12_24.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(5), "HM_c_freq"], "weibull", method = "mle")
HM.child.12_24.mcstoc <- mcstoc(rweibull, type="V", shape = HM.child.12_24.dist$estimate[[1]], scale = HM.child.12_24.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.child.nd.12_24.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(5), "HM_c_nd_freq"], "weibull", method = "mle")
HM.child.nd.12_24.mcstoc <- mcstoc(rweibull, type="V", shape = HM.child.d.12_24.dist$estimate[[1]], scale = HM.child.nd.12_24.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.child.d.12_24.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(5), "HM_c_d_freq"], "weibull", method = "mle")
HM.child.d.12_24.mcstoc <- mcstoc(rweibull, type="V", shape = HM.child.nd.12_24.dist$estimate[[1]], scale = HM.child.d.12_24.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.mom.nd.12_24.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(5), "HM_m_nd_freq"], "weibull", method = "mle")
HM.mom.nd.12_24.mcstoc <- mcstoc(rweibull, type="V", shape = HM.mom.nd.12_24.dist$estimate[[1]], scale = HM.mom.nd.12_24.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.mom.d.events.12_24.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(5), "m_feeding_freq"], "weibull", method = "mle")
HM.mom.d.events.12_24.mcstoc <- mcstoc(rweibull, type="V", shape = HM.mom.d.events.12_24.dist$estimate[[1]], scale = HM.mom.d.events.12_24.dist$estimate[[2]], rtrunc=TRUE, linf=0)

HM.child.24_36.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(6), "HM_c_freq"], "weibull", method = "mle")
HM.child.24_36.mcstoc <- mcstoc(rweibull, type="V", shape = HM.child.24_36.dist$estimate[[1]], scale = HM.child.24_36.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.child.nd.24_36.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(6), "HM_c_nd_freq"], "weibull", method = "mle")
HM.child.nd.24_36.mcstoc <- mcstoc(rweibull, type="V", shape = HM.child.d.24_36.dist$estimate[[1]], scale = HM.child.nd.24_36.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.child.d.24_36.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(6), "HM_c_d_freq"], "weibull", method = "mle")
HM.child.d.24_36.mcstoc <- mcstoc(rweibull, type="V", shape = HM.child.nd.24_36.dist$estimate[[1]], scale = HM.child.d.24_36.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.mom.nd.24_36.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(6), "HM_m_nd_freq"], "weibull", method = "mle")
HM.mom.nd.24_36.mcstoc <- mcstoc(rweibull, type="V", shape = HM.mom.nd.24_36.dist$estimate[[1]], scale = HM.mom.nd.24_36.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.mom.d.events.24_36.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(6), "m_feeding_freq"], "weibull", method = "mle")
HM.mom.d.events.24_36.mcstoc <- mcstoc(rweibull, type="V", shape = HM.mom.d.events.24_36.dist$estimate[[1]], scale = HM.mom.d.events.24_36.dist$estimate[[2]], rtrunc=TRUE, linf=0)

HM.child.36_48.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(7), "HM_c_freq"], "weibull", method = "mle")
HM.child.36_48.mcstoc <- mcstoc(rweibull, type="V", shape = HM.child.36_48.dist$estimate[[1]], scale = HM.child.36_48.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.child.nd.36_48.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(7), "HM_c_nd_freq"], "weibull", method = "mle")
HM.child.nd.36_48.mcstoc <- mcstoc(rweibull, type="V", shape = HM.child.d.36_48.dist$estimate[[1]], scale = HM.child.nd.36_48.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.child.d.36_48.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(7), "HM_c_d_freq"], "weibull", method = "mle")
HM.child.d.36_48.mcstoc <- mcstoc(rweibull, type="V", shape = HM.child.nd.36_48.dist$estimate[[1]], scale = HM.child.d.36_48.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.mom.nd.36_48.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(7), "HM_m_nd_freq"], "weibull", method = "mle")
HM.mom.nd.36_48.mcstoc <- mcstoc(rweibull, type="V", shape = HM.mom.nd.36_48.dist$estimate[[1]], scale = HM.mom.nd.36_48.dist$estimate[[2]], rtrunc=TRUE, linf=0)
HM.mom.d.events.36_48.dist <- fitdist(HM.SM[HM.SM$age.group %in% c(7), "m_feeding_freq"], "weibull", method = "mle")
HM.mom.d.events.36_48.mcstoc <- mcstoc(rweibull, type="V", shape = HM.mom.d.events.36_48.dist$estimate[[1]], scale = HM.mom.d.events.36_48.dist$estimate[[2]], rtrunc=TRUE, linf=0)



#################### Combine load_child/load_mom with HM_child/HM_mom with mouth SA dist, removal efficiency dist

#HF	Child own hand fracion mouthed/event	-	day	beta	3.7	25				Zartarian 2005 as presented in Ozkaynak 2011	Zartarian 2005
HF.child.mcstoc <- mcstoc(rbeta, type="V", shape1 = 3.7, shape2 = 25, rtrunc=TRUE, linf=0)
##########################################################
# I can do better than this!! Work to replace this value, which is based only on 20 children. 


























# Children will mouth a smaller portion of mom's hand than of own hand 
# Children are "mouth-limited" in how much surface area they can put in their mouth. 
# Children will be able to put more surface area of their own hands in their mouth 
# than the surface area of their mother's hands because they have more surface area per unit volume (smaller finger = can put more finger = more SA in mouth)

## Use age groups rather than age by month. There aren't very many observations for children of specific months, 
## and some kids who weren't observed may have been a bit bigger or smaller so seems okay to include the entire age group
## Could change this so XX is by age group, not by month....This is better because below I change it into age groups anyway...
C.hand.SA.WASHB.inMouth.f6.mcstoc <- C.hand.SA.WASHB.f6.mcstoc * HF.child.mcstoc # mean 12.3, med 11.4, min 0.373, max 50.5
C.hand.SA.WASHB.inMouth.6_12.mcstoc <- C.hand.SA.WASHB.6_12.mcstoc * HF.child.mcstoc # mean 13.4, med 12.4, min 0.236, max 58.5
C.hand.SA.WASHB.inMouth.12_24.mcstoc <- C.hand.SA.WASHB.12_24.mcstoc * HF.child.mcstoc # mean 16.2, med 15.0, min 0.278, max 65.6
C.hand.SA.WASHB.inMouth.24_36.mcstoc <- C.hand.SA.WASHB.24_36.mcstoc * HF.child.mcstoc # mean 17.5, med 16.2, min 0.304, max 67.7

# HF.ofmom.xx <- C.hand.SA.WASHB.inMouth.xx.mcstoc / median(M.hand.SA.WASHB.mcstoc)
HF.ofmom.f6 <- C.hand.SA.WASHB.inMouth.f6.mcstoc / median(M.hand.SA.WASHB.mcstoc) # use for age.group = 3 = <6 mo
HF.ofmom.6_12 <- C.hand.SA.WASHB.inMouth.6_12.mcstoc / median(M.hand.SA.WASHB.mcstoc) # use for age.group = 4 = 6-12 mo
HF.ofmom.12_24 <- C.hand.SA.WASHB.inMouth.12_24.mcstoc / median(M.hand.SA.WASHB.mcstoc) # use for age.group = 5 = 12_24 mo
HF.ofmom.24_36 <- C.hand.SA.WASHB.inMouth.24_36.mcstoc / median(M.hand.SA.WASHB.mcstoc) # use for age.group = 6 = 12_36 mo

## Don't use months because there isn't enough data for some of the months


####################### SEE Saliva extraction efficiency  = HMRE	Hand mouthing removal = transfer efficiency #######################

# For the many reasons listed in my soil paper, I will NOT use ths Ozkaynak 2011 HMRE per	day	beta	2	8			
# we estimated a triangular distribution with a lower bound of 24%, mode of 75%, and upper bound of 100% with the beta distribution using maximum likelihood estimation. The resulting distribution, beta(5.1, 2.6), has a 5th percentile of 38.1%, median of 68.2%, and a 95th percentile of 90.1%. 
fitdist(rtriangle(100000, 0.24, 0.75, (0.24+0.75)/2), "beta", method = "mle") # For a triangular dist, the values must be 0-1
hist(rbeta(10000, shape1 = 11.02, shape2 = 11.24))
SEE.mcstoc <- mcstoc(rbeta, type="V", shape1 = 11.02, shape2 = 11.24, rtrunc=TRUE, linf=0)




#####################################################################################
#####################################################################################

################## Soil ingestion from contacts with objects #########################

#####################################################################################
######################################################################################

################### Concentration of soil on objects ###################

# Ozkaynak uses Concentration of soil on obj = Load of dust/soil on surface * Object:Surface loading ratio
# For load of dust/soil on surface, Oz cites Adgate and uses lognormal(7.8, 2.9) g/m2
# lognormal(7.8, 2.9)  = meanlog 7.8 g/m2
# For Object:surface loading ratio, Oz cites Gurunuthan (2-70%) and decides on a uniform(0,20)
# Also says assumes that plastic has lower loading but higher transfer efficiency (though doesn't explicitly break out plastics and treat it differently)

floor.soil.concentration <- mcstoc(rlnorm, meanlog = log(7.8), sdlog = log(2.9), type = "V") * 1000/10000 # Conversion for mg/cm2
object.floor.loading.ratio <- mcstoc(runif, min = 0.001, max = 0.20, type = "V")
object.soil.concentration.mcstoc <- floor.soil.concentration * object.floor.loading.ratio
hist(object.soil.concentration.mcstoc)
hist(log10(object.soil.concentration.mcstoc))


# In Bangladesh, the load of dust/soil on surface is infinite, so we can't use this method
# Instead assume that loading on object is the SAME (upper bound) as loading on floor
# Floor loading in Giza, Egypt (where I recall SUPER bad air quality)(Khoder, 2010 n = 6 homes): indoor median = 1.45 g /m2 = 0.145 mg/cm2 , outdoor median = 8.22 g/m2
# This study is not clear on how these values were obtained....
# Floor loading in Delhi based on Kumar 2009 lead on floors and Banerjee 2003 lead in dust (n = 8 homes):  0.001576 mg/cm2

fitdist(rtriangle(100000, 0.00, 0.145, 0.001576), "beta", method = "mle") # For a triangular dist, the values must be 0-1
hist(rtriangle(100000, 0.00, 0.145, 0.001576))
hist(rbeta(10000, 1.47, 28.76))
summary(rbeta(10000, 1.47, 28.76))

object.soil.concentration.mcstoc <- mcstoc(rbeta, shape1 = 1.47, shape2 = 28.76, type = "V")


## The values used by Oz provide a greater range over which to test the importance of this factor, so for now I will use those


####################### contacts with objects #####################

# To examine ingestion of soil and dust from mouthing objects, I should consider on non-cloth objects,  
## Must parse soil ingestion into H-M, O-M, and direct, because for the fecal intake model I will use direct ingestion only
# FIB ingestion due to HM and OM contacts will be modeled using the EC on hands and obj. 
# This soil ingestion is for estimating ingestion of heavy metals, pesticides, etc. 

# fomites <- c("PlantMaterial","Metal","Paper","Cloth","OtherObject","Plastic","Wood/Bricks","Utensil")

# Consider only OM contacts with non-cloth and non-utensil objects, as these aren't expected to be contaminated with soil
# though the could (and are) contaminated with EC. 

## Structured observations
SO.OM <- SO.allObs %>%
  filter(location == 3) %>%
  mutate(round = "so.1") %>%
  select(participant_id,  round, age_SO_mo, age_SO_group, toy_m_tot_freq) %>% # Omit cloth_m_tot_freq (obj_m_tot_freq = cloth_m_tot_freq + toy_m_tot_freq)
  replace(is.na(.), 0) # Slightly risky because this will replace all na with 0, but I've checked many times and the only na values are in toy_m_tot_freq (two kids age 2 and 4 months), which are the ones I want to replace with 0
names(SO.OM) <- c("hh", "round", "age", "age.group", "OM")


## Video observations
vo123.obj.base <- read.csv("C:/Users/Tareq/Box Sync/VO R123/vo.11.obj.csv")
vo123.obj.base <- vo123.obj.base[!is.na(vo123.obj.base$round),]
VO.OM.sepfreqs <- vo123.obj.base %>%
  filter(actobj %in% c("Mouth_PlantMaterial", "Mouth_Metal", "Mouth_Paper", "Mouth_Plastic", "Mouth_OtherObject", "Mouth_Wood/Bricks")) %>% # Omit Mouth_Cloth & Mouth_Utensil
  select(hh, round, age.mo, age.group, actobj, freq) %>%
  spread(actobj, freq) %>%
  mutate(OM.freq = rowSums(.[5:10], na.rm = TRUE))

VO.OM <- VO.OM.sepfreqs %>%
  select(hh, round, age.mo, age.group, OM.freq) %>%
  arrange(age.mo, round, hh)

# There are only 50 (hh, vo.num) combinations in VO.OM.summary but there were 57 hh observations in vo.1 and vo.2 
# so in 7 hh there was NO mouthing of an object that was not cloth or utensils. 
# Seems unlikely but maybe these were older kids. 
# Find which are missing using setdiff() and add them to VO.OM.summary
hh.zero.OM <- setdiff(unique(vo123.obj.base[,c("hh","round", "age.mo", "age.group")]), unique(VO.OM[,c("hh","round", "age.mo", "age.group")]))
hh.zero.OM$OM.freq <- 0

VO.OM <- rbind(VO.OM, hh.zero.OM)

names(VO.OM) <- c("hh", "round", "age", "age.group", "OM")

OM <- rbind(SO.OM, VO.OM)
OM[OM$OM < 0.001, "OM"] <- 0.001
#### will need to remove this line when fix 2302 associated with NAs
OM <- OM[1:nrow(OM)-1,]

hist(OM[OM$age > 6 & OM$age <= 12, "OM"])
hist(OM[OM$age > 12 & OM$age <= 24, "OM"])

### Which distribution is suitable? Weibull, based on Xue meta-analysis, also the best fit for our dist based on code in 1120
# Doest seem to fit pretty well for the data here
bestFitDist(data = sample((OM[OM$age < 6, "OM"] + 10e-4), ndvar(), replace = TRUE), "Hand to mouth", "Hand to Mouth")
### all of the distributions fit terribly, normal sort of fits the best
OM.f6.dist <- fitdist(OM[OM$age < 6, "OM"], "weibull", method = "mle")
OM.f6.mcstoc <- mcstoc(rweibull, type="V", shape = OM.f6.dist$estimate[[1]], scale = OM.f6.dist$estimate[[2]], rtrunc=TRUE, linf=0)

bestFitDist(data = sample((OM[OM$age > 6 & OM$age <= 12, "OM"] + 10e-4), ndvar(), replace = TRUE), "Hand to mouth", "Hand to Mouth")
OM.6_12.dist <- fitdist(OM[OM$age > 6 & OM$age <= 12, "OM"], "weibull", method = "mle")
OM.6_12.mcstoc <- mcstoc(rweibull, type="V", shape = OM.f6.dist$estimate[[1]], scale = OM.f6.dist$estimate[[2]], rtrunc=TRUE, linf=0)

bestFitDist(data = sample((OM[OM$age > 12 & OM$age <= 24, "OM"] + 10e-4), ndvar(), replace = TRUE), "Hand to mouth", "Hand to Mouth")
OM.12_24.dist <- fitdist(OM[OM$age > 12 & OM$age <= 24, "OM"], "weibull", method = "mle")
OM.12_24.mcstoc <- mcstoc(rweibull, type="V", shape = OM.f6.dist$estimate[[1]], scale = OM.f6.dist$estimate[[2]], rtrunc=TRUE, linf=0)

bestFitDist(data = sample((OM[OM$age > 24, "OM"] + 10e-4), ndvar(), replace = TRUE), "Hand to mouth", "Hand to Mouth")
OM.24_36.dist <- fitdist(OM[OM$age > 24, "OM"], "weibull", method = "mle")
OM.24_36.mcstoc <- mcstoc(rweibull, type="V", shape = OM.f6.dist$estimate[[1]], scale = OM.f6.dist$estimate[[2]], rtrunc=TRUE, linf=0)

bestFitDist(data = sample((OM[OM$age > 24, "OM"] + 10e-4), ndvar(), replace = TRUE), "Hand to mouth", "Hand to Mouth")
OM.36_48.dist <- fitdist(OM[OM$age > 24, "OM"], "weibull", method = "mle")
OM.36_48.mcstoc <- mcstoc(rweibull, type="V", shape = OM.f6.dist$estimate[[1]], scale = OM.f6.dist$estimate[[2]], rtrunc=TRUE, linf=0)

################# Surface area of object that is mouthed ################
# Oz cites Leckie exponential(min = 1, mean = 10, max = 50)
# low <- 1
# high <- 50
# r <- 0.11 ## Just adjusted r by trial and error until I got the mean = 10 that Oz used
# C <- exp(-r*high); D <- exp(-r*low) 
# n <- 10000 
# U <- runif(n, min = C, max = D) 
# X <- (1/r)*log(1/U) 
# hist(X, breaks=100, xlim=c(0,50)) 
# summary(X)

obj.SA.mouthed.mcstoc <- mcstoc(rexp, type = "V", rate = 0.11, rtrunc = TRUE, linf = 1, lsup = 50)
hist(obj.SA.mouthed.mcstoc)


####### Object-to-mouth transfer efficiency ##########

# Oz again used beta(2, 8) which I think is too low, but it should be similar to SEE 
#Use the saliva extraction efficency (SEE)
OMRE.mcstoc <- SEE.mcstoc




#####################################################################################
#####################################################################################

################## Soil ingestion directly placing soil in the mouth per day #################

#####################################################################################
######################################################################################

####################### Amount of soil ingested per contact #######################
###### Lick contacts 
# Huang 2015 Shape dynamics
# How many licks does it take to get to the centre of a lollipop?
# For a lollipop for radius = 1 cm at 1 lick/cm2, it would take approx 1000 licks

# Hard candy density 
# (https://prezi.com/pii7hfuaq-vz/jolly-good-metrics/)
# Jolly rancher mass = 6.3 g, width = 1.5 cm, length = 2.5 cm, height = 1 cm, vol = 3.75 cm; density = 6.3g / 3.75 cm2 = 1.68 g/cm3
# 
# Mass of Huang's hard candy, r = 1 cm, vol = 4/3 * pi * r^3 = 4/3*pi = 4.18879 cm3, mass = density * vol = 1.68 g/ cm3 * 4.18879 cm3 = 7.037 g 
# 
# 7.037 g/ approx 1000 licks = 0.007037/licks = 7.037 mg/licks


##### Eating pieces of dirt 
# mass of a grain of rice = 1/64 g = 0.015625 g 
# Volume of a small grain of rice is 0.02 cm3, a large grain of rice is 0.07 cm3
# Assuming soil density of 1.0 g/cm3 --> each ingestion of a small piece of soil would be 0.02 g = 20 mg or 0.07 g = 70 mg of soil.

## BIg change made 26 May 2017: I had been using a triangular dist in mg improperly converted to g, currently all number are in the correct g amts
#triangular distribution(min = 0.007  g , max =  0.070  g , mode =  0.020  mg) --> beta distribution
fitdist(rtriangle(100000, 0.007, 0.070, 0.020), "beta", method = "mle") # For a triangular dist, the values must be 0-1
hist(rtriangle(10000, 0.007, 0.070, 0.020)) # amt in g

# Fitting of the distribution ' beta ' by maximum likelihood 
# Parameters:
#   estimate Std. Error
# shape1   5.295952 0.02297762
# shape2 158.621280 0.72073235
hist(rbeta(10000, 5.30, 158.62)) # amt in g
## How do I convert this to a beta dist for mg? maybe it's just 100*beta(shape1, shape2)
fitdist(rbeta(10000, 5.30, 158.62)*1000, "beta", method = "mle") # amt in mg? 
summary(rbeta(10000, 5.30, 158.62)*1000) # amt in mg
summary(rtriangle(100000, 0.007, 0.070, 0.020)*1000)

# amt directly consumed per direction ingestion (in mg)
SM.ingested.amt.mcstoc <- mcstoc(rbeta, type = "V", shape1 = 5.30, shape2 = 158.62, rtrunc=TRUE, linf=0)*1000




######################### Fraction of children every ingesting soil ####################
# HM.SM$SM !is.na()/n for each age group would give us the fraction of children that consumed soil in age each group
# HM.SM$SM is the freq of soil-to-mouth (direct ingestion) contacts for all the SO and VO observations
# More informative for the reader to report 1) % of children that ever placed soil in their mouths and 2) frequency of direct ingestion for consumers only. 
# With these numbers they could calc freq for direct ingestion for all kids, if they wanted to, but given only freq of direct ingestion for all kids, they could not back out the freq of consumers only

# # Rate for ALL children by age.group (na replaced by 0 so na.rm = TRUE has no effectmakes this consumers only)
# SM.allchildren <- HM.SM %>% 
#   replace(is.na(.), 0) %>% # replace(is.na(.), 0)
#   group_by(age.group) %>%
#   summarise(n = n(), SM.mean = mean(SM, na.rm = TRUE), SM.sd = sd(SM, na.rm = TRUE))
# # age.group     n    SM.mean     SM.sd
# # (int) (int)      (dbl)     (dbl)
# # 1         2     1 0.00000000       NaN
# # 2         3    23 0.04912988 0.1656385
# # 3         4   123 0.31192981 0.8484990
# # 4         5    54 0.23501688 0.4707441
# # 5         6     5 0.00000000 0.0000000

# # Rate for consumers only by age.group (the na.rm = TRUE makes this consumers only)
# SM.consumersOnly <- HM.SM %>% 
#   group_by(age.group) %>%
#   summarise(n = n(), NA.count = sum(is.na(SM)), consumers = n-NA.count, frac.consumers = consumers/n, 
#             SM.mean = mean(SM, na.rm = TRUE), SM.sd = sd(SM, na.rm = TRUE), SM.median = median(SM, na.rm = TRUE))
# 
# # age.group n NA.count consumers frac.consumers   SM.mean     SM.sd
# # (int) (int)    (int)     (int)          (dbl)     (dbl)     (dbl)
# # 1         2     1        1         0     0.00000000       NaN       NaN
# # 2         3    23       21         2     0.08695652 0.5649936 0.1437875
# # 3         4   123       90        33     0.26829268 1.1626475 1.3133233
# # 4         5    54       40        14     0.25925926 0.9064937 0.4979302
# # 5         6     5        5         0     0.00000000       NaN       NaN
# 
# ## Could use SM.consumersOnly$frac.consumers to get the fraction of children who ingested soil in each age group.
# ## Instead of determining the fraction of consumers based on my small observational dataset, use the WASHB survey results


# ## WASHB midline and endline survey had the question:
# "Has the [NAME] eaten any dirt or soil? Ask for each recall period: Today, yesterday, day before yesterday, in the past 7 days"
# Since the survey was conducted before the day was over, I would expect the results for "Has your child eaten any dirt or soil today" to be underestimates
# The results for yesterday and the day before yesterday should be similar (both reflect a 24-hr recall period) and higher than the results for "today"
# Use the mean (?) of the yesterday and two day results

soil.sevenday.survey.mid <- anth.soilMass.wide[anth.soilMass.wide$arm == "Control", c("survey.age.group.mid", "soil.today.mid", "soil.yest.mid", "soil.twoday.mid", "soil.sevenday.mid")]
names(soil.sevenday.survey.mid) <- c("age.group", "soil.today", "soil.yest", "soil.twoday", "soil.sevenday")

soil.sevenday.survey.end <- anth.soilMass.wide[anth.soilMass.wide$arm == "Control", c("survey.age.group.end", "soil.today.end", "soil.yest.end", "soil.twoday.end", "soil.sevenday.end")]
names(soil.sevenday.survey.end) <- c("age.group", "soil.today", "soil.yest", "soil.twoday", "soil.sevenday")

soil.survey <- rbind(soil.sevenday.survey.mid, soil.sevenday.survey.end)

soil.survey %>%
  group_by(age.group) %>%
  summarise(n = n(), today.mean = mean(soil.today, na.rm = TRUE), yest.mean = mean(soil.yest, na.rm = TRUE),
            twoday.mean = mean(soil.twoday, na.rm = TRUE), sevenday.mean = mean(soil.sevenday, na.rm = TRUE), 
            oneday.recall.diff = (twoday.mean - yest.mean), oneday.recall.mean = (yest.mean + twoday.mean)/2)

# age.group     n today.mean  yest.mean twoday.mean sevenday.mean oneday.recall.diff oneday.recall.mean
# (int) (int)      (dbl)      (dbl)       (dbl)         (dbl)              (dbl)              (dbl)
# 1         3    44 0.00000000 0.00000000  0.00000000     0.0000000        0.000000000         0.00000000
# 2         4   985 0.38414634 0.52182741  0.52899288     0.5524975        0.007165468         0.52541015
# 3         5    40 0.35000000 0.47500000  0.52500000     0.6000000        0.050000000         0.50000000
# 4         6   857 0.10630841 0.26635514  0.28271028     0.4228972        0.016355140         0.27453271
# 5         7   202 0.06467662 0.07462687  0.07462687     0.1044776        0.000000000         0.07462687 ## There are some children in age group SEVEN (? = 36-48 months old, tested because the target child wasn't around?) 
# 6        NA     4 0.00000000 0.00000000  0.00000000     0.0000000        0.000000000         0.00000000 

# Good, the diff between the survey results for "yesterday" and "twoday" are similar (max diff is 5%, which is not small, which is 10% of the total 50% of children reported to eat soil)

# ## Observed consumption vs surveyed
# <6 mo:     2/24 =  8.3% vs today: 0.0%,  yest: 0.0%,  twoday: 0.0%,  sevenday: 0.0%, oneday.recall: 0.0%     ##### Review the <6 contacts and ensure they are actually mouthing - why didn't mother's report?
# 6-12 mo: 33/123 = 26.8% vs today: 38.5%, yest: 52.2%, twoday: 52.9%, sevenday: 55.2%, oneday.recall: 52.5% 
# 12-24 mo: 14/54 = 25.9% vs today: 35.0%, yest: 47.5%, twoday: 52.5%, sevenday: 60.0%, oneday.recall: 50.0% 
# 24-36 mo:   0/5 =  0.0% vs today: 10.6%, yest: 26.6%, twoday: 28.2%, sevenday: 42.3%, oneday.recall: 27.4% 


# Can justify because our observation covered only about half of the child's waking hours. 
# If direct soil consumption rates are the same over the course of the day, this means that we missed observing approximately half of the events
# We observed this half of events among different children than the ones we observed, this would double the fraction of children who mouthed soil. 
# This supports applying the same ingestion frequency for consumer only to all of those who consumed. 
# --> i.e. if consumers-only rate is 0.5 events/hr , this applies to the 50% of kids who consumed rather than the 0.25 events/hr rate applying to every kid (which assumes they all consumed)

# Could use sevenday parental report as long-term average fraction of children consuming soil, BUT I am only modeling for ONE day so I should use only the one-day results. 
# Otherwise I will overestimate the daily (though be closer to the long-term average?) 

## The daily recall values are the fraction of ingesters per DAY. The fraction that ingests per HOUR is much smaller
# --> the daily rate divided by the number of hours awake, with the assumption that each child ingesting soil in that hour does not ingest soil at any other time of day 
# (this assumption seems bogus, but i'm not sure how toto do it otherwise)
# <6 mo:  oneday.recall: 0.0%  but we'll use observed 8.3%  --> 8.3/(24 - rnorm(ndvar(), 13.6, 2.1)) = 0.00610 # hr asleep mean = 13.6, sd = 2.1
# 6-12 mo: oneday.recall: 52.5%  --> 52.5/(24 - rnorm(ndvar(), 12.9, 1.3)) = 0.04069 # hr asleep mean = 12.9, sd = 1.3
# 12-24 mo: oneday.recall: 50.0% --> 50.0/(24 - rnorm(ndvar(), 12.6, 1.3)) = 0.03968 # hr asleep mean = 12.6, sd = 1.3
# 24-36 mo:   oneday.recall: 27.4% --> 27.4/(24 - rnorm(ndvar(), 12.0, 1.2)) = 0.022833  # hr asleep mean = 12.0, sd = 1.2

## Don't do soil consumption by hour, just do by day
SM.consumerfrac.p1.f6 <- (8.3)/100 
SM.consumerfrac.p1.6_12 <- (52.5)/100
SM.consumerfrac.p1.12_24 <- (50.0)/100
SM.consumerfrac.p1.24_36 <- (27.4)/100
SM.consumerfrac.p1.36_48 <- (27.4)/100 ## I don't have a value for children this age, so use the 24_36 soil ingetstion consumer fractions

# # percent consumers/day divided by the number of hours awake per day --> to get this into a fraction between 0 and 1, must divide by 100
# SM.consumerfrac.p1.f6 <- (8.3/(24 - rnorm(ndvar(), 13.6, 2.1)))/100 
# SM.consumerfrac.p1.6_12 <- (52.5/(24 - rnorm(ndvar(), 12.9, 1.3)))/100
# SM.consumerfrac.p1.12_24 <- (50.0/(24 - rnorm(ndvar(), 12.6, 1.3)))/100
# SM.consumerfrac.p1.24_36 <- (27.4/(24 - rnorm(ndvar(), 12.0, 1.2)))/100

######################### Number of direct soil-to-mouth contacts ####################
# Soil ingestion rate for each child within each age group (not summarized by age group)
## OF THOSE WHO CONSUMED SOIL - NOT CONSUMERS WILL BE REMOVD IN THE MCPROBTREE
## recall HM.SM$SM_freq is per HOUR - we need per DAY
SM.f6 <- HM.SM[HM.SM$age.group %in% c(2,3) & !(is.na(HM.SM$SM_freq)) & HM.SM$SM_freq > 0, "SM_freq"]
SM.f6.mcstoc <- mcstoc(rempiricalD, values = SM.f6, type="V", nsv = ndvar())
SM.6_12 <- HM.SM[HM.SM$age.group %in% c(4) & !(is.na(HM.SM$SM_freq)) & HM.SM$SM_freq > 0, "SM_freq"]
SM.6_12.mcstoc <- mcstoc(rempiricalD, values = SM.6_12, type="V", nsv = ndvar())
SM.12_24 <- HM.SM[HM.SM$age.group %in% c(5) & !(is.na(HM.SM$SM_freq)) & HM.SM$SM_freq > 0, "SM_freq"]
SM.12_24.mcstoc <- mcstoc(rempiricalD, values = SM.12_24, type="V", nsv = ndvar())
SM.24_36 <- HM.SM[HM.SM$age.group %in% c(6) & !(is.na(HM.SM$SM_freq)) & HM.SM$SM_freq > 0, "SM_freq"]
SM.24_36.mcstoc <- mcstoc(rempiricalD, values = SM.12_24, type="V", nsv = ndvar())
SM.36_48 <- HM.SM[HM.SM$age.group %in% c(7) & !(is.na(HM.SM$SM_freq)) & HM.SM$SM_freq > 0, "SM_freq"]
SM.36_48.mcstoc <- mcstoc(rempiricalD, values = SM.12_24, type="V", nsv = ndvar())

# In 2/24 (8.3%) observations of children <6 months old, children directly consumed soil with an average of 0.56 times/hr (sd = 0.14 events/hr). 
# Among 123 observations of children 6-12 months old, there were 33 (26.8%) observations of children putting soil into their mouths  (mean = 1.16 events/hr, sd = 1.31 events/hr) 
# and among 54 observations children 12-24 months old, there were 14 (25.9%) observations of soil consumption (mean = 0.90 events/hr, sd = 0.50 events/hr).
# None of the five children 24-36 months old directly consumed soil during the observation period. 


##########################################################################

######## Convert times per hour to times per day by selecting DIFF hour freq for each hour awake during the day then summing ###

####################### Duration awake each day #########################
# From Gallan, 2011
# By age group
# Don't use specific months of age because there isn't enough data for some of the months
## Seems a bit anamolous that hours slept is lower for 9-month-old children than 6 and 12 month old children so I'll use the 6- and 12-month-old values
awake.hr.f6 <- 24 - mcstoc(rnorm, type = "V", mean = 13.6, sd = 2.1, rtrunc=TRUE, linf=0) # use for f6
awake.hr.6_12 <- 24 - mcstoc(rnorm, type = "V", mean = 12.9, sd = 1.3, rtrunc=TRUE, linf=0) # use for 6_12
awake.hr.12_24 <- 24 - mcstoc(rnorm, type = "V", mean = 12.6, sd = 1.3, rtrunc=TRUE, linf=0) # use for 12-24
awake.hr.24_36 <- 24 - mcstoc(rnorm, type = "V", mean = 12.0, sd = 1.2, rtrunc=TRUE, linf=0) # use for 24-36
awake.hr.36_48 <- 24 - mcstoc(rnorm, type = "V", mean = 12.0, sd = 1.2, rtrunc=TRUE, linf=0) # use for 36-48


# types of mouthing (these are all in /hr and need to be converted to /day)
mouthing.types <- c(
  "HM.child.nd", # non-dietary         ## in the output mouthingTable, child != child.nd + child.d because they are calc from diff sample draws
  "HM.child.d", # dietary
  "HM.mom.nd", # caregiver non-dietary
  "HM.mom.d.events", # caregiver dietary
  "OM",
  "SM",
  "HM.child" # non-dietary + dietary
)

age.group.names <- c("f6", "6_12", "12_24", "24_36", "36_48")
for(i in 1:length(age.group.names)){
  age.group.name <- age.group.names[i]
  assign(paste("mouthing.", age.group.name, sep = ""), matrix(NA, ncol = length(mouthing.types), nrow = ndvar, dimnames = list(NULL, mouthing.types)))
}

## create a table for each age group f6, 6_12, 12_24, 24_36, 36_48 to store the daily
# The table is ndvar() long with the col each of the mouthing types 
allMouthingTable <- function(mouthing.type, age.group, i, j, mouthingTable){
  Aw.age.group <- paste("Aw.", age.group, sep = "")
  assign(Aw.age.group, sample(get(paste("awake.hr.", age.group, sep="")), 1))
  assign(paste(mouthing.type, ".", Aw.age.group, ".base", sep = ""), sample(get(paste(mouthing.type, ".", age.group, ".mcstoc", sep = "")), floor(get(Aw.age.group))))
  assign(paste(mouthing.type, ".", Aw.age.group, ".extra", sep = ""), sample(get(paste(mouthing.type, ".", age.group, ".mcstoc", sep = "")), 1) * get(Aw.age.group) %% 1)
  assign(paste(mouthing.type, ".", age.group, sep = ""), sum(get(paste(mouthing.type, ".", Aw.age.group, ".base", sep = "")), na.rm = TRUE) + get(paste(mouthing.type, ".", Aw.age.group, ".extra", sep = "")))
  mouthingTable[i, j] <- round(get(paste(mouthing.type, ".", age.group, sep = "")), digits = 1)
  mouthingTable[, "HM.child"] <- rowSums(mouthingTable[, c("HM.child.nd", "HM.child.d")])
  return(mouthingTable) ## yaya! Just like in Java. Otherwise I need to set mouthing.f6 to save in the global environment (I've passed a copy, not by reference into the funciton; don't know how to pass by ref in R)
}

for(k in 1:length(age.group.names)){
  age.group <- age.group.names[k]
  for(j in 1:length(mouthing.types)){
    mouthing.type <- mouthing.types[j]
    for(i in 1:ndvar){
      mouthingTable <- get(paste("mouthing.", age.group, sep = ""))
      assign(paste("mouthing.", age.group, sep = ""), allMouthingTable(mouthing.type, age.group, i, j, mouthingTable))
    }
  }
}


################ Back to soil consumption ##############
### For a child that is a consumer, multiply by the DAILY rate of consumption (in mouthingTable)
## Children are either consumers or not -> use a distribution of 0s and 1s, where 1s represent consumers
#Use mcprobtree instead of this: SM.consumerfrac.f6 <-    mcstoc(rempiricalD, values = c(rep(0, 10000), rep(1, 6)), type="V", nsv = ndvar()) # This is observed, not survey report, as the other values are. 

SM.freq.0 <- mcdata(0, type = "V")
SM.f6.day.mcnode <- mcstoc(rempiricalD, values = mouthing.f6[,"SM"], type = "V", nsv = ndvar())
SM.fracfreq.f6 <- mcprobtree(c(1 - SM.consumerfrac.p1.f6, SM.consumerfrac.p1.f6), list("0" = SM.freq.0, "1" = SM.f6.day.mcnode), type = "V")
SM.6_12.day.mcnode <- mcstoc(rempiricalD, values = mouthing.6_12[,"SM"], type = "V", nsv = ndvar())
SM.fracfreq.6_12 <- mcprobtree(c(1 - SM.consumerfrac.p1.6_12, SM.consumerfrac.p1.6_12), list("0" = SM.freq.0, "1" = SM.6_12.day.mcnode), type = "V")
SM.12_24.day.mcnode <- mcstoc(rempiricalD, values = mouthing.12_24[,"SM"], type = "V", nsv = ndvar())
SM.fracfreq.12_24 <- mcprobtree(c(1 - SM.consumerfrac.p1.12_24, SM.consumerfrac.p1.12_24), list("0" = SM.freq.0, "1" = SM.12_24.day.mcnode), type = "V")
SM.24_36.day.mcnode <- mcstoc(rempiricalD, values = mouthing.24_36[,"SM"], type = "V", nsv = ndvar())
SM.fracfreq.24_36 <- mcprobtree(c(1 - SM.consumerfrac.p1.24_36, SM.consumerfrac.p1.24_36), list("0" = SM.freq.0, "1" = SM.24_36.day.mcnode), type = "V")
SM.36_48.day.mcnode <- mcstoc(rempiricalD, values = mouthing.36_48[,"SM"], type = "V", nsv = ndvar())
SM.fracfreq.36_48 <- mcprobtree(c(1 - SM.consumerfrac.p1.36_48, SM.consumerfrac.p1.36_48), list("0" = SM.freq.0, "1" = SM.36_48.day.mcnode), type = "V")


### Based on rough estimates 
# Fraction of consumers * freq of direct ingestion/hr * est amt ingested/ingestion [mg] * hours awake
# > 3/24*0.565*(0.02*1000)*awake.hr.f6
#   node    mode  nsv nsu nva variate  min mean median  max Nas type outm
# 1    x numeric 1001   1   1       1 6.27 14.6   14.6 23.9   0    V each
# > 33/123*1.163*(0.02*1000)*awake.hr.6_12
#   node    mode  nsv nsu nva variate  min mean median max Nas type outm
# 1    x numeric 1001   1   1       1 42.1 69.4   69.3 101   0    V each
# > 14/54*0.906*(0.02*1000)*awake.hr.12_24
#   node    mode  nsv nsu nva variate  min mean median  max Nas type outm
# 1    x numeric 1001   1   1       1 31.5 53.5   53.6 70.6   0    V each

# would increase the amount of soil ingested by children <6 months old an average of 14.6 mg/day. 
# Total soil ingestion by children 6-12 months old would increase by 69.4 mg/day and 
# for children 12-24 months old, the average increase would be 53.5 mg/day. 
# As no direct soil consumption was observed among children 24-36 months old, 
# their estimated soil ingestion rate would not increase from the inclusion of direct soil ingest in the model. 

soil.direct.f6 <-  SM.fracfreq.f6 * SM.ingested.amt.mcstoc # mg # as of 25 May 2017 - mean 12 mg
soil.direct.6_12 <-  SM.fracfreq.6_12 * SM.ingested.amt.mcstoc # mg # as of 25 May 2017 - mean 210 mg
soil.direct.12_24 <- SM.fracfreq.12_24 * SM.ingested.amt.mcstoc # mg # as of 25 May 2017 - mean 170 mg
soil.direct.24_36 <-  SM.fracfreq.24_36 * SM.ingested.amt.mcstoc # mg # as of 25 May 2017 - mean 95 mg
soil.direct.36_48 <-  SM.fracfreq.36_48 * SM.ingested.amt.mcstoc # mg # as of 25 May 2017 - mean 95.3 mg








########## Eval the mc model ###########
# parameter.names <- c("child.hand.soil.concentration.mg.cm2", "child.hand.surface.area.cm2", 
#                      "child.hand.mouth.frequency.events.day", "child.hand.fraction.mouthed", 
#                      "mom.hand.soil.load.mg", "mom.hand.surface.area.cm2", "mom.hand.mouth.nonfood.frequency.events.day",
#                      "mom.feeding.frequency.events.day", "mom.hand.fraction.mouthed", 
#                      "saliva.removal.efficiency", "hours.awake")
# parameter.values_f6 <- c(soil.C.conc, C.hand.SA.WASHB.f6.mcstoc, HM.child.f6.mcstoc, HF.child.mcstoc, soil.M.load.mcstoc, 
#                       M.hand.SA.WASHB.mcstoc, HM.other.nd.f6.mcstoc, HM.other.d.f6.mcstoc, SEE.mcstoc, awake.day.f6)
# names(parameter.values_f6) <- parameter.names
# parameter.values_6_12 <- c(soil.C.conc, C.hand.SA.WASHB.6_12.mcstoc, HM.child.6_12.mcstoc, HF.child.mcstoc, soil.M.load.mcstoc, 
#                          M.hand.SA.WASHB.mcstoc, HM.other.nd.6_12.mcstoc, HM.other.d.6_12.mcstoc, SEE.mcstoc, awake.day.6_12)
# names(parameter.values_6_12) <- parameter.names
# parameter.values_12_24 <- c(soil.C.conc, C.hand.SA.WASHB.12_24.mcstoc, HM.child.12_24.mcstoc, HF.child.mcstoc, soil.M.load.mcstoc, 
#                            M.hand.SA.WASHB.mcstoc, HM.other.nd.12_24.mcstoc, HM.other.d.12_24.mcstoc, SEE.mcstoc, awake.day.12_24)
# names(parameter.values_12_24) <- parameter.names
# parameter.values_24_36 <- c(soil.C.conc, C.hand.SA.WASHB.24_36.mcstoc, HM.child.24_36.mcstoc, HF.child.mcstoc, soil.M.load.mcstoc, 
#                            M.hand.SA.WASHB.mcstoc, HM.other.nd.24_36.mcstoc, HM.other.d.24_36.mcstoc, SEE.mcstoc, awake.day.24_36)
# names(parameter.values_24_36) <- parameter.names

child.hand.soil.concentration.mg.cm2 <- soil.C.conc
# child.hand.mouth.frequency.events.hr
# child.hand.surface.area.cm2
child.hand.fraction.mouthed <- HF.child.mcstoc
mom.hand.soil.load.mg <- soil.M.load.mcstoc
# mom.hand.mouth.nonfood.frequency.events.hr
# mom.feeding.frequency.events.hr
mom.hand.surface.area.cm2 <- M.hand.SA.WASHB.mcstoc
# mom.hand.fraction.mouthed 
saliva.removal.efficiency <- SEE.mcstoc
obj.soil.concentration <- object.soil.concentration.mcstoc
# OM.mcstoc 
obj.SA.mouthed <- obj.SA.mouthed.mcstoc
# object.saliva.removal.efficiency <- OMRE.mcstoc
soil.directly.ingested.mg <- SM.ingested.amt.mcstoc
# hours.awake

# resampling frequency
# sample only once per simulation: child.hand.surface.area.cm2, hours.awake, mom.hand.surface.area.cm2
# sample once per hour: all other parameters 

#### Function to get results  ##############
# age.group = c(f6, 6_12, 12_24, 24_36)
# title = "\n Soil ingestion among children <6 months"

soil.ingestion.results <- function(age.group, title){ # title = "\n Soil ingestion among children <6 months"
  
  mouthingTable <- data.frame(get(paste("mouthing.", age.group,sep = "")))
  child.hand.surface.area.cm2 <- get(paste("C.hand.SA.WASHB.", age.group, ".mcstoc", sep = ""))
  child.hand.mouth.frequency.events.day <- mcdata(mouthingTable[,"HM.child"], type = "V", nsv = ndvar())  ## This data is already selected from a dist, so use the actual data with mcdata rather than creating a distribution using rempiricalD
  mom.hand.mouth.nonfood.frequency.events.day <- mcdata(mouthingTable[,"HM.mom.nd"], type = "V", nsv = ndvar())
  mom.feeding.frequency.events.day <- mcdata(mouthingTable[,"HM.mom.d.events"], type = "V", nsv = ndvar())
  mom.hand.fraction.mouthed <- get(paste("HF.ofmom.", age.group, sep = ""))
  child.obj.mouth.frequency.events.day  <- mcdata(mouthingTable[,"OM"], type = "V", nsv = ndvar())
  child.frequency.ingested.soil.events.day <- get(paste("SM.fracfreq.", age.group, sep = "")) # Includes the fraction of children ingesting and the frequency of those ingesting
  #hours.awake <- get(paste("awake.day.", age.group, sep = ""))
  
  soil.mcmodel.plot <- mcmodel({
    child.hand.soil.day <- (child.hand.soil.concentration.mg.cm2 * (child.hand.surface.area.cm2 * child.hand.fraction.mouthed)) * child.hand.mouth.frequency.events.day * saliva.removal.efficiency
    mom.hand.soil.day <- ((mom.hand.soil.load.mg / mom.hand.surface.area.cm2) * (mom.hand.surface.area.cm2 * mom.hand.fraction.mouthed) * (mom.hand.mouth.nonfood.frequency.events.day + mom.feeding.frequency.events.day)) * saliva.removal.efficiency
    child.obj.soil.day <- (obj.soil.concentration * obj.SA.mouthed) * child.obj.mouth.frequency.events.day  * saliva.removal.efficiency
    child.direct.soil.day <- soil.directly.ingested.mg * child.frequency.ingested.soil.events.day
    amt.soil.ingested.mg.day <- (child.hand.soil.day + mom.hand.soil.day + child.obj.soil.day + child.direct.soil.day)
    
    
    mc(child.hand.soil.concentration.mg.cm2, child.hand.surface.area.cm2, child.hand.mouth.frequency.events.day, child.hand.fraction.mouthed,
       mom.hand.soil.load.mg, mom.hand.surface.area.cm2, mom.hand.mouth.nonfood.frequency.events.day, mom.feeding.frequency.events.day, mom.hand.fraction.mouthed,
       obj.soil.concentration, child.obj.mouth.frequency.events.day, obj.SA.mouthed, saliva.removal.efficiency,
       soil.directly.ingested.mg, child.frequency.ingested.soil.events.day, #hours.awake, 
       child.hand.soil.day, mom.hand.soil.day, child.obj.soil.day, child.direct.soil.day, amt.soil.ingested.mg.day) #sleep.day.24_36,
    #   mc(child.hand.soil.day, mom.hand.soil.day, child.obj.soil.day, hours.awake, amt.soil.ingested.mg.day)
  })
  
  soil.evalmcmodel.plot <- evalmcmod(soil.mcmodel.plot, nsv=ndvar(), nsu=ndunc(), seed=seed)
  plot(soil.evalmcmodel.plot, prec = 0.001, stat = c("median", "mean"), lim = c(0.025, 0.25, 0.75, 0.975), na.rm =
         TRUE, griddim = NULL, xlab = NULL, ylab = "Fn(x)", main = "", draw = TRUE, paint = TRUE)
  hist(soil.evalmcmodel.plot)
  summary(soil.evalmcmodel.plot, probs = c(0, 0.5, 0.95, 1), lim = c(0.025, 0.975))
  
  
  # soil.mcmodel <- mcmodel({
  #   mc(child.hand.soil.concentration.mg.cm2, child.hand.surface.area.cm2, child.hand.mouth.frequency.events.hr, child.hand.fraction.mouthed,
  #      mom.hand.soil.load.mg, mom.hand.surface.area.cm2, mom.hand.mouth.nonfood.frequency.events.hr, mom.feeding.frequency.events.hr, mom.hand.fraction.mouthed,
  #      obj.soil.concentration, child.obj.mouth.frequency.events.hr, obj.SA.mouthed, saliva.removal.efficiency,
  #      soil.directly.ingested.mg, child.frequency.ingested.soil.events.hr,
  #      hours.awake, child.hand.soil.hr, mom.hand.soil.hr, child.obj.soil.hr, child.direct.soil.hr, amt.soil.ingested.mg.day) #sleep.hr.24_36,
  #   #   mc(child.hand.soil.hr, mom.hand.soil.hr, child.obj.soil.hr, hours.awake, amt.soil.ingested.mg.day)
  # })
  # soil.evalmcmodel <- evalmcmod(soil.mcmodel, nsv=ndvar(), nsu=ndunc(), seed=seed)
  # summary(soil.evalmcmodel, probs = c(0, 0.5, 0.95, 1), lim = c(0.025, 0.975))
  # 
  # # print(soil.evalmcmodel, digits = 2)
  # # write.csv(summary(soil.evalmcmodel), "C:/Users/Tareq/Box Sync/VO Soil/Mass of soil on children's hands/summary.csv")
  # 
  # #tornado(x, output=length(x), use=all.obs, method=c(spearman, kendall,pearson), lim=c(0.025,
  # #                                                                                     0.975))
  # #tornadounc(mc,output = length(mc), quant=c(0.5,0.75,0.975), use = all.obs,
  # #           method=c(spearman,kendall,pearson), ...)s
  # soil.mcmodel.tor <- mcmodel({
  #   mc(child.hand.soil.concentration.mg.cm2, child.hand.surface.area.cm2, child.hand.mouth.frequency.events.hr, child.hand.fraction.mouthed,
  #      mom.hand.soil.load.mg, mom.hand.surface.area.cm2, mom.hand.mouth.nonfood.frequency.events.hr, mom.feeding.frequency.events.hr, mom.hand.fraction.mouthed,
  #      obj.soil.concentration, child.obj.mouth.frequency.events.hr, obj.SA.mouthed, saliva.removal.efficiency,
  #      soil.directly.ingested.mg, child.frequency.ingested.soil.events.hr,
  #      hours.awake,  amt.soil.ingested.mg.day) #sleep.hr.24_36,
  #   #   mc(child.hand.soil.hr, mom.hand.soil.hr, child.obj.soil.hr, hours.awake, amt.soil.ingested.mg.day)
  # })
  # soil.evalmcmodel.tor <- evalmcmod(soil.mcmodel.tor, nsv=ndvar(), nsu=ndunc(), seed=seed)
  # tor.soil.evalmcmodel  <- tornado(soil.evalmcmodel.tor, use = "pairwise.complete.obs") #  I'd prefer to use complete.obs rather than pairwise (I think) but this doesn't work. method = "spearman", lim = c(0.025, 0.975)
  # plot(tor.soil.evalmcmodel )
  # #print(tor.soil.evalmcmodel )
  # 
  # # torunc.soil.evalmcmodel  <- tornadounc(soil.soil.evalmcmodel, output = "amt.soil.ingested.mg.day", quant = c(0.5, 0.75,0.975), use = "pairwise.complete.obs") #  I'd prefer to use complete.obs rather than pairwise (I think) but this doesn't work. method = "spearman", lim = c(0.025, 0.975)
  # # plot(torunc.soil.evalmcmodel )
  # # print(torunc.soil.evalmcmodel )
  # 
  # print(soil.evalmcmodel, digits = 2)
  # summary(soil.evalmcmodel, probs = c(0, 0.5, 0.95, 1), lim = c(0.025, 0.975))
  # #print(tor.soil.evalmcmodel)
  # 
  # #soil.eat.summary <- summary(amt.soil.ingested.mg.day, probs = c(0, 0.5, 0.95, 1), lim = c(0.025, 0.975))
  # #print(soil.eat.summary)
  # #soil.eat.summary.log <- summary(log10(amt.soil.ingested.mg.day), probs = c(0, 0.5, 0.95, 1), lim = c(0.025, 0.975))
  # #print(soil.eat.summary.log)
  # #hist(log10(amt.soil.ingested.mg.day), xlab = "log10(quantity of soil ingested) (mg)", ylab = "number of children", ylim = c(0, 20000), main = "\n Soil ingestion among children <6 months")
}

###################
# # Report the mean and median of the exposure means with a 95% confidence interval (CI95).
# mean.expo <- sapply(1:ndunc(), function(j) mean(eatSoil_f6$soil.test_f6[, j, ]))
# mean(mean.expo)
# quantile(mean.expo, probs = seq(0, 1, 0.025))[c("2.5%", "50%", "97.5%")]
# 
# # Generate an "ecdf" = empirical cumulative distribution function plot. This actually calls plot.mcnode().
# plot(soilevalmcmodel.plot, xlim = c(0, 100), ylim = c(0, 1), main = "Daily Soil Consumption by Children in Rural Bangladesh <6 months old", ylab = "Proportion of Population", xlab = "Daily Consumption of Soil")

soil.ingestion.results(age.group = "f6", title = "\n Soil ingestion among children <6 months") # median 117
soil.ingestion.results(age.group = "6_12", title = "\n Soil ingestion among children 6-12 months") # 377
soil.ingestion.results(age.group = "12_24", title = "\n Soil ingestion among children 12-24 months") #407
soil.ingestion.results(age.group = "24_36", title = "\n Soil ingestion among children 24-36 months") # 499
# but one val missing the 36_48 category 646


child.hand.soil <- (child.hand.soil.concentration.mg.cm2 * (child.hand.surface.area.cm2 * child.hand.fraction.mouthed)) * child.hand.mouth.frequency.events.hr * saliva.removal.efficiency * hours.awake
mom.hand.soil <- ((mom.hand.soil.load.mg / mom.hand.surface.area.cm2) * (mom.hand.surface.area.cm2 * mom.hand.fraction.mouthed) * (mom.hand.mouth.nonfood.frequency.events.hr + mom.feeding.frequency.events.hr)) * saliva.removal.efficiency * hours.awake
child.obj.soil <- (obj.soil.concentration * obj.SA.mouthed) * child.obj.mouth.frequency.events.hr  * saliva.removal.efficiency* hours.awake
child.direct.soil <- soil.directly.ingested.mg * child.frequency.ingested.soil.events.hr * hours.awake
amt.soil.ingested.mg.day <- (child.hand.soil.hr + mom.hand.soil.hr + child.obj.soil.hr + child.direct.soil.hr) * hours.awake

plot(child.hand.soil)
plot(mom.hand.soil, add = TRUE)
plot(child.direct.soil, add = TRUE)



#### Sensitivity Analysis comparing the p75/p50, p50/p25 and p75/p25  #####
#### Sensitivity Analysis comparing the p75/p50, p50/p25 and p75/p25  #####
age.groups <- c("f6", "6_12", "12_24", "24_36")

parameters1 <- c("child.hand.soil.concentration.mg.cm2", "child.hand.fraction.mouthed",
                 "mom.hand.soil.load.mg", "mom.hand.surface.area.cm2", "saliva.removal.efficiency",
                 "obj.soil.concentration", "obj.SA.mouthed", "soil.directly.ingested.mg")

parameters2 <- c("child.hand.surface.area.cm2", "child.hand.mouth.frequency.events.hr", 
                 "mom.hand.fraction.mouthed", "mom.hand.mouth.nonfood.frequency.events.hr", "mom.feeding.frequency.events.hr",
                 "child.obj.mouth.frequency.events.hr",
                 "child.frequency.ingested.soil.events.hr",  "hours.awake")

parameters <- c(parameters1, parameters2)

# Make sensitivity data.frame to hold results
sens.output <- data.frame(parameter = character(), age.group = character(), p25 = numeric(), p50 = numeric(), p75 = numeric(), p75.p50 = numeric(), p50.p25 = numeric(), p75.p25 = numeric())


for(age.group in age.groups){
  # Assign all parameters their distributions
  child.hand.soil.concentration.mg.cm2 <- soil.C.conc # should resample every hour = ndvar * hours
  child.hand.fraction.mouthed <- HF.child.mcstoc # should resample every hour = ndvar * hours
  mom.hand.soil.load.mg <- soil.M.load.mcstoc # should resample every hour = ndvar * hours
  mom.hand.surface.area.cm2 <- M.hand.SA.WASHB.mcstoc # should resample every hour = ndvar * hours
  saliva.removal.efficiency <- SEE.mcstoc 
  obj.soil.concentration <- object.soil.concentration.mcstoc # should resample every hour = ndvar * hours
  obj.SA.mouthed <- obj.SA.mouthed.mcstoc # should resample every hour = ndvar * hours
  soil.directly.ingested.mg <- SM.ingested.amt.mcstoc # should resample every hour = ndvar * hours
  
  child.hand.surface.area.cm2 <- get(paste("C.hand.SA.WASHB.", age.group, ".mcstoc", sep = ""))
  child.hand.mouth.frequency.events.day <- get(paste("HM.child.", age.group, ".mcstoc", sep = "")) # should resample every hour = ndvar * hours
  mom.hand.fraction.mouthed <- get(paste("HF.ofmom.", age.group, sep = "")) # should resample every hour = ndvar * hours
  mom.hand.mouth.nonfood.frequency.events.day <- get(paste("HM.other.nd.", age.group, ".mcstoc", sep = "")) # should resample every hour = ndvar * hours
  mom.feeding.frequency.events.day <- get(paste("HM.other.d.", age.group, ".mcstoc", sep = "")) # should resample every hour = ndvar * hours
  child.obj.mouth.frequency.events.day  <- get(paste("OM.", age.group, ".mcstoc", sep = "")) # should resample every hour = ndvar * hours
  child.frequency.ingested.soil.events.day <- get(paste("SM.fracfreq.", age.group, sep = "")) # Includes the fraction of children ingesting and the frequency of those ingesting, # should resample every hour = ndvar * hours
  #hours.awake <- get(paste("awake.day.", age.group, sep = ""))
  
  # Set all to p50
  for(parameter in parameters){
    assign(paste(parameter, ".p", sep = ""), as.numeric(unlist(summary(get(parameter)))[6])) #median(parameter, na.rm = TRUE))
  }
  
  # Calculate soil ingested based on median values for all parameters
  child.hand.soil.day.p50 <- child.hand.soil.concentration.mg.cm2.p * (child.hand.surface.area.cm2.p * child.hand.fraction.mouthed.p) * child.hand.mouth.frequency.events.day.p * saliva.removal.efficiency.p
  mom.hand.soil.day.p50 <- ((mom.hand.soil.load.mg.p / mom.hand.surface.area.cm2.p) * (mom.hand.surface.area.cm2.p * mom.hand.fraction.mouthed.p) * (mom.hand.mouth.nonfood.frequency.events.day.p + mom.feeding.frequency.events.day.p)) * saliva.removal.efficiency.p
  child.obj.soil.day.p50 <- (obj.soil.concentration.p * obj.SA.mouthed.p) * child.obj.mouth.frequency.events.day.p  * saliva.removal.efficiency.p
  child.direct.soil.day.p50 <- (soil.directly.ingested.mg.p * child.frequency.ingested.soil.events.day.p)
  amt.soil.ingested.mg.day.p50 <- (child.hand.soil.day.p50 + mom.hand.soil.day.p50 + child.obj.soil.day.p50 + child.direct.soil.day.p50)
  p50 <- amt.soil.ingested.mg.day.p50
  
  # Set ONE parameter to p25, calc amt ingested, then set to p75 and calc amt ingested
  for(parameter in parameters){
    #param.p50 <- as.numeric(unlist(summary(get(parameter)))[6])
    
    # 5th percentile
    param.p2.5 <- as.numeric(unlist(summary(get(parameter)))[4])
    assign(paste(parameter, ".p", sep = ""), param.p2.5)
    
    # Recalc outcomes will all at median except for the one paramter that was changed
    child.hand.soil.day.p2.5 <- child.hand.soil.concentration.mg.cm2.p * (child.hand.surface.area.cm2.p * child.hand.fraction.mouthed.p) * child.hand.mouth.frequency.events.day.p * saliva.removal.efficiency.p
    mom.hand.soil.day.p2.5 <- ((mom.hand.soil.load.mg.p / mom.hand.surface.area.cm2.p) * (mom.hand.surface.area.cm2.p * mom.hand.fraction.mouthed.p) * (mom.hand.mouth.nonfood.frequency.events.day.p + mom.feeding.frequency.events.day.p)) * saliva.removal.efficiency.p
    child.obj.soil.day.p2.5 <- (obj.soil.concentration.p * obj.SA.mouthed.p) * child.obj.mouth.frequency.events.day.p  * saliva.removal.efficiency.p
    child.direct.soil.day.p2.5 <- (soil.directly.ingested.mg.p * child.frequency.ingested.soil.events.day.p)
    amt.soil.ingested.mg.day.p2.5 <- (child.hand.soil.day.p2.5 + mom.hand.soil.day.p2.5 + child.obj.soil.day.p2.5 + child.direct.soil.day.p2.5) 
    p2.5 <- amt.soil.ingested.mg.day.p2.5
    
    # 25th percentile
    param.p25 <- as.numeric(unlist(summary(get(parameter)))[5])
    assign(paste(parameter, ".p", sep = ""), param.p25)
    
    # Recalc outcomes will all at median except for the one paramter that was changed
    child.hand.soil.day.p25 <- child.hand.soil.concentration.mg.cm2.p * (child.hand.surface.area.cm2.p * child.hand.fraction.mouthed.p) * child.hand.mouth.frequency.events.day.p * saliva.removal.efficiency.p
    mom.hand.soil.day.p25 <- ((mom.hand.soil.load.mg.p / mom.hand.surface.area.cm2.p) * (mom.hand.surface.area.cm2.p * mom.hand.fraction.mouthed.p) * (mom.hand.mouth.nonfood.frequency.events.day.p + mom.feeding.frequency.events.day.p)) * saliva.removal.efficiency.p
    child.obj.soil.day.p25 <- (obj.soil.concentration.p * obj.SA.mouthed.p) * child.obj.mouth.frequency.events.day.p  * saliva.removal.efficiency.p
    child.direct.soil.day.p25 <- (soil.directly.ingested.mg.p * child.frequency.ingested.soil.events.day.p)
    amt.soil.ingested.mg.day.p25 <- (child.hand.soil.day.p25 + mom.hand.soil.day.p25 + child.obj.soil.day.p25 + child.direct.soil.day.p25) 
    p25 <- amt.soil.ingested.mg.day.p25
    
    # Now set one value to the 75th percentile
    param.p75 <- as.numeric(unlist(summary(get(parameter)))[7])
    assign(paste(parameter, ".p", sep = ""), param.p75)
    
    # Recalc outcomes will all at median except for the one paramter that was changed
    child.hand.soil.day.p75 <- child.hand.soil.concentration.mg.cm2.p * (child.hand.surface.area.cm2.p * child.hand.fraction.mouthed.p) * child.hand.mouth.frequency.events.day.p * saliva.removal.efficiency.p
    mom.hand.soil.day.p75 <- ((mom.hand.soil.load.mg.p / mom.hand.surface.area.cm2.p) * (mom.hand.surface.area.cm2.p * mom.hand.fraction.mouthed.p) * (mom.hand.mouth.nonfood.frequency.events.day.p + mom.feeding.frequency.events.day.p)) * saliva.removal.efficiency.p
    child.obj.soil.day.p75 <- (obj.soil.concentration.p * obj.SA.mouthed.p) * child.obj.mouth.frequency.events.day.p  * saliva.removal.efficiency.p
    child.direct.soil.day.p75 <- (soil.directly.ingested.mg.p * child.frequency.ingested.soil.events.day.p)
    amt.soil.ingested.mg.day.p75 <- (child.hand.soil.day.p75 + mom.hand.soil.day.p75 + child.obj.soil.day.p75 + child.direct.soil.day.p75)
    p75 <- amt.soil.ingested.mg.day.p75
    
    # 97.5th percentile
    param.p97.5 <- as.numeric(unlist(summary(get(parameter)))[8])
    assign(paste(parameter, ".p", sep = ""), param.p97.5)
    
    # Recalc outcomes will all at median except for the one paramter that was changed
    child.hand.soil.day.p97.5 <- child.hand.soil.concentration.mg.cm2.p * (child.hand.surface.area.cm2.p * child.hand.fraction.mouthed.p) * child.hand.mouth.frequency.events.day.p * saliva.removal.efficiency.p
    mom.hand.soil.day.p97.5 <- ((mom.hand.soil.load.mg.p / mom.hand.surface.area.cm2.p) * (mom.hand.surface.area.cm2.p * mom.hand.fraction.mouthed.p) * (mom.hand.mouth.nonfood.frequency.events.day.p + mom.feeding.frequency.events.day.p)) * saliva.removal.efficiency.p
    child.obj.soil.day.p97.5 <- (obj.soil.concentration.p * obj.SA.mouthed.p) * child.obj.mouth.frequency.events.day.p  * saliva.removal.efficiency.p
    child.direct.soil.day.p97.5 <- (soil.directly.ingested.mg.p * child.frequency.ingested.soil.events.day.p)
    amt.soil.ingested.mg.day.p97.5 <- (child.hand.soil.day.p97.5 + mom.hand.soil.day.p97.5 + child.obj.soil.day.p97.5 + child.direct.soil.day.p97.5) 
    p97.5 <- amt.soil.ingested.mg.day.p97.5
    
    # Before the next parameter is tested, the parameter that was just tested has to be test back to it's median
    param.p50 <- as.numeric(unlist(summary(get(parameter)))[6])
    assign(paste(parameter, ".p", sep = ""), param.p50)
    
    # Calc the ratios
    p75.p50 <- p75/p50
    p50.p25 <- p50/p25
    p75.p25 <- p75/p25
    p97.5.p50 <- p97.5/p50
    p50.p2.5 <- p50/p2.5
    p97.5.p2.5 <- p97.5/p2.5
    
    # And add them to the dataframe
    new_row <- data.frame(parameter = parameter, age.group = age.group, 
                          param.p2.5 = param.p2.5, param.p25 = param.p25, param.p50 = param.p50, param.p75 = param.p75, param.p97.5 = param.p97.5,
                          p2.5 = p2.5, p25 = p25, p50 = p50, p75 = p75, p97.5 = p97.5,
                          p97.5.p50 = p97.5.p50, p75.p50 = p75.p50, p50.p25 = p50.p25, p50.p2.5 = p50.p2.5, p75.p25 = p75.p25, p97.5.p2.5 = p97.5.p2.5)
    sens.output <- rbind(sens.output, new_row)
    
    
  }
}

sens.output$age.group <- factor(sens.output$age.group, levels = c("", "f6", "6_12", "12_24", "24_36"),
                                labels = c("", "<6 months", "6 to <12 months", "12 to <24 months", "24 to <36 months"))

write.csv(sens.output, paste("C:/Users/Tareq/Box Sync/VO Soil/sensitivity.csv", sep = ""))




# # mcratio(mcnode)
# The mcratio function provides measures of variability, uncertainty, and both combined propose by [5] for an mc
# or an mcnode object - NOT FOR THE OUTCOME OF AMT - ONLY FOR THE SPECIFIC MCNODE. Given:
# A the median of uncertainty for the median of variability;
# B the median of uncertainty for the 97.5th percentile of variability;
# C the 97.5th percentile of uncertainty for the median percentile of variability;
# D the 97.5th percentile of uncertainty for the 97.5th percentile of variability.
# The following ratio are estimated:
#   - variability ratio B/A
#   - uncertainty ratio C/A
#   - overall uncertainty ratio D/A



# 
# # Processing 
# > tmp <- unmc(EC2, drop=TRUE) #unmc removes all attributes
# > dimu <- ncol(tmp$risk)
# > coef <- sapply(1:dimu, function(x) lm(tmp$risk[,x] ~ tmp$dose[,x])$coef)
# > apply(coef,1,summary)
# 
# # From mc2d documentation
# > library(fitdistrplus)
# > pzero <- sum(inca==0)/length(inca)
# > inca_non_0 <- inca[inca!=0]
# > descdist(inca_non_0)
# 
# > Adj_water <- fitdist(inca_non_0,"lnorm",method="mle")
# > meanlog <- Adj_water$est[1]
# > sdlog <- Adj_water$est[2]
# > summary(Adj_water)
# 
# # Consider uncertainty with bootstrapping
# > Boot <- bootdist(Adj_water, bootmethod="param", niter=ndunc())
# > Mean_conso <- mcdata(Boot$estim$meanlog, type="U")
# > Sd_conso <- mcdata(Boot$estim$sdlog, type="U")
# > conso1 <- mcstoc(rlnorm, type="VU", meanlog= Mean_conso, sdlog= Sd_conso)
# 
# #But for simplicity, we will not consider uncertainty around the estimates.
# #We will use the mcprobtree function to construct a mixture of 0 and non-0 distributions:
# > conso0 <- mcdata(0,type="V")
# > conso1 <- mcstoc(rlnorm, type="V", meanlog=meanlog, sdlog=sdlog)
# > v <- mcprobtree(c(pzero,1-pzero), list("0"=conso0,"1"=conso1), type = "V")
# > summary(v)



# Multiplot: par(mfrow=c(4,1))


# Diff between children <>6 months? --> NO
t.test(rlnorm(24,2.15,1.64), rlnorm(123,2.34, 1.58)) # p = 0.228

# Diff between children <>12 months? --> YES
t.test(c(rlnorm(24,2.15,1.64), rlnorm(123,2.34, 1.58)), c(rlnorm(54, 3.44, 1.77), rlnorm(5, 4.57, 1.34))) #p = 0.014

# Diff between children <>24 months? --> YES
t.test(c(rlnorm(24,2.15,1.64), rlnorm(123,2.34, 1.58), rlnorm(54, 3.44, 1.77)), c(rlnorm(5, 4.57, 1.34))) #p = 0.36




#### Why is the rate so high among children > 24 months? Due to high HM? #### What if we use the HM rate for younger kids?
### Children > 24 months old ###
child.hand.surface.area.cm2 <- C.hand.SA.WASHB.24_36.mcstoc
child.hand.mouth.frequency.events.day <- HM.child.24_36.mcstoc
mom.hand.mouth.nonfood.frequency.events.day <- HM.other.nd.24_36.mcstoc
mom.feeding.frequency.events.day <- HM.other.d.24_36.mcstoc
#mom.hand.fraction.mouthed <- HF.ofmom.24_36
#hours.awake <- awake.day.24_36

# child.hand.mouth.frequency.events.day <- HM.child.6_12.mcstoc
mom.hand.fraction.mouthed <- HF.ofmom.6_12


soil.test_24_36 <- mcmodel({ 
  amt.soil.ingested.mg.day <- ((((child.hand.soil.concentration.mg.cm2 * child.hand.surface.area.cm2) * child.hand.mouth.frequency.events.day * child.hand.fraction.mouthed) + 
                                  ((mom.hand.soil.load.mg / mom.hand.surface.area.cm2) * mom.hand.surface.area.cm2 * (mom.hand.mouth.nonfood.frequency.events.day + mom.feeding.frequency.events.day) * mom.hand.fraction.mouthed))  * saliva.removal.efficiency) 
  mc(child.hand.soil.concentration.mg.cm2, child.hand.surface.area.cm2, child.hand.mouth.frequency.events.day, child.hand.fraction.mouthed, mom.hand.soil.load.mg , mom.hand.surface.area.cm2,
     mom.hand.mouth.nonfood.frequency.events.day, mom.feeding.frequency.events.day, mom.hand.fraction.mouthed, saliva.removal.efficiency, hours.awake, amt.soil.ingested.mg.day) #sleep.day.24_36,
})
eatSoil_24_36 <- evalmcmod(soil.test_24_36, nsv=ndvar(), nsu=ndunc(), seed=seed)
summary(eatSoil_24_36)
tor.eatSoil_24_36  <- tornado(eatSoil_24_36, use = "pairwise.complete.obs") #  I'd prefer to use complete.obs rather than pairwise (I think) but this doesn't work. method = "spearman", lim = c(0.025, 0.975)
plot(tor.eatSoil_24_36 )
print(tor.eatSoil_24_36 )

# If we replace HM.child.24_36.mcstoc with HM.child.12_24.mcstoc
# amt.soil.ingested.mg.day :
#         mean  sd    Min  2.5%   25%  50%   75% 97.5%  Max  nsv Na's
# median 102.7 216 0.0893 0.846 11.37 37.4 105.6   587 2785 1001  426
# mean   103.9 238 0.1072 0.862 11.45 37.8 105.7   592 3428 1001  426
# 2.5%    86.6 150 0.0216 0.610  9.03 31.5  90.3   456 1294 1001  392
# 97.5%  125.6 394 0.2523 1.252 13.80 44.4 119.0   761 7315 1001  45

# If we replace HM.child.24_36.mcstoc with HM.child.6_12.mcstoc
# amt.soil.ingested.mg.day :
#        mean    sd    Min  2.5%  25%  50%  75% 97.5%  Max  nsv Na's
# median 47.0 113.7 0.0691 0.495 4.22 13.4 40.5   300 1533 1001  426
# mean   47.4 123.9 0.0795 0.504 4.23 13.2 41.1   307 1812 1001  426
# 2.5%   36.7  74.5 0.0168 0.358 3.49 11.0 33.5   205  694 1001  392
# 97.5%  59.4 248.9 0.1789 0.721 4.97 15.1 49.6   423 5268 1001  452


# What's the influence of the fraction of mom's hand mouthed?
HF.ofmom.f6
HF.ofmom.6_12
HF.ofmom.12_24
HF.ofmom.24_36

# Using HF.ofmom.24_36
# amt.soil.ingested.mg.day :
#   mean  sd   Min 2.5%  25%   50% 75% 97.5%   Max  nsv Na's
# median  186 307 1.059 6.28 38.2  92.0 206   912  3793 1001  426
# mean    186 324 1.156 6.26 38.2  91.1 207   915  4162 1001  426
# 2.5%    162 220 0.231 4.41 32.0  78.6 182   710  1651 1001  392
# 97.5%   213 561 2.639 8.74 43.4 102.0 233  1139 10373 1001  452

# If we replace HF.ofmom.24_36 with HF.ofmom.12_24
# amt.soil.ingested.mg.day :
#         mean  sd   Min 2.5%  25%   50% 75% 97.5%   Max  nsv Na's
# median  186 307 1.054 6.25 38.1  91.2 206   908  3793 1001  428
# mean    186 323 1.155 6.22 38.1  90.9 207   913  4161 1001  427
# 2.5%    162 220 0.224 4.38 31.6  78.8 183   711  1650 1001  392
# 97.5%   213 561 2.646 8.65 43.3 102.0 233  1140 10373 1001  453

# If we replace HF.ofmom.24_36 with HF.ofmom.6_12
# amt.soil.ingested.mg.day :
#         mean  sd   Min 2.5%  25%   50% 75% 97.5%   Max  nsv Na's
# median  185 307 1.050 6.14 37.7  91.0 205   906  3792 1001  428
# mean    185 323 1.132 6.13 37.8  90.3 206   912  4159 1001  428
# 2.5%    162 220 0.223 4.45 31.4  78.2 181   710  1646 1001  394
# 97.5%   213 561 2.587 8.66 43.1 100.6 232  1131 10372 1001  454



summary(anth.soilMass.wide$age.soil.child.7, na.rm = TRUE)

names(anth.soil.wide)
summary(anth.soilMass.wide$survey.age.mo.mid) # median 8.8
sd(anth.soilMass.wide$survey.age.mo.mid, na.rm = TRUE) # sd = 1.72
summary(anth.soilMass.wide$survey.age.mo.end) # median 22.5
sd(anth.soilMass.wide$survey.age.mo.end) # sd = 2.09 (how did the sd get bigger if the same children were invovled?.....)







# ##### Incorp the bootstrap #########
# ## https://cran.r-project.org/web/packages/mc2d/mc2d.pdf ## pg 36
# 
# ##Use of rempiricalD with nodes
# ##A bootstrap
# ndunc(5)
# ndvar(5)
# dataset <- c(1:9)
# (b <- mcstoc(rempiricalD, "U", nvariates=9, values=dataset))
# unclass(b)
# ##Then we build a VU node by sampling in each set of bootstrap
# (node <- mcstoc(rempiricalD, "VU", values=b))
# unclass(node)







