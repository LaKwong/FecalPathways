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
library(lubridate)
library(plyr)
library(dplyr)
#library(compositions)

ndvar <- 101 #1001, 10001
ndunc <- 10
seed <- 888

# Load soil mass data

soilMass.base <- read_excel("C:/Users/Laura Kwong/Box Sync/VO Soil/Mass of soil on children's hands/Soil mass_RO1_analysis.xlsx", sheet = "final")
soilMass.base <- soilMass.base[1:272,c("hh", "RO1.round", "survey.day", "survey.month", "rep.dup", "motherOrChild", "volTotal.ml", "volSample.ml", "filterMass.g", "filterSoilMass.g")]
soilMass.base$survey.yr <- 2016
soilMass.base$surveyDate.soil <- dmy(paste(as.numeric(soilMass.base$survey.day), as.numeric(soilMass.base$survey.month), soilMass.base$survey.yr, sep = "-"))
soilMass.base$soil.g <- soilMass.base$filterSoilMass.g - soilMass.base$filterMass.g
soilMass.base$soil.mg <- (soilMass.base$filterSoilMass.g - soilMass.base$filterMass.g)*1000
soilMass.base$soil.mg.halfLOD <- NA

#### Think about the detection limit! The precision of the scale is 5 mg = 0.005 g, so anything less than 0.005 g should be treated as...
## THIS ASSUMPTION WILL MAKE A BIG DIFFERENCE
## How many samples are below the detection limit of 0.005 g?
qplot(soilMass.base$soil.mg, geom="histogram", binwidth = 1)

sum(soilMass.base$soil.mg <= 5)/length(soilMass.base$soil.mg) # 57.0% of observations (including the reps) are less than the LOD

## Can try 1/2 the detection limit, though I remember hearing that there are better ways to handle this. 



##################################################################
###################################################################
## Conduct sensitivity analysis to det effect of choosing what happens to values < LOD
#belowLOD_replacementValue <- 2.5
#soilMass.base$soil.mg.halfLOD <- ifelse(soilMass.base$soil.mg <= 5, belowLOD_replacementValue, soilMass.base$soil.mg)

#belowLOD_replacementValue <- runif(1, 0, 5)
soilMass.base[soilMass.base$soil.mg <= 5, c("soil.mg.halfLOD")] <- runif(nrow(soilMass.base[soilMass.base$soil.mg <= 5, c("soil.mg.halfLOD")]), 0, 5)
soilMass.base[soilMass.base$soil.mg > 5, c("soil.mg.halfLOD")] <- soilMass.base[soilMass.base$soil.mg > 5, c("soil.mg.halfLOD")] 
soilMass.base$soil.mg.halfLOD
qplot(soilMass.base$soil.mg.halfLOD, geom="histogram", binwidth = 0.1)

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

#anth.soil.wide <- read.csv("C:/Users/Laura Kwong/Box Sync/VO Fecal Intake Model/anth.soil.wide.csv")
anth.soil.wide <- read.csv("C:/Users/Laura Kwong/Box Sync/VO Fecal Intake Model/anth.soil.wide_const_HandSAratio.csv")
anth.soil.wide$C.body.SA.end
anth.soil.wide <- anth.soil.wide %>%
  filter(arm == "Control")

anth.soilMass.wide <- left_join(anth.soil.wide, soilMass.wide, by = "hh")
# The dates had the correct format when saved but the format was not imported so need to change from factor to date
anth.soilMass.wide$unique_dob <- ymd(anth.soilMass.wide$unique_dob)
anth.soilMass.wide$surveyDate.mid <- ymd(anth.soilMass.wide$surveyDate.mid)
anth.soilMass.wide$surveyDate.end <- ymd(anth.soilMass.wide$surveyDate.end)

anth.soilMass.wide$age.soil.child.7 <- (anth.soilMass.wide$surveyDate.soil.child.7 - anth.soilMass.wide$unique_dob) / 365 * 12
anth.soilMass.wide$age.soil.child.8 <- (anth.soilMass.wide$surveyDate.soil.child.8 - anth.soilMass.wide$unique_dob) / 365 * 12
anth.soilMass.wide$round7DaysFromEnd <- anth.soilMass.wide$surveyDate.soil.child.7 - anth.soilMass.wide$surveyDate.end # at the time the soil round 7 sample was collected, children were 360 year older than the were at WASHB endline 
anth.soilMass.wide$round8DaysFromEnd <- anth.soilMass.wide$surveyDate.soil.child.8 - anth.soilMass.wide$surveyDate.end# at the time the soil round 8 sample was collected, children were 460 year older than the were at WASHB endline 


anth.soilMass.wide %>%
  select(age.soil.child.7) %>%
  filter(!is.na(age.soil.child.7)) / 365*12 # ~ most kids about 32-38 months old when sampled

anth.soilMass.wide %>%
  select(round8DaysFromEnd) %>%
  filter(!is.na(round8DaysFromEnd)) / 365*12
######################################################
######################################################
########## Bummer, the kids whos hand we sampled are WAY older than the kids in WASHB, so can't use the WASHB anthro data to est g/cm2
### How else can we get the SA of the hands - I think Ayse did not collect anthropometry data

### Maybe DHS data? - nope, birthweights only --> Ask Steve, kishor, mahbub?

# ######### For now, let's use the surface area of the oldest children in the WASHB dataset
# # Hand SA distribution for WASHB children OVER 30 months at endline
# ## There are only 7 Control arm kids > 30 months at WASHB Endline
# anth.C.over30.hand.SA.WASHB <- anth.soilMass.wide %>%
#   filter(survey.age.mo.end > 30) %>%
#   select(C.hand.SA.end)
# qplot(anth.C.over30.hand.SA.WASHB$C.hand.SA.end, geom = "histogram", binwidth = 1)
# # # Dist doesn't fit so use empirical
# # anth.C.over30.hand.SA.WASHB.dist <- fitdist(anth.C.over30.hand.SA.WASHB$C.hand.SA.end, "norm", method = "mle")
# # C.over30.hand.SA.WASHB.mcstoc <- mcstoc(rweibull, type="V", shape = anth.C.over30.hand.SA.WASHB.dist$estimate[[1]], scale = anth.C.over30.hand.SA.WASHB.dist$estimate[[2]], rtrunc=TRUE, linf=0)
# 
# anth.C.over30.hand.SA.WASHB.mcdata <- mcdata(sample(anth.C.over30.hand.SA.WASHB$C.hand.SA.end, size=ndvar(), replace=TRUE),type="V")


########## There's a pretty linear trend between hand size and age so try modling hand surface area based on age - just linea 
age.SA.Round7 <- anth.soilMass.wide[,c("survey.age.mo.mid", "C.hand.SA.mid")]
names(age.SA.Round7) <- c("age.mo", "C.hand.SA")
age.SA.Round8 <- anth.soilMass.wide[,c("survey.age.mo.end", "C.hand.SA.end")]
names(age.SA.Round8) <- c("age.mo", "C.hand.SA")

C.hand.SA.RO1 <- data.frame(rbind(age.SA.Round7, age.SA.Round8))

scatter.smooth(C.hand.SA.RO1$age.mo, C.hand.SA.RO1$C.hand.SA)
C.hand.SA.RO1.lm <- lm(C.hand.SA ~ age.mo, C.hand.SA.RO1)
C.hand.SA.RO1.lm$coefficients[1]
C.hand.SA.RO1.lm$coefficients[2]
# Call:
#   lm(formula = C.hand.SA ~ age.mo, data = C.hand.SA.RO1)
# 
# Coefficients:
#   (Intercept)       age.mo  
#      88.090        1.776 


anth.soilMass.wide$child.hand.SA.age.soil.7 <- as.numeric(C.hand.SA.RO1.lm$coefficients[1] + anth.soilMass.wide$age.soil.child.7 * C.hand.SA.RO1.lm$coefficients[2])
anth.soilMass.wide$child.hand.SA.age.soil.8 <- as.numeric(C.hand.SA.RO1.lm$coefficients[1] + anth.soilMass.wide$age.soil.child.8 * C.hand.SA.RO1.lm$coefficients[2])
C.hand.SA <- data.frame(c(anth.soilMass.wide[!is.na(anth.soilMass.wide$child.hand.SA.age.soil.7), "child.hand.SA.age.soil.7"], 
                                anth.soilMass.wide[!is.na(anth.soilMass.wide$child.hand.SA.age.soil.8), "child.hand.SA.age.soil.8"]))
names(C.hand.SA) <- c("C.hand.SA")
qplot(C.hand.SA$C.hand.SA, geom = "histogram", binwidth = 1)
# These hand surface areas are pretty close to what I'd expect! ~150

# Hand SA distribution for WASHB mothers at endline
anth.M.hand.SA.WASHB.end <- anth.soilMass.wide %>%
  select(M.hand.SA.end) %>%
  filter(M.hand.SA.end < 1000) ## Filter out the unreasonable hand SA values - there are not outliers but rather mistakes with height or weight measurement
qplot(anth.M.hand.SA.WASHB.end$M.hand.SA.end, geom = "histogram", binwidth = 1)
# In Yu (2009), XS slim women had a Measured Body SA of 15,220.56 cm2, assume ONE HAND is 2.26% = 343.9847 cm2 - my values are MUCH higher. Why?

# Dist DOES fit so use empirical
anth.M.hand.SA.WASHB.end.dist <- fitdist(anth.M.hand.SA.WASHB.end$M.hand.SA.end, "norm", method = "mle")
M.hand.SA.WASHB.end.mcstoc <- mcstoc(rnorm, type="V", mean = anth.M.hand.SA.WASHB.end.dist$estimate[[1]], sd = anth.M.hand.SA.WASHB.end.dist$estimate[[2]], rtrunc=TRUE, linf=0)
#### More parsimonious to use values from midline and endline combined # M.hand.SA.WASHB.mcstoc
# Anyway values are very close:
# > M.hand.SA.WASHB.end.mcstoc
# node    mode  nsv nsu nva variate min mean median max Nas type outm
# 1    x numeric 1001   1   1       1 234  322    322 404   0    V each
# > M.hand.SA.WASHB.mcstoc
# node    mode  nsv nsu nva variate min mean median max Nas type outm
# 1    x numeric 1001   1   1       1 232  321    321 388   0    V each
#

## This line also appears in line 320 below:
## Determine the all-ages hand surface area distribution for mothers
M.hand.SA.WASHB <- anth.soilMass.wide %>%
  select(hh, M.hand.SA.mid, M.hand.SA.end) %>%
  gather(key = round, value = M.hand.SA, M.hand.SA.mid, M.hand.SA.end) %>%
  select(hh, M.hand.SA) %>%
  filter(M.hand.SA < 1000 & !is.na(M.hand.SA))

qplot(M.hand.SA.WASHB$M.hand.SA, geom = "histogram", binwidth = 1)

# Dist DOES fit so use empirical
M.hand.SA.WASHB.dist <- fitdist(M.hand.SA.WASHB$M.hand.SA, "norm", method = "mle")

# Fitting of the distribution ' norm ' by maximum likelihood 
# Parameters:
#       estimate Std. Error
# mean 320.9020  0.5920519
# sd    26.2847  0.4186439
M.hand.SA.WASHB.mcstoc <- mcstoc(rnorm, type="V", mean = M.hand.SA.WASHB.dist$estimate[[1]], sd = M.hand.SA.WASHB.dist$estimate[[2]], rtrunc=TRUE, linf=0)
#M.hand.SA.WASHB.mcdata <- mcdata(sample(M.hand.SA.WASHB$M.hand.SA, size=ndvar(), replace=TRUE),type="V")



############## Est soil concentration on hands based on measured soil load for RO1 children and mothers and WASHB data for hand surface area ################
# soil.C.conc <- soil.C.load.mcdata / anth.C.over30.hand.SA.WASHB.mcdata
# Use the actual load divided by the est hand.SA based on age
anth.soilMass.wide[anth.soilMass.wide$hh %in% soilMass.wide$hh, c("child.hand.SA.age.soil.7", "soilOneHand.mg.child.7", "child.hand.SA.age.soil.8", "soilOneHand.mg.child.8")]
soil.C.conc.wide <- anth.soilMass.wide %>%
  filter(hh %in% soilMass.wide$hh) %>%
  mutate(soil.C.conc.7 = soilOneHand.mg.child.7 / child.hand.SA.age.soil.7, soil.C.conc.8 = soilOneHand.mg.child.8 / child.hand.SA.age.soil.8) %>%
  select(soil.C.conc.7, soil.C.conc.8)
soil.C.conc.long <- data.frame(combine(soil.C.conc.wide$soil.C.conc.7, soil.C.conc.wide$soil.C.conc.8))
names(soil.C.conc.long) <- "soil.C.conc"

##################################################################################################
##########################################################################################
################## should soil.C.conc be mcdata or a dist like soil.M.conc?
######### Should soil.C.conc or soil.M.conc be a logrithmic dist? 
qplot(soil.C.conc.long, geom = "histogram", binwidth = 0.01)

soil.C.conc <- mcdata(sample(soil.C.conc.long$soil.C.conc, size=ndvar(), replace=TRUE),type="V")
soil.M.conc <- soil.M.load.mcdata / M.hand.SA.WASHB.mcstoc #M.hand.SA.WASHB.end.mcstoc 
# use M.hand.SA.WASHB.mcstoc instead of M.hand.SA.WASHB.end.mcstoc to for simplicity of the methodology. Also .end has mean SA of 323 and the pooled sample has mean SA of 320

# Ozkaynak Soil adherence to hands mg/cm2 : mcstoc(rlnorm, type = "V", meanlog = 0.11, sdlog = 2.0)
# Produces values WAY higher than my estimates
# My est
# > soil.C.conc
# node    mode  nsv nsu nva variate     min   mean median  max Nas type outm
# 1    x numeric 1001   1   1       1 0.00935 0.0808 0.0433 0.86  37    V each

# Ozkaynak est
# > mcstoc(rlnorm, type = "V", meanlog = 0.11, sdlog = 2.0)
# node    mode  nsv nsu nva variate     min mean median max Nas type outm
# 1    x numeric 1001   1   1       1 0.00151 7.03   1.02 529   0    V each


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

scatter.smooth(C.hand.SA.WASHB$age.mo, C.hand.SA.WASHB$C.hand.SA)
C.hand.SA.WASHB.lm <- lm(C.hand.SA ~ age.mo, C.hand.SA.WASHB)
C.hand.SA.WASHB.lm$coefficients[1]
C.hand.SA.WASHB.lm$coefficients[2]

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

C.hand.SA.WASHB.f6 <- C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo < 6, "C.hand.SA"] # This variable isn't actually used
C.hand.SA.WASHB.f6.mcdata <- mcdata(sample(C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo < 6, "C.hand.SA"], size = ndvar(), replace = TRUE), type = "V")

C.hand.SA.WASHB.6_12 <- C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo >= 6 & C.hand.SA.WASHB$age.mo < 12, "C.hand.SA"] # This variable isn't actually used
C.hand.SA.WASHB.6_12.mcdata <- mcdata(sample(C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo < 6, "C.hand.SA"], size = ndvar(), replace = TRUE), type = "V")

C.hand.SA.WASHB.12_24 <- C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo >= 12 & C.hand.SA.WASHB$age.mo < 24, "C.hand.SA"] # This variable isn't actually used
C.hand.SA.WASHB.12_24.mcdata <- mcdata(sample(C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo < 6, "C.hand.SA"], size = ndvar(), replace = TRUE), type = "V")

C.hand.SA.WASHB.o24 <- C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo >= 24, "C.hand.SA"] # This variable isn't actually used
C.hand.SA.WASHB.o24.mcdata <- mcdata(sample(C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo > 24, "C.hand.SA"], size = ndvar(), replace = TRUE), type = "V")

qplot(C.hand.SA.WASHB.o24, geom = "histogram", binwidth = 1)

for(i in c(6:23)){ #max(HM.SM$age)
  assign(paste("C.hand.SA.WASHB.", i, sep = ""), C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo == i, "C.hand.SA"])
  assign(paste("C.hand.SA.WASHB.", i, ".mcdata", sep = ""), mcdata(sample(C.hand.SA.WASHB[C.hand.SA.WASHB$age.mo == i, "C.hand.SA"], size = ndvar(), replace = TRUE), type = "V"))
}

# C.hand.SA.WASHB.9.mcdata
# node    mode  nsv nsu nva variate  min mean median max Nas type outm
# 1    x numeric 1001   1   1       1 71.1 89.1   89.4 187   5    V each

# # Exposure Factors Handbook, 2011 for 6-12 mo old
# # 6 to < 12 mo	0.51	5.3	0.02703	270.3	0.024	0.027	120	135 - ONE HAND mean 120, 95th pc 135 - my children are WAY smaller 
# My est for a 6 mo old are about the same as for a US 1-3 mo old. 
# --> Check formula


## Determine the all-ages hand surface area distribution for mothers
M.hand.SA.WASHB <- anth.soilMass.wide %>%
  select(hh, M.hand.SA.mid, M.hand.SA.end) %>%
  gather(key = round, value = M.hand.SA, M.hand.SA.mid, M.hand.SA.end) %>%
  select(hh, M.hand.SA) %>%
  filter(M.hand.SA < 1000 & !is.na(M.hand.SA))

qplot(M.hand.SA.WASHB$M.hand.SA, geom = "histogram", binwidth = 1)

# Dist DOES fit so use empirical
M.hand.SA.WASHB.dist <- fitdist(M.hand.SA.WASHB$M.hand.SA, "norm", method = "mle")
M.hand.SA.WASHB.mcstoc <- mcstoc(rnorm, type="V", mean = M.hand.SA.WASHB.dist$estimate[[1]], sd = M.hand.SA.WASHB.dist$estimate[[2]], rtrunc=TRUE, linf=0)
#M.hand.SA.WASHB.mcdata <- mcdata(sample(M.hand.SA.WASHB$M.hand.SA, size=ndvar(), replace=TRUE),type="V")

################## Calc distribution of loading (mg) ###########
##### Est load on hands of all WASHB children using the est soil concentration on RO1 hands and the dist of hand SA of the WASHB hands
# Probablistic multiply by conc (g/cm2) by SA (cm2) of observed children to get est load g soil 

#Don't use the age month-specfic code for children below
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

# This file is pre-analysis in the R code for SO analysis that results in SO.HM.SM
#dem_obs_full <- read.csv("C:/Users/Laura Kwong/Box Sync/VO Longitudinal/count_mat_food_dem_fromSO_submitted to IJERPH 160313.csv")
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

## Be CAREFUL about FREQUENCIES vs COUNTS


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


# For all the kids that only have total HM values from the SO, need to est child, mom_d, mom_nd
# Est using the VO results for children 6-12 mo
# For the other children, use the actual child, mom_d, and mom_nd recorded in the VO


## For kids who only have SO values for total HM

hands_summary <- read.csv("C:/Users/Laura Kwong/Box Sync/VO R123/hands_summary.csv")
hands_summary$hh <- as.factor(hands_summary$hh)

### Must change to reflect that most mom HM contacts do not have recontam so really only want number of indep feeding events
feeding.mom.freq <- read.csv("C:/Users/Laura Kwong/Box Sync/VO R123/feeding.mom.freq.csv")
feeding.mom.freq <- feeding.mom.freq %>%
  select(-X, -X.1)
feeding.mom.freq$hh <- as.factor(feeding.mom.freq$hh)

HM.proportions <- left_join(hands_summary, feeding.mom.freq, by = c("vo.num", "hh", "age.vo")) %>%
  select(vo.num, hh, age.vo, age.vo.group, pf_child_hands_nd, pf_child_hands_d, pf_other_hands_nd, pc_other_hands_d, feeding.mom.freq)

scatter.smooth(feeding.mom.freq$age.vo, feeding.mom.freq$feeding.mom.freq)
# There is no trend by age, so use feeding freq from any age for a child of any age

# Make an array of vectors, where each vector is the child_nd, child_d, mom_nd, feeding.freq for each child in the age group
# Then create a random sample of the indicies
# Select the random index and get the child_nd, child_d, mom_nd, or feeding.freq col and save as mcdata

HM.proportions.4 <- HM.proportions %>% filter(age.vo.group == 4) %>% select(pf_child_hands_nd, pf_child_hands_d, pf_other_hands_nd, feeding.mom.freq) 
HM.proportions.index.4 <- sample(c(1:nrow(HM.proportions.4)), size = ndvar(), replace = TRUE)
HM.child.4.mcdata <- mcdata(rowSums(HM.proportions.4[HM.proportions.index.4, c(1,2)]), type = "V")
HM.child.nd.4.mcdata <- mcdata(HM.proportions.4[HM.proportions.index.4, 1], type = "V")
HM.child.nd.4.mcdata <- mcdata(HM.proportions.4[HM.proportions.index.4, 1], type = "V")
HM.child.d.4.mcdata <- mcdata(HM.proportions.4[HM.proportions.index.4, 2], type = "V")
HM.mom.nd.4.mcdata <- mcdata(HM.proportions.4[HM.proportions.index.4, 3], type = "V")
HM.mom.d.events.4.mcdata <- mcdata(HM.proportions.4[HM.proportions.index.4, 4], type = "V")

HM.proportions.5 <- HM.proportions %>% filter(age.vo.group == 5) %>% select(pf_child_hands_nd, pf_child_hands_d, pf_other_hands_nd, feeding.mom.freq) 
HM.proportions.index.5 <- sample(c(1:nrow(HM.proportions.5)), size = ndvar(), replace = TRUE)
HM.child.5.mcdata <- mcdata(rowSums(HM.proportions.5[HM.proportions.index.5, c(1,2)]), type = "V")
HM.child.nd.5.mcdata <- mcdata(HM.proportions.5[HM.proportions.index.5, 1], type = "V")
HM.child.d.5.mcdata <- mcdata(HM.proportions.5[HM.proportions.index.5, 2], type = "V")
HM.mom.nd.5.mcdata <- mcdata(HM.proportions.5[HM.proportions.index.5, 3], type = "V")
HM.mom.d.events.5.mcdata <- mcdata(HM.proportions.5[HM.proportions.index.5, 4], type = "V")

HM.proportions.6 <- HM.proportions %>% filter(age.vo.group == 6) %>% select(pf_child_hands_nd, pf_child_hands_d, pf_other_hands_nd, feeding.mom.freq) 
HM.proportions.index.6 <- sample(c(1:nrow(HM.proportions.6)), size = ndvar(), replace = TRUE)
HM.child.6.mcdata <- mcdata(rowSums(HM.proportions.6[HM.proportions.index.6, c(1,2)]), type = "V")
HM.child.nd.6.mcdata <- mcdata(HM.proportions.6[HM.proportions.index.6, 1], type = "V")
HM.child.d.6.mcdata <- mcdata(HM.proportions.6[HM.proportions.index.6, 2], type = "V")
HM.mom.nd.6.mcdata <- mcdata(HM.proportions.6[HM.proportions.index.6, 3], type = "V")
HM.mom.d.events.6.mcdata <- mcdata(HM.proportions.6[HM.proportions.index.6, 4], type = "V")



#######################################################################################
#######################################################################################
############## Use the data on total HM contacts from SO and VO to create distributions and use this distributions to create mcstoc
## For all HM_child
HM.SM.f6.dist <- fitdist(HM.SM[HM.SM$age < 6, "HM"], "weibull", method = "mle")
HM.SM.f6.mcstoc <- mcstoc(rweibull, type="V", shape = HM.SM.f6.dist$estimate[[1]], scale = HM.SM.f6.dist$estimate[[2]], rtrunc=TRUE, linf=0)

HM.SM.6_12.dist <- fitdist(HM.SM[HM.SM$age >= 6 & HM.SM$age < 12, "HM"], "weibull", method = "mle")
HM.SM.6_12.mcstoc <- mcstoc(rweibull, type="V", shape = HM.SM.f6.dist$estimate[[1]], scale = HM.SM.6_12.dist$estimate[[2]], rtrunc=TRUE, linf=0)

HM.SM.12_24.dist <- fitdist(HM.SM[HM.SM$age >=12 & HM.SM$age < 24, "HM"], "weibull", method = "mle")
HM.SM.12_24.mcstoc <- mcstoc(rweibull, type="V", shape = HM.SM.f6.dist$estimate[[1]], scale = HM.SM.12_24.dist$estimate[[2]], rtrunc=TRUE, linf=0)

HM.SM.o24.dist <- fitdist(HM.SM[HM.SM$age >= 24, "HM"], "weibull", method = "mle")
HM.SM.o24.mcstoc <- mcstoc(rweibull, type="V", shape = HM.SM.o24.dist$estimate[[1]], scale = HM.SM.o24.dist$estimate[[2]], rtrunc=TRUE, linf=0)

# # For the month-specific distributions, can I make weibull for each month?
#find.dist <- function(data, x, dist, who) {fitdist(data[data$age == x, who], dist, method = "mle")} #who is "HM_child" or "HM_mom"
find.dist <- function(data, x, dist) {fitdist(data[data$age == x, "HM"], dist, method = "mle")} #who is "HM_child" or "HM_mom"

for(i in 6:23){ #max(HM.SM$age)
  assign(paste("HM.SM.", i, ".dist", sep = ""), find.dist(HM.SM, i, "weibull"))
  assign(paste("HM.SM.", i, ".mcdata", sep = ""), mcdata(sample(HM.SM[HM.SM$age == i, "HM"], size = ndvar(), replace = TRUE), type = "V"))
}

####################################################################################################
############### Apportion HM to child (sum child_nd and child_d) and mom_d and feeding.freq
# age.group = 4 is 12-24 mo, age.group = 5 is 12-24 mo, age.group = 6 is 24-36 mo

# the percentages in hand_summary were in the form of XX% for the purpose of graphing boxplots, so need to divide by 100 to get to a percent between 0-1. 
HM.child.f6.mcdata <- HM.SM.f6.mcstoc * (HM.child.4.mcdata/100)
HM.other.nd.f6.mcdata <- HM.SM.f6.mcstoc * (HM.mom.nd.4.mcdata/100)
HM.other.d.f6.mcdata <- HM.mom.d.events.4.mcdata

HM.child.6_12.mcdata <- HM.SM.6_12.mcstoc * (HM.child.4.mcdata/100)
HM.other.nd.6_12.mcdata <- HM.SM.6_12.mcstoc * (HM.mom.nd.4.mcdata/100)
HM.other.d.6_12.mcdata <- HM.mom.d.events.4.mcdata

HM.child.12_24.mcdata <- HM.SM.12_24.mcstoc * (HM.child.5.mcdata/100)
HM.other.nd.12_24.mcdata <- HM.SM.12_24.mcstoc * (HM.mom.nd.5.mcdata/100)
HM.other.d.12_24.mcdata <- HM.mom.d.events.5.mcdata


HM.child.o24.mcdata <- HM.SM.o24.mcstoc * (HM.child.6.mcdata/100)
HM.other.nd.o24.mcdata <- HM.SM.o24.mcstoc * (HM.mom.nd.6.mcdata/100)
HM.other.d.o24.mcdata <- HM.mom.d.events.6.mcdata
# 
# HM.SM.o24.child.mcstoc
# HM.SM.f6.mom.mcstoc
# HM.SM.o24.mom.mcstoc
# 
# HM.child.nd.6.mcdata
# HM.child.d.6.mcdata
# HM.mom.nd.6.mcdata
# HM.mom.d.events.6.mcdata


#################### Combine load_child/load_mom with HM_child/HM_mom with mouth SA dist, removal efficiency dist

#HF	Child own hand fracion mouthed/event	-	day	beta	3.7	25				Zartarian 2005 as presented in Ozkaynak 2011	Zartarian 2005
HF.child.mcstoc <- mcstoc(rbeta, type="V", shape1 = 3.7, shape2 = 25, rtrunc=TRUE, linf=0)
##########################################################
# I can do better than this!! Work to replace this value, which is based only on 20 children. 


























# Children will mouth a smaller portion of mom's hand than of own hand 
# Children are "mouth-limited" in how much surface area they can put in their mouth. 
# Children will be able to put more surface area of their own hands in their mouth 
# than the surface area of their mother's hands because they have more surface area per unit volume (smaller finger = can put more finger = more SA in mouth)

for(i in c(4:40)){ #max(HM.SM$age)
  assign(paste("C.hand.SA.WASHB.inMouth.", i, ".mcdata", sep = ""), get(paste("C.hand.SA.WASHB.", i, ".mcdata", sep = "")) * HF.child.mcstoc)
  #assign(paste("C.hand.SA.WASHB.inMouth.", i, ".mcdata", sep = ""), mcdata(sample(get(paste("C.hand.SA.WASHB.inMouth.", i, sep = "")) , size = ndvar(), replace = TRUE), type = "V"))
}

C.hand.SA.WASHB.inMouth.6.mcdata
# HF.ofmom.xx <- C.hand.SA.WASHB.inMouth.xx.mcdata / median(M.hand.SA.WASHB.mcstoc)
HF.ofmom.6 <- 10.2/ median(M.hand.SA.WASHB.mcstoc)
HF.ofmom.9 <- 11.1/ median(M.hand.SA.WASHB.mcstoc) # use for age.group = 4 6-12 mo
HF.ofmom.12 <- 11.5/ median(M.hand.SA.WASHB.mcstoc)
HF.ofmom.18 <- mean(11.5, 14.1)/median(M.hand.SA.WASHB.mcstoc) # use for age.group = 5 12-24 mo
HF.ofmom.24 <- 14.1/ median(M.hand.SA.WASHB.mcstoc)
HF.ofmom.30 <- 16.1/ median(M.hand.SA.WASHB.mcstoc) # use for age.grop = 6 24-36 mo
#C.hand.SA.WASHB.inMouth.30.mcdata

# ## Could chane this so XX is by age group, not by month....This is better because below I change it into ae groups anyway...
# C.hand.SA.WASHB.inMouth.f6.mcdata <- C.hand.SA.WASHB.f6.mcdata * HF.child.mcstoc # mean 12.2, med 11.4
# C.hand.SA.WASHB.inMouth.6_12.mcdata <- C.hand.SA.WASHB.6_12.mcdata * HF.child.mcstoc # mean 12.2, med 11.3
# C.hand.SA.WASHB.inMouth.12_24.mcdata <- C.hand.SA.WASHB.12_24.mcdata * HF.child.mcstoc # mean 12.2, med 11.4
# C.hand.SA.WASHB.inMouth.o24.mcdata <- C.hand.SA.WASHB.o24.mcdata * HF.child.mcstoc # mean 18.4, med 16.8
# 
# # HF.ofmom.xx <- C.hand.SA.WASHB.inMouth.xx.mcdata / median(M.hand.SA.WASHB.mcstoc)
# HF.ofmom.f6 <- 12.2/ median(M.hand.SA.WASHB.mcstoc)
# HF.ofmom.6_12 <- 12.2/ median(M.hand.SA.WASHB.mcstoc) # use for age.group = 4 6-12 mo
# HF.ofmom.12_24 <- 12.2/ median(M.hand.SA.WASHB.mcstoc) # use for age.group = 4 6-12 mo
# HF.ofmom.o24 <- 18.4/ median(M.hand.SA.WASHB.mcstoc) # use for age.group = 4 6-12 mo


# HF.mom.mcstoc  <- mcstoc(runif, type="V", min = HF.ofmom.6, max = HF.ofmom.30, rtrunc=TRUE, linf=0) 
HF.mom.mcstoc.f6  <- mcstoc(runif, type="V", min = HF.ofmom.6, max = HF.ofmom.6, rtrunc=TRUE, linf=0) 
HF.mom.mcstoc.6_12  <- mcstoc(runif, type="V", min = HF.ofmom.6, max = HF.ofmom.12, rtrunc=TRUE, linf=0) 
HF.mom.mcstoc.12_24  <- mcstoc(runif, type="V", min = HF.ofmom.12, max = HF.ofmom.24, rtrunc=TRUE, linf=0) 
HF.mom.mcstoc.24_36  <- mcstoc(runif, type="V", min = HF.ofmom.24, max = HF.ofmom.30, rtrunc=TRUE, linf=0) 

# HMRE	Hand mouthing removal(transfer) eff.	-	day	beta	2	8				Cohen Hubal 2008 as presented in Ozkaynak 2011
HMRE.mcstoc <- mcstoc(rbeta, type="V", shape1 = 2, shape2 = 8, rtrunc=TRUE, linf=0)

### For now, for children under 6 mo, use the soil loading est for children 4 mo old
# C.soil.simple.f6 <- ((C.soil.load.4 * HM.child.f6.mcdata * HF.child.mcstoc) + (M.soil.load * (HM.other.nd.f6.mcdata + HM.other.d.f6.mcdata) * HF.mom.mcstoc))  * HMRE.mcstoc
# C.soil.simple.o24 <- ((C.soil.load.30 * HM.child.o24.mcdata * HF.child.mcstoc) + (M.soil.load * (HM.other.nd.o24.mcdata + HM.other.d.o24.mcdata) * HF.mom.mcstoc))  * HMRE.mcstoc
# 
# for(i in 6:23){
#   assign(paste("C.soil.simple.", i, sep = ""), ((get(paste("C.soil.load.", i, sep = "")) * get(paste("HM.SM.", i, ".child.mcdata", sep = ""))) + 
#                                                   (M.soil.load * get(paste("HM.SM.", i, ".mom.mcdata", sep = "")))) * HF.mcstoc * HMRE.mcstoc)
#   #C.soil.simple.10 <- ((C.soil.load.10 * HM.SM.10.child.mcdata) + (M.soil.load * HM.SM.10.mom.mcdata)) * HF.mcstoc * HMRE.mcstoc
# }

# From Gallan, 2011
awake.hr.3 <- 24 - mcstoc(rnorm, type = "V", mean = 13.6, sd = 2.1, rtrunc=TRUE, linf=0) # use for f6
awake.hr.6 <- 24 - mcstoc(rnorm, type = "V", mean = 12.9, sd = 2.1, rtrunc=TRUE, linf=0) # use for 6, 7
awake.hr.9 <- 24 - mcstoc(rnorm, type = "V", mean = 12.6, sd = 1.6, rtrunc=TRUE, linf=0) # use for 8, 9, 10
awake.hr.12 <- 24 - mcstoc(rnorm, type = "V", mean = 12.9, sd = 1.4, rtrunc=TRUE, linf=0) # use for 11, 12
awake.hr.12_24 <- 24 - mcstoc(rnorm, type = "V", mean = 12.6, sd = 1.3, rtrunc=TRUE, linf=0) # use for 12-24
awake.hr.24_36 <- 24 - mcstoc(rnorm, type = "V", mean = 12.0, sd = 1.2, rtrunc=TRUE, linf=0) # use for 24-36
awake.hr.36_48 <- 24 - mcstoc(rnorm, type = "V", mean = 12.0, sd = 1.2, rtrunc=TRUE, linf=0) # use for 36-48



# soil.C.conc <- mcdata(sample(soil.C.conc.long$soil.C.conc, size=ndvar(), replace=TRUE),type="V")
# 
# for(i in c(4:40)){
# }
# 
# M.soil.load <- (soil.M.load.mcdata / M.hand.SA.WASHB.end.mcstoc) * M.hand.SA.WASHB.mcstoc



########## Eval the mc model ###########
child.hand.soil.concentration.mg.cm2 <- soil.C.conc
child.hand.surface.area.cm2 <- C.hand.SA.WASHB.o24.mcdata  #C.hand.SA.WASHB.o24.mcdata
child.hand.mouth.frequency.events.hr <- HM.child.o24.mcdata  #HM.child.o24.mcdata
child.hand.fraction.mouthed <- HF.child.mcstoc
mom.hand.soil.load.mg <- soil.M.load.mcdata
#mom.soil.sample.hand.surface.area.cm2 <- M.hand.SA.WASHB.end.mcstoc
mom.hand.surface.area.cm2 <- M.hand.SA.WASHB.mcstoc
mom.hand.mouth.nonfood.frequency.events.hr <- HM.other.nd.o24.mcdata  #HM.other.nd.o24.mcdata
mom.feeding.frequency.events.hr <- HM.other.d.o24.mcdata  #HM.other.d.o24.mcdata
mom.hand.fraction.mouthed <- HF.mom.mcstoc.24_36  #HF.mom.mcstoc.24_30   #HF.mom.mcstoc.f6 is only one value so not good to use --> use HF.mom.mcstoc.6_12 instead, even though it is for the wrong age group
# > HF.mom.mcstoc.f6
# node    mode  nsv nsu nva variate    min   mean median    max Nas type outm
# 1    x numeric 1001   1   1       1 0.0317 0.0317 0.0317 0.0317   0    V each
# 
# > HF.mom.mcstoc.6_12
# node    mode  nsv nsu nva variate    min   mean median    max Nas type outm
# 1    x numeric 1001   1   1       1 0.0317 0.0337 0.0338 0.0358   0    V each
saliva.removal.efficiency <- HMRE.mcstoc
hours.awake <- awake.hr.24_36   # awake.hr.24_36
#amt.soil.ingested.mg.day <- C.soil.simple.o24

soil.test3 <- mcmodel({ 
  amt.soil.ingested.mg.day <- ((((child.hand.soil.concentration.mg.cm2 * child.hand.surface.area.cm2) * child.hand.mouth.frequency.events.hr * child.hand.fraction.mouthed) + 
                                  ((mom.hand.soil.load.mg / mom.hand.surface.area.cm2) * mom.hand.surface.area.cm2 * (mom.hand.mouth.nonfood.frequency.events.hr + mom.feeding.frequency.events.hr) * mom.hand.fraction.mouthed))  * saliva.removal.efficiency) * hours.awake 
  res.30 <- mc(child.hand.soil.concentration.mg.cm2, child.hand.surface.area.cm2, child.hand.mouth.frequency.events.hr, child.hand.fraction.mouthed, mom.hand.soil.load.mg , mom.hand.surface.area.cm2, mom.hand.mouth.nonfood.frequency.events.hr, mom.feeding.frequency.events.hr, mom.hand.fraction.mouthed, saliva.removal.efficiency, hours.awake, amt.soil.ingested.mg.day) #sleep.hr.24_36, 
})

res3 <- evalmcmod(soil.test3, nsv=ndvar(), nsu=ndunc(), seed=seed)
plot(res3, prec = 0.001, stat = c("median", "mean"), lim = c(0.025, 0.25, 0.75, 0.975), na.rm =
       TRUE, griddim = NULL, xlab = NULL, ylab = "Fn(x)", main = "", draw = TRUE, paint = TRUE)
tor.res3 <- tornado(res3, use = "complete.obs" )
plot(tor.res3)
summary(res3)
hist(res3)

soil.test1 <- mcmodel({ 
  C.soil.simple.f6 <- ((((soil.C.conc * C.hand.SA.WASHB.f6.mcdata) * HM.child.f6.mcdata * HF.child.mcstoc) + ((soil.M.load.mcdata / M.hand.SA.WASHB.mcstoc) * M.hand.SA.WASHB.mcstoc * (HM.other.nd.f6.mcdata + HM.other.d.f6.mcdata) * HF.mom.mcstoc.6_12))  * HMRE.mcstoc) * awake.hr.3
  # for(i in 6:23){
  #   assign(paste("C.soil.simple.", i, sep = ""), ((get(paste("C.soil.load.", i, sep = "")) * get(paste("HM.SM.", i, ".child.mcdata", sep = ""))) + 
  #                                                   (M.soil.load * get(paste("HM.SM.", i, ".mom.mcdata", sep = "")))) * HF.mcstoc * HMRE.mcstoc)
  #   #C.soil.simple.10 <- ((C.soil.load.10 * HM.SM.10.child.mcdata) + (M.soil.load * HM.SM.10.mom.mcdata)) * HF.mcstoc * HMRE.mcstoc
  # }
  #C.soil.simple.6 <- ((C.soil.load.6 * HM.SM.6.child.mcdata * HF.child.mcstoc) + (M.soil.load * HM.SM.6.mom.mcdata * HF.mom.mcstoc))  * HMRE.mcstoc * (24-sleep.hr.6)
  #C.soil.simple.7 <- ((C.soil.load.7 * HM.SM.7.child.mcdata * HF.child.mcstoc) + (M.soil.load * HM.SM.7.mom.mcdata * HF.mom.mcstoc))  * HMRE.mcstoc * (24-sleep.hr.6)
  
  res.4 <- mc(soil.C.conc, C.hand.SA.WASHB.f6.mcdata, HM.child.f6.mcdata, HF.child.mcstoc, soil.M.load.mcdata, M.hand.SA.WASHB.mcstoc, HM.other.nd.f6.mcdata, HM.other.d.f6.mcdata, HF.mom.mcstoc.6_12, HMRE.mcstoc, awake.hr.3, C.soil.simple.f6) # sleep.hr.3, 
  #res.6 <- mc(C.soil.load.6, HM.SM.6.child.mcdata, HF.child.mcstoc, M.soil.load, HM.SM.6.mom.mcdata, HF.mom.mcstoc, HMRE.mcstoc, sleep.hr.6, C.soil.simple.f6)
  #res.7 <- mc(C.soil.load.7, HM.SM.7.child.mcdata, HF.child.mcstoc, M.soil.load, HM.SM.7.mom.mcdata, HF.mom.mcstoc, HMRE.mcstoc, sleep.hr.6, C.soil.simple.f6)
})

soil.test2 <- mcmodel({ 
  C.soil.simple.o24 <- ((((soil.C.conc * C.hand.SA.WASHB.o24.mcdata) * HM.child.o24.mcdata * HF.child.mcstoc) + ((soil.M.load.mcdata / M.hand.SA.WASHB.mcstoc) * M.hand.SA.WASHB.mcstoc * (HM.other.nd.o24.mcdata + HM.other.d.o24.mcdata) * HF.mom.mcstoc.24_30))  * HMRE.mcstoc) * awake.hr.24_36
  res.30 <- mc(soil.C.conc, C.hand.SA.WASHB.o24.mcdata, HM.child.o24.mcdata, HF.child.mcstoc, soil.M.load.mcdata, M.hand.SA.WASHB.mcstoc, HM.other.nd.o24.mcdata, HM.other.d.o24.mcdata, HF.mom.mcstoc.24_30, HMRE.mcstoc, awake.hr.24_36, C.soil.simple.o24) #sleep.hr.24_36, 
})

res <- evalmcmod(soil.test1, nsv=ndvar(), nsu=ndunc(), seed=seed)
plot(res, prec = 0.001, stat = c("median", "mean"), lim = c(0.025, 0.25, 0.75, 0.975), na.rm =
       TRUE, griddim = NULL, xlab = NULL, ylab = "Fn(x)", main = "", draw = TRUE, paint = TRUE)
tor.res <- tornado(res, use = "complete.obs" )
plot(tor.res)
summary(res)
hist(res)

res2 <- evalmcmod(soil.test2, nsv=ndvar(), nsu=ndunc(), seed=seed)
plot(res2, prec = 0.001, stat = c("median", "mean"), lim = c(0.025, 0.25, 0.75, 0.975), na.rm =
       TRUE, griddim = NULL, xlab = NULL, ylab = "Fn(x)", main = "", draw = TRUE, paint = TRUE)
tor.res2 <- tornado(res2, use = "complete.obs" )
plot(tor.res2)
summary(res2)
hist(res2)

# C.soil.simple.o24 :
#   mean   sd  Min 2.5% 25% 50% 75% 97.5%   Max  nsv Nas
# NoUnc  731 1391 8.55   18 116 328 751  4600 15838 1001  723
# 


###### Translate from ingestion per hour to ingestion per day based on number of hours awake per day ###############
######## This makes the assumption that the frequency of HM and SM soil mouthing averaged over the duration of observation represents the average over the entire day
##### Check this assumption by.......
#### Rationalizing that we captured at least one (of how many?) feeding periods.



# Load hours awake per day (VO.R3 Status_161210)
obs.stat.v3.raw <- read_excel("C:/Users/Laura Kwong/Box Sync/VO R3/VO.R3 Status_170221.xlsx")

soil.survey.v3 <- obs.stat.v3.raw %>%
  select(soil.today, soil.yest, soil.sevenday, soil.typicalday, soil.reaction, soil.healthy, soil.parent.eat, soil.parent.eat.day, soil.parent.eat.wk, soil.parent.eat.preg)


############### Calc hours awake per day at diff age groups ########################
#### Ignores within-child correlation of awake hours at different ages
hr.awake.base <- obs.stat.v3.raw %>%
  select(pid, calc.age.mo, awake.day.hr, minus12.months, minus12.awake.time, minus24.months, minus24.awake.time)


# The average number of hours awake was collected by recall on during the third observation
# Recall error could cause 

names(hr.awake.base) <- c("hh", "age", "awake.hr", "age", "awake.hr", "age", "awake.hr")
hr.awake <- rbind(hr.awake.base[, c(1, 2:3)], hr.awake.base[, c(1, 3:4)], hr.awake.base[, c(1, 5:6)]) %>%
  mutate(age.mo = floor(age)) %>%
  select(hh, age, age.mo, awake.hr) %>%
  filter(!is.na(hh), age.mo > 0)

hr.awake.find.dist <- function(data, x, dist) {fitdist(data[data$age.mo == x, "awake.hr"], dist, method = "mle")}

max(hr.awake$age.mo, na.rm = TRUE)
for(i in 4:43){ #max(HM.SM$age)
  assign(paste("hr.awake.", i, ".dist", sep = ""), hr.awake.find.dist(hr.awake, i, "weibull"))
  assign(paste("hr.awake", i, ".mcdata", sep = ""), mcdata(sample(hr.awake[hr.awake$age.mo == i, "awake.hr"], size = ndvar(), replace = TRUE), type = "V"))
}

scatter.smooth(hr.awake$age, hr.awake$awake.hr)

