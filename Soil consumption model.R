
### THis is a a practice for GitHub # Hi Annanda!


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

#### Think about the detection limit! The precision of the scale is 5 mg = 0.005 g, so anything less than 0.005 g should be treated as...
## THIS ASSUMPTION WILL MAKE A BIG DIFFERENCE
## How many samples are below the detection limit of 0.005 g?
qplot(soilMass.base$soil.g, geom="histogram", binwidth = 0.001)

sum(soilMass.base$soil.g <= 0.005)/length(soilMass.base$soil.g) # 57.0% of observations (including the reps) are less than the LOD

## Can try 1/2 the detection limit, though I remember hearing that there are better ways to handle this. 

##################################################################
###################################################################
## For now I will use this simple method and dig more into it later. 
soilMass.base$soil.g.halfLOD <- ifelse(soilMass.base$soil.g <= 0.005, 0.0025, soilMass.base$soil.g)
qplot(soilMass.base$soil.g.halfLOD, geom="histogram", binwidth = 0.001)

#####################################################################
####################################################################


# AFTER replacing the values below the detection limit with 1/2 the detection limit, average the reps
soilMass.reps <- soilMass.base %>%
  filter(!is.na(rep.dup))

soilMass.rep.means <- soilMass.reps %>%
  group_by(hh, motherOrChild) %>%
  #group_by_(setdiff(names(soilMass.base), "rep.dup")) %>% # uses "group_by_" so I can use quoted columns; groups by all col except rep.dup
  summarise(soil.g.halfLOD = mean(soil.g.halfLOD))

# in soilMass.base, replace the reps with the rep mean
# keep only REP 1, set the rep.dup col to null, remove the soil.g col and merge with the soilMass.rep.means dataset by hh and motherOrChild to add the soil.g col based on the means of the reps
soilMass.reps %>%
  filter(rep.dup == "REP1") %>%
  select(-soil.g.halfLOD) %>%
  left_join(soilMass.rep.means[,c("hh", "motherOrChild", "soil.g")], by = c("hh", "motherOrChild"))
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
soilMass$soil.g.halfLOD <- ifelse(soilMass$soil.g < 0.005, 0.0025, soilMass$soil.g)
qplot(soilMass$soil.g.halfLOD, geom="histogram", binwidth = 0.001)

#####################################################################
####################################################################

# Calculate the mass of soil on one hand by finding the mass / ml of sample tested, multiplying by the total sample volume (250 mL for children and 350 mL for mothers) and dividing by 2 to get the mass on one hand
soilMass$soilOneHand.g <- ((soilMass$soil.g.halfLOD/soilMass$volSample.ml) * soilMass$volTotal.ml)/2
qplot(soilMass$soilOneHand.g, geom="histogram", binwidth = 0.001)

soilMass.narrow <- soilMass[,c("hh", "RO1.round", "surveyDate.soil", "motherOrChild", "soilOneHand.g")]


# Separate into mother and child mass
soilMass.child.Round7 <- soilMass.narrow %>%
  filter(motherOrChild == "CH", RO1.round == 7) %>%
  select(hh, RO1.round, surveyDate.soil, soilOneHand.g)

soilMass.child.Round8 <- soilMass.narrow %>%
  filter(motherOrChild == "CH", RO1.round == 8) %>%
  select(hh, RO1.round, surveyDate.soil, soilOneHand.g)

soilMass.mom.Round7 <- soilMass.narrow %>%
  filter(motherOrChild == "MH", RO1.round == 7) %>%
  select(hh, RO1.round, surveyDate.soil, soilOneHand.g)

soilMass.mom.Round8 <- soilMass.narrow %>%
  filter(motherOrChild == "MH", RO1.round == 8) %>%
  select(hh, RO1.round, surveyDate.soil, soilOneHand.g)

### nope

soilMass.wide <- soilMass.child.Round7 %>%
  left_join(soilMass.child.Round8, by = "hh") %>%
  left_join(soilMass.mom.Round7, by = "hh") %>%
  left_join(soilMass.mom.Round8, by = "hh")

names(soilMass.wide) <- c("hh", "RO1.c7", "surveyDate.soil.child.7", "soilOneHand.g.child.7",
                          "RO1.c8", "surveyDate.soil.child.8", "soilOneHand.g.child.8",
                          "RO1.m7", "surveyDate.soil.mom.7", "soilOneHand.g.mom.7",
                          "RO1.m8", "surveyDate.soil.mom.8", "soilOneHand.g.mom.8")
soilMass.wide <- soilMass.wide[, c("hh", "surveyDate.soil.child.7", "soilOneHand.g.child.7",
           "surveyDate.soil.child.8", "soilOneHand.g.child.8",
           "surveyDate.soil.mom.7", "soilOneHand.g.mom.7",
           "surveyDate.soil.mom.8", "soilOneHand.g.mom.8")]

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

anth.soilMass.wide <- left_join(anth.soil.wide, soilMass.wide, by = "hh")
# The dates had the correct format when saved but the format was not imported so need to change from factor to date
anth.soilMass.wide$unique_dob <- ymd(anth.soilMass.wide$unique_dob)
anth.soilMass.wide$surveyDate.mid <- ymd(anth.soilMass.wide$surveyDate.mid)
anth.soilMass.wide$surveyDate.end <- ymd(anth.soilMass.wide$surveyDate.end)

anth.soilMass.wide$age.soil <- anth.soilMass.wide$surveyDate.soil - anth.soilMass.wide$unique_dob
anth.soilMass.wide$daysFromMid <- anth.soilMass.wide$surveyDate.soil - anth.soilMass.wide$surveyDate.mid
anth.soilMass.wide$daysFromEnd <- anth.soilMass.wide$surveyDate.soil - anth.soilMass.wide$surveyDate.end

anth.soilMass.wide %>%
  select(daysFromMid) %>%
  filter(!is.na(daysFromMid))

# Calc grams soil/cm2 for emprical data



# Multiply by cm2 of observed children to get est g soil for children in this age group
## Or probabilistic multiply soil/cm2 * cm2 distribution
# -> If this method then need to # Determine a hand size distribution for children 3-6 mo old, 6-12, 12-24, 24-36 months old

#### Assumes loading is constant by age, which is likely NOT true

# Determine a mass of soil on hands distribution for children in diff age groups

# Load HM frequencies and soil ingestion frequencies for each individual, not hh
vo123.objclass.base <- read.csv("C:/Users/Laura Kwong/Box Sync/VO R123/vo.11.objclass.csv")
vo123.objclass.base %>%
  filter(actobj.class %in% c("Mouth_hands", "Mouth_hands_d", "Mouth_hands_nd", "Mouth_soil"))

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
