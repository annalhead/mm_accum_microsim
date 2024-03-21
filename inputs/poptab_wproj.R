# Creating baseline population

# This file creates a simulant population of 30-90 year-olds that is 1% of the ONS
# population mid-year estimates for 2019.
# Available from https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/lowersuperoutputareamidyearpopulationestimates

# For simulants who are aged 30 in years 2020-2049 (final year of the simulation),
# I have used ONS population projections.
# Available from https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationprojections/datasets/z1zippedpopulationprojectionsdatafilesuk
# As only projections only available by age and sex, have assigned region & IMD
# randomly based on distribution among 30 year-olds in 2019

# For start state (healthy, 1 condition, BMM, CMM), assigned randomly based on
# distribution of start states in CPRD data sample by sex, IMD, and 5-year age-group


library(data.table)
data.table::setDTthreads(5) #this is so that don't use all the processors
library(fst)
fst::threads_fst(5)
library(readxl)
library(ggplot2) # for cut width
library(hesim)


data_dir_lookup <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/Dictionaries etc/", x)


data_dir_CPRD <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/Data May 2020/", x)

create_patients <- function(pat, n_pats){
  patients <- pat[rep(1, n_pats), ]
  return(patients[, ])
}


lsoa_imd <- read_fst(data_dir_lookup("lsoa_to_locality_indx.fst"),
                     as.data.table = T)
#View(lsoa_imd)
lsoa_imd <- as.data.table(lsoa_imd)
lsoa_pop_males <- read_excel(data_dir_lookup("SAPE22DT2-mid-2019-lsoa-syoa-estimates-unformatted.xlsx"),
                             sheet = "Mid-2019 Males", skip = 3 )
#View(lsoa_pop_males)
lsoa_pop_males <- as.data.table(lsoa_pop_males)
lsoa_pop_females <- read_excel(data_dir_lookup("SAPE22DT2-mid-2019-lsoa-syoa-estimates-unformatted.xlsx"),
                               sheet = "Mid-2019 Females", skip = 3 )
#View(lsoa_pop_females)
lsoa_pop_females <- as.data.table(lsoa_pop_females)

lsoa_pop_males <- lsoa_pop_males[`LSOA Code` %like% "E"]
lsoa_pop_females <- lsoa_pop_females[`LSOA Code` %like% "E"]

setnames(lsoa_imd, "LSOA11CD","LSOA Code")
lsoa_pop_males[lsoa_imd, on = 'LSOA Code', `:=` (qimd = i.qimd, region = i.SHA11NM)]
lsoa_pop_females[lsoa_imd, on = 'LSOA Code', `:=` (qimd = i.qimd, region = i.SHA11NM)]
rm(lsoa_imd)

#want sum of people in each age, imd & region category

#Don't want people under 30 - will add these in
agecols <- c(as.character(30:89), "90+")
maletab <- lsoa_pop_males[, lapply(.SD, sum),
                          .SDcols = agecols,
                          keyby = .(qimd, region)]

femaletab <- lsoa_pop_females[, lapply(.SD, sum),
                              .SDcols = agecols,
                              keyby = .(qimd,  region)]
poptab <- rbind(maletab[, gender := "M"],
                femaletab[, gender := "F"])
rm(lsoa_pop_females, lsoa_pop_males)
cohtab <- melt(poptab, measure.vars = agecols,
               variable.name = "age", value.name = "N")
rm(agecols, poptab)

cohtab[, gender := factor(gender, levels = c("M", "F"))]

# swapping imd order so that 1 is least deprived & 5 most deprived
cohtab[, imd := factor(levels = c(1:5),
                       labels = c("1", "2", "3", "4", "5"))]
cohtab[qimd == "1 most deprived", imd := "5"]
cohtab[qimd == "2", imd := "4"]
cohtab[qimd == "3", imd := "3"]
cohtab[qimd == "4", imd := "2"]
cohtab[qimd == "5 least deprived", imd := "1"]

#Adding in yob
cohtab[, age := as.numeric(as.character(age))]
cohtab[is.na(age), age := 90]
cohtab[, yob := 2019-age]

cohtab[, cohortnum := c(1:nrow(cohtab))]

patients <- data.table()
for (i in 1:nrow(cohtab)) {
  subtab <- create_patients (cohtab[cohortnum == i, .(gender, age, imd, region)],
                             n_pats = round(cohtab[cohortnum == i, N/100]))
  subtab[, cohortnum := i]
  patients <- rbind(patients, subtab)
  rm(subtab)
}
rm(i)
patients[, `:=` (yob = 2019-age, init_year = 2019L)]






popproject <- read_excel("/mnt/alhead/UoL/CPRD2019mm/Dictionaries etc/ukpppopendata2020_popprojections.xls",
               sheet = "Population")
setDT(popproject)
popproject <- popproject[Age == 30]
popproject[, gender := factor(Sex,
                              levels = c(1,2),
                              labels = c("M", "F"))][
                                ,c("Sex", "Age"):= NULL]
popproject <- melt(popproject, id.vars = c("gender"), variable.name = "init_year", value.name = "N", variable.factor = FALSE)
popproject <- popproject[ init_year < 2050][, init_year := as.integer(init_year)]
popproject[, group := c(1:nrow(popproject))]

thirty <- cohtab[age == 30, .(gender, imd, region, N)]
thirtyall <- thirty[, sum(N)]
thirty[, pc := N/thirtyall]
setkey(thirty, "gender", "imd", "region")
thirty[, group2 := rep(c(1:50), 2)]




patients30 <- data.table()
for (i in 1:60) {
  subtab <- create_patients (popproject[group == i, .(gender, init_year)],
                             n_pats = round(popproject[group == i, N/100]))
  subtab[, cohortnum := i]
  patients30 <- rbind(patients30, subtab)
  rm(subtab)
}
rm(i)
n_pats <- patients30[, .N]
patients30$patient_id <- 1:n_pats

pc <- thirty[gender == "M", pc]
for (i in 1:30){
  npats <- patients30[cohortnum == i, .N]
  patients30[cohortnum == i, group2 := sample(1:50, npats, replace = T, pc)]
  rm(npats)
  patients30
}

pc <- thirty[gender == "F", pc]
for (i in 31:60){
  npats <- patients30[cohortnum == i, .N]
  patients30[cohortnum == i, group2 := sample(1:50, npats, replace = T, pc)]
  rm(npats)
  patients30
}

patients30[thirty, on = 'group2', `:=` (imd = i.imd, region = i.region) ]
patients30[, `:=` (yob = init_year - 30, ageenter = 30, age = (30-(init_year-2019)))][, c("patient_id", "group2") := NULL]

n_cohort <- cohtab[, max(cohortnum)]
patients30[, cohortnum := cohortnum + n_cohort]

cohtab30 <- patients30[patients30[, .I[1], by = cohortnum]$V1]
tmp <- patients30[, .N, keyby = .(cohortnum )]
cohtab30[tmp, on = 'cohortnum', N := i.N]
rm(tmp)

patients <- rbind(patients, patients30, fill = TRUE)
setkey(patients, cohortnum, init_year)

n_pats <- patients[, .N]
patients$patient_id <- 1:n_pats

# Want to add in start states based on the distribution of start states in my data
# Doing this based on sex, imd & agegroup at entry.
# myDT <-
#   read_fst(data_dir_CPRD("accumDT_atf.fst"), as.data.table = TRUE)[imd != "" & gender != "I"]
#
# myDT[, ageentergrp5 := cut_width(ageenter, width = 5, boundary = 20, closed = "left")]
# myDT[ageenter >=90 ,ageentergrp5 := ">=90"]
# startstateprob <- myDT[, .N, by = .(gender, imd, ageentergrp5, startstate)]
# setnames(startstateprob, "ageentergrp5", "agegrp5")

# Want to add in start states based on the distribution of start states in my data
# Doing this based on sex, imd & agegroup in 2019
myDT <- read_fst(data_dir_CPRD("combi_mm_detailed.fst"), as.data.table = T)[
  imd != "" & gender != "Indeterminate"]
myDT[, imd := factor(imd)]
myDT[, gender := factor(gender,
                        levels = c("Men", "Women"),
                        labels = c("M", "F"))]
myDT[, `:=` (agegrp5 = cut_width(age, width = 5, boundary = 20, closed = "left"))]
myDT[age >= 90, agegrp5 := ">=90"]
myDT[, startstate := ifelse(cmm != 0, 4,
                            ifelse(bmm != 0, 3,
                                   ifelse(onecond != 0, 2,
                                          1)))]
startstateprob <- myDT[, .N, by = .(gender, imd, agegrp5, startstate)]

cohtab <- rbind(cohtab, cohtab30, fill = TRUE)
cohtab[, `:=` (agegrp5 = cut_width(age, width = 5, boundary = 20, closed = "left"))]
cohtab[age == 90, agegrp5 := ">=90"]
cohtab2 <- cohtab[, .N, by = .(imd, gender, agegrp5, init_year)]
cohtab2[, `:=` (N = NULL, cohortnum2 = c(1:nrow(cohtab2)))]
cohtab <- cohtab[cohtab2, on = .(imd, gender, agegrp5)]

cohtab2[agegrp5 == "[85,90]", agegrp5 := "[85,90)" ]
#cohtab2[agegrp5 %in% c("[0,5)", "[5,10)", "[10,15)"), agegrp5 := "[15,20)"] #For the younger cohorts to feed in, using the same proportions as the 12-20 cohort
#As starting simulation at age 30 - using startstates based on those age 30
cohtab2[agegrp5 %in% c("[0,5)", "[5,10)", "[10,15)", "[15,20)", "[20,25)", "[25,30)"),
        agegrp5 := "[30,35)"] #For the younger cohorts to feed in, using the same proportions as the 30-35 cohort
startstateprob <- startstateprob[cohtab2, on = .(gender, imd, agegrp5)]
setkey(startstateprob, cohortnum2, startstate)
startstateprob[, `:=` (totalN = sum(N), prop = N/sum(N)), by = cohortnum2]
patients[cohtab, on = 'cohortnum', cohortnum2 := i.cohortnum2 ]
patients[, state_id := 0]


for (i in 1:nrow(cohtab2)){
  pc <- startstateprob[cohortnum2 == i, prop]
  npats <- patients[cohortnum2 == i, .N]
  patients[cohortnum2 == i, state_id := sample(1:4, npats, replace = T, pc)]
  rm(pc, npats)
  patients
}


cohtab[, birthcohort := cut_width(yob, width = 5, boundary = 1899)]
cohtab[yob <1930, birthcohort := "<1930"]
cohtab[, birthcohort :=
         factor (birthcohort,
                 levels = c("<1930", "[1929,1934]", "(1934,1939]",
                            "(1939,1944]", "(1944,1949]", "(1949,1954]",
                            "(1954,1959]", "(1959,1964]","(1964,1969]",
                            "(1969,1974]" ,"(1974,1979]", "(1979,1984]",
                            "(1984,1989]", "(1989,1994]", "(1994,1999]",
                            "(1999,2004]", "(2004,2009]","(2009,2014]", "(2014,2019]" ),
                 labels = c("<1930",  "(1929,1934]", #changing factor names to match those in the fitted model
                            "(1934,1939]",
                            "(1939,1944]", "(1944,1949]", "(1949,1954]",
                            "(1954,1959]", "(1959,1964]","(1964,1969]",
                            "(1969,1974]" ,"(1974,1979]", "(1979,1984]",
                            "(1984,1989]", "(1989,1994]", "(1994,1999]",
                            "(1999,2004]", "(2004,2009]","(2009,2014]", "(2014,2019]"))]



#renaming yorkshire
cohtab[ region == "Yorkshire and the Humber", region:= "Yorkshire And The Humber"]


patients[, birthcohort := cut_width(yob, width = 5, boundary = 1899)]
patients[yob <1930, birthcohort := "<1930"]
patients[, birthcohort :=
         factor (birthcohort,
                 levels = c("<1930", "[1929,1934]", "(1934,1939]",
                            "(1939,1944]", "(1944,1949]", "(1949,1954]",
                            "(1954,1959]", "(1959,1964]","(1964,1969]",
                            "(1969,1974]" ,"(1974,1979]", "(1979,1984]",
                            "(1984,1989]", "(1989,1994]", "(1994,1999]",
                            "(1999,2004]", "(2004,2009]","(2009,2014]", "(2014,2019]" ),
                 labels = c("<1930",  "(1929,1934]", #changing factor names to match those in the fitted model
                            "(1934,1939]",
                            "(1939,1944]", "(1944,1949]", "(1949,1954]",
                            "(1954,1959]", "(1959,1964]","(1964,1969]",
                            "(1969,1974]" ,"(1974,1979]", "(1979,1984]",
                            "(1984,1989]", "(1989,1994]", "(1994,1999]",
                            "(1999,2004]", "(2004,2009]","(2009,2014]", "(2014,2019]"))]



#renaming yorkshire
patients[ region == "Yorkshire and the Humber", region:= "Yorkshire And The Humber"]

tmp <- patients[,.N, keyby = cohortnum]
cohtab[tmp, on = 'cohortnum', N_1pc := i.N]


write_fst(patients, data_dir_lookup("cohtab_ONS2019_indage_1pc_30proj.fst"), 100)
write_fst(cohtab, data_dir_lookup("cohtab_ONS2019_indage_summary_1pc_30proj.fst"), 100)

