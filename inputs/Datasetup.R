# Data setup for parametric survival analysis models
# MM accumulation

# This file uses CPRD data from the patient and observation files to create data
# in counting process format for fitting parametric survival models to explore
# multimorbidity accumulation. Time is age, reset so that age 18 = time 0

# 4 states of multimorbidity accumulation:
#   1. Healthy (no chronic conditions of interest)
#   2. Initial chronic condition (1 chronic condition of interest)
#   3. Basic multimorbidity (2 or more chronic conditions of interest; BMM)
#   4. Complex multimorbidity (3 or more chronic conditions of intetest in 3 + body systems; CMM)

# Assuming all conditions are chronic (no recovery possible) this results in
# 7 unidirectional possible tranisitons:
#   T1:Healthy--> initial chronic condition
#   T2:Healthy--> death
#   T3:Initial chronic condition--> basic multimorbidity
#   T4:Initial chronic condition--> death
#   T5:Basic multimorbidity--> complex multimorbidity
#   T6:Basic multimorbidity--> death
#   T7:Complex multimorbidity--> death

# Input files:
#   1. 'Patient' file
#   2. Table of disease diagnosis dates:data from 'observation' file, cleaned so
#   that each individual has 1 line per 'diagnosis' of conditions of interest,
#   with the date of 'diagnosis'
# https://github.com/annalhead/CPRD_multimorbidity_trends contains code for preparing
# these files


# Set up ----

library(data.table)
data.table::setDTthreads(10) #this is so that don't use all the processors
library(fst)
library(flexsurv)
data.table::setDTthreads(10) #this is so that don't use all the processors
library(ggplot2)
library(ggthemes)
library(mstate)


data_dir_CPRD <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/Data May 2020/", x)

data_dir_Diseaselist <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/Disease_lists/", x)


data_dir_lookup <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/Dictionaries etc/", x)

data_dir_PM <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/forthesis/", x)


# Load patient file for patient characteristics
patient <- read_fst(data_dir_CPRD("patient.fst"),
                    columns = c("patid", "pracid", "gender", "imd",
                                "yob", "dob", "reg1yr", "censordate",
                                "censorreason", "region"),
                    as.data.table = TRUE)

# Load cleaned observation files - each row is 'diagnosis' of a condition
obstab <- read_fst(data_dir_CPRD("obstab14Dec.fst"),
                   columns =c("patid", "eventdate", "disease_num"),
                   as.data.table = T)


diseasesum <- fread(data_dir_Diseaselist("DiseaseSumm.csv"),
                    header = TRUE, sep = ",",
                    select = c("disease_num", "system_num"))

cprd_reg <- fread(file = data_dir_lookup("Region.txt"),
                  header = TRUE, sep = "\t")




# MM states setup ----

# Adding the MM states and calculating time in between states for each individual
# 1 row = 1 individual

if (file.exists(data_dir_PM("accumDT_atf.fst"))) {
  myDT <-
    read_fst(data_dir_PM("accumDT_atf.fst"), as.data.table = TRUE)
} else {

# Checking that have the same included individuals here as in the panel dataset
#  prepared in https://github.com/annalhead/CPRD_multimorbidity_trends
  combi_mm <-
    read_fst(data_dir_CPRD("combi_mm_detailed.fst"), as.data.table = T)

  study_pats <- combi_mm[, unique(patid)]

  patient <- patient[patid %in% study_pats]
  rm(combi_mm, study_pats)
  gc()


  patient[, enterstudy := reg1yr]
  patient[reg1yr < "2004-01-01", enterstudy := as.IDate("2004-01-01")]
  patient[enterstudy - dob < 6574, enterstudy := dob + 6574]
  myDT <-
    patient[enterstudy < censordate,] #getting rid of people who died before 1yr reg/1Jan2004
  #Rms 425

  #Adding in region name
  setnames(myDT, "region", "regionid")
  cprd_reg[, regionid := as.factor(regionid)]
  myDT[cprd_reg, on = 'regionid', region := i.Description]
  myDT[, region := factor(
    region,
    levels = c(
      "London",
      "South West",
      "South Central",
      "South East Coast",
      "West Midlands",
      "East Midlands",
      "East of England",
      "North West",
      "Yorkshire And The Humber",
      "North East"
    ),
    labels = c(
      "London",
      "South West",
      "South Central",
      "South East Coast",
      "West Midlands",
      "East Midlands",
      "East of England",
      "North West",
      "Yorkshire And The Humber",
      "North East"
    )
  )]


  # Let's make a table so that each patient has 1 line with all the key dates
  setkey(obstab, patid, eventdate)

  #adding in date & number the 1st condition
  tmp <- obstab[obstab[, .I[1], keyby = c("patid")]$V1]
  myDT[tmp, on = 'patid', `:=` (firstcond = i.disease_num, firstdate = i.eventdate)]

  #adding in date of 2nd cond (bmm)
  tmp <- obstab[obstab[, .I[2], keyby = c("patid")]$V1]
  myDT[tmp, on = 'patid', `:=` (bmmdate = i.eventdate)]

  #Adding in date of 3rd body system (cmm)
  obstab[diseasesum, on = 'disease_num', system_num := i.system_num]
  tmp <- obstab[,-c("disease_num")]
  #rm(diseasesum)
  setkey(tmp, patid, eventdate)
  tmp <-
    tmp[tmp[, .I[1], keyby = c("patid", "system_num")]$V1]  #takes the nth line - faster than .SD
  setkey(tmp, patid, eventdate)
  tmp <- tmp[tmp[, .I[3], keyby = c("patid")]$V1][!is.na(patid)]
  myDT[tmp, on = 'patid', `:=` (cmmdate = i.eventdate)]
  rm(obstab, tmp)


  myDT[, `:=` (
    firstyn = fifelse(is.na(firstdate), 0L, 1L),
    bmmyn = fifelse(is.na(bmmdate), 0L, 1L),
    cmmyn = fifelse(is.na(cmmdate), 0L, 1L)
  )]
  myDT[, `:=` (
    incfirst = fifelse(is.na(firstdate), 0L, fifelse(firstdate > enterstudy, 1L, 0L)),
    incbmm = fifelse(is.na(bmmdate), 0L, fifelse(bmmdate > enterstudy, 1L, 0L)),
    inccmm = fifelse(is.na(cmmdate), 0L, fifelse(cmmdate > enterstudy, 1L, 0L))
  )]

  myDT <-
    myDT[, .(
      patid,
      gender,
      imd,
      region,
      yob,
      dob,
      enterstudy,
      censordate,
      censorreason,
      firstyn,
      firstdate,
      firstcond,
      bmmyn,
      bmmdate,
      cmmyn,
      cmmdate
    )]


  #How many diagnoses happen on the same date
  myDT[firstdate == bmmdate , .N / nrow(myDT)] #4%
  myDT[bmmdate == cmmdate, .N / nrow(myDT)] #1%
  myDT[bmmdate == cmmdate & firstdate == cmmdate, .N / nrow(myDT)] #0.4%

  myDT[firstdate == bmmdate &
         (firstdate != cmmdate | is.na(cmmdate)), .N] # 38482
  myDT[firstdate != bmmdate & bmmdate == cmmdate, .N] #6759
  myDT[bmmdate == cmmdate & firstdate == cmmdate, .N] #3843


  #shifting by 1 day if 2 diagnoses happen on the same day
  myDT[firstdate == bmmdate, bmmdate := bmmdate + 1]
  myDT[bmmdate == cmmdate |
         bmmdate > cmmdate, #need bmmdate > cmmdate to account for prev line
       cmmdate := bmmdate + 1]
  myDT[censordate == firstdate, censordate := firstdate + 1]
  myDT[censordate == bmmdate, censordate := bmmdate + 1]
  myDT[censordate == cmmdate, censordate := cmmdate + 1]


  myDT[, startstate := 1]
  myDT[firstdate <= enterstudy, startstate := 2]
  myDT[bmmdate <= enterstudy, startstate := 3]
  myDT[cmmdate <= enterstudy, startstate := 4]

  # Adding days in between states
  #This is with days since reg as the timeframe
  myDT[, `:=` (
    yearenterstudy = year(enterstudy),
    firstcond_t = ifelse(firstyn == 1, firstdate - enterstudy, censordate - enterstudy),
    bmm_t = ifelse(bmmyn == 1, bmmdate - enterstudy, censordate - enterstudy),
    cmm_t = ifelse(cmmyn == 1, cmmdate - enterstudy, censordate - enterstudy),
    death_t = censordate - enterstudy,
    censor_t =  censordate - enterstudy,
    death = ifelse(censorreason == 1, 1, 0),
    censored = ifelse(censorreason == 0, 1, 0)
  )]
  myDT[startstate == 2, firstcond_t := 0]
  myDT[startstate == 3, c("firstcond_t", "bmm_t") := 0]
  myDT[startstate == 4, c("firstcond_t", "bmm_t", "cmm_t") := 0]
  myDT[, starttime := 0]


  # Doing time in age
  myDT[, ageenter := (enterstudy - dob) / 365.24]
  myDT[, agefirst := ifelse(firstdate > enterstudy,
                            (firstdate - dob) / 365.24, NA)]
  myDT[, agebmm := ifelse(bmmdate  > enterstudy,
                          (bmmdate - dob) / 365.24,  NA)]
  myDT[, agecmm := ifelse(cmmdate > enterstudy,
                          (cmmdate - dob) / 365.24, NA)]
  myDT[, agedeath :=  (censordate - dob) / 365.24]
  myDT[, agecensor :=  (censordate - dob) / 365.24]


  myDT[, birthcohort := cut_width(yob, width = 5, boundary = 1899)]


  #Where age~ is na, change to 0
  myDT[is.na(agefirst), agefirst := 0]
  myDT[is.na(agebmm), agebmm := 0]
  myDT[is.na(agecmm), agecmm := 0]

  myDT[, birthcohort := cut_width(yob, width = 5, boundary = 1899)]
  myDT[yob < 1920, birthcohort := "<1920"]
  myDT[, birthcohort :=
         factor (
           birthcohort,
           levels = c(
             "<1920",
             "(1919,1924]" ,
             "(1924,1929]",
             "(1929,1934]",
             "(1934,1939]",
             "(1939,1944]",
             "(1944,1949]",
             "(1949,1954]",
             "(1954,1959]",
             "(1959,1964]",
             "(1964,1969]",
             "(1969,1974]" ,
             "(1974,1979]",
             "(1979,1984]",
             "(1984,1989]",
             "(1989,1994]",
             "(1994,1999]",
             "(1999,2004]"
           ),
           labels = c(
             "<1920",
             "(1919,1924]",
             "(1924,1929]",
             "(1929,1934]",
             "(1934,1939]",
             "(1939,1944]",
             "(1944,1949]",
             "(1949,1954]",
             "(1954,1959]",
             "(1959,1964]",
             "(1964,1969]",
             "(1969,1974]" ,
             "(1974,1979]",
             "(1979,1984]",
             "(1984,1989]",
             "(1989,1994]",
             "(1994,1999]",
             "(1999,2004]"
           )
         )]

  #doing 5 year birth-cohorts for those that match the ONS 2019 population
  myDT[, birthcohortONS := birthcohort]
  myDT[yob < 1930, birthcohortONS := "<1930"]
  myDT[, birthcohortONS :=
         factor (
           birthcohortONS,
           levels = c(
             "<1930",
             "(1929,1934]",
             "(1934,1939]",
             "(1939,1944]",
             "(1944,1949]",
             "(1949,1954]",
             "(1954,1959]",
             "(1959,1964]",
             "(1964,1969]",
             "(1969,1974]" ,
             "(1974,1979]",
             "(1979,1984]",
             "(1984,1989]",
             "(1989,1994]",
             "(1994,1999]",
             "(1999,2004]"
           ),
           labels = c(
             "<1930",
             "(1929,1934]",
             "(1934,1939]",
             "(1939,1944]",
             "(1944,1949]",
             "(1949,1954]",
             "(1954,1959]",
             "(1959,1964]",
             "(1964,1969]",
             "(1969,1974]" ,
             "(1974,1979]",
             "(1979,1984]",
             "(1984,1989]",
             "(1989,1994]",
             "(1994,1999]",
             "(1999,2004]"
           )
         )]

  myDT[, birthcohort10 := cut_width(yob, width = 10, boundary = 1894)]
  myDT[yob < 1915, birthcohort10 := "<1915"]
  myDT[, birthcohort10 :=
         factor (
           birthcohort10,
           levels = c(
             "<1915",
             "(1.91e+03,1.92e+03]",
             "(1.92e+03,1.93e+03]",
             "(1.93e+03,1.94e+03]",
             "(1.94e+03,1.95e+03]",
             "(1.95e+03,1.96e+03]",
             "(1.96e+03,1.97e+03]",
             "(1.97e+03,1.98e+03]",
             "(1.98e+03,1.99e+03]",
             "(1.99e+03,2e+03]"
           ),
           labels = c(
             "<1915",
             "(1914,1924]",
             "(1924,1934]",
             "(1934,1944]",
             "(1944,1954]",
             "(1954,1964]",
             "(1964,1974]",
             "(1974,1984]",
             "(1984,1994]",
             "(1994,2004]"
           )
         )]

  myDT[, birthcohort10b := cut_width(yob, width = 10, boundary = 1899)]
  myDT[yob < 1920, birthcohort10b := "<1920"]
  myDT[yob >= 1990, birthcohort10b := ">=1990"]
  myDT[, birthcohort10b :=
         factor (
           birthcohort10b,
           levels = c(
             "<1920",
             "(1.92e+03,1.93e+03]",
             "(1.93e+03,1.94e+03]",
             "(1.94e+03,1.95e+03]",
             "(1.95e+03,1.96e+03]",
             "(1.96e+03,1.97e+03]",
             "(1.97e+03,1.98e+03]",
             "(1.98e+03,1.99e+03]",
             ">=1990"
           ),
           labels = c(
             "<1920",
             "(1920,1930]",
             "(1930,1940]",
             "(1940,1950]",
             "(1950,1960]",
             "(1960,1970]",
             "(1970,1980]",
             "(1980,1990]",
             ">=1990"
           )
         )]


  myDT[, ageentergrp5 := cut_width(ageenter,
                                   width = 5,
                                   boundary = 20,
                                   closed = "left")]
  myDT[ageenter >= 90 , ageentergrp5 := ">=90"]

  write_fst(myDT, data_dir_PM("accumDT_atf.fst"), 100)
}



# Counting process setup ----
# following this tutorial https://www.r-bloggers.com/2018/12/simulating-multi-state-models-with-r/
# for the mstate bit

#Making the transition matrix - don't have censoring as a state
tmat <- mstate::transMat(x = list(c(2, 5),
                                  c(3, 5),
                                  c(4, 5),
                                  c(5),
                                  c()),
                         names = c("Healthy", "IncCond", "BMM", "CMM", "Death"))
tmat

if (file.exists(data_dir_PM("accumDTmstate_age.fst"))) {
  myDTmstate_age <-
    read_fst(data_dir_PM("accumDTmstate_age.fst"), as.data.table = TRUE)
} else {

  myDTmstate_age <-  msprep(
    data = as.data.frame(myDT),
    trans = tmat,
    time = c(NA, "agefirst", "agebmm", "agecmm", "agedeath"),
    status  = c(NA, "firstyn", "bmmyn", "cmmyn", "death"),
    start = list(state = myDT[, startstate],
                 time = myDT[, ageenter]),
    keep = c(
      "patid",
      "gender",
      "ageenter",
      "yearenterstudy",
      "imd",
      "region",
      "firstcond",
      "yob",
      "birthcohort",
      "birthcohort10",
      "birthcohortONS"
    )
  )
  myDTmstate_age <- as.data.table(myDTmstate_age)

  #For now, taking out indeterminate gender, imd == "" & region == ""
  myDTmstate_age <-
    myDTmstate_age[imd != "" & region != "" & gender != "I"]
  myDTmstate_age <-
    myDTmstate_age[, gender := factor(gender,
                                      levels = c("M", "F"),
                                      labels = c("M", "F"))]
  myDTmstate_age <-
    myDTmstate_age[, imd := factor(imd,
                                   levels = c(1, 2, 3, 4, 5),
                                   labels = c(1, 2, 3, 4, 5))]

  #Creating new start/stop where age 18 == time 0
  #this is because models don't fit otherwise, and nothing happens in the first 18 years
  myDTmstate_age[, Tstart2 := Tstart - 18]
  myDTmstate_age[, Tstop2 := Tstop - 18]
  myDTmstate_age[Tstart < 18 , Tstart2 := 0]
  myDTmstate_age[Tstart < 18 , Tstop2 := Tstop - Tstart]
  write_fst(myDTmstate_age, data_dir_PM("accumDTmstate_age.fst"), 100)
}



# No left censoring set-up ----

# Fitting the models to data without left-censoring

#Take all the healthy people
healthy <- myDTmstate_age[from == 1, id]
myDTmstate_age_nolc <- myDTmstate_age[id %in% healthy]

#Take all the people where know the BMM start date
knowbmmdate <- myDTmstate_age[!id %in% healthy & from == 2 & to == 3, id]
myDTmstate_age_nolc <- rbind(myDTmstate_age_nolc, myDTmstate_age[id %in% knowbmmdate & from > 2])

#Take all the people where know the BMM start date
knowcmmdate <- myDTmstate_age[!id %in% healthy & !id %in% knowbmmdate & from == 3 & to == 4, id]
myDTmstate_age_nolc <- rbind(myDTmstate_age_nolc, myDTmstate_age[id %in% knowcmmdate & from > 3])
write_fst(myDTmstate_age_nolc, data_dir_PM("myDTmstate_age_nolc.fst"), 100)
