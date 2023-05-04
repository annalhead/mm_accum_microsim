# Parametric survival analysis model fitting

# This file takes counting process format data for 4 health states (+death)
# and 7 transitions of MM accumulation

# Time between states is age in years, reset so that age 18 is time 0
# Distribution for each transition was determined in the hazplots_leftcensored.Rmd file
# which plots the baseline hazard against the intercept-only parametric curves

# Input data is prepared using the 'Datasetup.R' script

# 4 states of multimorbidity accumulation:
#   1. Healthy (no chronic conditions of interest)
#   2. Initial chronic condition (1 chronic condition of interest)
#   3. Basic multimorbidity (2 or more chronic conditions of interest; BMM)
#   4. Complex multimorbidity (3 or more chronic conditions of intetest in 3 + body systems; CMM)

# 7 unidirectional possible tranisitons:
#   T1:Healthy--> initial chronic condition
#   T2:Healthy--> death
#   T3:Initial chronic condition--> basic multimorbidity
#   T4:Initial chronic condition--> death
#   T5:Basic multimorbidity--> complex multimorbidity
#   T6:Basic multimorbidity--> death
#   T7:Complex multimorbidity--> death


library(data.table)
data.table::setDTthreads(3) #this is so that don't use all the processors
library(fst)
fst::threads_fst(3)
library(qs)
library(flexsurv)
library(ggplot2)
library(future.apply)
library(future)
library(splines)
options(future.fork.enable = TRUE) # TODO remove for production
options(future.rng.onMisuse = "ignore") # Remove false warning
plan(multicore(workers = 3))

data_dir_PM <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/forthesis/", x)

data_dir_Diseaselist <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/Disease_lists/", x)


data_dir_lookup <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/Dictionaries etc/", x)

output_path <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/forthesis/SAmodelfits/", x)




myDTmstate_age <- read_fst(data_dir_PM("accumDTmstate_age.fst"), as.data.table = TRUE)[
  yob >= 1919 & yob <= 1986 & imd!= "" & region != ""]

t1 <- flexsurvreg(Surv(Tstart2, Tstop2, status) ~
                    gender + imd + region + birthcohortONS,
                  anc = list(sigma = ~ gender+ imd + region,
                             Q = ~ gender+ imd + region),
                  data = subset(myDTmstate_age, trans == 1),
                  dist = "gengamma")
print(t1)
AHnew_lc_bcONS <- list((t1))
qsave(AHnew_lc_bcONS , output_path("AHnew_lc_bcONS_ancDEV.q"), nthreads = 10)

t2 <- flexsurvreg(Surv(Tstart2, Tstop2, status) ~
                    gender + imd + region + birthcohortONS,
                  anc = list(shape = ~ gender+ imd + region),
                  data = subset(myDTmstate_age, trans == 2),
                  dist = "gompertz")
print(t2)
AHnew_lc_bcONS <- list((t1), (t2))
qsave(AHnew_lc_bcONS , output_path("AHnew_lc_bcONS_ancDEV.q"), nthreads = 10)

t3 <- flexsurvreg(Surv(Tstart2, Tstop2, status) ~
                    gender + imd + region + birthcohortONS,
                  anc = list(sigma = ~ gender+ imd + region,
                             Q = ~ gender+ imd + region),
                  data = subset(myDTmstate_age, trans == 3),
                  dist = "gengamma")
if(!exists("t3")){
  t3 <- flexsurvreg(Surv(Tstart2, Tstop2, status) ~
                      gender + imd + region + birthcohortONS,
                    anc = list(shape = ~ gender+ imd + region),
                    data = subset(myDTmstate_age, trans == 3),
                    dist = "gompertz") #for now gompertz as the gengamma wont fit
  print("t3 gomp")
}
print(t3)

AHnew_lc_bcONS <- list((t1), (t2), (t3))
qsave(AHnew_lc_bcONS , output_path("AHnew_lc_bcONS_ancDEV.q"), nthreads = 10)

t4 <- flexsurvreg(Surv(Tstart2, Tstop2, status) ~
                    gender + imd + region + birthcohortONS,
                  anc = list(shape = ~ gender + imd),
                  data = subset(myDTmstate_age, trans == 4),
                  dist = "gompertz")
print(t4)
AHnew_lc_bcONS <- list((t1), (t2), (t3), (t4))
qsave(AHnew_lc_bcONS , output_path("AHnew_lc_bcONS_ancDEV.q"), nthreads = 10)

t5 <- flexsurvreg(Surv(Tstart2, Tstop2, status) ~
                    gender + imd + region + birthcohortONS, #doesnt fit with the anc parameters
                  data = subset(myDTmstate_age, trans == 5),
                  dist = "gengamma")
if(!exists("t5")) {
  t5<- flexsurvreg(Surv(Tstart2, Tstop2, status) ~
                     gender + imd + region + birthcohortONS, #doesnt fit with the anc parameters
                   data = subset(myDTmstate_age, trans == 5),
                   dist = "gengamma")
  print("t5 no anc")
}
print(t5)
AHnew_lc_bcONS <- list((t1), (t2), (t3), (t4), (t5))
qsave(AHnew_lc_bcONS , output_path("AHnew_lc_bcONS_ancDEV.q"), nthreads = 10)

t6 <- flexsurvreg(Surv(Tstart2, Tstop2, status) ~
                    gender + imd + region + birthcohortONS,
                  anc = list(shape = ~ gender + imd),
                  data = subset(myDTmstate_age, trans == 6),
                  dist = "gompertz")
print(t6)
AHnew_lc_bcONS <- list((t1), (t2), (t3), (t4), (t5), (t6))
qsave(AHnew_lc_bcONS , output_path("AHnew_lc_bcONS_ancDEV.q"), nthreads = 10)

t7 <- flexsurvreg(Surv(Tstart2, Tstop2, status) ~
                    gender + imd + region + birthcohortONS,
                  anc = list(shape = ~ gender+ imd),
                  data = subset(myDTmstate_age, trans == 7),
                  dist = "gompertz")
print(t7)
AHnew_lc_bcONS <- list((t1), (t2), (t3), (t4), (t5), (t6), (t7))
qsave(AHnew_lc_bcONS , output_path("AHnew_lc_bcONS_ancDEV.q"), nthreads = 10)




# Saving the model fits as tables - 1 for Gengamma distribution & 1 for Gompertz distribution
gengam <- data.table(rn = c(  "mu", "sigma" ,                                "Q"       ,
                              "genderF"                          ,     "imd2"                        ,          "imd3"       ,
                              "imd4"                             ,     "imd5"                         ,         "regionSouth West"    ,
                              "regionSouth Central"             ,      "regionSouth East Coast"       ,         "regionWest Midlands"   ,
                              "regionEast Midlands"             ,      "regionEast of England"         ,        "regionNorth West"       ,
                              "regionYorkshire And The Humber"  ,      "regionNorth East"              ,        "birthcohortONS(1929,1934]" ,
                              "birthcohortONS(1934,1939]"       ,      "birthcohortONS(1939,1944]"     ,        "birthcohortONS(1944,1949]"  ,
                              "birthcohortONS(1949,1954]"       ,      "birthcohortONS(1954,1959]"     ,        "birthcohortONS(1959,1964]"   ,
                              "birthcohortONS(1964,1969]"       ,      "birthcohortONS(1969,1974]"     ,        "birthcohortONS(1974,1979]"    ,
                              "birthcohortONS(1979,1984]"       ,      "birthcohortONS(1984,1989]"     ,        "sigma(genderF)"                ,
                              "sigma(imd2)"                    ,       "sigma(imd3)"                   ,        "sigma(imd4)"                    ,
                              "sigma(imd5)"                   ,        "sigma(regionSouth West)"       ,        "sigma(regionSouth Central)"      ,
                              "sigma(regionSouth East Coast)"   ,      "sigma(regionWest Midlands)"    ,        "sigma(regionEast Midlands)"       ,
                              "sigma(regionEast of England)"   ,       "sigma(regionNorth West)"       ,        "sigma(regionYorkshire And The Humber)",
                              "sigma(regionNorth East)"         ,      "Q(genderF)"                    ,        "Q(imd2)"                 ,
                              "Q(imd3)"                         ,      "Q(imd4)"                       ,        "Q(imd5)"                  ,
                              "Q(regionSouth West)"             ,      "Q(regionSouth Central)"        ,        "Q(regionSouth East Coast)" ,
                              "Q(regionWest Midlands)"          ,      "Q(regionEast Midlands)"        ,        "Q(regionEast of England)"   ,
                              "Q(regionNorth West)"             ,      "Q(regionYorkshire And The Humber)",     "Q(regionNorth East)"))
gomp <- data.table(rn = c(  "shape", "rate" ,
                            "genderF"                          ,     "imd2"                        ,          "imd3"       ,
                            "imd4"                             ,     "imd5"                         ,         "regionSouth West"    ,
                            "regionSouth Central"             ,      "regionSouth East Coast"       ,         "regionWest Midlands"   ,
                            "regionEast Midlands"             ,      "regionEast of England"         ,        "regionNorth West"       ,
                            "regionYorkshire And The Humber"  ,      "regionNorth East"              ,        "birthcohortONS(1929,1934]" ,
                            "birthcohortONS(1934,1939]"       ,      "birthcohortONS(1939,1944]"     ,        "birthcohortONS(1944,1949]"  ,
                            "birthcohortONS(1949,1954]"       ,      "birthcohortONS(1954,1959]"     ,        "birthcohortONS(1959,1964]"   ,
                            "birthcohortONS(1964,1969]"       ,      "birthcohortONS(1969,1974]"     ,        "birthcohortONS(1974,1979]"    ,
                            "birthcohortONS(1979,1984]"       ,      "birthcohortONS(1984,1989]"     ,        "shape(genderF)"                ,
                            "shape(imd2)"                    ,       "shape(imd3)"                   ,        "shape(imd4)"                    ,
                            "shape(imd5)"                   ,        "shape(regionSouth West)"       ,        "shape(regionSouth Central)"      ,
                            "shape(regionSouth East Coast)"   ,      "shape(regionWest Midlands)"    ,        "shape(regionEast Midlands)"       ,
                            "shape(regionEast of England)"   ,       "shape(regionNorth West)"       ,        "shape(regionYorkshire And The Humber)",
                            "shape(regionNorth East)" ))
colnm_old <-  c("V1",  "2.5 %", "97.5 %")
colnm_new <-  c("AFT",  "LCI", "UCI")


for(i in 1:7){
  tmp <- as.data.table(round(cbind(exp(coef(AHnew_lc_bcONS[[i]])), exp(confint(AHnew_lc_bcONS[[i]]))), 3), keep.rownames = TRUE)
  setnames(tmp, colnm_old, paste0("t", i, "_", colnm_new))

  if(AHnew_lc_bcONS[[i]]$dlist$name == "gengamma"){
    gengam <- tmp[gengam, on = 'rn']}else{
      gomp <-  tmp[gomp, on = 'rn']
    }

}
fwrite(gengam, output_path("AHnew_lc_bcONS_anc_gengamtab.csv"), row.names = FALSE)
fwrite(gomp, output_path("AHnew_lc_bcONS_anc_gomptab.csv"), row.names = FALSE)

