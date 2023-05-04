# Validation and calibration 

# This file is for validation and calibration of the microsimulation against a 
# sub-set of CPRD data. 
# The resulting calibration factors are used in the final microsimulation version

# Validation and calibration is against a subset of CPRD data that has 10 years of follow-up

library(data.table)
data.table::setDTthreads(5) #this is so that don't use all the processors
library(fst)
fst::threads_fst(10)
library(Rcpp)
library(dqrng)
library(ggplot2)
library(qs)
library(flexsurv)
# library(future.apply)
# library(future)
# options(future.fork.enable = TRUE) # TODO remove for production
# options(future.rng.onMisuse = "ignore") # Remove false warning
#plan(multicore(workers = 14))

project_path <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/msmagetranstimeV2/", x)

input_path <- 
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/forthesis/", x)

data_dir_lookup <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/Dictionaries etc/", x)


data_dir_CPRD <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2019mm/Data May 2020/", x)

sourceCpp(project_path("/mmsim_mod_bc.cpp"), cacheDir = project_path("./.cache"))
source(project_path("/fn.R"))



# load look up table with durations
  lu <-   read_fst(project_path("state_duration_lookup_bcONS_DEV10se_cap110.fst"), as.data.table = T)
  lu[, birthcohortONS := droplevels(birthcohortONS)]
  setnames(lu, "birthcohortONS", "bc")
  setkey(lu,ageenter,gender, imd, region,  bc, quantile) # key on all cols except t1-7, s1-7 quantiles
  
  
lu[, rownum := .I]
directory <- lu[quantile == 0, rownum, keyby = key(lu)] # all quantile == 0 rows


cm <- readRDS(project_path("state_duration_corr.rds"))


#Calibration 
myDT <-
  read_fst(input_path("accumDT_atf.fst"), as.data.table = TRUE)

#With left censoring
myDT_lc <- myDT[ gender != "I" &
                     imd != "" &
                     region!= "" &
                     yob <= 1986 &
                     yob >= 1919,]
myDT_lc[, gender := droplevels(gender)]
myDT_lc[, imd := factor(imd)]
setnames(myDT_lc, "birthcohortONS", "bc")

myDT_lc[, startstate := factor(startstate,
                                 levels = c(1:4),
                                 labels = c("Healthy", "IncCond", "BMM", "CMM"))]

myDT_lc <- myDT_lc[(year(censordate) - yearenterstudy  >= 10 |
                          (year(censordate) - yearenterstudy  < 10 & death == 1)), ]
myDT_lc[year(censordate) - yearenterstudy < 10 & death == 1, state10y := "Death"]
myDT_lc[!state10y %in% "Death" & cmmyn == 1 & agecmm - ageenter < 10 , state10y := "CMM"]
myDT_lc[!state10y %in% c("Death", "CMM") &  bmmyn == 1  & agebmm - ageenter < 10 , state10y := "BMM"]
myDT_lc[!state10y %in% c("Death", "CMM", "BMM") & firstyn == 1 & agefirst - ageenter < 10 , state10y := "IncCond"]
myDT_lc[!state10y %in% c("Death", "CMM", "BMM", "IncCond") , state10y := "Healthy"]
myDT_lc[, state10y := factor(state10y,
                               levels = c("Healthy", "IncCond", "BMM", "CMM", "Death"),
                               labels = c("Healthy", "IncCond", "BMM", "CMM", "Death"))]
myDT_lcsim <- myDT_lc[,.(patid, gender ,imd , region  ,bc, yob,init_year = yearenterstudy, init_state = startstate, state10y)]
myDT_lcsim[, ageenter := init_year - yob]
setkey(myDT_lcsim,ageenter,gender, imd, region, bc, yob) # key on all cols of interest
myDT_lcsim[directory, on = setdiff(key(lu), "quantile"), rownum := i.rownum]

#grp_size <- 10
# myDT_lcsim <- rbindlist(rep(list(myDT_lcsim), grp_size), idcol = "mc")
truth <- myDT_lcsim[, .(prop_cprd =  .N / nrow(myDT_lcsim)), keyby =.(gender, state10y)]
n <- nrow(myDT_lcsim)
myDT_lcsim[, unique(init_state)]

rank_mtx <- generate_corr_unifs(n, cm, seed = 42L)
# rank_mtx[, 7] <- rank_mtx[, 7]/20 # Temp fix for t7 model

nam <-  c("healthy", "incCond", "bmm", "cmm", "last_state",
          "state10y_sim", "max_age", "max_year", "year_incCond")


#Without calibration
myDT_lcsim[, c("healthy", "incCond", "bmm", "cmm", "last_state") :=
             mmsim(.SD, lu, rank_mtx
             )]
postprocess(myDT_lcsim, TRUE)
bias <-
  myDT_lcsim[, .(prop_sim = .N / n), keyby = .(gender, state10y_sim)]
bias[truth, on = c("gender == gender", "state10y_sim ==state10y")][
  , absdiff := prop_sim - prop_cprd]


calib_thresh <- 0.0001
Mcf12 <- 1
Mcf34 <-  1
Mcf56 <-  1 
Mcf7  <-  1
Fcf12 <- 1# 
Fcf34 <- 1
Fcf56 <-  1
Fcf7  <-1
i <- 1L
bias <- 1
# calibrate
#while (i < 20 && any(abs(bias) > calib_thresh)) { #do these sequentially 
while (i < 50 && any(abs(bias[c(1,5)]) > calib_thresh)) {
  i <- i + 1L
  for (j in nam) {
    if (j %in% names(myDT_lcsim)) {
      set(myDT_lcsim, NULL, j, NULL)
    }
  }
  myDT_lcsim[, c("healthy", "incCond", "bmm", "cmm", "last_state") :=
                 mmsim(.SD, lu, rank_mtx, Mcalib_t12 = Mcf12, Mcalib_t34 = Mcf34,
                       Mcalib_t56 = Mcf56, Mcalib_t7 = Mcf7, Fcalib_t12 = Fcf12,
                       Fcalib_t34 = Fcf34, Fcalib_t56 = Fcf56, Fcalib_t7 = Fcf7
                 )]
  postprocess(myDT_lcsim, TRUE)
  #bias <-
  #  myDT_lcsim[, .N / n, keyby = .(gender, state10y_sim)]$V1[c(1:4, 6:9)] - truth$V1[c(1:4, 6:9)]
  bias <-
    myDT_lcsim[, .(prop_sim = .N / n), keyby = .(gender, state10y_sim)]
  bias <- bias[truth, on = c("gender == gender", "state10y_sim ==state10y")][c(1:4,6:9), prop_sim - prop_cprd]
  
  print(bias)
    Mcf12 <- fifelse(abs(bias[1]) > calib_thresh, Mcf12 * (1 - bias[1]), Mcf12)
   # Mcf34 <- fifelse(abs(bias[2]) > calib_thresh, Mcf34 * (1 - bias[2]), Mcf34)
   # Mcf56 <- fifelse(abs(bias[3]) > calib_thresh, Mcf56 * (1 - bias[3]), Mcf56)
   # Mcf7  <- fifelse(abs(bias[4]) > calib_thresh, Mcf7  * (1 - bias[4]), Mcf7 )
    Fcf12 <- fifelse(abs(bias[5]) > calib_thresh, Fcf12 * (1 - bias[5]), Fcf12)
   # Fcf34 <- fifelse(abs(bias[6]) > calib_thresh, Fcf34 * (1 - bias[6]), Fcf34)
   # Fcf56 <- fifelse(abs(bias[7]) > calib_thresh, Fcf56 * (1 - bias[7]), Fcf56)
   # Fcf7  <- fifelse(abs(bias[8]) > calib_thresh, Fcf7  * (1 - bias[8]), Fcf7 )
}

print(c(Mcf12, Mcf34, Mcf56, Mcf7, Fcf12, Fcf34, Fcf56, Fcf7 ))
# [1] 0.9819753 0.9553712 0.9491299 0.8168817 0.9720173 0.9750439 0.9505665 0.8628891

postprocess(myDT_lcsim, TRUE)
bias <-
  myDT_lcsim[, .(prop_sim = .N / n), keyby = .(gender, state10y_sim)]
bias <- bias[truth, on = c("gender == gender", "state10y_sim ==state10y")][
  , absdiff := prop_sim - prop_cprd]



# Validation
for (j in nam) {
  if (j %in% names(myDT_lcsim)) {
    set(myDT_lcsim, NULL, j, NULL)
  }
}
myDT_lcsim[, c("healthy", "incCond", "bmm", "cmm", "last_state") :=
               mmsim(.SD, lu, rank_mtx, Mcalib_t12 = Mcf12, Mcalib_t34 = Mcf34,
                     Mcalib_t56 = Mcf56, Mcalib_t7 = Mcf7, Fcalib_t12 = Fcf12,
                     Fcalib_t34 = Fcf34, Fcalib_t56 = Fcf56, Fcalib_t7 = Fcf7)]
postprocess(myDT_lcsim, TRUE)
postprocess(myDT_lcsim, FALSE)
myDT_lcsim[, hist(max_age)]
hist(myDT_lcsim[, .N / n, keyby = .(gender, state10y_sim)]$V1 -
       myDT_lcsim[, .N / n, keyby = .(gender, state10y)]$V1)
hist(myDT_lcsim[, .N / n, keyby = .(imd, state10y_sim)]$V1 -
       myDT_lcsim[, .N / n, keyby = .(imd, state10y)]$V1)
hist(myDT_lcsim[, .N / n, keyby = .(region, state10y_sim)]$V1 -
       myDT_lcsim[, .N / n, keyby = .(region, state10y)]$V1)
hist(myDT_lcsim[, .N / n, keyby = .(yob, state10y_sim)]$V1 -
       myDT_lcsim[, .N / n, keyby = .(yob, state10y)]$V1)



