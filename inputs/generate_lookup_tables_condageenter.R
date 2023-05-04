# Creating lookup tables conditional on age entry

# This file takes fitted survival analysis models, and uses the predict() function
# to predict time to transition for all combinations of covariates conditional
# on age at entry to state, and at deciles of the survival of the distribution



library(data.table)
data.table::setDTthreads(5) #this is so that don't use all the processors
library(fst)
fst::threads_fst(5)
library(qs)
library(flexsurv)
library(future)
library(splines)
options(future.fork.enable = TRUE) # TODO remove for production
options(future.rng.onMisuse = "ignore") # Remove false warning
plan(multicore(workers = 5))

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
           "/UoL/CPRD2019mm/forthesis/", x)


myDT <- read_fst(data_dir_PM("accumDTmstate_age.fst"), as.data.table = TRUE)[
  yob >= 1919 & yob <= 1986 & imd!= "" & region != ""]




lu <- CJ(gender = myDT$gender,
         imd = myDT$imd,
         region = myDT$region,
         birthcohortONS = myDT$birthcohortONS,
         unique = TRUE)  #
setkeyv(lu, sort(names(lu)))
lu[, ln := .I]


# generate distribution of duration in each state
gen_lu_tbl <-
  function(i,
           mdl_list, # a list of flexsurv models
           newdata, # unique combinations of simulant attributes
           j, #conditional on survival to age2 j
           granularity = gran, # of quantiles
           max_time = 110, # will replace Inf with this
           plot = FALSE) {
    for(j in c(0:82)){
      tt <- predict(
        mdl_list[[i]],
        newdata = newdata,
        type = "quantile",
        start = j, #conditional on survival to age2 j
        se.fit = TRUE,
        p = seq(0, 1, granularity)
      )
      assign(paste0("t", i), rbindlist(tt$.pred, idcol = "ln"))
      tt <- rbindlist(tt$.pred, idcol = "ln")
      tt <- tt[newdata, on = "ln"][, ln := NULL][, ageenter := j + 18]
      tt[.pred_quantile > max_time, `:=` (.pred_quantile = max_time, .std_error = 0)]
      if (plot)
        print(tt[, hist(.pred_quantile)])
      setnames(tt,
               c(".quantile", ".pred_quantile", ".std_error"),
               c("quantile", paste0("t", i), paste0("se", i)))
      setkeyv(tt, c(key(newdata), "quantile"))
      #invisible(tt) # the t# contains the duration in state # before jump. se# is the standard error of the duration
      bigtt <- rbind(bigtt,tt)
    }
    invisible(bigtt)
  }





# Load the survival analysis fits (from modelfitscript.R)
mdl_list <- qread(data_dir_PM("SAmodelfits/AHnew_lc_bcONS_ancDEV.q"), nthreads = 10)

gran <- 1e-1


bigtt <- data.table()
print("tt1 start")
Sys.time()
tt1 %<-% gen_lu_tbl(1, mdl_list, lu, j)
qsave(tt1, output_path("t1_DEV10se.qs") )
print("tt1 stop")
Sys.time()

print("tt2 start")
Sys.time()
tt2 %<-% gen_lu_tbl(2, mdl_list, lu, j)
qsave(tt2, output_path("t2_DEV10se.qs") )
print("tt2 stop")
Sys.time()

print("tt3 start")
Sys.time()
tt3 %<-% gen_lu_tbl(3, mdl_list, lu, j)
qsave(tt3, output_path("t3_DEV10se.qs") )
print("tt3 stop")
Sys.time()

print("tt4 start")
Sys.time()
tt4 %<-% gen_lu_tbl(4, mdl_list, lu, j)
qsave(tt4, output_path("t4_DEV10se.qs") )
print("tt4 stop")
Sys.time()

print("tt5 start")
Sys.time()
tt5 %<-% gen_lu_tbl(5, mdl_list, lu, j)
qsave(tt5, output_path("t5_DEV10se.qs") )
print("tt5 stop")
Sys.time()

print("tt6 start")
Sys.time()
tt6 %<-% gen_lu_tbl(6, mdl_list, lu, j)
qsave(tt6, output_path("t6_DEV10se.qs") )
print("tt6 stop")
Sys.time()

print("tt7 start")
Sys.time()
tt7 %<-% gen_lu_tbl(7, mdl_list, lu, j)
qsave(tt7, output_path("t7_DEV10se.qs") )
print("tt7 stop")
Sys.time()


ll <- tt1[tt2, on = c(key(lu), "quantile", "ageenter")][
  tt3, on = c(key(lu), "quantile", "ageenter")][
    tt4, on = c(key(lu), "quantile", "ageenter")][
      tt5, on = c(key(lu), "quantile", "ageenter")][
        tt6, on = c(key(lu), "quantile", "ageenter")][
          tt7, on = c(key(lu), "quantile", "ageenter")]
setcolorder(ll, c(key(lu), "ageenter", "quantile", paste0("t", 1:7), paste0("se", 1:7)))
write_fst(ll, output_path("state_duration_lookup_bcONS_DEV10se.fst") , compress = 80)
print("ll saved")



# Setting a max age: 110 - ONS cap
llcap <- copy(ll)

max_ttime <- 92 #this is 110 in age

ts <- names(llcap)[7:13]
for(i in ts){
  llcap[, paste0(i) :=  ifelse(get(i) > max_ttime, max_ttime, get(i))]
}

##need to change these into time to next transition
for(i in ts){
  setnames(llcap, paste0(i), "tt")
  llcap[, tt:= tt - ageenter + 18][
    tt < 0, tt := 0] #make sure no negative times
  setnames(llcap, "tt", paste0(i))
}
write_fst(llcap, output_path("state_duration_lookup_bcONS_DEV10se_cap110.fst"), compress = 80)


