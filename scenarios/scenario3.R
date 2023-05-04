#scenarios
library(data.table)
data.table::setDTthreads(10) #this is so that don't use all the processors
library(fst)
fst::threads_fst(10)
library(Rcpp)
library(dqrng)
library(ggplot2)
library(qs)
library(flexsurv)
library(RColorBrewer)
library(ggthemes)
library(viridis)
library(CKutils)

scenario <- "scenario3"
# redistributive policy
# reduce gap between 2/3 to 1 by 50% and then 4/5 to 3 by 50%

n_sim <- 500


sim_path <-
  function(x = character(0))
    paste0("./simulation/", x)

outpt_pth <-
  function(x = character(0))
    paste0("./output/",scenario,"/", x)

# Load all the extra functions etc.
source(sim_path("/fn.R"))

# Load the model
sourceCpp(sim_path("/mmsim_mod_bc.cpp"), cacheDir = sim_path("./.cache"))


if(!dir.exists(outpt_pth())) dir.create(outpt_pth(), recursive = TRUE)
if(!dir.exists(outpt_pth("lifeyears/"))) dir.create(outpt_pth("lifeyears/"), recursive = TRUE)
if(!dir.exists(outpt_pth("cumlyears/"))) dir.create(outpt_pth("cumlyears/"), recursive = TRUE)


# load look up table with durations
#ll <- read_fst( input_pth("state_duration_lookup_bcONS_DEV10_cap130.fst"), as.data.table = T)
ll <- read_ll()

tcols <- names(ll)[7:13]



# reduce gap between 2/3 to 1 by 50% and then 4/5 to 3 by 50%
ll_l <- melt.data.table(ll,
                        id.vars = c( "bc", "gender", "imd", "region", "quantile", "ageenter"),
                        measure.vars = tcols,
                        variable.name = "transition",
                        value.name = "time"
                        )
#moving 2 & 3 towards 1
ll_l_imd1 <- ll_l[imd == "1"]
ll_l <- ll_l[ll_l_imd1, on = c("bc", "gender", "region", "quantile", "ageenter", "transition"), imd1 := i.time]
ll_l[, gap := ifelse(imd1 > time, imd1- time, 0)] #only want to reduce the gap where it exists
ll_l[, new_time := gap * 0.5 + time]
ll_l[imd %in% c(2,3), time := new_time][
  , c("imd1", "gap", "new_time") := NULL]
#moving 4 & 5 towards 3
ll_l_imd3 <- ll_l[imd == "3"]
ll_l <- ll_l[ll_l_imd3, on = c("bc", "gender", "region", "quantile", "ageenter", "transition"), imd3 := i.time]
ll_l[, gap := ifelse(imd3 > time, imd3- time, 0)] #only want to reduce the gap where it exists
ll_l[, new_time := gap * 0.5 + time]
ll_l[imd %in% c(4,5), time := new_time][
  , c("imd3", "gap", "new_time") := NULL]

# ll_w <- dcast(ll_l, #if wanted to look at different imds side by side
#               formula = bc + gender + region + quantile + ageenter + transition ~ imd,
#               value.var = "time"
#               )
ll <- dcast(ll_l,
            formula = bc + gender + region + quantile + ageenter + imd ~ transition,
            value.var = "time"
            )
rm(ll_l, ll_l_imd1, ll_l_imd3)




# Pop to be simulated
# Pop to be simulated
patients <- read_patients()

# The simulation
outstrata <- c("mc", "year", "gender", "age","imd", "region")


tt <- data.table( # Auxiliary table
  y = c(2019L:2049L),
  year = 2019L:2049L)

# The simulation
for (i in 1:n_sim){
  sim_output <- run_sim(patients, ll, 41L + i)
  sim_output[, mc := i]
  print(paste0("completed iteration #",i))
  # fwrite(sim_output, outpt_pth(paste0("RAW/", scenario, "_sim_output_mc",i,".csv.gz")))

  lifeyears <- sim_output[ageenter %in% c(30,50,65) & init_year %in% c(2019, 2029, 2039, 2049),
                          .(patient_id, gender, imd, region, yob, ageenter, init_year, init_state,
                            healthy, incCond, bmm, cmm, max_age, last_state)][, mc := i]
  fwrite(lifeyears, outpt_pth(paste0("lifeyears/", scenario,"lifeyears_mc",i,".csv.gz")))
  print(paste0("iteration #",i, "lifeyears saved"))


  pats <- mk_pnl(sim_output, tt)
  rm(sim_output)
  gc()
  #fwrite(pats, outpt_pth(paste0("panel/",scenario,"_panel_mc",i,".csv.gz")))
  #print(paste0("iteration #",i, "panel saved"))

  prevtab <- mk_prevtab(pats, outstrata)
  fn_state(prevtab)
  # agegroup5yfn(prevtabtmp) #will add this in later as it's just extra info to store
  gc()
  fwrite_safe(prevtab, outpt_pth(paste0(scenario,"_prvl.csv.gz")), append = TRUE)
  print(paste0("iteration #",i, "prevtab saved"))
  rm(prevtab)
  gc()

  inctab <- mk_inctab(pats, outstrata)
  fn_state(inctab)
  #agegroup5yfn(inctabtmp)
  gc()
  fwrite_safe(inctab, outpt_pth(paste0(scenario,"_incd.csv.gz")), append = TRUE)
  print(paste0("iteration #",i, "inctab saved"))
  rm(inctab)
  gc()

  cftab <- mk_cftab(pats, outstrata)
  fn_state_cf(cftab)
  gc()
  #agegroup5yfn(cftabtmp)

  cftab[is.na(N), N := 0]
  cftab[is.na(atrisk), atrisk := 0]
  gc()
  fwrite_safe(cftab, outpt_pth(paste0(scenario,"_cf.csv.gz")), append = T)
  print(paste0("iteration #",i, "cftab saved"))
  rm(cftab)
  gc()



  #cumulative years in state - making these so can just easily sum them
  #halving for incident years. Halving again for year die
  pats[, `:=` (healthy = ifelse(state_incCond == 0 & state_BMM == 0 & state_CMM == 0, 1,
                                ifelse(state_incCond == 1, 0.5, 0)),
               state_incCond = state_incCond/2,
               state_BMM = state_BMM/2,
               state_CMM = state_CMM/2)]
  pats[year == max_year, `:=` (healthy = healthy/2,
                               state_incCond = state_incCond/2,
                               state_BMM = state_BMM/2,
                               state_CMM = state_CMM/2)]
  yrs_in_state_sum <- pats[, .(N = .N,
                               Healthy = sum(healthy ) ,
                               IncCond = sum(state_incCond),
                               BMM = sum(state_BMM),
                               CMM = sum(state_CMM ) ),
                           keyby = .(year, gender, imd, region)]
  rm(pats)
  gc()
  setkey(yrs_in_state_sum, year, gender, imd, region)
  yrs_in_state_sum[, `:=`( Healthy_cml = cumsum(Healthy),
                           IncCond_cml = cumsum(IncCond),
                           BMM_cml = cumsum(BMM),
                           CMM_cml = cumsum(CMM),
                           N_cml = cumsum(N),
                           Healthy_cml_pp = cumsum(Healthy/N),
                           IncCond_cml_pp = cumsum(IncCond/N),
                           BMM_cml_pp = cumsum(BMM/N),
                           CMM_cml_pp = cumsum(CMM/N)
  ), keyby = .(gender, imd, region)][, mc:= i]


  fwrite(yrs_in_state_sum, outpt_pth(paste0("cumlyears/", scenario,"cumlyears_mc",i,".csv.gz")))
  print(paste0("iteration #",i, "yrs_in_state_sum saved"))
  rm(yrs_in_state_sum)
  gc()

}
