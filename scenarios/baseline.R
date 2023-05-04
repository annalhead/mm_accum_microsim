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
library(scales)
library(cowplot)

#uncertainty intervals for output
prbl <- c( 0.5, 0.25, 0.75, 0.025, 0.975)

# blue palette
mypalette <- c("#A8DDB5" ,"#7BCCC4" ,"#4EB3D3" ,"#2B8CBE", "#08589E")
mypalette2 <- c("#CCEBC5", "#A8DDB5" ,"#7BCCC4" ,"#4EB3D3" ,"#2B8CBE", "#0868AC" , "#084081")

scenario <- "baseline"

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



runmodel <- FALSE

if(runmodel){

  # load look up table with durations
  #read in the lookup table
  ll <- read_ll()

  # Pop to be simulated
  patients <- read_patients()


  # The simulation & post-processing

  n_sim <- 500
  outstrata <- c("mc", "year", "gender", "age", "imd", "region")
  prevtab <- data.table()
  inctab <- data.table()
  cftab <- data.table()

  tt <- data.table(# Auxiliary table
    y = c(2019L:2049L),
    year = 2019L:2049L)

  for (i in 251:n_sim) {
    sim_output <- run_sim(patients, ll, 41L + i)
    sim_output[, mc := i]
    print(paste0("completed iteration #", i))
    # fwrite(sim_output, outpt_pth(paste0("RAW/", scenario, "_sim_output_mc",i,".csv.gz")))

    lifeyears <-
      sim_output[ageenter %in% c(30, 50, 65) &
                   init_year %in% c(2019, 2029, 2039, 2049),
                 .(
                   patient_id,
                   gender,
                   imd,
                   region,
                   yob,
                   ageenter,
                   init_year,
                   init_state,
                   healthy,
                   incCond,
                   bmm,
                   cmm,
                   max_age,
                   last_state
                 )][, mc := i]
    fwrite(lifeyears, outpt_pth(paste0(
      "lifeyears/", scenario, "lifeyears_mc", i, ".csv.gz"
    )))
    print(paste0("iteration #", i, "lifeyears saved"))


    pats <- mk_pnl(sim_output, tt)
    rm(sim_output)
    gc()
    #fwrite(pats, outpt_pth(paste0("panel/",scenario,"_panel_mc",i,".csv.gz")))
    #print(paste0("iteration #",i, "panel saved"))

    prevtab <- mk_prevtab(pats, outstrata)
    fn_state(prevtab)
    # agegroup5yfn(prevtabtmp) #will add this in later as it's just extra info to store
    gc()
    fwrite_safe(prevtab, outpt_pth(paste0(scenario, "_prvl.csv.gz")), append = TRUE)
    print(paste0("iteration #", i, "prevtab saved"))
    rm(prevtab)
    gc()

    inctab <- mk_inctab(pats, outstrata)
    fn_state(inctab)
    #agegroup5yfn(inctabtmp)
    gc()
    fwrite_safe(inctab, outpt_pth(paste0(scenario, "_incd.csv.gz")), append = TRUE)
    print(paste0("iteration #", i, "inctab saved"))
    rm(inctab)
    gc()

    cftab <- mk_cftab(pats, outstrata)
    fn_state_cf(cftab)
    gc()
    #agegroup5yfn(cftabtmp)

    cftab[is.na(N), N := 0]
    cftab[is.na(atrisk), atrisk := 0]
    gc()
    fwrite_safe(cftab, outpt_pth(paste0(scenario, "_cf.csv.gz")), append = T)
    print(paste0("iteration #", i, "cftab saved"))
    rm(cftab)
    gc()



    #cumulative years in state - making these so can just easily sum them
    #halving for incident years. Halving again for year die
    pats[, `:=` (
      healthy = ifelse(
        state_incCond == 0 & state_BMM == 0 & state_CMM == 0,
        1,
        ifelse(state_incCond == 1, 0.5, 0)
      ),
      state_incCond = state_incCond / 2,
      state_BMM = state_BMM / 2,
      state_CMM = state_CMM / 2
    )]
    pats[year == max_year, `:=` (
      healthy = healthy / 2,
      state_incCond = state_incCond / 2,
      state_BMM = state_BMM / 2,
      state_CMM = state_CMM / 2
    )]
    yrs_in_state_sum <- pats[, .(
      N = .N,
      Healthy = sum(healthy) ,
      IncCond = sum(state_incCond),
      BMM = sum(state_BMM),
      CMM = sum(state_CMM)
    ),
    keyby = .(year, gender, imd, region)]
    rm(pats)
    gc()
    setkey(yrs_in_state_sum, year, gender, imd, region)
    yrs_in_state_sum[, `:=`(
      Healthy_cml = cumsum(Healthy),
      IncCond_cml = cumsum(IncCond),
      BMM_cml = cumsum(BMM),
      CMM_cml = cumsum(CMM),
      N_cml = cumsum(N),
      Healthy_cml_pp = cumsum(Healthy / N),
      IncCond_cml_pp = cumsum(IncCond / N),
      BMM_cml_pp = cumsum(BMM / N),
      CMM_cml_pp = cumsum(CMM / N)
    ), keyby = .(gender, imd, region)][, mc := i]

    # #states aren't exclusive though ...
    # yrs_in_state_sum[, `:=`( Healthy_cml_prop = cumsum(Healthy_cml/N_cml),
    #                          IncCond_cml_prop = cumsum(IncCond_cml/N_cml),
    #                          BMM_cml_prop = cumsum(BMM_cml/N_cml),
    #                          CMM_cml_prop = cumsum(CMM_cml/N_cml)
    # ), keyby = .(gender, imd, region)]
    fwrite(yrs_in_state_sum, outpt_pth(paste0(
      "cumlyears/", scenario, "cumlyears_mc", i, ".csv.gz"
    )))
    print(paste0("iteration #", i, "yrs_in_state_sum saved"))
    rm(yrs_in_state_sum)
    gc()
  }



}
