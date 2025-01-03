#scenarios
library(data.table)
data.table::setDTthreads(3) #this is so that don't use all the processors
library(fst)
fst::threads_fst(3)
library(Rcpp)
library(dqrng)
library(ggplot2)
library(qs)
library(flexsurv)
library(RColorBrewer)
library(ggthemes)
library(viridis)
library(CKutils)
library(foreach)
library(doParallel)
registerDoParallel(cores=3)
library(scales)

scenario <- "scenario1"
#targeeted intervention on worst off
# reduce the gap between imd3 & imd4/5 by x% amount to prevent 13,000 deaths

n_sim <- 200

set_scenario <- FALSE

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



tt <- data.table( # Auxiliary table
  y = c(2020L:2050L),
  year = 2020L:2050L)

# Setting up the reduction amounts
prvl <- fread("output/baseline/baseline_prvl.csv.gz", header = T )[state == "BMM" & year == 2019]
agegroup5yfn(prvl)
prvl[, imd := factor(imd, c(1:5))]
esp <- mk_esp()
strata <- c("mc", "year", "gender", "imd", "agegrp5")
prvl_st <- fn_stnd(prvl, strata)
outstrata <- c("mc", "year","imd")
d <- prvl_st[, .(prvl = sum(N * wt_esp) / sum(popsize_wtd)*100), keyby = eval(outstrata)] #doing it this way keeps the variable name as popsize
d <- fn_mc(d, outstrata, prbl, "outcome", "prvl_") #this then gives the median and various uncertainty intervals as defined by prbl above
rm(prvl, prvl_st)
gc()
r1 <- 1
r2 <- 1
r3 <- 1
r4 <- d[imd == 4, `prvl_50%`]/d[imd == 3, `prvl_50%`]
# 1.03448
r5 <- d[imd == 5, `prvl_50%`]/d[imd == 3, `prvl_50%`]
# 1.072873

chngtab <- data.table(imd = factor(c(1:5)),
                      ratio = c(r1, r2, r3, r4, r5))

chngtab[, ratio2 := (ratio - 1)/5 + 1]


if(set_scenario){
  # Pop to be simulated - for this scenario only change imd 4&5,
  # so for finding the reduction value don't need to simulate everyone
  patients <- read_patients()

  patients <- patients[imd %in% c(4,5)]

  # The simulation setup
  outstrata <- c("mc", "year", "imd")
  outstrata_cml <- outstrata[!outstrata %in% ("year")]




  scenario_nm <- "baseline"
  baseline <- data.table()
  for (scenario in scenario_nm){
    baseline <- fread(outpt_pth(paste0(scenario,"_cf.csv.gz")))
    baseline[, scenario := scenario]
  }

  baseline[, imd := factor(imd)]
  baseline[, state := factor(state,
                             levels = c("Healthy", "IncCond", "BMM", "CMM", "Overall"))]
  setkey(baseline, year, imd)

  outstrata <- c("mc")
  d <- baseline[state == "Overall" & mc <= n_sim, .(N = sum(N)), keyby = outstrata]
  setkey(d, mc)

  deaths_overall <- d[, median(N)]

  outstrata <- c("mc")
  d <- baseline[state == "Overall" & mc <= n_sim & imd %in% c(4, 5), .(N = sum(N)), keyby = outstrata]
  setkey(d, mc)

  deaths_imd45 <- d[, median(N)]
  rm(baseline)
  gc()


  deaths_overall
  # 231354
  deaths_imd45
  # 91833

  # target is reducing overall deaths by 5%, which is
  deaths_overall * 0.03 #6938.82

  target <- deaths_imd45 - round(deaths_overall * 0.03) # 3% of the baseline deaths



  reduc <- 0.75
  j <- 0
  threshold <- 10
  diff <- 20000

  # Want to prevent *at least* 100,000 cases, but not more than 101,000 cases
  while (j < 10 && any(abs(diff) > threshold | diff > 0)) {
    j <- j + 1L

    ll <- read_ll()
    tcols <- names(ll)[7:13]


    #For imd 4 & 5 reducing the gap with imd3 by 50%
    ll_l <- melt.data.table(ll,
                            id.vars = c( "bc", "gender", "imd", "region", "quantile", "ageenter"),
                            measure.vars = tcols,
                            variable.name = "transition",
                            value.name = "time"
    )
    ll_l_imd3 <- ll_l[imd == "3"]
    ll_l <- ll_l[ll_l_imd3, on = c("bc", "gender", "region", "quantile", "ageenter", "transition"), imd3 := i.time]
    ll_l[, gap := ifelse(imd3 > time, imd3-time, 0)] #only want to reduce the gap where it exists
    ll_l[chngtab, on = 'imd', ratio := i.ratio2] #Add in the ratios for the gradient
    ll_l[imd == 4, new_time := gap * (reduc * ratio ) + time]
    ll_l[imd == 5, new_time := gap * (reduc * ratio ) + time]

    ll_l[imd %in% c(4,5), time := new_time][
      , c("imd3", "gap", "new_time", "ratio") := NULL]
     ll <- dcast(ll_l,
                formula = bc + gender + region + quantile + ageenter + imd ~ transition,
                value.var = "time"
    )
    rm(ll_l, ll_l_imd3)


    scenario_setup <- function(i){
      sim_output <- run_sim(patients, ll, 41L + i)
      sim_output[, mc := i]
      sim_output
    }
    output <- foreach(i = 1:n_sim, .combine = rbind) %dopar% scenario_setup(i)


    #only interested in comparing the ones have changed
    deaths <- output[imd %in% c(4:5)][year_max_year <= 2050, .(deaths = .N), keyby = mc][, median(deaths)]
    diff <- deaths - target

    if(abs(diff) > 20){
      if(diff <0){
        reduc <- reduc - 0.005 }
      else{
        reduc <- reduc + 0.005 }
    } else{
      if(abs(diff) > 10 | diff > 0){
        reduc <- reduc * (cases/target)^100
      }
    }


    print(paste0("diff ",diff, "cases, new reduc ratio ", reduc))
    "diff -0.5cases, new reduc ratio 0.35754709271749"
  }
}


reduc <- 0.735

# Run the whole thing for results ----
scenario <- "scenario1"
patients <- read_patients()
ll <- read_ll()
tcols <- names(ll)[7:13]


# For imd 4 & 5 reducing the gap with imd3
ll_l <- melt.data.table(ll,
                        id.vars = c( "bc", "gender", "imd", "region", "quantile", "ageenter"),
                        measure.vars = tcols,
                        variable.name = "transition",
                        value.name = "time"
                        )
ll_l_imd3 <- ll_l[imd == "3"]
ll_l <- ll_l[ll_l_imd3, on = c("bc", "gender", "region", "quantile", "ageenter", "transition"), imd3 := i.time]
ll_l[, gap := ifelse(imd3 > time, imd3-time, 0)] #only want to reduce the gap where it exists
ll_l[chngtab, on = 'imd', ratio := i.ratio2] #Add in the ratios for the gradient
ll_l[imd == 4 , new_time := gap * (reduc * ratio ) + time]
ll_l[imd == 5 , new_time := gap * (reduc * ratio ) + time]
ll_l[imd %in% c(4,5) , time := new_time][
  , c("imd3", "gap", "new_time", "ratio") := NULL]
ll <- dcast(ll_l,
            formula = bc + gender + region + quantile + ageenter + imd ~ transition,
            value.var = "time"
)
rm(ll_l, ll_l_imd3)




# Pop to be simulated
patients <- read_patients()

# The simulation
outstrata <- c("mc", "year", "gender", "age","imd", "region")

outstrata <- c("mc", "year", "gender", "age","imd", "region")

foreach(i = 1:n_sim, .combine = rbind) %dopar% process_results(i)
