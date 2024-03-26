# Scenarios - baseline with parameter uncertainty

# This version of the baseline uses the standard errors for each estimated transition and
# for each iteration, a normal distribution is used to inject parameter uncertainty into the
# transition time estimates

library(data.table)
data.table::setDTthreads(5) #this is so that don't use all the processors
library(fst)
fst::threads_fst(5)
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

scenario <- "baseline_param"

#uncertainty intervals for output
prbl <- c( 0.5, 0.25, 0.75, 0.025, 0.975)

# blue palette
mypalette <- c("#A8DDB5" ,"#7BCCC4" ,"#4EB3D3" ,"#2B8CBE", "#08589E")
mypalette2 <- c("#CCEBC5", "#A8DDB5" ,"#7BCCC4" ,"#4EB3D3" ,"#2B8CBE", "#0868AC" , "#084081")


theme_set(new = theme_few())
theme_update(axis.text = element_text(size = 10),
             axis.title = element_text(size = 10,
                                       margin = margin(t = 10, r = 10, b = 10, l = 10)),
             legend.position = "bottom",
             legend.text=element_text(size=10),
             legend.title=element_text(size=10)
             #plot.title = element_text(hjust = 0.5),
)



sim_path <-
  function(x = character(0))
    paste0("./simulation/", x)

outpt_pth <-
  function(x = character(0))
    paste0("./output/",scenario,"/", x)

outpt_pth_summary <-
  function(x = character(0))
    paste0("./output/",scenario,"/summary/", x)

# Load all the extra functions etc.
source(sim_path("/fn.R"))

# Load the model
sourceCpp(sim_path("/mmsim_mod_bc.cpp"), cacheDir = sim_path("./.cache"))


if(!dir.exists(outpt_pth())) dir.create(outpt_pth(), recursive = TRUE)
if(!dir.exists(outpt_pth("lifeyears/"))) dir.create(outpt_pth("lifeyears/"), recursive = TRUE)
if(!dir.exists(outpt_pth("cumlyears/"))) dir.create(outpt_pth("cumlyears/"), recursive = TRUE)
if(!dir.exists(outpt_pth("summary/"))) dir.create(outpt_pth("summary/"), recursive = TRUE)


# Flags for different parts
runmodel <- FALSE
runtables <- FALSE
runplots <- FALSE
plotsave <- FALSE
runtables <- FALSE
runONScalib <- TRUE
runONScalib_alt <- FALSE
ineq <- FALSE

if(runmodel){
# load look up table with durations
  #read in the lookup table
  lu <- read_ll()
  ts <- names(lu)[7:13]

  # Pop to be simulated
  patients <- read_patients()


# The simulation & post-processing

n_sim <-100
outstrata <- c("mc", "year", "gender", "age","imd", "region")
prevtab <- data.table()
inctab <- data.table()
cftab <- data.table()

tt <- data.table( # Auxiliary table
  y = c(2019L:2050L),
  year = 2019L:2050L)

for (i in 401:n_sim){
  ll <- copy(lu)
  ll[, `:=` (t1 = rnorm(.N, t1, se1),
             t2 = rnorm(.N, t2, se2),
             t3 = rnorm(.N, t3, se3),
             t4 = rnorm(.N, t4, se4),
             t5 = rnorm(.N, t5, se5),
             t6 = rnorm(.N, t6, se6),
             t7 = rnorm(.N, t7, se7))][,
               c("se1", "se2", "se3", "se4", "se5", "se6", "se7") := NULL]
  for(j in ts){
    setnames(ll, paste0(j), "tt")
    ll[ageenter + tt > 110, tt:= 110 - ageenter][
      tt < 0, tt := 0] #make sure no negative times
    setnames(ll, "tt", paste0(j))
  }

  sim_output <- run_sim(patients, ll, 41L + i)
  rm(ll)
  gc()
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
yrs_in_state_sum <- pats[new_pid != 1, # Don't want to include people in the first year before start the simulation
                          .(N = .N,
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

# #states aren't exclusive though ...
# yrs_in_state_sum[, `:=`( Healthy_cml_prop = cumsum(Healthy_cml/N_cml),
#                          IncCond_cml_prop = cumsum(IncCond_cml/N_cml),
#                          BMM_cml_prop = cumsum(BMM_cml/N_cml),
#                          CMM_cml_prop = cumsum(CMM_cml/N_cml)
# ), keyby = .(gender, imd, region)]
fwrite(yrs_in_state_sum, outpt_pth(paste0("cumlyears/", scenario,"cumlyears_mc",i,".csv.gz")))
print(paste0("iteration #",i, "yrs_in_state_sum saved"))
rm(yrs_in_state_sum)
gc()
}



}

# Calibrating to the ONS population projections
if(runONScalib){
  prevtab <- fread(outpt_pth(paste0(scenario,"_prvl.csv.gz")), header = TRUE)
  # Make the ONScalib table
  tot_num <- prevtab[state == "IncCond", #only need 1 state as the atrisk is all the same
                     .(N = sum(atrisk)), keyby = .(mc, year, gender, age)]
  tot_num[, age2 := ifelse(age >= 105, "105+", age)][, age := NULL]

  # Want to calculate the weights based on the average total popsize,
  # otherwise all the scenarios will have the same total pop
  tot_num <- tot_num[, .(N = mean(N)), keyby = .(year, gender, age2)]

  popproject <- fread(project_path("ONSpopproj.csv"), header = T)

  tot_num[popproject, on = c("year", "gender", "age2"), ons_N := i.N]
  tot_num[, N := N * 100]
  tot_num[, ONSwt := ons_N/N]
  # tot_num[age2 == 30, ONSwt := 1] # not weighting the 30 year olds as have their numbers already
  fwrite(tot_num, outpt_pth(paste0(scenario,"_ONSwt.csv.gz")))

  prevtab[, age2 := ifelse(age >= 105, "105+", age)]
  prevtab[tot_num, on = c("year", "gender", "age2"), ONSwt := i.ONSwt]
  #  prevtab[, `:=` (N_wtd = N* ONSwt, atrisk_wtd = atrisk * ONSwt )]
  fwrite(prevtab, outpt_pth(paste0(scenario,"_prvl.csv.gz")), row.names  = FALSE)
  rm(prevtab)

  inctab <- fread(outpt_pth(paste0(scenario,"_incd.csv.gz")), header = TRUE)
  inctab[, age2 := ifelse(age >= 105, "105+", age)]
  inctab[tot_num, on = c("year", "gender", "age2"), ONSwt := i.ONSwt]
  #  prevtab[, `:=` (N_wtd = N* ONSwt, atrisk_wtd = atrisk * ONSwt )]
  fwrite(inctab, outpt_pth(paste0(scenario,"_incd.csv.gz")), row.names  = FALSE)
  rm(inctab)

  cftab <- fread(outpt_pth(paste0(scenario,"_cf.csv.gz")), header = TRUE)
  cftab[, age2 := ifelse(age >= 105, "105+", age)]
  cftab[tot_num, on = c("year", "gender", "age2"), ONSwt := i.ONSwt]
  #  prevtab[, `:=` (N_wtd = N* ONSwt, atrisk_wtd = atrisk * ONSwt )]
  fwrite(cftab, outpt_pth(paste0(scenario,"_cf.csv.gz")), row.names  = FALSE)
  rm(cftab)
}

# Calibrating to alternative ONS projections (low and high migration)
if(runONScalib_alt){
  prevtab <- fread(outpt_pth(paste0(scenario,"_prvl.csv.gz")), header = TRUE)
  # Make the ONScalib table
  agegroup5yfn(prevtab)
  prevtab[, Ages := ifelse(age < 90, substr(agegrp5, 1,5),
                           ifelse(age >= 90 & age < 95, "90-94",
                                  ifelse(age >= 95 & age < 100, "95-99", "100 & over")))]
  tot_num <- prevtab[state == "IncCond", #only need 1 state as the atrisk is all the same
                     .(N = sum(atrisk)), keyby = .(mc, year, gender, Ages)]
  # tot_num[, age2 := ifelse(age >= 105, "105+", age)][, age := NULL]

  # Want to calculate the weights based on the average total popsize,
  # otherwise all the scenarios will have the same total pop
  tot_num <- tot_num[, .(N = mean(N)), keyby = .(year, gender, Ages)]

  popproject <- fread(project_path("ONSpopproj_highmig.csv"), header = T)

  tot_num[popproject, on = c("year", "gender", "Ages"), ons_N := i.N]
  tot_num[, N := N * 100]
  tot_num[, ONSwt := ons_N/N]
  # tot_num[age2 == 30, ONSwt := 1] # not weighting the 30 year olds as have their numbers already
  fwrite(tot_num, outpt_pth(paste0(scenario,"_ONSwt_highmig.csv.gz")))

  prevtab[tot_num, on = c("year", "gender", "Ages"), ONSwt_highmig := i.ONSwt]
  #  prevtab[, `:=` (N_wtd = N* ONSwt, atrisk_wtd = atrisk * ONSwt )]
  #fwrite(prevtab, outpt_pth(paste0(scenario,"_prvl.csv.gz")), row.names  = FALSE)

  inctab <- fread(outpt_pth(paste0(scenario,"_incd.csv.gz")), header = TRUE)
  agegroup5yfn(inctab)
  inctab[, Ages := ifelse(age < 90, substr(agegrp5, 1,5),
                          ifelse(age >= 90 & age < 95, "90-94",
                                 ifelse(age >= 95 & age < 100, "95-99", "100 & over")))]
  inctab[tot_num, on = c("year", "gender", "Ages"), ONSwt_highmig := i.ONSwt]
  # fwrite(inctab, outpt_pth(paste0(scenario,"_incd.csv.gz")), row.names  = FALSE)




  # low migration
  tot_num <- prevtab[state == "IncCond", #only need 1 state as the atrisk is all the same
                     .(N = sum(atrisk)), keyby = .(mc, year, gender, Ages)]
  # tot_num[, age2 := ifelse(age >= 105, "105+", age)][, age := NULL]

  # Want to calculate the weights based on the average total popsize,
  # otherwise all the scenarios will have the same total pop
  tot_num <- tot_num[, .(N = mean(N)), keyby = .(year, gender, Ages)]

  popproject <- fread(project_path("ONSpopproj_lowmig.csv"), header = T)
  tot_num[popproject, on = c("year", "gender", "Ages"), ons_N := i.N]
  tot_num[, N := N * 100]
  tot_num[, ONSwt := ons_N/N]
  # tot_num[age2 == 30, ONSwt := 1] # not weighting the 30 year olds as have their numbers already
  fwrite(tot_num, outpt_pth(paste0(scenario,"_ONSwt_lowmig.csv.gz")))
  prevtab[tot_num, on = c("year", "gender", "Ages"), ONSwt_lowhmig := i.ONSwt]
  #  prevtab[, `:=` (N_wtd = N* ONSwt, atrisk_wtd = atrisk * ONSwt )]
  fwrite(prevtab, outpt_pth(paste0(scenario,"_prvl.csv.gz")), row.names  = FALSE)
  rm(prevtab)

  inctab[tot_num, on = c("year", "gender", "Ages"), ONSwt_lowmig := i.ONSwt]
  fwrite(inctab, outpt_pth(paste0(scenario,"_incd.csv.gz")), row.names  = FALSE)


  #
  # cftab <- fread(outpt_pth(paste0(scenario,"_cf.csv.gz")), header = TRUE)
  # cftab[, age2 := ifelse(age >= 105, "105+", age)]
  # cftab[tot_num, on = c("year", "gender", "age2"), ONSwt := i.ONSwt]
  # #  prevtab[, `:=` (N_wtd = N* ONSwt, atrisk_wtd = atrisk * ONSwt )]
  # fwrite(cftab, outpt_pth(paste0(scenario,"_cf.csv.gz")), row.names  = FALSE)
  # rm(cftab)
}


# Making big summary table for select years
if(runtables){
  yrs <- c(2019,2029,2039,2049)
  bigtab <- fread(outpt_pth(paste0(scenario,"_prvl.csv.gz")), header = TRUE)[
    year %in% yrs]


  bigtab[, imd := factor(imd)]
  healthy <- bigtab[state == "IncCond"]
  healthy[, `:=` (state = "Healthy", N = atrisk - N)]
  bigtab <- rbind(bigtab, healthy)

  bigtab[, agegrp_big := factor(ifelse(age < 65, "wrk_age", "over65"))]

  bigtab[, state := factor(state,
                           levels = c("Healthy", "IncCond", "BMM", "CMM"))]



  #adding in 10yr agegroups - there is definitely a better way to do this
  agegroup5yfn(bigtab)
  bigtab[, agegrp10 := agegrp5]
  bigtab[agegrp5 %in% c("30-34", "35-39"), agegrp10 := "30-39"]
  bigtab[agegrp5 %in% c("40-44", "45-49"), agegrp10 := "40-49"]
  bigtab[agegrp5 %in% c("50-54", "55-59"), agegrp10 := "50-59"]
  bigtab[agegrp5 %in% c("60-64", "65-69"), agegrp10 := "60-69"]
  bigtab[agegrp5 %in% c("70-74", "75-79"), agegrp10 := "70-79"]
  bigtab[agegrp5 %in% c("80-84", "85-89"), agegrp10 := "80-89"]
  bigtab[, agegrp10 := factor(agegrp10)]

  #bigtab[, outcome := N/atrisk *100]
  fwrite(bigtab, outpt_pth_summary("prvltab.csv.gz"), row.names = FALSE)

  # Standardising to esp
  esp <- mk_esp()
  strata <- c("mc", "year", "gender", "imd", "agegrp_big", "agegrp5", "agegrp10", "state")

  #prev
  bigtabsum <- fn_stnd(bigtab, strata)


  #Proportion by state - unstandardised
  outstrata <- c("mc", "year", "state")
  d <- bigtabsum[,
                 .(prop = sum(N)/sum(atrisk)), keyby = outstrata]
  d <- fn_mc(d, outstrata, prbl, "outcome", "prop_") #this then gives the median and various uncertainty intervals as defined by prbl above
  #fwrite(d, outpt_pth_summary("prvl_inclsv_states.csv"), row.names = FALSE)


  outstrata <- c("mc", "year", "imd","state")
  d <- bigtabsum[,
                 .(prop = sum(N)/sum(atrisk)), keyby = outstrata]
  d <- fn_mc(d, outstrata, prbl, "outcome", "prop_") #this then gives the median and various uncertainty intervals as defined by prbl above
  d[, state := factor(state,
                      levels = c( "Healthy", "IncCond", "BMM", "CMM" ),
                      labels = c("Healthy", "1 condition", "Basic\nmultimorbidity",
                                  "Complex\nmultimorbidity"))]
  d[, year := factor(year)]
  ggplot(d[state != "Healthy"],
         aes(x = year, y = `prop_50%` * 100, fill = imd)) +
    geom_col(position = "dodge") +
    facet_grid(rows = vars(state)) +
    ylim(0,90) +
    ylab("Crude prevalence (%)") +
    xlab("Year") +
    labs(fill = "IMD quintile") +
    guides(fill = guide_legend(reverse = FALSE)) +
    scale_fill_manual(values = mypalette) +
    theme(legend.position = "bottom")

  if(plotsave){
    #ggsave(filename = outpt_pth("simplots/prvl_overall_st_excl.png"))
    ggsave2(outpt_pth_summary("prvl_overall_crude_incl.svg") , scale = 0.7, width = 12, height = 8)

  }

  extraplot <- FALSE
  if(extraplot){
    bigtab2 <- fread(outpt_pth(paste0(scenario,"_prvl.csv.gz")), header = TRUE)[
      , .(mc, year, gender, imd, state, N, atrisk)
    ]
    bigtab2[, imd := factor(imd)]

    setkey(bigtab2, mc, year, gender, imd, state)


    outstrata <- c("mc", "year", "gender","imd","state")
    d <- bigtab2[,
                   .(prop = sum(N)/sum(atrisk)), keyby = outstrata]
    d <- fn_mc(d, outstrata, prbl, "outcome", "prop_") #this then gives the median and various uncertainty intervals as defined by prbl above
    d[, state := factor(state,
                        levels = c( "Healthy", "IncCond", "BMM", "CMM" ),
                        labels = c("Healthy", "1 condition", "Basic\nmultimorbidity",
                                   "Complex\nmultimorbidity"))]
    fn_sex_lbls(d)

    ggplot(d[state != "Healthy"],
           aes(x = year, y = `prop_50%` * 100, col = imd)) +
      geom_smooth(linewidth = 0.5, se = FALSE) +
      facet_grid(state ~ gender ) +
      ylim(0,100) +
      ylab("Crude prevalence (%)") +
      xlab("Year") +
      labs(col = "IMD quintile") +
      guides(col = guide_legend(reverse = FALSE)) +
      scale_color_manual(values = mypalette) +
      theme(legend.position = "bottom")

    if(plotsave){
      #ggsave(filename = outpt_pth("simplots/prvl_overall_st_excl.png"))
      ggsave2(outpt_pth_summary("prvl_overall_crude_incl_allyears.svg") , scale = 0.7, width = 12, height = 8)

    }

    #Looking with the different migration assumptions


    outstrata <- c("mc", "year","state")
    main <- bigtab2[,
                    .(prop = sum(N * ONSwt)/sum(atrisk * ONSwt)), keyby = outstrata]
    main <- fn_mc(main, outstrata, prbl, "outcome", "prop_") #this then gives the median and various uncertainty intervals as defined by prbl above

    low <- bigtab2[,
                   .(prop = sum(N * ONSwt_lowhmig)/sum(atrisk * ONSwt_lowhmig)), keyby = outstrata]
    low <- fn_mc(low, outstrata, prbl, "outcome", "prop_") #this then gives the median and various uncertainty intervals as defined by prbl above

    high <- bigtab2[,
                    .(prop = sum(N * ONSwt_highmig)/sum(atrisk * ONSwt_highmig)), keyby = outstrata]
    high <- fn_mc(high, outstrata, prbl, "outcome", "prop_") #this then gives the median and various uncertainty intervals as defined by prbl above

    main <- rbind(main[, assump := "main"], low[, assump := "low"], high[, assump := "high"])

    main[, state := factor(state,
                           levels = c( "Healthy", "IncCond", "BMM", "CMM" ),
                           labels = c("Healthy", "1 condition", "Basic\nmultimorbidity",
                                      "Complex\nmultimorbidity"))]

    main[, assump := factor(assump,
                            levels = c("main", "low", "high"),
                            labels = c("Main (2020)", "Low migration (2018)", "High migration (2018)"))]



    ggplot(main[state != "Healthy"],
           aes(x = year, y = `prop_50%` * 100, col = assump)) +
      geom_smooth(linewidth = 0.5, se = FALSE) +
      facet_grid( ~state  ) +
      ylim(0,100) +
      ylab("Crude prevalence (%)") +
      xlab("Year") +
      labs(col = "ONS projection variant") +
      guides(col = guide_legend(reverse = FALSE)) +
      #      scale_color_manual(values = mypalette) +
      theme(legend.position = "bottom")
    if(plotsave){
      #ggsave(filename = outpt_pth("simplots/prvl_overall_st_excl.png"))
      ggsave2(outpt_pth_summary("prvl_overall_crude_excl_migassump.svg") , scale = 0.7, width = 14, height = 8)

    }








    outstrata <- c("mc", "year","state")
    main <- bigtab2[,
                    .(prvl = sum(N * ONSwt)), keyby = outstrata]
    main <- fn_mc(main, outstrata, prbl, "outcome", "prvl_") #this then gives the median and various uncertainty intervals as defined by prbl above

    low <- bigtab2[,
                   .(prvl = sum(N * ONSwt_lowhmig)), keyby = outstrata]
    low <- fn_mc(low, outstrata, prbl, "outcome", "prvl_") #this then gives the median and various uncertainty intervals as defined by prbl above

    high <- bigtab2[,
                    .(prvl = sum(N * ONSwt_highmig)), keyby = outstrata]
    high <- fn_mc(high, outstrata, prbl, "outcome", "prvl_") #this then gives the median and various uncertainty intervals as defined by prbl above

    main <- rbind(main[, assump := "main"], low[, assump := "low"], high[, assump := "high"])

    main[, state := factor(state,
                           levels = c( "Healthy", "IncCond", "BMM", "CMM" ),
                           labels = c("Healthy", "1 condition", "Basic\nmultimorbidity",
                                      "Complex\nmultimorbidity"))]



    ggplot(main[state != "Healthy"],
           aes(x = year, y = `prvl_50%` * 100, col = assump)) +
      geom_smooth(linewidth = 0.5, se = FALSE) +
      facet_grid( ~state  ) +
      ylim(0,100) +
      ylab("Crude prevalence (%)") +
      xlab("Year") +
      labs(col = "IMD quintile") +
      guides(col = guide_legend(reverse = FALSE)) +
      scale_color_manual(values = mypalette) +
      theme(legend.position = "bottom")



  }



  # mutually exclusive health states

  prvlbystate <- bigtabsum[, .(N = sum(N), atrisk = sum(atrisk)),
                           by = .(mc, year, gender, imd, agegrp5 ,
                                  agegrp10,  agegrp_big, state)]

  prvlbystate <- dcast(prvlbystate,
                       mc + year + gender + imd + agegrp5 + agegrp10 + agegrp_big  + atrisk ~ state,
                       value.var = "N")

  prvlbystate[, Healthy := atrisk - IncCond][
    , IncCond := IncCond - BMM][
      , BMM := BMM - CMM]

  prvlbystate <- melt(prvlbystate,
                      id.vars = c("mc","year", "gender", "imd", "agegrp5", "agegrp10", "agegrp_big", "atrisk"),
                      measure.vars = c("Healthy", "IncCond", "BMM", "CMM"),
                      value.factor = TRUE,
                      value.name = "N",
                      variable.name = "state")
  prvlbystate[, imd := factor(imd)]

  # reverse the levels for plotting
  prvlbystate[, state := factor(state,
                                levels = c("CMM", "BMM", "IncCond", "Healthy"))]


  # Making a baseline only stacked plot - crude
  outstrata <- c("mc", "year", "imd", "state")
  d <- prvlbystate[,
                   .(prvl_st = sum(N)/ sum(atrisk)), keyby = outstrata]
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above
  d[, state := factor(state,
                      levels = c("CMM", "BMM", "IncCond", "Healthy"),
                      labels = c("Complex multimorbidity", "Basic multimorbidity",
                                 "1 condition", "Healthy"))]
  ggplot(d[year %in% c(2019, 2029, 2039,2049) ],
         aes(x = imd, y = `outcome_50%` * 100, fill = state)) +
    geom_col(position = "stack") +
    facet_grid( ~ year) +
    ylab("Crude prevalence (%)") +
    xlab("IMD quintile (1 = least deprived)") +
    labs(fill = "State") +
    guides(fill = guide_legend(reverse = TRUE)) +
    scale_fill_manual(values = mypalette) +
    theme(legend.position = "bottom")

  if(plotsave){
    #ggsave(filename = outpt_pth("simplots/prvl_overall_st_excl.png"))
    ggsave2(outpt_pth_summary("prvl_overall_crude_excl.svg") , scale = 0.7, width = 14, height = 8)

  }




  # Making a baseline only stacked plot - age-sex standardised
  prvlbystate <- fn_stnd(prvlbystate, strata)

  outstrata <- c("mc", "year", "imd", "state")
  d <- prvlbystate[,
                   .(prvl_st = sum(N * wt_esp)/ sum(popsize_wtd)), keyby = outstrata]
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above
  d[, state := factor(state,
                      levels = c("CMM", "BMM", "IncCond", "Healthy"),
                      labels = c("Complex multimorbidity", "Basic multimorbidity",
                                 "1 condition", "Healthy"))]
  ggplot(d[year %in% c(2019, 2029, 2039,2049) ],
         aes(x = imd, y = `outcome_50%` * 100, fill = state)) +
    geom_col(position = "stack") +
    facet_grid( ~ year) +
    ylab("Standardised prevalence (%)") +
    xlab("IMD quintile (1 = least deprived)") +
    labs(fill = "State") +
    guides(fill = guide_legend(reverse = TRUE)) +
    scale_fill_manual(values = mypalette) +
    theme(legend.position = "bottom")

  if(plotsave){
    #ggsave(filename = outpt_pth("simplots/prvl_overall_st_excl.png"))
    ggsave2(outpt_pth_summary("prvl_overall_st_excl.svg") , scale = 0.7, width = 14, height = 8)

  }


  #Making this a table
  outstrata <- c("mc", "year", "imd", "state")

  d <- prvlbystate[,
                   .(N = sum(N) * 100), keyby = outstrata]
  d_w <- dcast(d, mc  + state + imd ~ year, value.var = "N")
  d <- fn_mc(d, outstrata, prbl, "outcome", "N_")
  fwrite(d, outpt_pth_summary("N_exclsv_states_imd.csv"), row.names = FALSE)
  d_w[, change := (`2049`/`2019`) * 100]
  outstrata <- c("mc", "imd", "state")
  d_w <- fn_mc(d_w[, .(mc, imd, state, change)], outstrata, prbl, "outcome", "change_")
  fwrite(d_w, outpt_pth_summary("N_exclsv_states_imd_change.csv"), row.names = FALSE)



  #Making this a table - not by imd
  outstrata <- c("mc", "year", "state")

  d <- prvlbystate[,
                   .(N = sum(N) * 100), keyby = outstrata]
  d_w <- dcast(d, mc  + state ~ year, value.var = "N")
  d <- fn_mc(d, outstrata, prbl, "outcome", "N_")
  fwrite(d, outpt_pth_summary("N_exclsv_states.csv"), row.names = FALSE)
  d_w[, change := (`2049`/`2019`) * 100]
  outstrata <- c("mc", "state")
  d_w <- fn_mc(d_w[, .(mc, state, change)], outstrata, prbl, "outcome", "change_")
  fwrite(d_w, outpt_pth_summary("N_exclsv_states_change.csv"), row.names = FALSE)




  #Making this a table - by agegroup
  outstrata <- c("mc", "year", "agegrp_big", "state")

  d <- prvlbystate[,
                   .(N = sum(N) * 100), keyby = outstrata]
  d_w <- dcast(d, mc + agegrp_big + state ~ year, value.var = "N")
  d <- fn_mc(d, outstrata, prbl, "outcome", "N_")
  fwrite(d, outpt_pth_summary("N_exclsv_states_agegrp.csv"), row.names = FALSE)
  d_w[, change := (`2049`/`2019`) * 100]
  outstrata <- c("mc", "agegrp_big", "state")
  d_w <- fn_mc(d_w[, .(mc, agegrp_big, state, change)], outstrata, prbl, "outcome", "change_")
  fwrite(d_w, outpt_pth_summary("N_exclsv_states_agegrp_change.csv"), row.names = FALSE)





  #Making this a tableby sex
  outstrata <- c("mc", "year", "gender", "state")

  d <- prvlbystate[,
                   .(N = sum(N) * 100), keyby = outstrata]
  d_w <- dcast(d, mc + gender + state ~ year, value.var = "N")
  d <- fn_mc(d, outstrata, prbl, "outcome", "N_")
  fwrite(d, outpt_pth_summary("N_exclsv_states_sex.csv"), row.names = FALSE)
  d_w[, change := (`2049`/`2019`) * 100]
  outstrata <- c("mc", "gender", "state")
  d_w <- fn_mc(d_w[, .(mc, gender, state, change)], outstrata, prbl, "outcome", "change_")
  fwrite(d_w, outpt_pth_summary("N_exclsv_states_sex_change.csv"), row.names = FALSE)



  #Making this a table of crude, exclusive prevalence
  outstrata <- c("mc", "year", "imd", "state")

  d <- prvlbystate[,
                   .(prvl = (sum(N)/sum(atrisk)) * 100), keyby = outstrata]
  d_w <- dcast(d, mc + imd + state ~ year, value.var = "prvl")
  d <- fn_mc(d, outstrata, prbl, "outcome", "prvl_")
  fwrite(d, outpt_pth_summary("prvl_exclsv_states_imd.csv"), row.names = FALSE)
  d_w[, change := (`2049`/`2019`) * 100]
  outstrata <- c("mc", "imd", "state")
  d_w <- fn_mc(d_w[, .(mc, imd, state, change)], outstrata, prbl, "outcome", "change_")
  fwrite(d_w, outpt_pth_summary("prvl_exclsv_states_imd_change.csv"), row.names = FALSE)

  #Making this a table - not by imd
  outstrata <- c("mc", "year", "state")

  d <- prvlbystate[,
                   .(prvl = (sum(N)/sum(atrisk)) * 100), keyby = outstrata]
  d_w <- dcast(d, mc + state ~ year, value.var = "prvl")
    d <- fn_mc(d, outstrata, prbl, "outcome", "N_")
  fwrite(d, outpt_pth_summary("prvl_exclsv_states.csv"), row.names = FALSE)
  d_w[, change := (`2049`/`2019`) * 100]
  outstrata <- c("mc", "state")
  d_w <- fn_mc(d_w[, .(mc, state, change)], outstrata, prbl, "outcome", "change_")
  fwrite(d_w, outpt_pth_summary("prvl_exclsv_states_change.csv"), row.names = FALSE)




  #Making this a table of crude, exclusive prevalence by sex
  outstrata <- c("mc", "year", "gender", "state")

  d <- prvlbystate[,
                   .(prvl = (sum(N)/sum(atrisk)) * 100), keyby = outstrata]
  d_w <- dcast(d, mc + gender + state ~ year, value.var = "prvl")
  d <- fn_mc(d, outstrata, prbl, "outcome", "prvl_")
  fwrite(d, outpt_pth_summary("prvl_exclsv_states_sex.csv"), row.names = FALSE)
  d_w[, change := (`2049`/`2019`) * 100]
  outstrata <- c("mc", "gender", "state")
  d_w <- fn_mc(d_w[, .(mc, gender, state, change)], outstrata, prbl, "outcome", "change_")
  fwrite(d_w, outpt_pth_summary("prvl_exclsv_states_gender_change.csv"), row.names = FALSE)




  #Making this a table of crude, exclusive prevalence by agegetp
  outstrata <- c("mc", "year", "agegrp_big", "state")

  d <- prvlbystate[,
                   .(prvl = (sum(N)/sum(atrisk)) * 100), keyby = outstrata]
  d_w <- dcast(d, mc + agegrp_big + state ~ year, value.var = "prvl")
  d <- fn_mc(d, outstrata, prbl, "outcome", "prvl_")
  fwrite(d, outpt_pth_summary("prvl_exclsv_states_agegrp.csv"), row.names = FALSE)
  d_w[, change := (`2049`/`2019`) * 100]
  outstrata <- c("mc", "agegrp_big", "state")
  d_w <- fn_mc(d_w[, .(mc, agegrp_big, state, change)], outstrata, prbl, "outcome", "change_")
  fwrite(d_w, outpt_pth_summary("prvl_exclsv_states_agegrp_big_change.csv"), row.names = FALSE)




  ### MM or no MM
  outstrata <- c("mc", "year", "state2")
  prvlbystate[, state2 := ifelse(state %in% c("BMM", "CMM"), "MM", "no MM")]

  d <- prvlbystate[,
                   .(N = sum(N) * 100), keyby = outstrata]
  d_w <- dcast(d, mc  + state2 ~ year, value.var = "N")
  d <- fn_mc(d, outstrata, prbl, "outcome", "N_")
  fwrite(d, outpt_pth_summary("N_exclsv_statesMMnoMM.csv"), row.names = FALSE)
  d_w[, change := (`2049`/`2019`) * 100]
  outstrata <- c("mc", "state2")
  d_w <- fn_mc(d_w[, .(mc, state2, change)], outstrata, prbl, "outcome", "change_")
  fwrite(d_w, outpt_pth_summary("N_exclsv_states_change.csv"), row.names = FALSE)







  ## Mean life expectancies

    all.files <- list.files(path = outpt_pth(paste0("lifeyears")) ,pattern = ".csv.gz", full.names = T)
    lifecourse <- rbindlist(lapply(all.files, fread))



  lifecourse[, imd := factor(imd)]

  outstrata <- c("mc", "imd" , "ageenter", "init_year")

  expect <- data.table()

  # LE
  d <- lifecourse[, .(le = mean(max_age)), keyby = outstrata]
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome") #this then gives the median and various uncertainty intervals as defined by prbl above
  #d <- dcast(d, scenario ~ imd, value.var = "le50%")
  expect <- rbind(expect, d[, type := "le"])


  #MM free LE
  d <- lifecourse[, .(le = mean(healthy + incCond)), keyby = outstrata]
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome") #this then gives the median and various uncertainty intervals as defined by prbl above
  #d <- dcast(d, scenario ~ imd, value.var = "mmfreele50%")
  expect <- rbind(expect, d[, type := "mmfree_le"])


  #CMM free LE
  d <- lifecourse[, .(le = mean(healthy + incCond + bmm)), keyby = outstrata]
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome") #this then gives the median and various uncertainty intervals as defined by prbl above
  #d <- dcast(d, scenario ~ imd, value.var = "cmmfreele50%")
  expect <- rbind(expect, d[, type := "cmmfree_le"])







  outstrata <- c("mc", "ageenter","init_year" )

  expect_overall <- data.table()
  # LE
  d <- lifecourse[, .(le = mean(max_age)), keyby = outstrata]
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome") #this then gives the median and various uncertainty intervals as defined by prbl above
  #d <- d[, .(scenario, overall = `le50%`)]
  expect_overall <- rbind(expect_overall, d[, `:=`(type = "le")])


  #MM free LE
  d <- lifecourse[, .(le = mean(healthy + incCond)), keyby = outstrata]
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome") #this then gives the median and various uncertainty intervals as defined by prbl above
  #d <- d[, .(scenario, overall = `mmfreele50%`)]
  expect_overall <- rbind(expect_overall, d[, `:=`(type = "mmfree_le")])


  #CMM free LE
  d <- lifecourse[, .(le = mean(healthy + incCond + bmm)), keyby = outstrata]
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome") #this then gives the median and various uncertainty intervals as defined by prbl above
  #d <- d[, .(scenario, overall = `cmmfreele50%`)]
  expect_overall <- rbind(expect_overall, d[, `:=`(type ="cmmfree_le")])


  expect <- rbind(expect_overall[, imd := "all"], expect)
  fwrite(expect, outpt_pth_summary("le.csv"), row.names = FALSE)

  #IMD Diff
  expect_diff <- data.table()

  instrata <- c("mc", "ageenter", "imd" ,"init_year" )
  outstrata <- c("mc", "ageenter","init_year" )

  # LE
  d <- lifecourse[, .(le = mean(max_age)), keyby = instrata]
  imd <- dcast(d[imd %in% c(1,5)], mc  + ageenter + init_year ~ imd, value.var = "le")
  imd[, `:=` (abs = `5`-`1`, rel = `5` / `1` *100)][, c("5", "1") := NULL]
  imd <- fn_mc(imd, outstrata, prbl, "outcome", "outcome")
  expect_diff <- rbind(expect_diff, imd[, `:=`(type ="le_diff")])

  # MMfree LE
  d <- lifecourse[, .(le = mean(healthy + incCond)), keyby = instrata]
  imd <- dcast(d[imd %in% c(1,5)], mc  + ageenter + init_year ~ imd, value.var = "le")
  imd[, `:=` (abs = `5`-`1`, rel = `5` / `1` *100)][, c("5", "1") := NULL]
  imd <- fn_mc(imd, outstrata, prbl, "outcome", "outcome")
  expect_diff <- rbind(expect_diff, imd[, `:=`(type ="mmfreele_diff")])

  # CMMfree LE
  d <- lifecourse[, .(le = mean(healthy + incCond + bmm)), keyby = instrata]
  imd <- dcast(d[imd %in% c(1,5)], mc  + ageenter + init_year ~ imd, value.var = "le")
  imd[, `:=` (abs = `5`-`1`, rel = `5` / `1` *100)][, c("5", "1") := NULL]
  imd <- fn_mc(imd, outstrata, prbl, "outcome", "outcome")
  expect_diff <- rbind(expect_diff, imd[, `:=`(type ="cmmfreele_diff")])

  fwrite(expect_diff, outpt_pth_summary("le_imddiff.csv"), row.names = FALSE)



  ## Mean life expectancies  by sex
  outstrata <- c("mc", "imd" , "ageenter", "init_year", "gender")
  expect_sex <- data.table()

  # LE
  d <- lifecourse[, .(le = mean(max_age)), keyby = outstrata]
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome") #this then gives the median and various uncertainty intervals as defined by prbl above
  #d <- dcast(d, scenario + gender ~ imd , value.var = "le50%")
  expect_sex <- rbind(expect_sex, d[, type := "le"])


  #MM free LE
  d <- lifecourse[, .(le = mean(healthy + incCond)), keyby = outstrata]
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome") #this then gives the median and various uncertainty intervals as defined by prbl above
  #d <- dcast(d, scenario + gender~ imd, value.var = "mmfreele50%")
  expect_sex <- rbind(expect_sex, d[, type := "mmfree_le"])


  #CMM free LE
  d <- lifecourse[, .(le = mean(healthy + incCond + bmm)), keyby = outstrata]
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome") #this then gives the median and various uncertainty intervals as defined by prbl above
  #d <- dcast(d, scenario + gender~ imd, value.var = "cmmfreele50%")
  expect_sex <- rbind(expect_sex, d[, type := "cmmfree_le"])



  outstrata <- c("mc" , "ageenter","init_year", "gender")

  expect_overall_sex <- data.table()
  # LE
  d <- lifecourse[, .(le = mean(max_age)), keyby = outstrata]
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome") #this then gives the median and various uncertainty intervals as defined by prbl above
  #d <- d[, .(scenario, gender, overall = `le50%`)]
  expect_overall_sex <- rbind(expect_overall_sex, d[, `:=`(type = "le")])


  #MM free LE
  d <- lifecourse[, .(le = mean(healthy + incCond)), keyby = outstrata]
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome") #this then gives the median and various uncertainty intervals as defined by prbl above
  #d <- d[, .(scenario, gender,overall = `mmfreele50%`)]
  expect_overall_sex <- rbind(expect_overall_sex, d[, `:=`(type = "mmfree_le")])


  #CMM free LE
  d <- lifecourse[, .(le = mean(healthy + incCond + bmm)), keyby = outstrata]
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome") #this then gives the median and various uncertainty intervals as defined by prbl above
  #d <- d[, .(scenario,gender, overall = `cmmfreele50%`)]
  expect_overall_sex <- rbind(expect_overall_sex, d[, `:=`(type ="cmmfree_le")])

  expect_sex <- rbind(expect_overall_sex[, imd := "all"], expect_sex )
  #expect_sex <- expect_overall_sex[expect_sex, on = c("scenario","gender","type")]
  fwrite(expect_sex, outpt_pth_summary("le_sex.csv"), row.names = FALSE)

  # expect_sex_l <- melt(expect_sex, id.vars = c("gender","type"),
  #                      measure.vars = c("overall", "1", "2", "3", "4", "5"),
  #                      variable.name = "imd",
  #                      value.name = "years"
  # )


  ## Mean life expectancies  by sex
  outstrata <- c("mc" , "ageenter", "init_year", "gender")
  expect_sex_diff <- data.table()

  # LE
  d <- lifecourse[, .(le = mean(max_age)), keyby = outstrata]
  d <- dcast(d, mc + ageenter + init_year ~ gender)
  d[, `:=` (absdiff = `F`-`M`, reldiff = `F`/`M`)]
  outstrata <- c("mc" , "ageenter", "init_year")
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome") #this then gives the median and various uncertainty intervals as defined by prbl above
  #d <- dcast(d, scenario + gender ~ imd , value.var = "le50%")
  expect_sex_diff <- rbind(expect_sex_diff, d[, type := "le"])


  #MM free LE
  outstrata <- c("mc" , "ageenter", "init_year", "gender")
    d <- lifecourse[, .(le = mean(healthy + incCond)), keyby = outstrata]
  d <- dcast(d, mc + ageenter + init_year ~ gender)
  d[, `:=` (absdiff = `F`-`M`, reldiff = `F`/`M`)]
  outstrata <- c("mc" , "ageenter", "init_year")
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome") #this then gives the median and various uncertainty intervals as defined by prbl above
  #d <- dcast(d, scenario + gender~ imd, value.var = "mmfreele50%")
  expect_sex_diff <- rbind(expect_sex_diff, d[, type := "mmfree_le"])


  #CMM free LE
  outstrata <- c("mc" , "ageenter", "init_year", "gender")
  d <- lifecourse[, .(le = mean(healthy + incCond + bmm)), keyby = outstrata]
  d <- dcast(d, mc + ageenter + init_year ~ gender)
  d[, `:=` (absdiff = `F`-`M`, reldiff = `F`/`M`)]
  outstrata <- c("mc" , "ageenter", "init_year")
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome") #this then gives the median and various uncertainty intervals as defined by prbl above
  #d <- dcast(d, scenario + gender~ imd, value.var = "cmmfreele50%")
  expect_sex_diff <- rbind(expect_sex_diff, d[, type := "cmmfree_le"])
  fwrite(expect_sex_diff, outpt_pth_summary("le_sex_diff.csv"), row.names = FALSE)

  #IMD Diff by sex
  expect_diff <- data.table()

  instrata <- c("mc", "ageenter", "gender","imd" ,"init_year" )
  outstrata <- c("mc", "ageenter", "gender", "init_year" )

  # LE
  d <- lifecourse[, .(le = mean(max_age)), keyby = instrata]
  imd <- dcast(d[imd %in% c(1,5)], mc  + ageenter + gender + init_year ~ imd, value.var = "le")
  imd[, `:=` (abs = `5`-`1`, rel = `5` / `1` *100)][, c("5", "1") := NULL]
  imd <- fn_mc(imd, outstrata, prbl, "outcome", "outcome")
  expect_diff <- rbind(expect_diff, imd[, `:=`(type ="le_diff")])

  # MMfree LE
  d <- lifecourse[, .(le = mean(healthy + incCond)), keyby = instrata]
  imd <- dcast(d[imd %in% c(1,5)], mc  + ageenter + gender + init_year ~ imd, value.var = "le")
  imd[, `:=` (abs = `5`-`1`, rel = `5` / `1` *100)][, c("5", "1") := NULL]
  imd <- fn_mc(imd, outstrata, prbl, "outcome", "outcome")
  expect_diff <- rbind(expect_diff, imd[, `:=`(type ="mmfreele_diff")])

  # CMMfree LE
  d <- lifecourse[, .(le = mean(healthy + incCond + bmm)), keyby = instrata]
  imd <- dcast(d[imd %in% c(1,5)], mc  + ageenter + gender + init_year ~ imd, value.var = "le")
  imd[, `:=` (abs = `5`-`1`, rel = `5` / `1` *100)][, c("5", "1") := NULL]
  imd <- fn_mc(imd, outstrata, prbl, "outcome", "outcome")
  expect_diff <- rbind(expect_diff, imd[, `:=`(type ="cmmfreele_diff")])

  fwrite(expect_diff, outpt_pth_summary("le_sex_imddiff.csv"), row.names = FALSE)







  ### What about cumulative incident cases?


  bigtab_inc <- data.table()
    tmptab <- fread(outpt_pth(paste0(scenario,"_incd.csv.gz")), header = TRUE)
    setkey(tmptab, mc, year, gender, age, region, state)
    agegroup5yfn(tmptab)
    agegroup10yfn_b(tmptab)
    tmptab[, agegrp_big := ifelse(age <= 65, "wrk_age","over65") ]
    tmptab <- tmptab[, .(N = sum(N), atrisk = sum(atrisk)),
                     keyby = .(mc, year, gender, imd, agegrp5, agegrp10, agegrp_big, state)]
    setkey(tmptab, mc, year, gender, imd, agegrp5, agegrp10, agegrp_big, state)
    tmptab[, cumN := cumsum(N), keyby = .(mc, gender, imd, agegrp5, agegrp10, agegrp_big, state)]
   # tmptab <- tmptab[year == 2049]
    tmptab[, scenario := scenario]
    bigtab_inc <- rbind(bigtab_inc, tmptab[, scenario := scenario])
    rm(tmptab)
    gc()

  fwrite(bigtab_inc, outpt_pth_summary(paste0("incd_sum.csv.gz")), )
  bigtab_inc[, imd := factor(imd)]



  #Standardised
  # # Standardising to esp
  esp <- mk_esp()
  strata <- c("mc","year", "gender", "imd", "agegrp5", "state")


  #inc
  inctabsum <- fn_stnd(bigtab_inc, strata)

  #inc
  outstrata <- c("mc", "year", "gender","imd", "state")
  fn_sex_lbls(inctabsum)

    d <- inctabsum[, .(incd = sum(N * wt_esp) / sum(popsize_wtd)*10000), keyby = eval(outstrata)] #doing it this way keeps the variable name as popsize
  d <- fn_mc(d, outstrata, prbl, "outcome", "incd_") #this then gives the median and various uncertainty intervals as defined by prbl above

  d[, state := factor(state,
                       levels = c("1 condition", "Basic multimorbidity", "Complex multimorbidity"),
                       labels = c("1 condition", "Basic\nmultimorbidity", "Complex\nmultimorbidity"))]
  ggplot(d, aes(x = year, y = `incd_50%`, colour = imd)) +
    geom_smooth(
      linewidth = 0.5,
      se = FALSE
    ) +
    facet_grid(state ~ gender) +
    expand_limits(y = 0) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Age-standardised incidence (per 10,000)") +
    labs(colour = "IMD quintile")  +
    guides(colour = guide_legend(reverse = FALSE), nrow = 1) +
    scale_colour_manual(values = mypalette) +
    theme(legend.position = "bottom")
  if(plotsave){
    ggsave2(filename = outpt_pth_summary("st_inc_imd.svg"),  scale = 0.7, width = 12, height = 8)
    }



  bigtab_inc <- bigtab_inc[year == 2049]

  setkey(bigtab_inc, mc, gender, imd, agegrp_big, state)

  #By IMD
  bigtab_inc_sum <- bigtab_inc[, .(cumN = sum(cumN)*100), keyby = .(mc, imd, state)]
  setkey(bigtab_inc_sum, mc, imd, state)
  bigtab_inc_sum_w <- dcast(bigtab_inc_sum, mc + state ~ imd, value.var = "cumN")
  bigtab_inc_sum_w[, `:=` (absdiff = `5` - `1`, reldiff = `5`/`1`)]
  outstrata <- c("mc", "state")
  bigtab_inc_sum_w <- fn_mc( bigtab_inc_sum_w  , outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above
  bigtab_inc_sum_w <- bigtab_inc_sum_w[outcome %in% c("absdiff", "reldiff")]
  #fwrite(bigtab_inc_sum_w, outpt_pth_summary("cumlcases_imd_diff.csv"), row.names = FALSE)

  outstrata <- c("mc", "imd", "state")
  d <- fn_mc(bigtab_inc_sum, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above

  tmp <- bigtab_inc[, .(cumN = sum(cumN)*100), keyby = .(mc, state)]
  setkey(tmp, mc, state)
  outstrata <- c("mc", "state")
  tmp <- fn_mc(tmp, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above
  d <- rbind(d, tmp[, imd := "all"])
    #fwrite(d, outpt_pth_summary("cumlcases_imd.csv"), row.names = FALSE)








  #By IMD & agegrp
  d <- bigtab_inc[, .(cumN = sum(cumN)*100), keyby = .(mc, imd, agegrp_big, state)]
  setkey(d, mc, imd, agegrp_big, state)

  d_w <- dcast(d, mc + agegrp_big + state ~ imd, value.var = "cumN")
  d_w <- d_w[, `:=` (absdiff = `5` - `1`, reldiff = `5`/`1`)][, .(mc, agegrp_big, state, absdiff, reldiff)]
  outstrata <- c("mc", "agegrp_big","state")
  d_w <- fn_mc( d_w  , outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above
  #fwrite(d_w, outpt_pth_summary("cumlcases_imd_age_diff.csv"), row.names = FALSE)

  outstrata <- c("mc", "imd","agegrp_big", "state")
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above

  tmp <- bigtab_inc[, .(cumN = sum(cumN)*100), keyby = .(mc, agegrp_big, state)]
  setkey(tmp, mc, agegrp_big, state)
  outstrata <- c("mc", "agegrp_big","state")
  tmp <- fn_mc(tmp, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above
  d <- rbind(d, tmp[, imd := "all"])

  #fwrite(d, outpt_pth_summary("cumlcases_imd_agegrp.csv"), row.names = FALSE)


  #By IMD & agegrp
  bigtab_inc_sum <- bigtab_inc[, .(cumN = sum(cumN)*100), keyby = .(mc, imd, state)]
  setkey(bigtab_inc_sum, mc, imd, state)


  outstrata <- c("mc", "imd", "state")
  d <- fn_mc(bigtab_inc_sum, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above

  tmp <- bigtab_inc[, .(cumN = sum(cumN)*100), keyby = .(mc, state)]
  setkey(tmp, mc, state)
  outstrata <- c("mc", "state")
  tmp <- fn_mc(tmp, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above
  d <- rbind(d, tmp[, imd := "all"])

  #fwrite(d, outpt_pth_summary("cumlcases_imd.csv"), row.names = FALSE)


  ### deaths

  bigtab_cf <- data.table()
  tmptab <- fread(outpt_pth(paste0(scenario,"_cf.csv.gz")), header = TRUE)
  setkey(tmptab, mc, year, gender, age, region, state)
  agegroup5yfn(tmptab)
  agegroup10yfn_b(tmptab)
  tmptab[, agegrp_big := ifelse(age <= 65, "wrk_age","over65") ]
  tmptab <- tmptab[, .(N = sum(N), atrisk = sum(atrisk)),
                   keyby = .(mc, year, gender, imd, agegrp5, agegrp10, agegrp_big, state)]
  setkey(tmptab, mc, year, gender, imd, agegrp5, agegrp10, agegrp_big, state)
  tmptab[, cumN := cumsum(N), keyby = .(mc, gender, imd, agegrp5, agegrp10, agegrp_big, state)]
  # tmptab <- tmptab[year == 2049]
  tmptab[, scenario := scenario]
  bigtab_cf <- rbind(bigtab_cf, tmptab[, scenario := scenario])
  rm(tmptab)
  gc()

  fwrite(bigtab_cf, outpt_pth_summary(paste0("cf_sum.csv.gz")), )
  bigtab_cf[, imd := factor(imd)]



  #Standardised
  # # Standardising to esp
  esp <- mk_esp()
  strata <- c("mc","year", "gender", "imd", "agegrp5", "state")


  cftabsum <- fn_stnd(bigtab_cf, strata)


  outstrata <- c("mc", "year", "gender","imd", "state")
  fn_sex_lbls(cftabsum)
  d <- cftabsum[, .(cf = sum(N * wt_esp) / sum(popsize_wtd)*10000), keyby = eval(outstrata)] #doing it this way keeps the variable name as popsize
  d <- fn_mc(d, outstrata, prbl, "outcome", "cf_") #this then gives the median and various uncertainty intervals as defined by prbl above
  d[, state := factor(state,
                      levels = c("Overall", "Healthy","IncCond", "BMM", "CMM"),
                      labels = c("Overall", "Healthy","1 condition", "Basic\nmultimorbidity", "Complex\nmultimorbidity") )]

  ggplot(d, aes(x = year, y = `cf_50%`, colour = imd)) +
    geom_smooth(
      linewidth = 0.5,
      se = FALSE
    ) +
    facet_grid(state ~ gender) +
    expand_limits(y = 0) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Age-standardised case-fatality") +
    labs(colour = "IMD quintile")  +
    guides(colour = guide_legend(reverse = FALSE), nrow = 1) +
    scale_colour_manual(values = mypalette) +
    theme(legend.position = "bottom")
  if(plotsave){
    ggsave2(filename = outpt_pth_summary("st_cf_imd.svg"),  scale = 0.7, width = 12, height = 10)
  }


  bigtab_cf <- bigtab_cf[year == 2049]

  setkey(bigtab_cf, mc, gender, imd, agegrp_big)

  #By IMD
  bigtab_cf_sum <- bigtab_cf[state == "Overall", .(cumN = sum(cumN)*100), keyby = .(mc, imd)]
  setkey(bigtab_cf_sum, mc, imd)
  bigtab_cf_sum_w <- dcast(bigtab_cf_sum, mc  ~ imd, value.var = "cumN")
  bigtab_cf_sum_w[, `:=` (absdiff = `5` - `1`, reldiff = `5`/`1`)]
  outstrata <- c("mc")
  bigtab_cf_sum_w <- fn_mc( bigtab_cf_sum_w  , outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above
  bigtab_cf_sum_w <- bigtab_cf_sum_w[outcome %in% c("absdiff", "reldiff")]
  #fwrite(bigtab_cf_sum_w, outpt_pth_summary("cumldeaths_imd_diff.csv"), row.names = FALSE)

  outstrata <- c("mc", "imd")
  d <- fn_mc(bigtab_cf_sum, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above

  tmp <- bigtab_cf[, .(cumN = sum(cumN)*100), keyby = .(mc)]
  setkey(tmp, mc)
  outstrata <- c("mc")
  tmp <- fn_mc(tmp, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above
  d <- rbind(d, tmp[, imd := "all"])
  #fwrite(d, outpt_pth_summary("cumldeaths_imd.csv"), row.names = FALSE)


  #By IMD & agegrp
  d <- bigtab_cf[state == "Overall", .(cumN = sum(cumN)*100), keyby = .(mc, imd, agegrp_big)]
  setkey(d, mc, imd, agegrp_big)

  d_w <- dcast(d, mc + agegrp_big  ~ imd, value.var = "cumN")
  d_w <- d_w[, `:=` (absdiff = `5` - `1`, reldiff = `5`/`1`)][, .(mc, agegrp_big, absdiff, reldiff)]
  outstrata <- c("mc", "agegrp_big")
  d_w <- fn_mc( d_w  , outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above
  #fwrite(d_w, outpt_pth_summary("cumldeaths_imd_age_diff.csv"), row.names = FALSE)

  outstrata <- c("mc", "imd","agegrp_big")
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above

  tmp <- bigtab_cf[state == "Overall", .(cumN = sum(cumN)*100), keyby = .(mc, agegrp_big)]
  setkey(tmp, mc, agegrp_big)
  outstrata <- c("mc", "agegrp_big")
  tmp <- fn_mc(tmp, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above
  d <- rbind(d, tmp[, imd := "all"])

  #fwrite(d, outpt_pth_summary("cumldeaths_imd_agegrp.csv"), row.names = FALSE)


  #By IMD & agegrp
  bigtab_cf_sum <- bigtab_cf[state == "Overall", .(cumN = sum(cumN)*100), keyby = .(mc, imd)]
  setkey(bigtab_cf_sum, mc, imd)


  outstrata <- c("mc", "imd")
  d <- fn_mc(bigtab_cf_sum, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above

  tmp <- bigtab_cf[state == "Overall", .(cumN = sum(cumN)*100), keyby = .(mc)]
  setkey(tmp, mc)
  outstrata <- c("mc")
  tmp <- fn_mc(tmp, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above
  d <- rbind(d, tmp[, imd := "all"])

  #fwrite(d, outpt_pth_summary("cumldeaths_imd.csv"), row.names = FALSE)



  #Prop deaths under/over 65
  outstrata <- c("mc", "agegrp_big")
  d <- bigtab_cf[state == "Overall", .(cumN = sum(cumN)*100), keyby = outstrata]
  d <- dcast(d,mc ~ agegrp_big, value.var = "cumN")
  d[, total := wrk_age + over65][, `:=` (over65 = over65/total * 100,
                                         wrk_age = wrk_age/total * 100)][,
                                                                         total := NULL]
  outstrata <- c("mc")
  d <- fn_mc(d, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above

   outstrata <- c("mc", "imd","agegrp_big")
  e <- bigtab_cf[state == "Overall", .(cumN = sum(cumN)*100), keyby = outstrata]
  e <- dcast(e,mc + imd ~ agegrp_big, value.var = "cumN")
  e[, total := wrk_age + over65][, `:=` (over65 = over65/total * 100,
                                         wrk_age = wrk_age/total * 100)][,
                                                                         total := NULL]
  outstrata <- c("mc", "imd")
  e <- fn_mc(e, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above
  d <- rbind(e, d[, imd := "all"])
  #fwrite(d, outpt_pth_summary("cumldeaths_propbyage.csv"), row.names = FALSE)




}


#### Needs checking below here ####

if(ineq){

  plotsave <- FALSE

  bigtab <- fread(outpt_pth_summary("prvltab.csv.gz"), header = TRUE)[
    year == 2049 & state != "Healthy"]



  library(PHEindicatormethods) #Need to do this with ridit scores so don't have to use the old package
  require(tidyverse)
  #the agegrp10 col keeps going weird
  bigtab[, agegrp10 := agegrp5]
  bigtab[agegrp5 %in% c("30-34 years", "35-39 years"), agegrp10 := "30-39 years"]
  bigtab[agegrp5 %in% c("40-44 years  ", "45-49 years"), agegrp10 := "40-49 years"]
  bigtab[agegrp5 %in% c("50-54 years", "55-59 years"), agegrp10 := "50-59 years"]
  bigtab[agegrp5 %in% c("60-64 years", "65-69 years"), agegrp10 := "60-69 years"]
  bigtab[agegrp5 %in% c("70-74 years", "75-79 years"), agegrp10 := "70-79 years"]
  bigtab[agegrp5 %in% c("80-84 years", "85-89 years"), agegrp10 := "80-89 years"]
  bigtab[, agegrp10 := factor(agegrp10)]
  #Don't want healthy in this becasue it confuses things - higher % is better
  data <- bigtab[year %in% c(2049) &
                   state != "Healthy", .(
                     N = sum(N),
                     atrisk = sum(atrisk),
                     outcome = sum(N) / sum(atrisk) * 100
                   ),
                 by = c("mc", "imd", "agegrp10", "state")][
                   , se := sqrt((outcome * (100 - outcome)) / atrisk)]
  prvl_sii <-
    as.data.table(
      phe_sii(
        group_by(data, mc, agegrp10, state),
        imd,
        atrisk,
        value_type = 1,
        # for rates
        value = outcome,
        se = se,
        confidence = 0.95,
        rii = TRUE,
        type = "standard"
      )
    )
  fwrite(prvl_sii, outpt_pth_summary("prvl_sii.csv"), row.names = FALSE)


  outstrata <- c("mc", "agegrp10", "state")
  d <- fn_mc(prvl_sii, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above

  fn_state_lbls(d)

  mypalette <- brewer.pal(9, name = "GnBu")[2:9]

  ggplot(d[ outcome == "sii"], aes(x = state, y = `outcome_50%`, fill = agegrp10)) +
    geom_col(position = "dodge") +
    ylab("Absolute difference in prevalence") +
    ylim(-10,15) +
    xlab("State") +
    geom_hline(aes(yintercept = 0)) +#add x-axis at 0
    # ggtitle("Absolute inequalities in prevalence in 2049 by age-group") +
    labs(fill = "10-year\n age-group") +
    scale_fill_manual(values = mypalette) +
    theme(legend.position = "bottom", axis.title.x = element_text(size = 12)) +
    guides(fill = guide_legend(nrow = 2))



  if(plotsave){
    ggsave2(outpt_pth_summary("prvl_sii_2049.svg") , scale = 0.6, width = 14, height = 9)

  }

}


