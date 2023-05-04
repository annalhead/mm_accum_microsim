# Functions for the simulation model and post-processing

outpt_pth_plt <-
  function(x = character(0))
    paste0("./output/",scenario,"/simplots/", x)


input_path <-
  function(x = character(0))
    paste0("./inputs/", x)


library(ggplot2)
library(ggthemes)
theme_set(theme_few())
theme_update(axis.text.x = element_text(size = 8),
             legend.background = element_rect(
               fill = '#FBFBFB',
               linewidth = 0.5,
               linetype = "solid"
             ),
             legend.title = element_text(size = 10))


# Ensures that when fwrite appends file colnames of file to be written, match
# those already in the file
fwrite_safe <- function(x,
                        file = "",
                        append = TRUE,
                        ...) {
  if (append) {
    if (file.exists(file)) {
      col_names_disk <- names(fread(file, nrows = 0))
      col_names_file <- names(x)
      col_names <- outersect(col_names_disk, col_names_file)
      if (length(col_names) > 0)
        x[, (col_names) := NA]
      setcolorder(x, col_names_disk)
    }
  }
  fwrite(x, file, append, ...)
}

# Reading in & preparing the simulant population
read_patients <- function(){
  patients <- read_fst(input_path("cohtab_ONS2019_indage_1pc_30proj.fst"), as.data.table = T)
  patients[,  c( "cohortnum", "cohortnum2") := NULL]
  patients <- patients[yob >= 1919]
  patients[, init_state := factor(state_id,
                                  levels = c(1:4),
                                  labels = c("Healthy", "IncCond", "BMM", "CMM"))]
  patients[, birthcohort := droplevels(birthcohort)]
  patients[yob > 1989, `:=` (birthcohort = "(1984,1989]",#set all the new cohorts to the latest birthcohort have info on
                             age = 30)] #and don't start them until they are 30
  patients <- patients[,.(patient_id, gender ,imd , region  ,yob, ageenter = age, bc = birthcohort, init_year , init_state)]
  patients[, unique(init_state)]
  patients[, gender := factor(gender,
                              levels = c("M", "F"),
                              labels = c("M", "F"))]
  patients[, region := droplevels(region)]
}

# Reading in the lookup table
read_ll <- function(){
  ll <- read_fst(input_path("state_duration_lookup_bcONS_DEV10se_cap110.fst"), as.data.table = T)
  ll[, birthcohortONS := droplevels(birthcohortONS)]
  setnames(ll, "birthcohortONS", "bc")
  setkey(ll,ageenter,gender, imd, region,  bc, quantile) # key on all cols except t1-7, s1-7 quantiles
}


# Making a matrix of random numbers
generate_corr_unifs <- function(n, M, seed, verbose = FALSE) {
  # from http://comisef.wikidot.com/tutorial:correlateduniformvariates
  stopifnot(is.matrix(M))
  # Check that matrix is semi-positive definite
  # stopifnot(min(eigen(M, only.values = TRUE)$values) >= 0)

  if (verbose) M_original <- M

  # adjust correlations for uniforms
  for (i in seq_len(dim(M)[[1L]])) {
    for (j in seq_len(dim(M)[[2L]])) {
      if (i != j) {
        M[i, j] <- 2 * sin(pi * M[i, j] / 6)
        M[j, i] <- 2 * sin(pi * M[j, i] / 6)
      }
    }
  }

  # X <- matrix(dqrnorm(n * dim(M)[[2]]), n)
  X <- parallel_random_matrix(n, dim(M)[[2]], seed, 2L)
  # colnames(X) <- colnames(M)

  # induce correlation, check correlations
  Y <- pnorm(X %*% chol(M))

  if (verbose) {
    message(paste0("Mean square error is: ", signif(sum((cor(Y) - M_original) ^
                                                          2), 3)))
  }

  return(Y)
}




# Post-process of the model outputs
postprocess <- function(pop, y10 = FALSE) {
  #Need to remove the -1s otherwise it messes up things
  for(j in c("incCond", "bmm", "cmm")){
    set(pop, NULL, j, fifelse(pop[[j]] < 0L, 0L, pop[[j]]))
  }
  rm(j)
  if (y10) {
    pop[, state10y_sim := factor(

      fifelse( # init_state != "Healthy"
        healthy >= 10,
        "Healthy",
        fifelse(
          healthy + incCond >= 10,
          "IncCond",
          fifelse(
            healthy + incCond + bmm >= 10,
            "BMM",
            fifelse(healthy + incCond + bmm + cmm >= 10, "CMM", "Death")
          )
        )
      )
      ,
      levels = c("Healthy", "IncCond", "BMM", "CMM", "Death"),
      labels = c("Healthy", "IncCond", "BMM", "CMM", "Death")
    )]
  } else {
    pop[, `:=` (
      rownum = NULL,
      max_age =
        ageenter + healthy + incCond + bmm + cmm
    )]
    pop[, max_year := yob + max_age][
      , year_max_year := as.integer(max_year)]#[, max_year := NULL]
    pop[init_state == "Healthy" & #got to start healthy
          last_state != "Healthy", #and not die healthy
        year_incCond := ageenter + healthy + yob]
    pop[init_state %in% c("Healthy", "IncCond") & #got to start before BMM
          !last_state %in% c("Healthy", "IncCond"),
        year_BMM := ageenter + healthy + yob + incCond]
    pop[init_state %in% c("Healthy", "IncCond", "BMM") &
          !last_state %in% c("Healthy", "IncCond", "BMM"), #got to start before CMM
        year_CMM :=
          ageenter + healthy + yob + incCond + bmm
    ]
  }
}




# The simulation
# Other model inputs.
# Based on this ll: input_path("state_duration_lookup_bcONS_DEV10.fst")
# load correlation matrix of durations
cm <- readRDS(input_path("state_duration_corr.rds"))



# Calibration parameters ----
# Calibration parameters derived from the vldtn_clbrtn.R file
Mcf12 <- 0.9819753
Mcf34 <-  0.9553712
Mcf56 <-  0.9491299
Mcf7  <-  0.8168817
Fcf12 <- 0.9720173
Fcf34 <- 0.9750439
Fcf56 <-  0.9505665
Fcf7  <-0.8628891


run_sim <- function(myDT, lu_tab, startseed){
  setkey(lu_tab,ageenter,gender, imd, region,  bc, quantile) # key on all cols except t1-7, s1-7 quantiles
  lu_tab[, rownum := .I]
  directory <- lu_tab[quantile == 0, rownum, keyby = key(lu_tab)] # the first row in ll the attributes combinations appears
  n <- nrow(myDT)
  simtab <- copy(myDT)
  simtab[directory, on = setdiff(key(lu_tab), "quantile"), rownum := i.rownum]
  rank_mtx <- generate_corr_unifs(n, cm, seed = startseed) #will need to remove the seed


  simtab[, c("healthy", "incCond", "bmm", "cmm", "last_state") :=
           mmsim(.SD, lu_tab, rank_mtx, Mcalib_t12 = Mcf12, Mcalib_t34 = Mcf34,
                 Mcalib_t56 = Mcf56, Mcalib_t7 = Mcf7, Fcalib_t12 = Fcf12,
                 Fcalib_t34 = Fcf34, Fcalib_t56 = Fcf56, Fcalib_t7 = Fcf7)]
  postprocess(simtab, TRUE)
  postprocess(simtab, FALSE)
}




# Make a panel
mk_pnl <- function(data, tt){
  data[, c("ageenter", "healthy",  "incCond", "bmm", "cmm",
           "last_state", "state10y_sim",  "max_age",
            "bc") := NULL]
  gc()
  yrs <- c("year_incCond", "year_BMM", "year_CMM", "max_year")
  for(col in yrs)
    set(data, j = col, value = as.integer(data[[col]]))
  rm(yrs)

  setkey(data, mc, patient_id)

  #THIS IS SLOW
  pats <- data[tt,
               on = c("init_year  <= y", "year_max_year >= y")]
  na.omit(pats, cols = "patient_id")
  pats[, c("init_year", "year_max_year") := NULL] #this year_max_year gets capped at 2049
  gc()

#this takes some time, but isn't horrible
  setkey(pats, mc, patient_id, year)
  pats[, `:=` (
    age = year - yob,
    state_incCond = ifelse(year_incCond < year | init_state != "Healthy", 2L,
      ifelse(year_incCond == year, 1L,
             0L)),
    state_BMM = ifelse(year_BMM < year | init_state %in% c("BMM", "CMM"), 2L,
      ifelse(year_BMM == year, 1L,
             0L)),
    state_CMM = ifelse(year_CMM < year | init_state == "CMM", 2L,
                             ifelse(year_CMM == year,1L, 0L)))]
  pats[is.na(state_incCond), state_incCond := 0L ]
  pats[is.na(state_BMM), state_BMM := 0L ]
  pats[is.na(state_CMM), state_CMM := 0L ]
  pats[, c("init_state", "year_incCond", "year_BMM", "year_CMM") := NULL]
  gc()
  pats
}


#Make a prevalence summary table
mk_prevtab <- function(data, outstrata){
  rbind(data[,
             .(N = sum(state_incCond  != 0L), atrisk = .N),
               #outcome = sum(state_incCond  != 0L) /.N *100),
             keyby = outstrata][
               , state := "IncCond"],
        data[,
             .(N = sum(state_BMM  != 0L), atrisk = .N),
              # outcome = sum(state_BMM  != 0L) /.N *100),
             keyby = outstrata][
               , state := "BMM"],
        data[,
             .(N = sum(state_CMM  != 0L), atrisk = .N),
              # outcome = sum(state_CMM  != 0L) /.N *100),
             keyby = outstrata][
               , state := "CMM"])
}


# Make an incidence summary table
mk_inctab <- function(data, outstrata){
  rbind(data[,
                     .(N = sum(state_incCond  == 1L), atrisk = sum(state_incCond  != 2L) ),
                     # outcome = sum(state_incCond  == 1L) /sum(state_incCond  != 2L) *10000),
                     keyby = outstrata][
                       , state := "IncCond"],
        data[,
                     .(N = sum(state_BMM  == 1L), atrisk = sum(state_BMM  != 2L)),
                     #outcome = sum(state_BMM  == 1L) /sum(state_BMM  != 2L) *10000),
                     keyby = outstrata][
                       , state := "BMM"],
        data[,
                     .(N = sum(state_CMM  == 1L), atrisk = sum(state_CMM  != 2L)),
                     #outcome = sum(state_CMM  == 1L) /sum(state_CMM  != 2L) *10000),
                     keyby = outstrata][
                       , state := "CMM"])
  }

# Make a case fatality summary table
mk_cftab <- function(data, outstrata){
  rbind(data[,
                    .(N = sum(max_year  == year), atrisk = .N ),
                    keyby = outstrata][
                      , state := "Overall"],
        data[,
                    .(N = sum(max_year  == year & state_incCond == 0), atrisk = sum(state_incCond == 0L) ),
                    keyby = outstrata][
                      , state := "Healthy"],
        data[,
                    .(N = sum(max_year == year & state_incCond  != 0L), atrisk = sum(state_incCond  != 0L ) ),
                    keyby = outstrata][
                      , state := "IncCond"],
        data[,
                    .(N = sum(max_year  == year & state_BMM  != 0L), atrisk = sum(state_BMM  != 0L)),
                    keyby = outstrata][
                      , state := "BMM"],
        data[,
                    .(N = sum(max_year  == year & state_CMM  != 0L), atrisk = sum(state_CMM  != 0L)),
                    keyby = outstrata][
                      , state := "CMM"])
}


fn_state <- function(data){
  data[, state := factor(state,
                            levels = c("IncCond", "BMM", "CMM"))]
}


fn_state_cf <- function(data){
  data[, state := factor(state,
                          levels = c("Overall", "Healthy","IncCond", "BMM", "CMM"),
                          labels = c("Overall", "Healthy","IncCond", "BMM", "CMM") )]
}






agegroup5yfn <- function(x) {
  lev <- c("15-19 years", "20-24 years", "25-29 years", "30-34 years",
           "35-39 years", "40-44 years", "45-49 years", "50-54 years", "55-59 years",
           "60-64 years", "65-69 years", "70-74 years", "75-79 years", "80-84 years",
           "85-89 years", "90plus years")
  x[, agegrp5 := fcase(
    (age) >= 90L,           factor(lev[16], levels = lev),
    between(age, 85, 89), factor(lev[15], levels = lev),
    between(age, 80, 84), factor(lev[14], levels = lev),
    between(age, 75, 79), factor(lev[13], levels = lev),
    between(age, 70, 74), factor(lev[12], levels = lev),
    between(age, 65, 69), factor(lev[11], levels = lev),
    between(age, 60, 64), factor(lev[10], levels = lev),
    between(age, 55, 59), factor(lev[9], levels = lev),
    between(age, 50, 54), factor(lev[8], levels = lev),
    between(age, 45, 49), factor(lev[7], levels = lev),
    between(age, 40, 44), factor(lev[6], levels = lev),
    between(age, 35, 39), factor(lev[5], levels = lev),
    between(age, 30, 34), factor(lev[4], levels = lev),
    between(age, 25, 29), factor(lev[3], levels = lev),
    between(age, 20, 24), factor(lev[2], levels = lev),
    between(age, 15, 19), factor(lev[1], levels = lev)
  )
  ]
}
agegroup10yfn <- function(x) {
  lev <- c("15-19 years", "20-29 years",
           "30-39 years", "40-49 years",
           "50-59 years", "60-69 years",
           "70-79 years", "80-89 years",
           "90plus years")
  x[, agegrp10 := fcase(
    (age) >= 90L,           factor(lev[9], levels = lev),
    between(age, 80, 89), factor(lev[8], levels = lev),
    between(age, 70, 79), factor(lev[7], levels = lev),
    between(age, 60, 69), factor(lev[6], levels = lev),
    between(age, 50, 59), factor(lev[5], levels = lev),
    between(age, 40, 49), factor(lev[4], levels = lev),
    between(age, 30, 39), factor(lev[3], levels = lev),
    between(age, 20, 29), factor(lev[2], levels = lev),
    between(age, 15, 19), factor(lev[1], levels = lev)
  )
  ]
}

agegroup10yfn_b <- function(x) {
  lev <- c("30-39 years", "40-49 years",
           "50-59 years", "60-69 years",
           "70-79 years", "80plus years")
  x[, agegrp10 := fcase(
    (age) >= 80L,           factor(lev[6], levels = lev),
    between(age, 70, 79), factor(lev[5], levels = lev),
    between(age, 60, 69), factor(lev[4], levels = lev),
    between(age, 50, 59), factor(lev[3], levels = lev),
    between(age, 40, 49), factor(lev[2], levels = lev),
    between(age, 30, 39), factor(lev[1], levels = lev)
  )
  ]
}



# European standardised population 2013 (esp) weights
mk_esp <- function(){
tt <- data.table(agegrp5 = agegrp_name(0, 99),
                 wt_esp  = c(1000, 4000, 5500, 5500, 5500, 6000, 6000, 6500,
                             7000, 7000, 7000, 7000, 6500, 6000, 5500, 5000,
                             4000, 2500, 1500, 800, 200))
esp <- CJ(agegrp5 = agegrp_name(0, 99, ),
          gender = c("M", "F"),
          imd = as.factor(c(1:5))
)
absorb_dt(esp, tt)


esp[agegrp5  %in% c("90-94", "95-99"), agegrp5 := "90plus"]
esp <- esp[, .(wt_esp = sum(wt_esp)), keyby = .(gender, imd, agegrp5)]
}


fn_stnd <- function(data, strata){
  data[, agegrp5 := gsub(" years", "", agegrp5)]
  sumtab <- data[, lapply(.SD, sum), .SDcols = c("N", "atrisk"), keyby = eval(strata)]
  absorb_dt(sumtab, esp)
  sumtab[, popsize_wtd := wt_esp]
  sumtab[, wt_esp := wt_esp / atrisk]
  sumtab
}



fn_mc <- function(dt, outstrata, prbl, outcome_nm, prefix){
dt <- melt(dt, id.vars = outstrata) #need to make it long for the fquantile func for uncertainty intervals
#setkey(dt, "variable") #this is if you were doing it for more than one variable
dt <- dt[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(dt, c(setdiff(outstrata, "mc"), outcome_nm, percent(prbl, prefix = prefix)))
setkeyv(dt, setdiff(outstrata, "mc")) #setkeyv is because the col names are in ""; setdiff removes the "mc"
}



process_results <- function(){
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



fn_sex_lbls <- function(dt){
  dt[, gender := factor(gender,
                        levels = c("M", "F"),
                        labels = c("Men", "Women"))]
}

fn_state_lbls <- function(dt){
  dt[, state := factor(state,
                       levels = c("IncCond", "BMM", "CMM"),
                       labels = c("1 condition", "Basic multimorbidity", "Complex multimorbidity"))]
}

fn_statehealthy_lbls <- function(dt){
  dt[, state := factor(state,
                       levels = c("Healthy","IncCond", "BMM", "CMM"),
                       labels = c("Healthy", "1 condition", "Basic multimorbidity", "Complex multimorbidity"))]
}

fn_statecf_lbls <- function(dt){
  dt[, state := factor(state,
                       levels = c("Overall", "Healthy","IncCond", "BMM", "CMM"),
                       labels = c("Overall", "Healthy", "1 condition", "Basic multimorbidity", "Complex multimorbidity"))]
}

fn_scenario_lbls <- function(dt){
  dt[, scenario := factor(scenario,
                         levels = c("scenario1","scenario2","scenario3","scenario4","scenario5"),
                         labels = c("1 Targeted", "2 Universal + focus on gap",
                                    "3 Redistributive", "4 Proportionate universalism",
                                    "5 Inequalities removal"))]
}
fn_scenario_lbls2 <- function(dt){
  dt[, scenario := factor(scenario,
                          levels = c("baseline","scenario1","scenario2","scenario3","scenario4","scenario5"),
                          labels = c("0 Baseline","1 Targeted", "2 Universal + focus on gap",
                                     "3 Redistributive", "4 Proportionate universalism",
                                     "5 Inequalities removal"))]
}
