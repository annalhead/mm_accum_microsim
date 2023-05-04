#Compare the different scenarios


#scenarios
library(data.table)
data.table::setDTthreads(2) #this is so that don't use all the processors
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
library(cowplot)
library(scales)
library(fst)


input_path <-
  function(x = character(0))
    paste0("./inputs/", x)

sim_path <-
  function(x = character(0))
    paste0("./simulation/", x)

outpt_pth <-
  function(x = character(0))
    paste0("./output/",scenario,"/", x)


outpt_pth_cmpr_plt <-
  function(x = character(0))
    paste0("./output/comparisons/", x)

plotsave <- FALSE
ineq <- FALSE

# blue palette
mypalette <- c("#A8DDB5" ,"#7BCCC4" ,"#4EB3D3" ,"#2B8CBE", "#08589E")
mypalette2 <- c("#CCEBC5", "#A8DDB5" ,"#7BCCC4" ,"#4EB3D3" ,"#2B8CBE", "#0868AC" , "#084081")

# Load all the extra functions etc.
source(sim_path("/fn.R"))

# Load the model
sourceCpp(sim_path("/mmsim_mod_bc.cpp"), cacheDir = sim_path("./.cache"))


scenario_nm <- c("baseline", "scenario1",
               "scenario2", "scenario3", "scenario4", "scenario5"
               )


#uncertainty intervals
prbl <- c( 0.5, 0.25, 0.75, 0.025, 0.975)



theme_set(new = theme_few())
theme_update(axis.text = element_text(size = 10),
             axis.title = element_text(size = 10,
                                       margin = margin(t = 10, r = 10, b = 10, l = 10)),
             legend.position = "bottom",
             legend.text=element_text(size=10),
             legend.title=element_text(size=10)
             #plot.title = element_text(hjust = 0.5),
)

bigtab <- data.table()
for (scenario in scenario_nm){
  tmptab <- fread(outpt_pth(paste0(scenario,"_prvl.csv.gz")), header = TRUE)[year == 2049]
  tmptab[, scenario := scenario]
  bigtab <- rbind(bigtab, tmptab[, scenario := scenario])
}

bigtab[, imd := factor(imd)]
healthy <- bigtab[scenario %in% sc_nm & state == "IncCond"]
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

bigtab[, outcome := N/atrisk *100]
#fwrite(bigtab, outpt_pth_cmpr_plt("prvltab.csv.gz"), row.names = FALSE)

# Standardising to esp
esp <- mk_esp()
strata <- c("mc", "year", "gender", "imd", "agegrp_big", "agegrp5", "agegrp10", "state", "scenario")

#prev
bigtabsum <- fn_stnd(bigtab, strata)


#Proportion by state - unstandardised
outstrata <- c("mc", "year", "scenario", "state")
d <- bigtabsum[,
               .(prop = sum(N)/sum(atrisk)), keyby = outstrata]
d <- fn_mc(d, outstrata, prbl, "outcome", "prop_") #this then gives the median and various uncertainty intervals as defined by prbl above
fwrite(d, outpt_pth_cmpr_plt("prvl_inclsv_states.csv"), row.names = FALSE)



#Differences in prev by agegroup
outstrata <- c("mc", "scenario","state", "agegrp_big")
d <- bigtabsum[year == 2049, .(N = sum(N)), keyby = outstrata]
tmp <- d[scenario == "baseline"]
d <- d[scenario != "baseline"]
d[tmp, on = c("mc", "state", "agegrp_big"), bl := i.N]
d[, diff := N-bl]
d <- fn_mc(d, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above
d[state %in% c("BMM", "CMM") & outcome == "diff"]
fwrite(d, outpt_pth_cmpr_plt("prvldiff_age.csv"), row.names = FALSE)


#Differences in prev by imd
outstrata <- c("mc", "scenario","state", "imd")
d <- bigtabsum[year == 2049, .(N = sum(N)), keyby = outstrata]
tmp <- d[scenario == "baseline"]
d <- d[scenario != "baseline"]
d[tmp, on = c("mc", "state", "imd"), bl := i.N]
d[, diff := N-bl]
d <- fn_mc(d, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above
d[state %in% c("BMM", "CMM") & outcome == "diff"]
fwrite(d, outpt_pth_cmpr_plt("prvldiff_imd.csv"), row.names = FALSE)


#mutually exclusive health states

prvlbystate <- bigtabsum[, .(mc, year, gender, imd, agegrp5 ,agegrp10,  agegrp_big, state,     scenario  ,  N ,atrisk)]
prvlbystate <- dcast(prvlbystate,
                    mc + year + gender + imd + agegrp5 + agegrp10 + agegrp_big + scenario + atrisk ~ state,
                    value.var = "N")

prvlbystate[, Healthy := atrisk - IncCond][
  , IncCond := IncCond - BMM][
    , BMM := BMM - CMM]
prvlbystate <- melt(prvlbystate,
                    id.vars = c("mc","year", "gender", "imd", "agegrp5", "agegrp10", "agegrp_big", "scenario", "atrisk"),
                     measure.vars = c("Healthy", "IncCond", "BMM", "CMM"),
                    value.factor = TRUE,
                    value.name = "N",
                    variable.name = "state")
prvlbystate[, imd := factor(imd)]

# reverse the levels for plotting
prvlbystate[, state := factor(state,
                              levels = c("CMM", "BMM", "IncCond", "Healthy"))]



#Making this a table
outstrata <- c("mc", "year", "imd", "scenario", "state")

d <- prvlbystate[,
                 .(N = sum(N) * 100), keyby = outstrata]
d <- fn_mc(d, outstrata, prbl, "outcome", "N_")
fwrite(d, outpt_pth_cmpr_plt("N_exclsv_states_imd.csv"), row.names = FALSE)

#Making this a table - not by imd
outstrata <- c("mc", "year", "scenario", "state")

d <- prvlbystate[,
                 .(N = sum(N) * 100), keyby = outstrata]
d <- fn_mc(d, outstrata, prbl, "outcome", "N_")
fwrite(d, outpt_pth_cmpr_plt("N_exclsv_states.csv"), row.names = FALSE)


d_w <- dcast(d[year == 2049], state ~ scenario, value.var = "N_50%")


#Making this a table - by agegroup
outstrata <- c("mc", "year", "agegrp_big","scenario", "state")

d <- prvlbystate[,
                 .(N = sum(N) * 100), keyby = outstrata]
d <- fn_mc(d, outstrata, prbl, "outcome", "N_")
fwrite(d, outpt_pth_cmpr_plt("N_exclsv_states_agegrp.csv"), row.names = FALSE)


d_w2 <- dcast(d[year == 2049], agegrp_big + state ~ scenario, value.var = "N_50%")
d_w <- rbind(d_w[, agegrp_big := "all"], d_w2)


#Making this a tableby sex
outstrata <- c("mc", "year", "gender", "scenario", "state")

d <- prvlbystate[,
                 .(N = sum(N) * 100), keyby = outstrata]
d <- fn_mc(d, outstrata, prbl, "outcome", "N_")
fwrite(d, outpt_pth_cmpr_plt("N_exclsv_states_sex.csv"), row.names = FALSE)




#Making this a table of crude, exclusive prevalence
outstrata <- c("mc", "year", "imd", "scenario", "state")

d <- prvlbystate[,
                 .(prvl = (sum(N)/sum(atrisk)) * 100), keyby = outstrata]
d <- fn_mc(d, outstrata, prbl, "outcome", "prvl_")
fwrite(d, outpt_pth_cmpr_plt("prvl_exclsv_states_imd.csv"), row.names = FALSE)

#Making this a table - not by imd
outstrata <- c("mc", "year", "scenario", "state")

d <- prvlbystate[,
                 .(prvl = (sum(N)/sum(atrisk)) * 100), keyby = outstrata]
d <- fn_mc(d, outstrata, prbl, "outcome", "N_")
fwrite(d, outpt_pth_cmpr_plt("prvl_exclsv_states.csv"), row.names = FALSE)

#Making this a table of crude, exclusive prevalence by sex
outstrata <- c("mc", "year", "gender", "scenario", "state")

d <- prvlbystate[,
                 .(prvl = (sum(N)/sum(atrisk)) * 100), keyby = outstrata]
d <- fn_mc(d, outstrata, prbl, "outcome", "prvl_")
fwrite(d, outpt_pth_cmpr_plt("prvl_exclsv_states_sex.csv"), row.names = FALSE)

#Making this a table of crude, exclusive prevalence by agegetp
outstrata <- c("mc", "year", "agegrp_big", "scenario", "state")

d <- prvlbystate[,
                 .(prvl = (sum(N)/sum(atrisk)) * 100), keyby = outstrata]
d <- fn_mc(d, outstrata, prbl, "outcome", "prvl_")
fwrite(d, outpt_pth_cmpr_plt("prvl_exclsv_states_agegrp.csv"), row.names = FALSE)




#absolute numbers
outstrata <- c("mc", "year", "imd", "scenario", "state")
d <- prvlbystate[,
                 .(N = sum(N)), keyby = outstrata]
d <- fn_mc(d, outstrata, prbl, "outcome", "N_") #this then gives the median and various uncertainty intervals as defined by prbl above



#Just plotting bog-standard prbl
outstrata <- c("mc", "year", "scenario", "state", "gender")
d <- bigtabsum[,
               .(prvl_cr = sum(N)/ sum(atrisk)), keyby = outstrata]
d <- fn_mc(d, outstrata, prbl, "outcome", "prvl_cr_") #this then gives the median and various uncertainty intervals as defined by prbl above

ggplot(d,
       aes(x = year, y = `prvl_cr_50%` , col = scenario)) +
  geom_smooth(se = FALSE, size = 0.5) +
  facet_grid(state ~ gender) +
  expand_limits(y = 0) +
  ylab("Crude prevalence") +
  ggtitle("Prevalence by scenario, state and gender")
if(plotsave){
  ggsave(filename = outpt_pth_cmpr_plt("prvl_gender.png"))
}
# Standardised
outstrata <- c("mc", "year", "scenario", "state")
d <- bigtabsum[,
               .(prvl_st = sum(N * wt_esp)/ sum(popsize_wtd)), keyby = outstrata]
d <- fn_mc(d, outstrata, prbl, "outcome", "prvl_st_") #this then gives the median and various uncertainty intervals as defined by prbl above
ggplot(d,
       aes(x = year, y = `prvl_st_50%`, col = scenario)) +
  geom_smooth(se = FALSE, size = 0.5) +
  facet_grid(rows = "state", scales = "free" ) +
  expand_limits(y = 0) +
  ylab("Standardised prevalence")

if(plotsave){
  ggsave(filename = outpt_pth_cmpr_plt("prvl_overall_st.png"))
}

outstrata <- c("mc", "year", "imd", "scenario", "state")
d <- bigtabsum[,
                 .(prvl_cr = sum(N)/ sum(atrisk)), keyby = outstrata]
d <- fn_mc(d, outstrata, prbl, "outcome", "prvl_cr_") #this then gives the median and various uncertainty intervals as defined by prbl above

ggplot(d,
       aes(x = year, y = `prvl_cr_50%` , col = imd)) +
  geom_smooth(se = FALSE, size = 0.5) +
  facet_grid(state ~ scenario) +
  expand_limits(y = 0) +
  ylab("Crude prevalence") +
  ggtitle("Prevalence by imd")

# Standardised
outstrata <- c("mc", "year", "imd", "scenario", "state")
d <- bigtabsum[,
               .(prvl_st = sum(N * wt_esp)/ sum(popsize_wtd)), keyby = outstrata]
d <- fn_mc(d, outstrata, prbl, "outcome", "prvl_st_") #this then gives the median and various uncertainty intervals as defined by prbl above
ggplot(d,
       aes(x = year, y = `prvl_st_50%`, col = imd)) +
  geom_smooth(se = FALSE, size = 0.5) +
  facet_grid(state ~ scenario) +
  expand_limits(y = 0) +
  ylab("Standardised prevalence")






# Making a baseline only stacked plot - age-sex standardised
scenario <- "baseline"
outstrata <- c("mc", "year", "imd", "state")
d <- prvlbystate[scenario %in% scenario,
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
  ggsave(filename = outpt_pth("simplots/prvl_overall_st_excl.png"))
}





if(ineq){

# RII & SII for prevalence - basically in crude rates no absolute or relative diff
library(PHEindicatormethods) #Need to do this with ridit scores so don't have to use the old package
require(tidyverse)
#the agegrp10 col keeps going weird
bigtab[, agegrp10 := agegrp5]
bigtab[agegrp5 %in% c("30-34", "35-39"), agegrp10 := "30-39"]
bigtab[agegrp5 %in% c("40-44", "45-49"), agegrp10 := "40-49"]
bigtab[agegrp5 %in% c("50-54", "55-59"), agegrp10 := "50-59"]
bigtab[agegrp5 %in% c("60-64", "65-69"), agegrp10 := "60-69"]
bigtab[agegrp5 %in% c("70-74", "75-79"), agegrp10 := "70-79"]
bigtab[agegrp5 %in% c("80-84", "85-89"), agegrp10 := "80-89"]
bigtab[, agegrp10 := factor(agegrp10)]
#Don't want healthy in this becasue it confuses things - higher % is better
data <- bigtab[year %in% c(2049)  &
                    state != "Healthy", .(
                      N = sum(N),
                      atrisk = sum(atrisk),
                      outcome = sum(N) / sum(atrisk) * 100
                    ),
                    by = c("mc", "imd", "agegrp10", "scenario", "state")][
                      , se := sqrt((outcome * (100 - outcome)) / atrisk)]
prvl_sii <-
  as.data.table(
    phe_sii(
      group_by(data, mc, agegrp10, scenario, state),
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


fn_state(prvl_sii)

outstrata <- c("mc", "year", "scenario", "state")
d <- prvl_sii[,c(1:6)]
d <- fn_mc(d, outstrata, prbl, "outcome", "ineq_") #this then gives the median and various uncertainty intervals as defined by prbl above

ggplot(d[year == 2049 & outcome == "sii"], aes(x = state, y = `ineq_50%`, fill = scenario)) +
  geom_col(position = "dodge") +
  ylab("Absolute difference in prevalence") + #between most & least deprived
  ggtitle("SII in 2049")
if(plotsave){
  ggsave(filename = outpt_pth_cmpr_plt("prvl_sii.png"))
}

ggplot(d[year == 2049 & outcome == "rii"], aes(x = state, y = `ineq_50%` - 1, fill = scenario)) +
  geom_col(position = "dodge") +
  scale_y_continuous(
    limits = c(-0.08, 0.01),
    breaks = c(-0.06, -0.04, -0.02, 0, 0.02),
    labels =  c(0.94, 0.96, 0.98, 1, 1.02),
    name = "Relative Inequalities") +
  geom_hline(aes(yintercept = 0)) +#add x-axis at 0
  ggtitle("RII in 2049")


#By agegroup - there are inequalities if you do it by agegroup


data <- bigtabsum[year %in% c(2019, 2029, 2039, 2049) &
                  state != "Healthy", .(
  N = sum(N),
  atrisk = sum(atrisk),
  outcome = sum(N) / sum(atrisk) * 100
),
by = c("mc", "year", "imd", "agegrp10", "scenario", "state")][
  , se := sqrt((outcome * (100 - outcome)) / atrisk)]

prvl_sii_age10 <-
  as.data.table(
    phe_sii(
      group_by(data, mc, year, agegrp10, scenario, state),
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

fn_state(prvl_sii_age10)


outstrata <- c("mc", "year", "scenario", "agegrp10", "state")
d <- prvl_sii_age10[,c(1:7)]
d <- fn_mc(d, outstrata, prbl, "outcome", "ineq_") #this then gives the median and various uncertainty intervals as defined by prbl above

ggplot(d[year == 2049 & outcome == "sii"], aes(x = agegrp10 , y = `ineq_50%`, fill = scenario)) +
  geom_col(position = "dodge") +
  facet_grid(~ state) +
  ylab("Absolute difference in prevalence") #between most & least deprived

#With error bars
ggplot(d[year == 2049 & outcome == "sii"], aes(x = agegrp10 , y = `ineq_50%`, fill = scenario)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin=`ineq_2%`, ymax=`ineq_98%`),
                width=.2, position=position_dodge(.9)) +
  facet_grid(~ state) +
  ylab("Absolute difference in prevalence") #between most & least deprived
if(plotsave){
  ggsave(filename = outpt_pth_cmpr_plt("prvl_sii_agegrp.png"))
}

ggplot(d[year == 2049 & outcome == "rii"], aes(x = agegrp10 , y = `ineq_50%` - 1, fill = scenario)) +
  geom_col(position = "dodge") +
  facet_grid(~ state) +
  scale_y_continuous(
    limits = c(-0.1, 0.6),
    breaks = c( 0, 0.2, 0.4, 0.6),
    labels =  c(1, 1.2, 1.4, 1.6),
    name = "Relative Inequalities")
}



### What about cumulative incident cases?

bigtab_inc <- data.table()
for (scenario in scenario_nm){
  tmptab <- fread(outpt_pth(paste0(scenario,"_incd.csv.gz")), header = TRUE)
  setkey(tmptab, mc, year, gender, age, region, state)
  tmptab[, agegrp_big := ifelse(age <= 65, "wrk_age","over65") ]
  tmptab <- tmptab[, .(N = sum(N), atrisk = sum(atrisk)), keyby = .(mc, year, gender, imd, agegrp_big, state)]
  setkey(tmptab, mc, year, gender, imd, agegrp_big, state)
  tmptab[, cumN := cumsum(N), keyby = .(mc, gender, imd, agegrp_big, state)]
  tmptab <- tmptab[year == 2049]
  tmptab[, scenario := scenario]
  bigtab_inc <- rbind(bigtab_inc, tmptab[, scenario := scenario])
  rm(tmptab)
  gc()
}
fwrite(bigtab_inc, outpt_pth_cmpr_plt(paste0("incd_sum.csv.gz")), )
bigtab_inc[, imd := factor(imd)]

setkey(bigtab_inc, mc, scenario, gender, imd, agegrp_big, state)
bigtab_inc_sum <- bigtab_inc[, .(cumN = sum(cumN)), keyby = .(mc, scenario, imd, state)]
setkey(bigtab_inc_sum, mc, scenario, imd, state)

#Cases prevented or postponed, but this doesn't take into account people living longer...
# which is why cpp are higher in 2039 than 2049
bl <- bigtab_inc_sum[scenario == "baseline"]
bigtab_inc_sum2 <- bigtab_inc_sum[bl, on = c("mc", "imd", "state"), bl := i.cumN]
bigtab_inc_sum2 <- bigtab_inc_sum2[scenario != "baseline"]
bigtab_inc_sum2[, ppcases := bl-cumN]
bigtab_inc_sum2[, relred := cumN/bl *100]
fn_state(bigtab_inc_sum2)


 outstrata <- c("mc", "scenario", "imd", "state")
 d <- fn_mc(bigtab_inc_sum2, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above

rslts_tab <- d[outcome == "cumN", .(scenario, state, imd, cumN = `outcome_50%` * 100)]
tmp <-  d[ outcome == "relred", .(scenario, state, imd, relred = `outcome_50%` )]
rslts_tab <- rslts_tab[tmp, on = c("scenario", "state", "imd")]
ggplot(tmp, aes(x = scenario, y = relred, fill = imd)) +
  geom_col(position = "dodge") +
  facet_grid(rows = vars(state)) +
  ylab("Cases prevented or postponed") +
  ggtitle("Cumulative incident cases prevented or postponed by imd")




#overall
outstrata <- c("mc", "scenario", "state")
d <- bigtab_inc_sum2[ , .(cumN = sum(cumN), bl = sum(bl)), keyby = .(mc, scenario, state)]
d[, relred := cumN/bl * 100][,ppcases := bl-cumN]
d <- fn_mc(d, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above
#fwrite(d, outpt_pth_cmpr_plt("ppcase_overall.csv"), row.names = FALSE)

tmp <- d[outcome == "cumN", .(scenario, state, imd = "Overall", cumN = `outcome_50%` * 100)]
tmp2 <- d[outcome == "relred", .(scenario, state, imd = "Overall", relred = `outcome_50%` )]
tmp[tmp2, on = c("scenario", "state"), relred := i.relred]
#rslts_tab <- rbind(rslts_tab, d)
rslts_tab <- rbind(rslts_tab, tmp)
d <- bigtab_inc_sum2[, .(cumN = sum(cumN)), keyby = .(mc, scenario, state)]
d <- fn_mc(d, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above
d <- d[, .(scenario, state, imd = "all", cumN = `outcome_50%` * 100)]

rslts_tab[, imd := factor(imd, levels = c("Overall", "1", "2", "3", "4", "5"))]
fn_state_lbls(rslts_tab)
mypalette3 <- c("#666666", mypalette)
fn_scenario_lbls(rslts_tab)
ggplot(rslts_tab[state != "1 condition"], aes(x = imd, y = relred - 100, fill = imd)) +
  geom_col(position = "dodge") +
  facet_grid(state ~scenario ) +
  xlab("IMD quintile") +
  ylab("Cumulative incident cases compared to baseline") +
  #ggtitle("Cumulative incident cases in 2019 as a proportion of the baseline scenario") +
  labs(fill = "IMD quintile") +
  scale_fill_manual(values = mypalette3) +
  scale_y_continuous(limits = c(-15,12),
                     breaks = c(-15, -10, -5, 0, 5, 10),
                     labels = c("85%","90%", "95%","100%","105%", "110%")) +
  theme(legend.position="bottom", axis.title.y = element_text(size = 12)) +
  guides(fill = guide_legend(nrow = 1))
if(plotsave){
  #ggsave(filename = outpt_pth_cmpr_plt("cpp_prop.png"))
  ggsave2(outpt_pth_cmpr_plt("cpp_prop.svg") , scale = 0.7, width = 16, height = 10)

}

if(ineq){
# Equity summary chart - Not 100% sure I have got this right...
#Also, how to calculate the equity curve
dt <- bigtab_inc_sum2[year %in% c(2019, 2029, 2039, 2049),
                          .(ppcases = sum(ppcases), atrisk = sum(atrisk)),
                            keyby = .(mc, imd, year, scenario, state)][
                              ,se := 0 # need to add this in...
                            ]


# SHould be able to do this by calculating rigit scores
require(tidyverse)
sii_inc <-
  as.data.table(
    phe_sii(
      group_by(dt, mc, year, scenario, state),
      imd,
      atrisk,
      value_type = 0,
      # for rates
      value = ppcases,
      se = se,
      confidence = 0.95,
      rii = TRUE,
      type = "standard"
    )
  )

sii_inc[, c(7:10) := NULL] # for now don't need the CIs

outstrata <- c("mc", "year", "scenario", "state")

sii_inc <- sii_inc[dt[, .(cpp =as.numeric(sum(ppcases))),
              keyby = outstrata],
           on = outstrata]
d <- fn_mc(sii_inc, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above

fn_state(d)

d <- dcast(d, year + scenario + state ~ `outcome`, value.var = "outcome_50%"
              )

ggplot(d[year == 2049 ],
       aes(x = cpp,
           y = sii, col = scenario)) +
  geom_point() +
  facet_grid(~ state) +
  expand_limits(x = 0) +
  xlab("Cases prevented or postponed") +
  ylab("Absolute inequality reduction (cases)") +
  labs(colour = "Scenario") +
  ggtitle("Reduction in cases and absolute inequality: 2049") +
  geom_vline(aes(xintercept = 0)) + #add y-axis at 0
  geom_hline(aes(yintercept = 0)) #add x-axis at 0

}

outstrata <- c("mc", "scenario", "state")
bigtab_inc_sum2w <- dcast(bigtab_inc_sum2, mc + scenario + state ~ imd, value.var = "ppcases")
d <- fn_mc(bigtab_inc_sum2w, outstrata, prbl, "absdiff", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above


#What about splitting the chart by over and under 65s
bigtab_inc_sum_65split <- bigtab_inc[, .(N = sum(N),
                                         atrisk = sum(atrisk),
                                         cumN = sum(cumN)),
                                     keyby = .(scenario, mc, agegrp_big, imd, state)]
setkey(bigtab_inc_sum_65split, mc, scenario, imd, agegrp_big, state)

outstrata <- c("mc", "scenario", "imd", "state", "agegrp_big")
d <- fn_mc(bigtab_inc_sum_65split, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above




#Cases prevented or postponed, but this doesn't take into account people living longer...
# which is why cpp are higher in 2039 than 2049
bl <- bigtab_inc_sum_65split[scenario == "baseline"]
bigtab_inc_sum_65split2 <- bigtab_inc_sum_65split[bl, on = c("mc", "agegrp_big","imd", "state"), bl := i.cumN]
bigtab_inc_sum_65split2 <- bigtab_inc_sum_65split2[scenario != "baseline"]
bigtab_inc_sum_65split2[, ppcases := bl-cumN]
bigtab_inc_sum_65split2[, relred := cumN/bl * 100]
fn_state(bigtab_inc_sum_65split2)

outstrata <- c("mc", "scenario", "imd", "agegrp_big","state")
d <- fn_mc(bigtab_inc_sum_65split2, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above


ggplot(d[outcome == "ppcases" & year %in% c(2019, 2029, 2039, 2049) &  agegrp_big == "wrk_age"],
       aes(x = year, y = `outcome_50%`, fill = imd)) +
  geom_col(position = "dodge") +
  facet_grid(state ~ scenario) +
  ylab("Cases prevented or postponed") +
  ggtitle("Cumulative incident cases prevented or postponed by imd - 65 and under")

ggplot(d[outcome == "ppcases" & year %in% c(2019, 2029, 2039, 2049) & agegrp_big == "over65"],
       aes(x = year, y = `outcome_50%`, fill = imd)) +
  geom_col(position = "dodge") +
  facet_grid(state ~ scenario) +
  ylab("Cases prevented or postponed") +
  ggtitle("Cumulative incident cases prevented or postponed by imd - over 65s")

d[outcome == "ppcases" & agegrp_big == "wrk_age" & scenario == "scenario4"]


rslts_tab <- d[outcome == "cumN", .(scenario, state, agegrp_big, imd, cumN = `outcome_50%` * 100)]
tmp <-  d[ outcome == "relred", .(scenario, state, agegrp_big, imd, relred = `outcome_50%` )]
rslts_tab <- rslts_tab[tmp, on = c("scenario", "state", "agegrp_big", "imd")]



#overall
outstrata <- c("mc", "scenario", "state", "agegrp_big")
d <- bigtab_inc_sum_65split2[, .(cumN = sum(cumN), bl = sum(bl)), keyby = .(mc, scenario, state, agegrp_big)]
d[, relred := cumN/bl * 100]
d[, ppcases := bl-cumN]

d <- fn_mc(d, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above
tmp <- d[outcome == "cumN", .(scenario, state, agegrp_big, imd = "Overall", cumN = `outcome_50%` * 100)]
tmp2 <- d[outcome == "relred", .(scenario, state, agegrp_big, imd = "Overall", relred = `outcome_50%` )]
tmp[tmp2, on = c("scenario", "state", "agegrp_big"), relred := i.relred]
#rslts_tab <- rbind(rslts_tab, d)
rslts_tab <- rbind(rslts_tab, tmp)

rslts_tab[, imd := factor(imd, levels = c("Overall", "1", "2", "3", "4", "5"))]
fn_state_lbls(rslts_tab)
mypalette3 <- c("#666666", mypalette)
fn_scenario_lbls(rslts_tab)
ggplot(rslts_tab[agegrp_big == "wrk_age" & state != "1 condition"], aes(x = imd, y = relred - 100, fill = imd)) +
  geom_col(position = "dodge") +
  facet_grid(state ~scenario ) +
  xlab("IMD quintile") +
  ylab("Cumulative incident cases compared to baseline") +
  #ggtitle("Cumulative incident cases in 2019 as a proportion of the baseline scenario") +
  labs(fill = "IMD quintile") +
  scale_fill_manual(values = mypalette3) +
  scale_y_continuous(limits = c(-15,15),
                     breaks = c(-15, -10, -5, 0, 5, 10, 15),
                     labels = c("85%","90%", "95%","100%","105%", "110%", "115%")) +
  theme(legend.position="bottom", axis.title.y = element_text(size = 12)) +
  guides(fill = guide_legend(nrow = 1))
if(plotsave){
#  ggsave(filename = outpt_pth_cmpr_plt("cpp_prop_wrkage.png"))
  ggsave2(outpt_pth_cmpr_plt("cpp_prop_wrkage.svg") , scale = 0.7, width = 16, height = 10)
}

ggplot(rslts_tab[agegrp_big == "over65" & state != "1 condition"], aes(x = imd, y = relred - 100, fill = imd)) +
  geom_col(position = "dodge") +
  facet_grid(state ~scenario ) +
  xlab("IMD quintile") +
  ylab("Cumulative incident cases compared to baseline") +
  #ggtitle("Cumulative incident cases in 2019 as a proportion of the baseline scenario") +
  labs(fill = "IMD quintile") +
  scale_fill_manual(values = mypalette3) +
  scale_y_continuous(limits = c(-15,15),
                     breaks = c(-15, -10, -5, 0, 5, 10, 15),
                     labels = c("85%","90%", "95%","100%","105%", "110%", "115%")) +
  theme(legend.position="bottom", axis.title.y = element_text(size = 12)) +
  guides(fill = guide_legend(nrow = 1))
if(plotsave){
  #ggsave(filename = outpt_pth_cmpr_plt("cpp_prop_65pl.png"))
  ggsave2(outpt_pth_cmpr_plt("cpp_prop_65pl.svg") , scale = 0.7, width = 16, height = 10)

}













#Just plotting bog-standard incidence
fn_state(bigtab_inc)

#overall by year
outstrata <- c("mc", "year", "scenario", "state", "gender")
d <- bigtab_inc[,
               .(incd_cr = sum(N)/ sum(atrisk)), keyby = outstrata]
d <- fn_mc(d, outstrata, prbl, "outcome", "incd_cr_") #this then gives the median and various uncertainty intervals as defined by prbl above


ggplot(d,
       aes(x = year, y = `incd_cr_50%` * 10000 , col = scenario)) +
  geom_smooth(se = FALSE, size = 0.5) +
  facet_grid(state ~ gender) +
  expand_limits(y = 0) +
  ylab("Crude incidence") +
  ggtitle("Incidence by scenario, state and gender")
if(plotsave){
  ggsave(filename = outpt_pth_cmpr_plt("incd_gender.png"))
}
# Standardised
# Standardising to esp
esp <- mk_esp()
strata <- c("mc","year", "gender", "imd", "agegrp5", "state", "scenario")
bigtab_inc <- fn_stnd(bigtab_inc, strata)

outstrata <- c("mc", "year", "scenario", "state")
d <- bigtab_inc[,
               .(incd_st = sum(N * wt_esp)/ sum(popsize_wtd)), keyby = outstrata]
d <- fn_mc(d, outstrata, prbl, "outcome", "incd_st_") #this then gives the median and various uncertainty intervals as defined by prbl above
ggplot(d,
       aes(x = year, y = `incd_st_50%` * 10000, col = scenario)) +
  geom_smooth(se = FALSE, size = 0.5) +
  facet_grid(rows = "state", scales = "free" ) +
  expand_limits(y = 0) +
  ylab("Standardised incidence")

if(plotsave){
  ggsave(filename = outpt_pth_cmpr_plt("incd_overall_st.png"))
}

ggplot(bigtab_inc[, sum(N)/ sum(atrisk), keyby = .(year, scenario, state, imd)],
       aes(x = year, y = V1 * 10000, col = imd)) +
  geom_smooth(se = FALSE, linewidth = 0.5) +
  facet_grid(state ~ scenario) +
  expand_limits(y = 0) +
  ylab("Crude incidence (per 10,000 persons") +
  ggtitle("Crude incidence by imd")

#standardised
ggplot(bigtab_inc[, sum(N * wt_esp)/ sum(popsize_wtd), keyby = .(year, scenario, state, imd)],
       aes(x = year, y = V1 * 10000, col = imd)) +
  geom_smooth(se = FALSE, linewidth = 0.5) +
  facet_grid(state ~ scenario) +
  expand_limits(y = 0) +
  ylab("Standardised incidence (per 10,000 persons") +
  ggtitle("Standardised incidence by imd")





### What about cumulative cf?
#these are not from mutually exclusive health states
bigtab_cf <- data.table()
for (scenario in scenario_nm){
  tmptab <- fread(outpt_pth(paste0(scenario,"_cf.fst")), as.data.table = TRUE)
  tmptab[, scenario := scenario]
  bigtab_cf <- rbind(bigtab_cf, tmptab[, scenario := scenario])
}

bigtab_cf[, imd := factor(imd)]
bigtab_cf[, state := factor(state,
                            levels = c("Healthy", "IncCond", "BMM", "CMM", "Overall"))]
setkey(bigtab_cf, scenario, year, gender, age, region, state)
bigtab_cf_sum <- bigtab_cf[, .(N = sum(N), atrisk = sum(atrisk)), keyby = .(mc, scenario, year, imd, state)]
setkey(bigtab_cf_sum, mc, scenario, imd, state, year)
bigtab_cf_sum[, cumN := cumsum(N), keyby = .(mc, scenario, imd, state)]


outstrata <- c("mc", "year", "scenario", "imd", "state")
d <- fn_mc(bigtab_cf_sum, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above

ggplot(d[outcome %in% "cumN"], aes(x = year, y = `outcome_50%`, col = imd)) +
    geom_line() +
  facet_grid(state ~ scenario, scales = "free")  +
  ylab("Deaths") +
  ggtitle("Cumulative deaths by imd")

ggplot(d[outcome %in% "cumN" & year %in% c(2019, 2029, 2039, 2049)],
       aes(x = year, y = `outcome_50%`, fill = imd)) +
  geom_col(position = "dodge") +
  facet_grid(state ~ scenario) +
  ylab("Deaths") +
  ggtitle("Cumulative deaths by imd")


#deaths prevented or postponed, but this doesn't take into account people living longer...
bl <- bigtab_cf_sum[scenario == "baseline"]
bigtab_cf_sum2 <- bigtab_cf_sum[bl, on = c("mc","year", "imd", "state"), bl := i.cumN]
bigtab_cf_sum2 <- bigtab_cf_sum2[scenario != "baseline"]
bigtab_cf_sum2[, ppdeaths := bl-cumN]

outstrata <- c("mc", "year", "scenario", "imd", "state")
d <- fn_mc(bigtab_cf_sum2, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above

ggplot(d[outcome == "ppdeaths" & year %in% c(2019, 2029, 2039, 2049)], aes(x = year, y = `outcome_50%`, fill = imd)) +
  geom_col(position = "dodge") +
  facet_grid(state ~ scenario) +
  ylab("Deaths prevented or postponed") +
  ggtitle("Cumulative deaths prevented or postponed by imd")

ggplot(d[outcome == "ppdeaths" & year %in% c(2019, 2029, 2039, 2049)],
       aes(x = year, y = `outcome_50%`, fill = scenario)) +
  geom_col(position = "dodge") +
  facet_grid(state ~ imd, scales = "free") +
  ylab("Deaths prevented or postponed") +
  ggtitle("Cumulative deaths prevented or postponed by imd")


#Making the equity plot
dt <- bigtab_cf_sum2[year %in% c(2019, 2029, 2039, 2049) & state == "Overall",
                      .(ppdeaths = sum(ppdeaths), atrisk = sum(atrisk)),
                      keyby = .(mc, imd, year, scenario)][
                        ,se := 0 # need to add this in...
                      ]

if(ineq){
# SHould be able to do this by calculating rigit scores
require(tidyverse)
sii_cf <-
  as.data.table(
    phe_sii(
      group_by(dt, mc, year, scenario),
      imd,
      atrisk,
      value_type = 0,
      # for rates
      value = ppdeaths,
      se = se,
      confidence = 0.95,
      rii = TRUE,
      type = "standard"
    )
  )

sii_cf[, c(5:8) := NULL] # for now don't need the CIs

outstrata <- c("mc", "year", "scenario")

sii_cf <- sii_cf[dt[, .(dpp =as.numeric(sum(ppdeaths))),
              keyby = outstrata],
           on = outstrata]
d <- fn_mc(sii_cf, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above


d <- dcast(d, year + scenario  ~ `outcome`, value.var = "outcome_50%"
)

ggplot(d[year == 2049], aes(x = dpp, y = sii, col = scenario)) +
  geom_point() +
  expand_limits(x = 0) +
  xlab("Deaths prevented or postponed") +
  ylab("Absolute inequality reduction (deaths)") +
  labs(colour = "Scenario") +
  ggtitle("Reduction in deaths and absolute inequality: 2049") +
  geom_vline(aes(xintercept = 0)) + #add y-axis at 0
  geom_hline(aes(yintercept = 0)) #add x-axis at 0
if(plotsave){
  ggsave(filename = outpt_pth_cmpr_plt("dpp_ineq.png"))
}
}

#setahs in uner/over 65s
bigtab_cf[, agegrp_big := ifelse(age <= 65, "wrk_age", "over65")]
bigtab_cf_sum_65split <- bigtab_cf[, .(N = sum(N), atrisk = sum(atrisk)),
                           keyby = .(mc, scenario, year, agegrp_big, imd, state)]
setkey(bigtab_cf_sum_65split, mc,scenario, agegrp_big, imd, state, year)
bigtab_cf_sum_65split[, cumN := cumsum(N), keyby = .(mc, scenario, agegrp_big, imd, state)]


outstrata <- c("mc", "year", "scenario", "imd", "state", "agegrp_big")
d <- fn_mc(bigtab_cf_sum_65split, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above

ggplot(d[outcome %in% "cumN"], aes(x = year, y = `outcome_50%`, col = imd, linetype = agegrp_big)) +
    geom_line() +
  facet_grid(state ~ scenario, scales = "free")  +
  ylab("Deaths") +
  ggtitle("Cumulative deaths by imd & under/over 65")

#deaths prevented or postponed, by agegrp
setkey(bigtab_cf_sum_65split, mc,  scenario, agegrp_big, imd, state, year)
bl <- bigtab_cf_sum_65split[scenario == "baseline"]
bigtab_cf_sum_65split2 <- bigtab_cf_sum_65split[bl, on = c("mc", "year","agegrp_big", "imd", "state"), bl := i.cumN]
bigtab_cf_sum_65split2 <- bigtab_cf_sum_65split2[scenario != "baseline"]
bigtab_cf_sum_65split2[, ppdeaths := bl-cumN]


outstrata <- c("mc", "year", "scenario", "imd", "agegrp_big","state")
d <- fn_mc(bigtab_cf_sum_65split2, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above

ggplot(d[outcome == "ppdeaths" & year %in% c(2019, 2029, 2039, 2049) & agegrp_big == "wrk_age"],
       aes(x = year, y = `outcome_50%`, fill = imd)) +
  geom_col(position = "dodge") +
  facet_grid(state ~ scenario) +
  ylab("Deaths prevented or postponed") +
  ggtitle("Cumulative deaths prevented or postponed by imd in ages 65 and under")

ggplot(d[outcome == "ppdeaths"  & year %in% c(2019, 2029, 2039, 2049) & agegrp_big == "wrk_age"],
       aes(x = year, y = `outcome_50%`, fill = scenario)) +
  geom_col(position = "dodge") +
  facet_grid(state ~ imd, scales = "free") +
  ylab("Deaths prevented or postponed") +
  ggtitle("Cumulative deaths prevented or postponed by imd in over 65s")


#Making the equity plot
dt <- bigtab_cf_sum_65split2[year %in% c(2019, 2029, 2039, 2049) & state == "Overall",
                     .(ppdeaths = sum(ppdeaths), atrisk = sum(atrisk)),
                     keyby = .(imd, mc, year, scenario, agegrp_big)][
                       ,se := 0 # need to add this in...
                     ]

if(ineq){
# SHould be able to do this by calculating rigit scores
require(tidyverse)
sii_cf_age <-
  as.data.table(
    phe_sii(
      group_by(dt, mc, year, scenario, agegrp_big),
      imd,
      atrisk,
      value_type = 0,
      # for rates
      value = ppdeaths,
      se = se,
      confidence = 0.95,
      rii = TRUE,
      type = "standard"
    )
  )

sii_cf_age[, c(7:10) := NULL] # for now don't need the CIs
sii_cf_age <- sii_cf_age[dt[, .(dpp =sum(ppdeaths)), keyby = .(mc, year, scenario, agegrp_big)],
                 on = c("mc","year", "scenario", "agegrp_big")]

outstrata <- c("year", "scenario", "agegrp_big")
d <- fn_mc(sii_cf_age, outstrata, prbl, "outcome", "outcome_") #this then gives the median and various uncertainty intervals as defined by prbl above


d <- dcast(d, year + scenario  + agegrp_big ~ `outcome`, value.var = "outcome_50%"
)


d[, agegrp_big := factor(agegrp_big,
                         levels = c("wrk_age", "over65"),
                         labels = c("Working age", "Over 65"))]
#Multiplying by 100 to scale back up
d[, `:=` (dpp = dpp * 100,
          sii = sii * 100)]

d[, ]

ggplot(d[year == 2049], aes(x = dpp, y = sii, col = scenario, shape = agegrp_big)) +
  geom_point(size = 2) +
  expand_limits(x = 0) +
  xlab("Deaths prevented or postponed") +
  ylab("Absolute inequality reduction (deaths)") +
  labs(colour = "Scenario", shape = "Age group") +
  ggtitle("Reduction in deaths and absolute inequality and under/over 65s: 2049") +
  geom_vline(aes(xintercept = 0)) + #add y-axis at 0
  geom_hline(aes(yintercept = 0)) + #add x-axis at 0
  scale_y_continuous( labels = label_comma()) +
  scale_x_continuous( labels = label_comma()) +
  scale_color_brewer(type = "qual", palette = 6) +
  theme(legend.position="bottom")
if(plotsave){
  ggsave(filename = outpt_pth_cmpr_plt("dpp_ineq_65split.png"))
}

}
# bog standard cf
#overall by year
outstrata <- c("mc", "year", "scenario", "state", "gender")
d <- bigtab_cf[,
                .(cf_cr = sum(N)/ sum(atrisk)), keyby = outstrata]
d <- fn_mc(d, outstrata, prbl, "outcome", "cf_cr_") #this then gives the median and various uncertainty intervals as defined by prbl above


ggplot(d,
       aes(x = year, y = `cf_cr_50%` * 10000 , col = scenario)) +
  geom_smooth(se = FALSE, size = 0.5) +
  facet_grid(state ~ gender) +
  expand_limits(y = 0) +
  ylab("Crude case fatality") +
  ggtitle("Case fatality by scenario, state and gender")
if(plotsave){
  ggsave(filename = outpt_pth_cmpr_plt("cf_gender.png"))
}
# Standardised
# Standardising to esp
esp <- mk_esp()
strata <- c("mc","year", "gender", "imd", "agegrp5", "state", "scenario")
bigtab_cf <- fn_stnd(bigtab_cf, strata)

outstrata <- c("mc", "year", "scenario", "state")
d <- bigtab_cf[,
                .(cf_st = sum(N * wt_esp)/ sum(popsize_wtd)), keyby = outstrata]
d <- fn_mc(d, outstrata, prbl, "outcome", "cf_st_") #this then gives the median and various uncertainty intervals as defined by prbl above
ggplot(d,
       aes(x = year, y = `cf_st_50%` * 10000, col = scenario)) +
  geom_smooth(se = FALSE, size = 0.5) +
  facet_grid(rows = "state", scales = "free" ) +
  expand_limits(y = 0) +
  ylab("Standardised fatality")

if(plotsave){
  ggsave(filename = outpt_pth_cmpr_plt("cf_overall_st.png"))
}







# years in each state

yrs_in_state <- data.table()
for (scenario in sc_nm){
  tmptab <- read_fst(outpt_pth(paste0(scenario,"_panel.fst")),
                     as.data.table = T,
                     columns = c("gender", "bc", "mc", "year",
                                 "state_incCond", "state_BMM", "state_CMM"))[
                                   mc == 1]  # doing only for 1 iteration for now
  tmptab[, scenario := scenario]
  yrs_in_state <- rbind(yrs_in_state, tmptab[, scenario := scenario])
}

#counting incident years as half years
yrs_in_state_sum <- yrs_in_state[mc == 1, .(N = .N,
                                     Healthy = sum(state_incCond == 0 & state_BMM == 0 & state_CMM == 0 ) +
                                       sum(state_incCond == 1)/2,
                                     IncCond = sum(state_incCond == 2) + sum(state_incCond == 1)/2,
                                     BMM = sum(state_BMM == 2) + sum(state_BMM == 1)/2,
                                     CMM = sum(state_CMM == 2) + sum(state_CMM == 1)/2),
                                 keyby = .(year, scenario)]

setkey(yrs_in_state_sum, scenario, year)
yrs_in_state_sum[, `:=`( Healthy_cml = cumsum(Healthy),
                         IncCond_cml = cumsum(IncCond),
                         BMM_cml = cumsum(BMM),
                         CMM_cml = cumsum(CMM),
                         N_cml = cumsum(N)
                         ), keyby = .(scenario)]
yrs_in_state_sum[, `:=`( Healthy_cml_pp = cumsum(Healthy/N),
                         IncCond_cml_pp = cumsum(IncCond/N),
                         BMM_cml_pp = cumsum(BMM/N),
                         CMM_cml_pp = cumsum(CMM/N)
), keyby = .(scenario)]

#states aren't exclusive though ...
yrs_in_state_sum[, `:=`( Healthy_cml_prop = cumsum(Healthy_cml/N_cml),
                         IncCond_cml_prop = cumsum(IncCond_cml/N_cml),
                         BMM_cml_prop = cumsum(BMM_cml/N_cml),
                         CMM_cml_prop = cumsum(CMM_cml/N_cml)
), keyby = .(scenario)]

yrs_in_state_sum[year == 2049, -c("N", "Healthy", "IncCond", "BMM", "CMM")]




bigtab_inc_sum2_all <- bigtab_inc_sum[bl, on = c("mc", "year", "state"), bl := i.cumN]
bigtab_inc_sum2_all <- bigtab_inc_sum2_all[scenario != "baseline"]

bigtab_inc_sum2_all[, ppcases := bl-cumN]




## Mean life expectancies

lifecourse <- data.table()
for (scenario in scenario_nm){
  all.files <- list.files(path = outpt_pth(paste0("lifeyears")) ,pattern = ".csv.gz", full.names = T)
  tmptab <- rbindlist(lapply(all.files, fread))
  tmptab[, scenario := scenario]
  lifecourse <- rbind(lifecourse, tmptab[, scenario := scenario])
}


lifecourse[, imd := factor(imd)]

outstrata <- c("mc", "scenario", "imd" , "ageenter", "init_year")

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







outstrata <- c("mc", "scenario", "ageenter","init_year" )

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
fwrite(expect, outpt_pth_cmpr_plt("le.csv"), row.names = FALSE)

#IMD Diff
expect_diff <- data.table()

instrata <- c("mc", "scenario", "ageenter", "imd" ,"init_year" )
outstrata <- c("mc", "scenario", "ageenter","init_year" )

# LE
d <- lifecourse[, .(le = mean(max_age)), keyby = instrata]
imd <- dcast(d[imd %in% c(1,5)], mc + scenario + ageenter + init_year ~ imd, value.var = "le")
imd[, `:=` (abs = `5`-`1`, rel = `5` / `1` *100)][, c("5", "1") := NULL]
imd <- fn_mc(imd, outstrata, prbl, "outcome", "outcome")
expect_diff <- rbind(expect_diff, imd[, `:=`(type ="le_diff")])

# MMfree LE
d <- lifecourse[, .(le = mean(healthy + incCond)), keyby = instrata]
imd <- dcast(d[imd %in% c(1,5)], mc + scenario + ageenter + init_year ~ imd, value.var = "le")
imd[, `:=` (abs = `5`-`1`, rel = `5` / `1` *100)][, c("5", "1") := NULL]
imd <- fn_mc(imd, outstrata, prbl, "outcome", "outcome")
expect_diff <- rbind(expect_diff, imd[, `:=`(type ="mmfreele_diff")])

# CMMfree LE
d <- lifecourse[, .(le = mean(healthy + incCond + bmm)), keyby = instrata]
imd <- dcast(d[imd %in% c(1,5)], mc + scenario + ageenter + init_year ~ imd, value.var = "le")
imd[, `:=` (abs = `5`-`1`, rel = `5` / `1` *100)][, c("5", "1") := NULL]
imd <- fn_mc(imd, outstrata, prbl, "outcome", "outcome")
expect_diff <- rbind(expect_diff, imd[, `:=`(type ="cmmfreele_diff")])

fwrite(expect_diff, outpt_pth_cmpr_plt("le_imddiff.csv"), row.names = FALSE)


#BL Diff
expect_bldiff <- data.table()

instrata <- c("mc", "scenario", "ageenter" ,"init_year" )
outstrata <- c("mc", "scenario", "ageenter","init_year" )

# LE
d <- lifecourse[, .(le = mean(max_age)), keyby = instrata]
tmp <- d[scenario == "baseline"]
d[tmp, on = c("mc", "ageenter", "init_year"), bl := i.le]
d[, `:=` (abs = le - bl, rel = le / bl *100)][, c("le", "bl") := NULL]
d <- fn_mc(d[scenario != "baseline"], outstrata, prbl, "outcome", "outcome")
expect_bldiff <- rbind(expect_bldiff, d[, `:=`(type ="le_diff")])

# MMfree LE
d <- lifecourse[, .(le = mean(healthy + incCond)), keyby = instrata]
tmp <- d[scenario == "baseline"]
d[tmp, on = c("mc", "ageenter", "init_year"), bl := i.le]
d[, `:=` (abs = le - bl, rel = le / bl *100)][, c("le", "bl") := NULL]
d <- fn_mc(d[scenario != "baseline"], outstrata, prbl, "outcome", "outcome")
expect_bldiff <- rbind(expect_bldiff, d[, `:=`(type ="mmfreele_diff")])

# CMMfree LE
d <- lifecourse[, .(le = mean(healthy + incCond + bmm)), keyby = instrata]
tmp <- d[scenario == "baseline"]
d[tmp, on = c("mc", "ageenter", "init_year"), bl := i.le]
d[, `:=` (abs = le - bl, rel = le / bl *100)][, c("le", "bl") := NULL]
d <- fn_mc(d[scenario != "baseline"], outstrata, prbl, "outcome", "outcome")
expect_bldiff <- rbind(expect_bldiff, d[, `:=`(type ="cmmfreele_diff")])

fwrite(expect_bldiff, outpt_pth_cmpr_plt("le_bldiff.csv"), row.names = FALSE)

#BL Diff by imd
expect_bldiff_imd <- data.table()

instrata <- c("mc", "scenario", "ageenter", "imd" ,"init_year" )
outstrata <- c("mc", "scenario", "ageenter", "imd" , "init_year" )

# LE
d <- lifecourse[, .(le = mean(max_age)), keyby = instrata]
tmp <- d[scenario == "baseline"]
d[tmp, on = c("mc", "imd" , "ageenter", "init_year"), bl := i.le]
d[, `:=` (abs = le - bl, rel = le / bl *100)][, c("le", "bl") := NULL]
d <- fn_mc(d[scenario != "baseline"], outstrata, prbl, "outcome", "outcome")
expect_bldiff_imd <- rbind(expect_bldiff_imd, d[, `:=`(type ="le_diff")])

# MMfree LE
d <- lifecourse[, .(le = mean(healthy + incCond)), keyby = instrata]
tmp <- d[scenario == "baseline"]
d[tmp, on = c("mc", "imd" , "ageenter", "init_year"), bl := i.le]
d[, `:=` (abs = le - bl, rel = le / bl *100)][, c("le", "bl") := NULL]
d <- fn_mc(d[scenario != "baseline"], outstrata, prbl, "outcome", "outcome")
expect_bldiff_imd <- rbind(expect_bldiff_imd, d[, `:=`(type ="mmfreele_diff")])

# CMMfree LE
d <- lifecourse[, .(le = mean(healthy + incCond + bmm)), keyby = instrata]
tmp <- d[scenario == "baseline"]
d[tmp, on = c("mc", "imd" , "ageenter", "init_year"), bl := i.le]
d[, `:=` (abs = le - bl, rel = le / bl *100)][, c("le", "bl") := NULL]
d <- fn_mc(d[scenario != "baseline"], outstrata, prbl, "outcome", "outcome")
expect_bldiff_imd <- rbind(expect_bldiff_imd, d[, `:=`(type ="cmmfreele_diff")])

#fwrite(expect_bldiff_imd, outpt_pth_cmpr_plt("le_bldiff_imd.csv"), row.names = FALSE)




#For the plot
tmp <- expect[type == "mmfree_le" & init_year == 2019 & ageenter == 30,
              .(scenario, imd, type, years = `outcome50%`)]
tmp2 <- expect[type == "cmmfree_le" & init_year == 2019 & ageenter == 30,
              .(scenario, imd, type, years = `outcome50%`)]
tmp2[tmp, on = c("scenario", "imd"), mmfree := i.years][, years := years - mmfree][, mmfree := NULL]
tmp <- rbind(tmp,tmp2)

tmp[, scenario := factor(scenario,
                         levels = c("baseline","scenario1","scenario2","scenario3","scenario4","scenario5"),
                         labels = c("0 Baseline","1 Targeted", "2 Universal +\nfocus on gap",
                                    "3 Redistributive", "4 Proportionate\nuniversalism",
                                    "5 Inequalities\nremoval"))]
tmp[, type := factor(type,
                     levels = c( "cmmfree_le", "mmfree_le"),
                     labels = c( "Complex", "Basic"))]
ggplot(tmp[imd != "all"],
       aes(x = imd, y = years, fill = type)) +
  geom_col(position = "stack") +
  facet_grid(cols =vars( scenario)) +
  scale_fill_brewer(type = "qual", palette = "Paired") +
  ylab("Years lived without multimorbidity") +
  xlab("IMD quintile") +
  labs(fill = "Multimorbidity type") +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10)) +
  guides(fill = guide_legend(reverse = TRUE))

if(plotsave){
#  ggsave(filename = outpt_pth_cmpr_plt("mmfreeyears.png"))
  ggsave2(outpt_pth_cmpr_plt("mmfreeyears.svg") , scale = 0.7, width = 14, height = 8)

}

#


## Mean life expectancies  by sex
outstrata <- c("mc", "scenario", "imd" , "ageenter", "init_year", "gender")
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



outstrata <- c("mc", "scenario" , "ageenter","init_year", "gender")

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
fwrite(expect_sex, outpt_pth_cmpr_plt("le_sex.csv"), row.names = FALSE)

expect_sex_l <- melt(expect_sex, id.vars = c("scenario", "gender","type"),
                 measure.vars = c("overall", "1", "2", "3", "4", "5"),
                 variable.name = "imd",
                 value.name = "years"
)

tmp <- expect_sex_l[type == "cmmfree_le"]
tmp[expect_sex_l[type == "mmfree_le"], on = c("scenario", "gender", "imd"), mmfree_le := i.years]
tmp[, years := years- mmfree_le]
tmp <- rbind(tmp[, .(scenario, gender, type, imd, years)],
             tmp[, .(scenario, gender, type = "mmfree_le" , imd, years = mmfree_le)])
fn_scenario_lbls2(tmp)
tmp[, type := factor(type,
                     levels = c( "cmmfree_le", "mmfree_le"),
                     labels = c( "Complex", "Basic"))]
ggplot(tmp[imd != "overall"],
       aes(x = imd, y = years, fill = type)) +
  geom_col(position = "stack") +
  facet_grid(scenario ~ gender) +
  scale_fill_brewer(type = "qual", palette = "Paired") +
  ylab("Years lived without multimorbidity") +
  xlab("IMD quintile") +
  labs(fill = "Multimorbidity type") +
  theme(legend.position = "bottom")

if(plotsave){
  ggsave(filename = outpt_pth_cmpr_plt("mmfreeyears_sex.png"))
}


