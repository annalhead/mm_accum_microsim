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
library(foreach)
library(doParallel)
registerDoParallel(cores=3)

# Number of simulations
n_sim <- 200

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


# Pop to be simulated
patients <- read_patients()
ll <- read_ll()
outstrata <- c("mc", "year", "gender", "age","imd", "region")
tt <- data.table( # Auxiliary table
  y = c(2020L:2050L),
  year = 2020L:2050L)

foreach(i = 1:n_sim, .combine = rbind) %dopar% process_results(i)
