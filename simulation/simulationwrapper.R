# Basic simulation wrapper

# This script is a basic R wrapper for the c++ simulation code
# Required functions are in the fn.R file

# Set up ----
library(data.table)
data.table::setDTthreads(5) #this is so that don't use all the processors
library(fst)
fst::threads_fst(10)
library(Rcpp)
library(dqrng)
library(qs)

sim_path <-
  function(x = character(0))
    paste0("./simulation/", x)



source(sim_path("/fn.R"))
sourceCpp(sim_path("/mmsim_mod_bc.cpp"), cacheDir = project_path("./.cache")) #this is w deciles

# Read in the inputs ----
#read in the lookup table
ll <- read_ll()

# Pop to be simulated
patients <- read_patients()



# The simulation ----

# the function run_sim takes three inputs:
# 1. the population to be simulated
# 2. the lookup table of distributions of transition times
# 3. a seed for drawing the random number

n_sim <- 10 #number of iterations
sim_output <- data.table()
for (i in 1:n_sim){
  tmptab <- run_sim(patients, ll, 41L + i)
  sim_output <- rbind(sim_output, tmptab[, mc := i])
  print(paste0("completed iteration #",i))
}
