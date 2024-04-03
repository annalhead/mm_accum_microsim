# mm_accum_microsim

This repository contains the code for the development and application of a dynamic microsilmulation to explore socioeconomic inequalities in multimorbidity. 
This version is the code used for the paper	 *Socioeconomic inequalities in accumulation of multimorbidity in England from 2019 to 2049: a microsimulation projection study* Head A, Birkett M, Fleming K, Kypridemos C, O'Flaherty M. Lancet Public Health 2024; 9: e231–39. Doi: 10.1016/S2468-2667(24)00028-8

This work forms part of my Public Health PhD: *Socioeconomic inequalities in multimorbidity: an epidemiological and microsimulation study*
This work was developed with support from Dr Chris Kypridemos, Dr Max Birkett, Dr Kate Fleming, and Prof Martin O'Flaherty. 

In short, the microsimulation uses primary care data from CPRD Aurum to estimate transition times between 4 health states and death, and then applies these to ONS 2019 population estimates, and projects them forwards for 30 years. 

To look at the accumulation of multimorbidity, we defined a multistate model with four sequential health states: 
1) healthy, 
2) one long-term condition, 
3) basic multimorbidity (2 or more chronic conditions)
4) complex multimorbidity (3 or more chronic conditions across 3 or more body systems)

Individuals could die from each of the four health states, resulting in a total of seven unidirectional transitions. 
 - T1:Healthy--> initial chronic condition
 - T2:Healthy--> death
 - T3:Initial chronic condition--> basic multimorbidity
 - T4:Initial chronic condition--> death
 - T5:Basic multimorbidity--> complex multimorbidity
 - T6:Basic multimorbidity--> death
 - T7:Complex multimorbidity--> death
 
This project builds on from the data cleaning and phenotyping outlined in the repositories: https://github.com/annalhead/CPRD_multimorbidity_codelists and https://github.com/annalhead/CPRD_multimorbidity_trends, and the work presented here: https://www.thelancet.com/journals/lanhl/article/PIIS2666-7568(21)00146-X/fulltext 

Currently, this repository only contains the analysis and simulation code, and not the model input files. If you have any questions, please email me at anna.head2@liverpool.ac.uk 

## Repository summary

There are 4 main folders: inputs, simulation, scenarios, auxil 

### Inputs
Together, these files prepare the inputs needed for the microsimulation model 
1. Datasetup.R: This takes cleaned CPRD on dates of diagnosis, and rearranges it into counting process format so that each line is 1 transition for each individual. 
2. Hazplots_leftcensored.Rmd: This compares parametric distributions to the baseline hazard plots, and is used for determining which distribution to use for each transition in parametric survival analysis
3. modelfitscript.R: This fits parametric survival models for each of the 7 transitions 
4. generate_lookup_tables_condageenter.R: This creates a look-up table of predicted transition times for each transition and combination of covariates and single year of age at entry 
5. poptab_wproj.R: This creates the simulant population, based on ONS population estimates and projections 

### Simulation 
This folder contains the c++ code for the simulation model alongside the supporting R code. The simulation uses the output of the files from the 'inputs' folder
1. mmsim_mod_bc.cpp: the simulation c++ code
2. fn.R: supplementary functions for loading the required data, running the model, and processing and saving the results
3. simulationwrapper.R: a simple R wrapper for running the model 

### Scenarios 
This folder contains the R scripts for running the various theoretical scenarios I included in my PhD thesis, plus the validation/calibration code 
1. vldtn_clbrtn.R: uses a subset of CPRD data with 10 years of follow-up to validate and calibrate the model. Final calibration factors are saved in the fn.R file above
2. baseline.R: applies the microsimulation model to the ONS 2019 estimates for 30-90 year olds from 2019-2049 
3. baseline_param.R: as baseline.R, but with added parameter uncertainty for the transition times between states 

### Auxil 
Additional files 
/

