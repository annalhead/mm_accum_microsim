`.sourceCpp_1_DLLInfo` <- dyn.load('/mnt/alhead/UoL/CPRD2019mm/mm_accum_microsim/simulation/.cache/sourceCpp-x86_64-pc-linux-gnu-1.0.10/sourcecpp_3f40fb16fb6d88/sourceCpp_2.so')

parallel_random_matrix <- Rcpp:::sourceCppFunction(function(n, m, seed, ncores) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_parallel_random_matrix')
mmsim <- Rcpp:::sourceCppFunction(function(pop, dfTransitionTimeQuantilesBySimulantType, rank, Mcalib_t12 = 1, Mcalib_t34 = 1, Mcalib_t56 = 1, Mcalib_t7 = 1, Fcalib_t12 = 1, Fcalib_t34 = 1, Fcalib_t56 = 1, Fcalib_t7 = 1, Mnudge_t5 = 1) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_mmsim')

rm(`.sourceCpp_1_DLLInfo`)
