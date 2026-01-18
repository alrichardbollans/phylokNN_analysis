library(doParallel)
library(foreach)

num_cores <- 12#set to 32 as this function seems to find all threads on cluster node detectCores()  # Use all cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

number_of_simulation_iterations = 10
foreach(iter = 1:number_of_simulation_iterations) %dopar% {
# for(iter in 1:number_of_simulation_iterations){
  repo_path = Sys.getenv('KEWSCRATCHPATH')
  source(file.path(repo_path, 'phyloKNN_analysis', 'imputation','R_continuous_imputation_helper_functions.R'))
  
  missingness_types = c('mcar')
  cases = c('ultrametric')
  binary_ev_models = c('large_BiSSE')
  continuous_ev_models = c('large_BMT')
  print(iter)
  for (missing_type in missingness_types) {
    for(simulation_ev_model in binary_ev_models){# Keep the inner loop sequential

            for(case in cases){
                run_corHMM_models(case, simulation_ev_model, iter, missing_type, 'binary')
             }


    }

     for(simulation_ev_model in continuous_ev_models){# Keep the inner loop sequential

            for(case in cases){
                run_phylopars_models(case, simulation_ev_model, iter, missing_type, 'continuous')
             }



    }
    
  }
}
