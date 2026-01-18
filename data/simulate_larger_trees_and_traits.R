source('helper_simulation_methods.R')
source('helpful_phyl_methods.R')
library(phytools)

get_tree <- function(include.extinct, birth,death,number_of_extant_taxa){
  
  SimTree <- phytools::pbtree(b = birth, d = death, n = number_of_extant_taxa, extant.only = !include.extinct)
  if (!is.null(SimTree) && class(SimTree) == "phylo" && length(phytools::getExtant(SimTree))==number_of_extant_taxa && 
      length(SimTree$tip.label)<1500) {
    
    return(list("tree" = SimTree, "birth" = birth))
  } else {
    # To avoid recursive stack overflow..
    get_tree(include.extinct, birth*1.1,death,number_of_extant_taxa)
  }
}

get_bhisse_sample <- function(hidden.traits,include.extinct, number_of_extant_taxa){
  # BISSE and HiSSE internally generate trees as the state influences the birth and death rates.
  
  if(hidden.traits==1){
    # Heterogeneous Transition Rate Models
    # https://revbayes.github.io/tutorials/sse/hisse
    # a HiSSE model with 1 hidden binary trait (2 hidden states) and 1 observed binary trait (2 observed states), totaling 4 states. 
    # This allows for interactions between hidden and observed traits in diversification rates.
    # For 4 states (e.g., 0A, 0B, 1A, 1B), define:
    death_rates = runif(4, min = 0, max = 1)
    birth_rates = runif(4, min = 0, max = 1)
    
  }
  if(hidden.traits==0){
    # BiSSE (Binary-State Speciation and Extinction)
    # https://revbayes.github.io/tutorials/sse/bisse-intro.html#bisse_theory
    # The BiSSE model (in {diversitree}) links a binary trait (0 or 1) to different birth/death rates.
    
    # Define birth/death rates depending on trait state
    death_rates = runif(2, min = 0, max = 1)
    birth_rates = runif(2, min = 0, max = 1)
  }
  
  turnover.rates <- death_rates+birth_rates  # λ + μ for each of 4 states
  eps.values <- death_rates/birth_rates     # μ/λ ratios for each state
  
  # Get indices for transition rates. Allow transition among hidden categories to vary.
  transition.rates <- hisse::TransMatMakerHiSSE(hidden.traits =hidden.traits, cat.trans.vary = TRUE)
  for(i in 1:6){
    transition.rates[transition.rates==i]<-runif(1, min = 0, max = 1)
  }
  simulated.result <- hisse::SimulateHisse(turnover.rates, eps.values, 
                                           transition.rates, max.taxa=number_of_extant_taxa, x0=0)
  hisse_tree = hisse::SimToPhylo(simulated.result, include.extinct=include.extinct, drop.stem=TRUE)
  # plot(hisse_tree)
  
  # # Define colors for binary states
  # trait_colors <- ifelse(traits == 1, "red", "blue")
  # 
  # # Plot tree with colored tip labels
  # plot(hisse_tree, tip.color = trait_colors, cex = 1.2)
  
  if (!is.null(hisse_tree) && class(hisse_tree) == "phylo" && length(phytools::getExtant(hisse_tree))==number_of_extant_taxa
      && length(hisse_tree$tip.label)<5*number_of_extant_taxa) {
    
    if(hidden.traits==1){
      # Convert states back into observed binary character
      # Extract tip states (0, 1, 2, 3)
      tip_states <- hisse_tree$tip.state
      
      # Map to observed binary trait (0 or 1)
      observed_traits <- tip_states %% 2  # 0,2 → 0, 1,3 → 1
      
      # Overwrite tip names (careful!)
      hisse_tree$tip.state <- observed_traits
    }
    
    traits = hisse_tree$tip.state
    traits = traits[match(hisse_tree$tip.label, names(traits))]
    ground_truth = data.frame(traits)
    param_dataframe = data.frame(turnover.rates=c(turnover.rates),eps.values=c(eps.values), transition.rates=c(transition.rates))
    
    return(list(tree=hisse_tree, FinalData= ground_truth, Dataframe=param_dataframe))
  } else {
    get_bhisse_sample(hidden.traits, include.extinct, number_of_extant_taxa)
  }
  
  
}

get_BM_T_sample <- function(simulation_path, with_trend){
  tree = ape::read.tree(file.path(simulation_path, 'tree.tre'))
  if(with_trend){
  # Brownian Motion with a Trend (BM + Trend)
  #Similar to Theta selection in TDIP (https://github.com/Matgend/TDIP/blob/62c6655f7da66b0f89a48554a8eba7e697ea36eb/R/data_simulator.R#L262)

    mu = runif(1, min=-10, max=10)
  }else{
    mu=0

  }

  trait_BM_trend <- fastBM(tree, sig2=1, a=0, mu=mu)  # sig2 = BM variance, mu = trend strength
  trait_BM_trend_scaled = scale(trait_BM_trend)
  names(trait_BM_trend_scaled) <- names(trait_BM_trend)
  # plot(trait_BM_trend_scaled, ylab="Trait Value", xlab="Species", main="BM with a Trend")
  # phenogram(tree, trait_BM_trend_scaled, fsize=0.8, main="Trait Evolution under BM with a Trend")

  ground_truth = data.frame(trait_BM_trend_scaled)

  param_dataframe = data.frame(mu=c(mu))

  min = min(ground_truth$trait_BM_trend_scaled)
  max = max(ground_truth$trait_BM_trend_scaled)
#   print('########## mu, min max')
#   print(mu)
#   print(min)
#   print(max)
#   print('##########')
  out = list(tree=tree, FinalData= ground_truth, Dataframe=param_dataframe)

  return(out)
}
number_of_taxa = 1000
for(i in 1:10){
  print(i)

    ## birth death tree
    this_sim_path = file.path("simulations", 'ultrametric', 'large_trees', i)
    death_rate = runif(1, min = 0, max = 1)
    birth_rate = runif(1, min = 0, max = 1)


    a = get_tree(FALSE, birth_rate,death_rate, number_of_taxa)
    testit::assert(ape::is.ultrametric(a$tree))
    testit::assert(length(a$tree$tip.label)==1000)

    tree = a$tree
    birth = a$birth
    output_tree(this_sim_path, tree, birth, death_rate)

       ## bisse tree
    ultra_bisse_sample = get_bhisse_sample(0,FALSE,number_of_taxa)
      sim_path= file.path("simulations", 'ultrametric', 'BiSSE_large_trees', i)
      output_tree(sim_path, ultra_bisse_sample$tree)
      output_simulation(sim_path, ultra_bisse_sample,  'large_BiSSE')
}

# Add simulated traits
for(i in 1:10){
  print(i)

    sim_path = file.path("simulations", 'ultrametric', 'large_trees', i)
    bmt_sample = get_BM_T_sample(sim_path,TRUE)
    output_simulation(sim_path, bmt_sample, 'large_BMT')

}

### Add PEMs
for(iteration in 1:10){
  print(iteration)
    missingness_types = c('mcar')
    for (missing_type in missingness_types) {
    setup_ = set_up('ultrametric', 'large_BiSSE', iteration, missing_type)
    output_folder = get_input_data_paths('ultrametric', 'large_BiSSE', iteration)$value_path
    get_PEMS(setup_, output_folder, missing_type)

    setup_ = set_up('ultrametric', 'large_BMT', iteration, missing_type)
    output_folder = get_input_data_paths('ultrametric', 'large_BMT', iteration)$value_path
    get_PEMS(setup_, output_folder, missing_type)
}
}
