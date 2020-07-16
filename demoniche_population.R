demoniche_population <-
  function(Matrix_projection,
           Matrix_projection_var,
           n,
           populationmax,
           K = NULL,
           Kweight = BEMDEM$Kweight,
           onepopulation_Niche,
           sumweight,
           noise,
           prob_scenario,
           prev_mx,
           transition_affected_demogr,
           transition_affected_niche,
           transition_affected_env,
           env_stochas_type,
           yx_tx)
  {
    prob_scenario_noise <-
      c(prob_scenario[prev_mx[yx_tx]] * noise
        , 1 - (prob_scenario[prev_mx[yx_tx]] * noise))
    
    rand_mxs     <-
      sample(1:2, 1, prob = prob_scenario_noise, replace = TRUE)
    one_mxs      <-
      Matrix_projection[, rand_mxs]    # select one matrix rand_mxs = 1
    prev_mx[yx_tx + 1]    <-
      rand_mxs # saves the choice of matrix
    
    ## Modify the chosen matrix with habitat suitability values and demographic stochasticity
    if (Matrix_projection_var[1] != FALSE)
    {
      one_mxs_var <- one_mxs * (Matrix_projection_var[, rand_mxs])
      
      if (is.numeric(transition_affected_niche)) {
        ## Niche values
        one_mxs[transition_affected_niche] <-
          one_mxs[transition_affected_niche] * onepopulation_Niche
      }
      
      if (is.numeric(transition_affected_env)) {
        ## environmental stochasticity
        switch(
          EXPR = env_stochas_type,
          normal =
            one_mxs[transition_affected_env] <-
            # normal distribution
            rnorm(length(one_mxs[transition_affected_env]), mean = one_mxs[transition_affected_env],
                  sd = one_mxs_var[transition_affected_env])
          ,
          lognormal = one_mxs[transition_affected_env] <-
            # lognormal distribution
            rlnorm(
              length(one_mxs[transition_affected_env]),
              meanlog = one_mxs[transition_affected_env],
              sdlog = one_mxs_var[transition_affected_env]
            )
        )
      }
    }
    
    
    one_mxs[one_mxs < 0] <-
      0 # Matrix values cannot be negative
    #  one_mxs[one_mxs == NA]
    
    A <-
      matrix(one_mxs,
             ncol = length(n),
             nrow = length(n),
             byrow = FALSE)
    
    # To check if surivial and persistence columns sum to more than one.
    Atest <- A
    Atest[1, ][-1] <- 0 # Remove number of recruits.
    #  colSums(Atest)
    if (sum(colSums(Atest) > 1)) {
      for (zerox in which(colSums(Atest) > 1))
        # zerox = 5
      {
        Atest[, zerox] <- Atest[, zerox] / sum(Atest[, zerox])
      }
      A[-1, ] <- Atest[-1, ] # the new A, <1
    } # end if
    
    
    n <-
      as.vector(A %*% n)                          # Matrix multiplication!
    
    n <-
      floor(n) # If the number of individuals is less than one, replace with zero
    
    ################################################################################
    # Allee effect here
    # Sample number of offsping from poisson distribution!
    
    # if(is.numeric(transition_affected_demogr)) { # demographic stochasticity
    
    #           one_mxs[transition_affected_demogr[transition_affected_demogr %in% 2:length(BEMDEM$stages)]]
    
    #         fecundity_index <- transition_affected_demogr[transition_affected_demogr %in% 2:length(BEMDEM$stages)]
    #rpois(5, lambda = one_mxs[transition_affected_demogr[transition_affected_demogr %in% 2:length(BEMDEM$stages)]])# Poisson for fecundities
    #          rpois(length(fecundity_index), lambda = one_mxs[fecundity_index]) # Poisson for fecundities
    
    #  rnorm(5, mean = one_mxs[transition_affected_demogr[transition_affected_demogr %in% 2:length(BEMDEM$stages)]],
    #         sd = one_mxs_sd[transition_affected_demogr[transition_affected_demogr %in% 2:length(BEMDEM$stages)]])
    
    #one_mxs[transition_affected_demogr] <-
    # rpois(one_mxs[transition_affected_demogr], lambda = ) * onepopulation_Niche
    #   survival_index <- transition_affected_demogr[transition_affected_demogr = fecundity_index]
    #    survival_index <- diag(matrix(one_mxs, ncol = length(n), byrow = FALSE))
    
    #    switch(EXPR = K,
    #             Ricker = r_mx*(sum(n)*((K - sum(n))/K))
    #             Scramble =
    ################################################################################
    
    
    #  Simple density dependence
    #  if population size exceeds populationmax, reduce population to populationmax
    if (sum(n) > 0) {
      if (is.numeric(populationmax)) {
        if (sum(n * Kweight) > populationmax)
        {
          n <- n * (populationmax / sum(n * sumweight))
        }
      }
    }
    return(n)
    
  }
