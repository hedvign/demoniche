demoniche_dispersal <-
  function(seeds_per_population,
           fraction_LDD,
           fraction_SDD,
           dispersal_probabilities,
           dist_latlong,
           neigh_index)
  {
    seeds_per_population_migrate_LDD <-
      seeds_per_population * fraction_LDD # how many seeds each patch lost.
    
    seeds_per_population_migrate_SDD <-
      seeds_per_population * fraction_SDD
    
    # seeds_per_population_stay <- seeds_per_population - seeds_per_population_migrate
    
    seeds_per_population_new_SDD <-
      seeds_per_population_new_LDD <- rep(0, length(seeds_per_population))
    
    ################################################################################
    if (fraction_SDD > 0)
      # Dispersal to contigous cells
      ################################################################################
    {
      # Create a subset of only occupied patches.
      source_patches <- which(seeds_per_population_migrate_SDD > 0)
      
      for (px_orig in source_patches)
        # for each origin patch    px_orig = 2
      {
        for (pxdisp_new in 1:length(seeds_per_population_migrate_SDD))
          # for each new possible patch    pxdisp_new = 2   px_orig =  1
        {
          #  if(px_orig == pxdisp_new) break
          
          # print(paste( px_orig, pxdisp_new))
          if (dist_latlong[pxdisp_new, px_orig] == neigh_index[1])
          {
            seeds_per_population_new_SDD[pxdisp_new] <-
              # where do they go?
              seeds_per_population_new_SDD[pxdisp_new] + (seeds_per_population_migrate_SDD[px_orig] * 0.2) # these are the dispersing seeds          #  addseeds[pxdisp_new] +
            #   print(paste("dispersal contingous: px_orig", px_orig, "to pxdisp_new",pxdisp_new, sep = " "))
          }
          if (length(neigh_index) == 2) {
            if (dist_latlong[pxdisp_new, px_orig] == neigh_index[2])
            {
              seeds_per_population_new_SDD[pxdisp_new] <-
                seeds_per_population_new_SDD[pxdisp_new] + (seeds_per_population_migrate_SDD[px_orig] * 0.05)        #addseeds[pxdisp_new] +
              #     print(paste("dispersal sideways: px_orig", px_orig, "to pxdisp_new", pxdisp_new, sep = " "))
            }
          }
        } # end pxdisp_new
      } # end px_orig
      
      # return
      
      #  plot( seeds_per_population-ret)
    }
    ################################################################################
    ################################################################################
    
    if (fraction_LDD > 0) {
      # If colonization of new patches
      #  dist_populations <-   dist_populations[1:10,1:10]
      
      # the probability that seeds reach the other patch
      # add <- (1- colSums(probs)) / (nrow(probs) - 1)   # each column has to sum one because no seeds 'die'
      # for(ix in 1:ncol(probs)){ probs[ix,] <- probs[ix,] + add }
      
      #        seeds_per_population_migrate2 <-       seeds_per_population_migrate #[891:900]
      #        probs[,1]
      
      #    dispersal_probabilities2 <-    dispersal_probabilities *1000000000000000000000000000
      
      seeds_per_population_new_LDD <-
        as.vector(dispersal_probabilities %*% seeds_per_population_migrate_LDD)
      #         plot(seeds_per_population_migrate3)
      
      
      ################################################################################
    } # end if numeric LDD
    seeds_stay <-
      (
        seeds_per_population - seeds_per_population_migrate_SDD -  seeds_per_population_migrate_LDD
      )
    
    return(seeds_stay + seeds_per_population_new_SDD + seeds_per_population_new_LDD)
    
  }   # end fcn
