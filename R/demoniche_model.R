#' @title Run spatially explicit demographic model
#' @export demoniche_model
#'
#' @param modelname Name of \code{list} object containing model parameters
#' @param Niche TRUE/FALSE that specifies whether to run the simulation accounting for effects of Niche values
#' @param Dispersal TRUE/FALSE that specifies whether to allow long-distance dispersal
#' @param repetitions Number of replicate simulations to carry out (default = 1000)
#' @param foldername Path name of folder where results of simulations will be written to.
#'
#' @return Something very special
#'
demoniche_model <-
  function(modelname,
           Niche,
           Dispersal,
           repetitions = 1000,
           foldername)
  {
    BEMDEM <-
      get(modelname, envir = .GlobalEnv)   #assign("BEMDEM", get(modelname))

    #load(paste(modelname, ".rda", sep = ""))
    #  rm(Hmontana)
    # require(sp)
    #require(popbio)
    #require(lattice)

    #### Create vectors ############################################################
    Projection <- array(
      0,
      # Projection is where all the results go
      dim =
        c(
          BEMDEM$no_yrs,
          length(BEMDEM$stages),
          nrow(BEMDEM$Niche_ID),
          length(BEMDEM$years_projections)
        )
      ,
      dimnames =
        list(
          paste("timesliceyear", 1:BEMDEM$no_yrs, sep = "_"),
          c(paste(BEMDEM$stages)),
          BEMDEM$Niche_ID[, "Niche_ID"],
          paste(BEMDEM$years_projections)
        )
    )

    eigen_results <-
      vector(mode = "list" , length(BEMDEM$list_names_matrices))
    names(eigen_results) <- unlist(BEMDEM$list_names_matrices)

    yrs_total <- BEMDEM$no_yrs * length(BEMDEM$years_projections)

    population_sizes <- array(
      NA,
      # per matrix!
      dim = c(yrs_total, length(BEMDEM$list_names_matrices), repetitions),
      dimnames = list(
        paste("year", 1:yrs_total, sep = ""),
        BEMDEM$list_names_matrices,
        paste("rep", 1:repetitions, sep = "_")
      )
    )

    population_results  <- array(
      1:200,
      # for all matrices!
      dim = c(yrs_total, 4, length(BEMDEM$list_names_matrices)),
      dimnames = list(
        paste("year", 1:yrs_total, sep = ""),
        c("Meanpop", "SD", "Max", "Min"),
        paste(BEMDEM$list_names_matrices)
      )
    )

    metapop_results <- array(
      NA,
      # per matrix!
      dim = c(yrs_total, length(BEMDEM$list_names_matrices), repetitions),
      dimnames = list(
        paste("year", 1:yrs_total, sep = ""),
        BEMDEM$list_names_matrices,
        paste("rep", 1:repetitions, sep = "_")
      )
    )

    simulation_results <- array(
      NA,
      dim = c(
        length(BEMDEM$list_names_matrices),
        7 + length(BEMDEM$years_projections)
      ),
      dimnames =
        list(
          BEMDEM$list_names_matrices,
          c(
            "lambda",
            "stoch_lambda",
            "mean_perc_ext_final",
            "initial_population_area",
            "initial_population",
            "mean_final_pop",
            "mean_no_patches_final",
            paste("EMA", BEMDEM$years_projections)
          )
        )
    )

    EMA <- array(
      0,
      dim = c(
        repetitions,
        length(BEMDEM$list_names_matrices),
        length(BEMDEM$years_projections),
        2
      ),
      dimnames =
        list(
          paste("rep", 1:repetitions, sep = "_"),
          BEMDEM$list_names_matrices,
          BEMDEM$years_projections,
          c("EMA", "No_populations")
        )
    )

    population_Niche <- rep(1, nrow(BEMDEM$Niche_ID))

    simulation_results[, "initial_population_area"] <-
      sum(BEMDEM$Orig_Populations[, "area_population"])
    # original population size
    simulation_results[, "initial_population"] <-
      round(sum(colSums(BEMDEM$n0_all) * BEMDEM$sumweight), 0)


    dir.create(paste(getwd(), "/", foldername, sep = ""), showWarnings = FALSE)




    #### Start repetitions #########################################################

    for (rx in 1:repetitions)
      # tx = 1   rx = 1
    {
      print(paste("Starting projections for repetition:", rx), quote = FALSE)

      for (mx in 1:length(BEMDEM$list_names_matrices))
        # selects two matrices for the simulations   mx = 1
      {
        print(paste(
          "Projecting for scenario/matrix:",
          (BEMDEM$list_names_matrices)[mx]
        ), quote = FALSE)

        yx_tx <- 0 # restart counter

        # this selects two matrices, the first basic matrix and the projection one mx = 1
        Matrix_projection <-
          cbind(BEMDEM$matrices[, 1], (BEMDEM$matrices[, mx]))

        # this selects two standard deviation matrices, the first basic matrix and the projection one (column mx)

        if (BEMDEM$matrices_var[1] != FALSE)
          # If matrices are included or not.
        {
          if (ncol(BEMDEM$matrices_var) > 1) {
            Matrix_projection_var <-
              cbind(BEMDEM$matrices_var[, 1], (BEMDEM$matrices_var[, mx]))
          } else {
            Matrix_projection_var <-
              cbind(BEMDEM$matrices_var[, 1], (BEMDEM$matrices_var[, 1]))
          }
        } else
        {
          Matrix_projection_var <- FALSE
        }

        prev_mx <- rep(1, times = yrs_total + 1)

        for (tx in 1:length(BEMDEM$years_projections))
          # selects which time-slice of the simulation tx = 1
        {
          if (Niche == TRUE) {
            # Niche values
            population_Niche <-
              BEMDEM$Niche_values[, tx]
          }
          for (yx in 1:BEMDEM$no_yrs)
          {
            yx_tx <- yx_tx + 1 # add one year to counter

            ####################################### PATCH PROJECTION when tx = 1 ###########     tx=1  px=1   tx = 2
            ################################################################################     yx=2  yx = 1

            if (tx == 1 && yx == 1)
            {
              n0s <-
                BEMDEM$n0_all[rowSums(BEMDEM$n0_all) > 0, ] # take stage vector from intial population, where over zero
              n0s_ID <-
                which(rowSums(BEMDEM$n0_all) > 0) # which populations/rows are above zero
            } else {
              # Second year running the model
              if (tx != 1 &&
                  yx == 1) {
                # when going to next time step, yx == 1

                n0s <-
                  t(Projection[BEMDEM$no_yrs, , colSums(Projection[BEMDEM$no_yrs, , , tx -
                                                                     1]) > 0, tx - 1])
                n0s_ID <-
                  which(colSums(Projection[BEMDEM$no_yrs, , , tx - 1]) > 0)
                BEMDEM$Niche_ID[colSums(Projection[BEMDEM$no_yrs, , , tx -
                                                     1]) > 0, 2]

              } else  {
                # take population from previous year, same time period
                n0s <-
                  t(Projection[yx - 1, , colSums(Projection[yx - 1, , , tx]) > 0, tx])
                n0s_ID <-
                  which(colSums(Projection[yx - 1, , , tx]) > 0)

              }
            }
            # n <- c(10,5,5,5,5,5)

            population_Niche_short <-
              population_Niche[n0s_ID]

            # run fcn for each population px=1
            if (nrow(n0s) > 0) {
              for (px in 1:nrow(n0s))  {
                # px = 1 tx = 1
                n <- as.vector(n0s[px, ])

                #  source('demoniche_population.r')
                # if (sum(n) > 0) {

                # selects the right K for this pop and timeperiod
                populationmax <-
                  BEMDEM$populationmax_all[n0s_ID[px], tx]

                Projection[yx, , n0s_ID[px], tx]    <-
                  demoniche_population(
                    Matrix_projection = Matrix_projection,
                    Matrix_projection_var = Matrix_projection_var,
                    n = n,
                    populationmax = populationmax,
                    onepopulation_Niche = population_Niche_short[px],
                    sumweight = BEMDEM$sumweight,
                    Kweight = BEMDEM$Kweight,
                    prob_scenario = BEMDEM$prob_scenario,
                    noise = BEMDEM$noise,
                    prev_mx = prev_mx,
                    transition_affected_demogr = BEMDEM$transition_affected_demogr,
                    transition_affected_niche = BEMDEM$transition_affected_niche,
                    transition_affected_env = BEMDEM$transition_affected_env,
                    env_stochas_type = BEMDEM$env_stochas_type,
                    yx_tx = yx_tx
                  )

                #     } # end if n > 0
              } # end if n0s > 0
            } # end px for

            #  Which many patches persist since last timestep (including seeds)?
            metapop_results[yx_tx, mx, rx] <-
              length(intersect(which(colSums(
                Projection[yx, , , tx]
              ) > 1), n0s_ID))
            ################################################################################

            ##################### DISPERSAL FOR ALL PATCHES ####################
            # Calculate dispersal for one repetition (rx), one matrix (mx), one time period(tx), and one year
            # and all populationes.
            # load('/noCC_nodispersal/Projection1_meanmatrix.rda')
            # print(paste(Projection[yx,1,,tx]))
            if (sum(Projection[yx, 1, , tx]) > 0)
              # break
            {
              if (Dispersal == TRUE)
              {
                if (Niche == TRUE) {
                  population_Niche <- BEMDEM$Niche_values[, tx]
                }

                # source("demoniche_dispersal.r")
                disp <-
                  demoniche_dispersal(
                    seeds_per_population = Projection[yx, 1, , tx],
                    fraction_LDD = BEMDEM$fraction_LDD,
                    dispersal_probabilities = BEMDEM$dispersal_probabilities,
                    dist_latlong = BEMDEM$dist_latlong,
                    neigh_index = BEMDEM$neigh_index,
                    fraction_SDD = BEMDEM$fraction_SDD
                  )

                Projection[yx, 1, , tx] <- disp

              } # end Dispersal if
            } # end sum if
            ################################################################################

            # save the results from one run
            population_sizes[yx_tx, mx, rx] <-
              sum(rowSums(Projection[yx, , , tx]) * BEMDEM$sumweight)

          } # end yx

          EMA[rx, mx, tx, 1] <- # the minimum population each repetition
            min(apply((Projection[, , , tx] * BEMDEM$sumweight), 1, sum))

          EMA[rx, mx, tx, 2] <-
            # the number of populations that exist each repetition
            sum(colSums(Projection[yx, , , tx]) > 1)

          simulation_results[mx, 7 + tx] <- mean(EMA[, mx, tx, 1])

        } # end tx loop

        save(
          Projection,
          file = paste(
            getwd(),
            "/",
            foldername,
            "/Projection_rep",
            rx,
            "_",
            BEMDEM$list_names_matrices[mx],
            ".rda",
            sep = ""
          )
        )

        save(
          population_sizes,
          file = paste(
            getwd(),
            "/",
            foldername,
            "/population_sizes",
            ".rda",
            sep = ""
          )
        )

        # PLOTS of spatial occupancy from the last repetition  yx = 1
        pop <- data.frame(cbind(BEMDEM$Niche_ID[, 2:3],
                                (
                                  colSums(Projection[yx, , , ] * BEMDEM$sumweight)
                                )))

        form <-
          as.formula(paste(paste(colnames(pop)[-c(1:2)], collapse = "+"), "X+Y", sep =
                             "~"))

        jpeg(
          filename = paste(
            getwd(),
            "/",
            foldername,
            "/map_",
            BEMDEM$list_names_matrices[mx],
            ".jpeg",
            sep = ""
          )
        )
        print(lattice::levelplot(
          form,
          pop,
          col.regions = rev(heat.colors(100)),
          allow.multiple = TRUE,
          main = paste(
            "Distribution",
            foldername,
            BEMDEM$list_names_matrices[mx],
            sep = "_"
          )
        ))
        dev.off()


      } # end mx loop

    } # end rx loop

    rm(Projection)

    #### END OF MODELLING #################################################
    ####################################################################################

    # extra loop to calculate mean, eigenvalues, etc
    print("Calculating summary values", quote = FALSE)

    for (mx in 1:length(BEMDEM$list_names_matrices))
      # mx = 1
    {
      # for 'eigen_results' object for each matrix.
      eigen_results[[mx]] <-
        c(popbio::eigen.analysis(matrix(
          BEMDEM$matrices[, mx],
          ncol = length(BEMDEM$stages),
          byrow = FALSE
        )),
        LTRE = list(matrix(
          popbio::LTRE(
            matrix(
              BEMDEM$matrices[, mx],
              ncol = length(BEMDEM$stages),
              byrow = FALSE
            ),
            matrix(
              BEMDEM$matrices[, 1],
              ncol = length(BEMDEM$stages),
              byrow = FALSE
            )
          ),
          nrow = length(BEMDEM$stages)
        )))

      simulation_results[mx, "lambda"] <-
        eigen_results[[mx]]$lambda1

      simulation_results[mx, "stoch_lambda"] <-
        popbio::stoch.growth.rate(
          list(
            matrix(BEMDEM$matrices[, 1],
                   ncol = length(BEMDEM$stages)),
            matrix(BEMDEM$matrices[, mx],
                   ncol = length(BEMDEM$stages))
          ),
          prob = NULL,
          maxt = 50000,
          verbose = FALSE
        )$sim

      simulation_results[mx, "mean_final_pop"] <-
        mean(population_sizes[yrs_total, mx, ])

      # Final mean percentage of patches extinct
      #   simulation_results[mx, "mean_perc_ext_final"] <- NA
      #                      (metapop_results[1,mx,] - mean(metapop_results[yrs_total,mx,]))/metapop_results[1,mx,]*100

      simulation_results[mx, "mean_no_patches_final"] <-
        mean(EMA[, mx, length(BEMDEM$years_projections), 2])


      for (yx_tx in 1:yrs_total)
      {
        # mean of all repetitions
        population_results[yx_tx, "Meanpop", mx] <-
          mean(population_sizes[yx_tx, mx, ])
        # sd of all repetitions
        population_results[yx_tx, "SD", mx] <-
          sd(population_sizes[yx_tx, mx, ])
        # Min, minimum abundance of all simualtions each year.
        population_results[yx_tx, "Min", mx] <-
          min(population_sizes[yx_tx, mx, ])
        # Max, maximum abundance of all simualtions each year.
        population_results[yx_tx, "Max", mx] <-
          max(population_sizes[yx_tx, mx, ])
      }

      #  jpeg(file = paste(getwd(), "/", foldername, "/mean_",
      #       BEMDEM$list_names_matrices[mx], ".jpeg", sep = ""))
      #   plot(population_results[,"Meanpop",mx], ylim = c(0, max(population_results[,"Meanpop",mx])),
      #       main = paste("meanpopulation", BEMDEM$list_names_matrices[mx], foldername, sep = "_"),type = "l",
      #        xlab = "Year", ylab = "Population")
      #   dev.off()

    } # end mx loop number 2, for eigenvalues.

    ######## PLOTS ################
    # plot EMAS
    jpeg(
      filename = paste(getwd(), "/", foldername, "/EMAs.jpeg", sep = ""),
      width = 580,
      height = 480
    )
    #  matplot(t(simulation_results[, 8:(7+length(BEMDEM$years_projections))]), pch = 15, type = "l")
    for (mx in 1:length(BEMDEM$list_names_matrices))
    {
      par(mar = c(7, 7, 4, 2) + 0.1, cex = 1.5)
      plot(
        simulation_results[mx, 8:(7 + length(BEMDEM$years_projections))],
        type = 'b',
        ylim =
          range(simulation_results[, 8:c(7 + length(BEMDEM$years_projections))])
        ,
        col = mx,
        xlab = "",
        ylab = "",
        axes = FALSE
      )
      par(new = TRUE)
    } # end mx loop number 3, for plotting EMAs
    axis(1,
         at = 1:length(BEMDEM$years_projections),
         labels = FALSE)

    text(
      1:length(BEMDEM$years_projections),
      y = par("usr")[3] - 5,
      srt = 45,
      adj = 1,
      labels = (BEMDEM$years_projections),
      xpd = TRUE
    )
    mtext(1,
          text = "Time",
          line = 6,
          cex = 1.7)
    mtext(2,
          text = "EMA (number of individuals)",
          line = 6,
          cex = 1.7)
    axis(2, xpd = TRUE, las = 1)

    legend(
      "topright",
      legend = BEMDEM$list_names_matrices,
      col = 1:mx,
      fill = 1:mx
    )
    title("EMA for different matrices", cex = 0.9)

    dev.off()

    #### plot population results
    jpeg(
      filename = paste(getwd(), "/", foldername, "/population_results.jpeg", sep = ""),
      width = 580,
      height = 480
    )
    for (mx in 1:length(BEMDEM$list_names_matrices))
    {
      par(mar = c(7, 4, 4, 2) + 0.1)
      plot(
        population_results[, "Meanpop", mx],
        type = 'l',
        ylim = range(population_results[, 1, ]),
        col = mx,
        xlab = "",
        ylab = "EMA",
        axes = TRUE,
        lwd = 1.2
      )
      par(new = TRUE)
      plot(
        c(population_results[, "Meanpop", mx] + population_results[, "SD", mx]),
        type = 'l',
        lty = 2,
        ylim = range(population_results[, 1, ]),
        col = mx,
        xlab = "",
        ylab = "EMA",
        axes = FALSE
      )
      par(new = TRUE)
      plot(
        c(population_results[, "Meanpop", mx] - population_results[, "SD", mx]),
        type = 'l',
        lty = 2,
        ylim = range(population_results[, 1, ]),
        col = mx,
        xlab = "",
        ylab = "EMA",
        axes = FALSE
      )
      par(new = TRUE)
    } # end mx loop number 3, for plotting EMAs
    legend(
      "topright",
      legend = BEMDEM$list_names_matrices,
      col = 1:mx,
      fill = 1:mx
    )
    title("Population sizes (+- 1 SD) for different scenarios")
    dev.off()

    #########################################################
    write.table(
      simulation_results,
      sep = ";",
      file =
        paste(getwd(), "/", foldername, "/population_results.csv", sep = ""),
      col.names = NA
    )

    save(
      simulation_results,
      file =
        paste(getwd(), "/", foldername, "/population_sizes.rda", sep = "")
    )

    save(metapop_results,
         file =
           paste(getwd(), "/", foldername, "/metapop_results.rda", sep = ""))

    save(eigen_results,
         file =
           paste(getwd(), "/", foldername, "/eigen_results.rda", sep = ""))

    save(EMA, file =
           paste(getwd(), "/", foldername, "/EMA.rda", sep = ""))

    save(
      population_results,
      file = paste(getwd(), "/",
                   foldername, "/population_results.rda", sep = "")
    )

    print("    All repetitions completed!", quote = FALSE)

    return(population_results)

  }     # end fcn
