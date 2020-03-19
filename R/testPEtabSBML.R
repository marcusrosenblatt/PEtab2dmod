#' Test PEtabSBML import
#'
#' @description This function imports, evaluates and tests the PEtab model.
#'
#' @param models model names to test
#'
#' @return evaluation data.frame
#'
#' @author Marcus Rosenblatt and Svenja Kemmer
#'
#' @export
#'
testPEtabSBML <- function(models = c(
                            #"Boehm_JProteomeRes2014"
                            # "Fujita_SciSignal2010",
                            # "Borghans_BiophysChem1997",
                            # "Elowitz_Nature2000",
                            # "Sneyd_PNAS2002",
                            # "Crauste_CellSystems2017",
                            # "Schwen_PONE2014",
                            # "Raia_CancerResearch2011",
                            # "Zheng_PNAS2012",
                            # "Beer_MolBioSystems2014",
                            # "Brannmark_JBC2010",
                            # "Bruno_JExpBio2016",
                            # "Chen_MSB2009",
                            # "Fiedler_BMC2016",
                            # "Weber_BMC2015",
                            # "Swameye_PNAS2003"
                            # "Bachmann_MSB2011"
                            # "Lucarelli_CellSystems2018",
                             "0001",
                             "0002",
                             "0003",
                             "0004",
                             "0005",
                             "0006",
                             "0007",
                             "0008",
                             "0009",
                             "0010",
                             "0011",
                             "0012",
                             "0013",
                             "0014",
                             "0015",
                             "0016"
                          ), testFit = TRUE, timelimit = 5000, testCases = FALSE) {
  cat(green("Start test function...\n"))
  mywd <- getwd()
  teststarttime <- Sys.time()
  output <- NULL
  predictions <- NULL
  for (model in models) {
    setwd(mywd)
    importtest <- F
    plottest <- F
    bestfit <- NA
    cat(blue(paste0("Testing ", model, "\n")))
    fgh <- try_with_time_limit(
      {
        test <- try(importPEtabSBML(model, compile = T, testCases = testCases), silent = T)
        if (inherits(test, "try-error")) "import error" else test
      },
      timelimit
    )
    if (fgh == "import error") {
      cat(yellow("Import error or time limit exceeded for", model, "\n\n\n"))
      output <- rbind(output, data.frame(
        modelname = model, import = importtest,
        fitting_time = NA, plot = plottest, chi2 = NA, LL = NA, bestfit = NA, difference = NA
      ))
    } else {
      importtest <- T
      testobj <- try(obj(pouter))
      if (inherits(testobj, "try-error")) {
        cat(red("Warning: Error in calculation of objective function.\n"))
        output <- rbind(output, data.frame(
          modelname = model, import = importtest,
          fitting_time = NA, plot = plottest, chi2 = NA, LL = NA, bestfit = NA, difference = NA
        ))
      } else {
        if (is.numeric(testobj$value)) {
          cat(green("Calculation of objective function successful.\n"))
          if (testCases){
            # calculate predictions for trajectory comparison
            mysimulations <- read.csv(paste0("PEtabTests/", model, "/_simulations.tsv"), sep = "\t")
            simu_time <- unique(mysimulations$time)
            prediction <- (g*x*p0)(simu_time, pouter)
            predictions <- rbind(predictions, data.frame(
              modelname = model, pred = prediction, obs.transformation = NA
            ))
            # append observable scale to predictions
            for (i in 1:length(observables)) {
              scale <- attr(observables, "obsscales")[i]
              predictions <- predictions %>% mutate(obs.transformation = ifelse(modelname == model & pred.name == names(observables)[i], scale, obs.transformation))
            }
          }
        } else {
          cat(red("Warning: obj(pouter) is not numeric.\n"))
        }
        # objLL <- mynormL2(mydata, g * x * p0, outputLL = T)
        # testLL <- try(-0.5 * objLL(pouter)$value)
        # if (inherits(testLL, "try-error")) testLL <- NA
        if (testFit) {
          fitstarttime <- Sys.time()
          myframe <- fitModelPEtabSBML(nrfits = 20)
          fitendtime <- Sys.time()
          if (is.parframe(myframe) & !is.null(myframe)) {
            if (is.numeric(obj(myframe[1, ])$value)) {
              cat(green("Fit test successful.\n"))
              bestfit <- obj(myframe[1, ])$value
            } else {
              cat(red("Warning: obj(myframe) is not numeric.\n"))
            }
          } else {
            cat(red("Warning: Fit test not successful..\n"))
          }
          mytimediff <- as.numeric(difftime(fitendtime, fitstarttime, unit = "secs"))
          if (mytimediff > 3600) {
            cat(green(paste0("Fitting done in ", as.character(format(as.numeric(difftime(fitendtime, fitstarttime, unit = "hours")), digits = 3)), " hours.\n")))
          } else
          if (mytimediff > 60) {
            cat(green(paste0("Fitting done in ", as.character(format(as.numeric(difftime(fitendtime, fitstarttime, unit = "mins")), digits = 3)), " minutes.\n")))
          } else {
            cat(green(paste0("Fitting done in ", as.character(format(as.numeric(difftime(fitendtime, fitstarttime, unit = "secs")), digits = 3)), " seconds.\n")))
          }
        }
        pdf(file = paste0("Test/", model, "_plotAll.pdf"))
        plotPEtabSBML()
        dev.off()
        pdf(file = paste0("Test/", model, "_plotTargetsObserved.pdf"))
        plotPEtabSBML(name %in% names(observables))
        dev.off()
        pdf(file = paste0("Test/", model, "_plotConditionsObserved.pdf"))
        plotPEtabSBML(condition %in% names(mydata))
        dev.off()
        plottest <- T
        cat(green("Import and plot test for ", fgh, " successful!\n\n\n"))

        output <- rbind(output, data.frame(
          modelname = model, import = importtest,
          # fitting_time = format(as.numeric(difftime(fitendtime, fitstarttime, unit = "mins")), digits = 3),
          plot = plottest, chi2 = attr(testobj,"chisquare"), LL = -0.5*testobj$value
          # , bestfit = bestfit, difference = bestfit - testobj$value
        ))
      }
    }
  }
  if (testCases) {
    simu_output <- output[1]
    output <- cbind(output, chi2_sol = NA, tol_chi2_sol = NA, LL_sol = NA, tol_LL_sol = NA)
    for (model in models) {
      mysolution <- read_yaml(paste0("PEtabTests/", model, "/_", model, "_solution.yaml"))
      output[which(output$modelname == model), "chi2_sol"] <- mysolution$chi2
      output[which(output$modelname == model), "tol_chi2_sol"] <- mysolution$tol_chi2
      output[which(output$modelname == model), "LL_sol"] <- mysolution$llh
      output[which(output$modelname == model), "tol_LL_sol"] <- mysolution$tol_llh
      
      # extract simulation values
      simu_output[which(simu_output$modelname == model), "tol_simus_sol"] <- mysolution$tol_simulations
      mysimulations <- read.csv(paste0("PEtabTests/", model, "/_simulations.tsv"), sep = "\t")
      simu_prediction <- subset(predictions, modelname == model)
      
      # iterate through simulation points
      for (nrow in 1:nrow(mysimulations)) {
        simu_row <- mysimulations[nrow,]
        simu_time <- simu_row$time
        simu_obs <- simu_row$observableId %>% as.character()
        simu_condi <- simu_row$simulationConditionId %>% as.character()
        simu_obspars <- simu_row$observableParameters
        
        if(!is.null(simu_obspars) & length(unique(simu_prediction$pred.condition)) > 1){
          simu_condi <- paste0(simu_condi, "_", simu_obspars)
        }
        pred_row <- subset(simu_prediction, pred.time == simu_time & pred.name == simu_obs & pred.condition == simu_condi)
        # retransform simulation value according to observable transformation
        if(pred_row$obs.transformation == "log10") pred_row$pred.value <- 10**pred_row$pred.value
        if(pred_row$obs.transformation == "log") pred_row$pred.value <- exp(pred_row$pred.value)
        simu_value <- pred_row$pred.value

        simu_output[which(simu_output$modelname == model), paste0("simu_", nrow)] <- simu_value
        simu_output[which(simu_output$modelname == model), paste0("simu_", nrow,"_sol")] <- mysimulations$simulation[nrow]
      }
    }
  }

  testendtime <- Sys.time()
  mytimediff <- as.numeric(difftime(testendtime, teststarttime, unit = "secs"))
  if (mytimediff > 3600) {
    cat(green(paste0("Test done in ", as.character(format(as.numeric(difftime(testendtime, teststarttime, unit = "hours")), digits = 3)), " hours.\n")))
  } else
  if (mytimediff > 60) {
    cat(green(paste0("Test done in ", as.character(format(as.numeric(difftime(testendtime, teststarttime, unit = "mins")), digits = 3)), " minutes.\n")))
  } else {
    cat(green(paste0("Test done in ", as.character(format(as.numeric(difftime(testendtime, teststarttime, unit = "secs")), digits = 3)), " seconds.\n")))
  }

  if (testCases) {
    
    # check simulations
    for (model in models) {
      correctORnot <- NULL
      modelrow <- subset(simu_output, modelname == model)
      modelrow_woNA <- modelrow[colSums(!is.na(modelrow)) > 0]
      for (ncol in seq(3,(ncol(modelrow_woNA)),2)) {
        simuCompare <- abs(modelrow_woNA[[ncol]]-modelrow_woNA[[ncol+1]]) < modelrow_woNA$tol_simus_sol
        correctORnot <- c(correctORnot, simuCompare)
      }
      if(length(unique(correctORnot)) == 1){
        SimuPassed <- unique(correctORnot)
      } else SimuPassed <- FALSE
      simu_output[which(simu_output$modelname == model), "Passed"] <- SimuPassed
    }

    output <- cbind(output,
                    X2Passed = (abs(output$chi2 - output$chi2_sol) < output$tol_chi2_sol),
                    LLPassed = (abs(output$LL - output$LL_sol) < output$tol_LL_sol)
    )
  }
  
  if (!testCases) simu_output <- NULL
  return(list(output = output,simu_output = simu_output))
}



try_with_time_limit <- function(expr, cpu = Inf, elapsed = Inf) {
  y <- try(
    {
      setTimeLimit(cpu, elapsed)
      expr
    },
    silent = TRUE
  )
  if (inherits(y, "try-error")) NULL else y
}
