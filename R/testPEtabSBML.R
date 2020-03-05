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
                            # "Boehm_JProteomeRes2014",
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
                            "0008"
                            # "0009",
                            # "0010",
                            # "0011",
                            # "0012",
                            # "0013",
                            # "0014"
                          ), testFit = TRUE, timelimit = 5000, tests = FALSE) {
  cat(green("Start test function...\n"))
  mywd <- getwd()
  teststarttime <- Sys.time()
  output <- NULL
  for (model in models) {
    setwd(mywd)
    importtest <- F
    plottest <- F
    bestfit <- NA
    cat(blue(paste0("Testing ", model, "\n")))
    fgh <- try_with_time_limit(
      {
        test <- try(importPEtabSBML(model, compile = T, TestCases = tests), silent = T)
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
        } else {
          cat(red("Warning: obj(pouter) is not numeric.\n"))
        }
        objLL <- mynormL2(mydata, g * x * p0, outputLL = T)
        testLL <- try(-0.5 * objLL(pouter)$value)
        if (inherits(testLL, "try-error")) testLL <- NA
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
          fitting_time = format(as.numeric(difftime(fitendtime, fitstarttime, unit = "mins")), digits = 3),
          plot = plottest, chi2 = testobj$value, LL = testLL, bestfit = bestfit, difference = bestfit - testobj$value
        ))
      }
    }
  }
  if (tests) {
    output <- cbind(output, chi2_sol = NA, tol_chi2_sol = NA, LL_sol = NA, tol_LL_sol = NA)
    for (model in models) {
      mysolution <- read_yaml(paste0("PEtabTests/", model, "/_", model, "_solution.yaml"))
      output[which(output$modelname == model), "chi2_sol"] <- mysolution$chi2
      output[which(output$modelname == model), "tol_chi2_sol"] <- mysolution$tol_chi2
      output[which(output$modelname == model), "LL_sol"] <- mysolution$llh
      output[which(output$modelname == model), "tol_LL_sol"] <- mysolution$tol_llh
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

  if (tests) {
    output <- cbind(output,
      Passed = (abs(output$chi2 - output$chi2_sol) < output$tol_chi2_sol) &
        (abs(output$LL - output$LL_sol) < output$tol_LL_sol)
    )
  }
  return(output)
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