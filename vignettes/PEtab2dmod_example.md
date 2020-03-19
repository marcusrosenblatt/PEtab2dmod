## PEtab

[PEtab](https://github.com/PEtab-dev/PEtab/) is a format thought for the straight forward exchange of models 
between tools with regard to parameter estimation. It covers SBML models along with their data, observables, 
parameters and condition specifications. dMod supports the import of PEtab files and provides some basic functions 
for visualization and further testing. 

## Installation

Usage of the PEtab import requires libSBML. The libSBML R source package can be downloaded as tar.gz file provided 
by [SBML](https://sourceforge.net/projects/sbml/files/libsbml/5.18.0/stable/R%20interface/) and installed as R package. 

```r

install.packages("~/Downloads/libSBML_5.18.0.tar.gz", repos = NULL, type = "source")

```

## Model import

In the following we illustrate the import and usage of PEtab models with a model example from the 
[benchmark collection](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab). We use the STAT5 activation model 
[Boehm_JProteomeRes2014](https://pubs.acs.org/doi/abs/10.1021/pr5006923) that is also provided with dMod in /BenchmarkModels/?.

A PEtab model consits of the subsequent files, exemplarily shown for Boehm_JProteomeRes2014 here, which have to meet the [requirements](https://github.com/PEtab-dev/PEtab/blob/master/doc/documentation_data_format.md) 
specified by PEtab.

```r

# SBML model file: model_Boehm_JProteomeRes2014.xml
# observable file: observables_Boehm_JProteomeRes2014.tsv
# condition file: experimentalCondition_Boehm_JProteomeRes2014.tsv
# data file: measurementData_Boehm_JProteomeRes2014.tsv
# parameter file: parameters_Boehm_JProteomeRes2014.tsv

```

As a first step, we load dMod, which includes all necessary functionalities, and import the model. 

```r

library(dMod)
importPEtabSBML(modelname = "Boehm_JProteomeRes2014", path2model = "../BenchmarkModels/") 

```

Now we can access model objects such as equations and observables loaded in the global environment.

```r

> reactions
                                    Conserved quantities: 2
1                1*STAT5A+2*pApA+1*pApB+2*nucpApA+1*nucpApB
2 -0.5*STAT5A-1*pApA+0.5*STAT5B+1*pBpB+-1*nucpApA+1*nucpBpB

  Check           Educt ->         Product                                                             Rate Description
1              2*STAT5A ->            pApA   (1.25e-7*exp(-1*Epo_degradation_BaF3*time))*(STAT5A**2)*k_phos            
2       STAT5A + STAT5B ->            pApB (1.25e-7*exp(-1*Epo_degradation_BaF3*time))*STAT5A*STAT5B*k_phos            
3              2*STAT5B ->            pBpB   (1.25e-7*exp(-1*Epo_degradation_BaF3*time))*(STAT5B**2)*k_phos            
4                  pApA ->         nucpApA                                                  k_imp_homo*pApA            
5                  pApB ->         nucpApB                                                k_imp_hetero*pApB            
6                  pBpB ->         nucpBpB                                                  k_imp_homo*pBpB            
7               nucpApA ->        2*STAT5A                                               k_exp_homo*nucpApA            
8               nucpApB -> STAT5A + STAT5B                                             k_exp_hetero*nucpApB            
9               nucpBpB ->        2*STAT5B                                               k_exp_homo*nucpBpB 

```

```r

> observables
Idx       Inner <- Outer
  1 pSTAT5A_rel <- (100*pApB+200*pApA*specC17)/(pApB+STAT5A*specC17+2*pApA*specC17)
  2 pSTAT5B_rel <- -(100*pApB-200*pBpB*(specC17-1))/((STAT5B*(specC17-1)-pApB)+2*pBpB*(specC17-1))
  3 rSTAT5A_rel <- (100*pApB+100*STAT5A*specC17+200*pApA*specC17)/(2*pApB+STAT5A*specC17+2*pApA*specC17-
                   STAT5B*(specC17-1)-2*pBpB*(specC17-1))

```

In addition, observation and prediction functions as well as the **objective function** are generated during the import and provided for further use.

## Parameter optimization

Published parameter values are provided in *pouter*. However, you might want to fit your model again in dMod. The default fitting function performs a **multistart optimization** with the number of fits defined in *nrfits*.

```r

myframe <- fitModelPEtabSBML(nrfits = 20)
bestfit <- as.parvec(myframe, 1)

```  

## Model visualization

A default plotting function can be used to visualize data and fits as well as model states.

```r

plotPEtabSBML(pouter1 = bestfit, name%in%names(observables))

```
![](README_files/figure-html/prediction-1.png)<!-- -->

## Feature Tests

To verify the support of PEtab by different modeling tools, a [test suite](https://github.com/PEtab-dev/petab_test_suite) of minimal example models was designed. It can be found in /TestModels/?. The function **testPEtabSBML()** is used to perform the tests, checking for the agreement of likelihood, chi2 and model simulations.



