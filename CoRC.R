install.packages("remotes")
remotes::install_github("jpahle/CoRC")
CoRC::getCopasi()

# key   name  reaction                   rate_law                                flux number_flux
# <chr> <chr> <chr>                      <chr>                                  <dbl>       <dbl>
#   1 (v_0) v_0   2 * STAT5A -> pApA         FunctionDB.Functions[Function for v_0]     0           0
# 2 (v_1) v_1   STAT5A + STAT5B -> pApB    FunctionDB.Functions[Function for v_1]     0           0
# 3 (v_2) v_2   2 * STAT5B -> pBpB         FunctionDB.Functions[Function for v_2]     0           0
# 4 (v_3) v_3   pApA -> nucpApA            FunctionDB.Functions[Function for v_3]     0           0
# 5 (v_4) v_4   pApB -> nucpApB            FunctionDB.Functions[Function for v_4]     0           0
# 6 (v_5) v_5   pBpB -> nucpBpB            FunctionDB.Functions[Function for v_5]     0           0
# 7 (v_6) v_6   nucpApA -> 2 * STAT5A      FunctionDB.Functions[Function for v_6]     0           0
# 8 (v_7) v_7   nucpApB -> STAT5A + STAT5B FunctionDB.Functions[Function for v_7]     0           0
# 9 (v_8) v_8   nucpBpB -> 2 * STAT5B      FunctionDB.Functions[Function for v_8]     0           0
# > getFunction("Function for v_0")