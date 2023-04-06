species = {[0]='H', [1]='H2', [2]='O', [3]='O2', [4]='H2O2', [5]='H2O', [6]='HO2', [7]='OH', [8]='Ar', [9]='He', [10]='N2', }
config = {
  tempLimits = {lower=300.000000, upper=50000.000000},
  odeStep = {method='rkf', errTol=1.000000e-03},
  tightTempCoupling = false,
  maxSubcycles = 10000,
  maxAttempts = 4
}

reaction = {}
