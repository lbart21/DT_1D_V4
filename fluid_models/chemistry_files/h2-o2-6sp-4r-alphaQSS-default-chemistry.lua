species = {[0]='H2', [1]='O2', [2]='OH', [3]='H2O', [4]='H', [5]='O', }
config = {
  tempLimits = {lower=300.000000, upper=50000.000000},
  odeStep = {method='alpha-qss', eps1= 1.000000e-03, eps2= 5.000000e-04, delta= 1.000000e-10, maxIters=10},
  tightTempCoupling = true,
  maxSubcycles = 10000,
  maxAttempts = 4
}

reaction = {}
reaction[1] = {
  equation = "H + O2 <=> O + OH",
  type = "elementary",
  frc = {model='Arrhenius', A=3.550000000000e+09, n=-0.410000, C=8.353406699140e+03, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 1, 4,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 2, 5,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[2] = {
  equation = "O + H2 <=> H + OH",
  type = "elementary",
  frc = {model='Arrhenius', A=5.080000000000e-02, n=2.670000, C=3.165236634795e+03, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 0, 5,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 2, 4,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[3] = {
  equation = "H2 + OH <=> H2O + H",
  type = "elementary",
  frc = {model='Arrhenius', A=2.160000000000e+02, n=1.510000, C=1.726035239642e+03, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 0, 2,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 3, 4,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[4] = {
  equation = "O + H2O <=> OH + OH",
  type = "elementary",
  frc = {model='Arrhenius', A=2.970000000000e+00, n=2.020000, C=6.743111431836e+03, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 3, 5,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 2,},
  prodCoeffs = { 2.000000e+00,},
}

