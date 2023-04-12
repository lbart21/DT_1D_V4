species = {[0]='O', [1]='O2', [2]='N2', [3]='H', [4]='H2', [5]='H2O', [6]='HO2', [7]='OH', [8]='H2O2', }
config = {
  tempLimits = {lower=300.000000, upper=50000.000000},
  odeStep = {method='alpha-qss', eps1= 1.000000e-03, eps2= 5.000000e-04, delta= 1.000000e-10, maxIters=10},
  tightTempCoupling = true,
  maxSubcycles = 10000,
  maxAttempts = 4
}

reaction = {}
reaction[1] = {
  equation = "H2 + OH <=> H2O + H",
  type = "elementary",
  frc = {model='Arrhenius', A=2.100000000000e+07, n=0.000000, C=2.566683442375e+03, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 4, 7,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 3, 5,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
  label = "r1",
}

reaction[2] = {
  equation = "H + O2 <=> OH + O",
  type = "elementary",
  frc = {model='Arrhenius', A=1.250000000000e+08, n=0.000000, C=8.203321590337e+03, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 1, 3,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 0, 7,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
  label = "r2",
}

reaction[3] = {
  equation = "O + H2 <=> OH + H",
  type = "elementary",
  frc = {model='Arrhenius', A=2.950000000000e+07, n=1.000000, C=4.932058379467e+03, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 0, 4,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 3, 7,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
  label = "r3",
}

reaction[4] = {
  equation = "H + O2 + M <=> HO2 + M",
  type = "anonymous_collider",
  frc = {model='Arrhenius', A=1.590000000000e+03, n=0.000000, C=-5.032712632109e+02, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 1, 3,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 6,},
  prodCoeffs = { 1.000000e+00,},
  efficiencies = {
    [0]=1.000000e+00,
    [1]=2.000000e+00,
    [2]=2.000000e+00,
    [3]=1.000000e+00,
    [4]=5.000000e+00,
    [5]=3.250000e+01,
    [6]=1.000000e+00,
    [7]=1.000000e+00,
    [8]=1.000000e+00,
  },
  label = "r4",
}

reaction[5] = {
  equation = "H + H + M <=> H2 + M",
  type = "anonymous_collider",
  frc = {model='Arrhenius', A=1.000000000000e+06, n=-1.000000, C=0.000000000000e+00, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 3,},
  reacCoeffs = { 2.000000e+00,},
  prodIdx = { 4,},
  prodCoeffs = { 1.000000e+00,},
  efficiencies = {
    [0]=1.000000e+00,
    [1]=2.000000e+00,
    [2]=2.000000e+00,
    [3]=1.000000e+00,
    [4]=5.000000e+00,
    [5]=1.500000e+01,
    [6]=1.000000e+00,
    [7]=1.000000e+00,
    [8]=1.000000e+00,
  },
  label = "r5",
}

reaction[6] = {
  equation = "H2 + HO2 <=> H2O2 + H",
  type = "elementary",
  frc = {model='Arrhenius', A=9.600000000000e+06, n=0.000000, C=1.207851031706e+04, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 4, 6,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 3, 8,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
  label = "r6",
}

reaction[7] = {
  equation = "M + H2O2 <=> OH + OH + M",
  type = "anonymous_collider",
  frc = {model='Arrhenius', A=1.170000000000e+11, n=0.000000, C=2.289884247609e+04, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 8,},
  reacCoeffs = { 1.000000e+00,},
  prodIdx = { 7,},
  prodCoeffs = { 2.000000e+00,},
  efficiencies = {
    [0]=1.000000e+00,
    [1]=7.800000e-01,
    [2]=1.000000e+00,
    [3]=1.000000e+00,
    [4]=2.300000e+00,
    [5]=6.000000e+00,
    [6]=1.000000e+00,
    [7]=1.000000e+00,
    [8]=6.600000e+00,
  },
  label = "r7",
}

reaction[8] = {
  equation = "HO2 + H <=> OH + OH",
  type = "elementary",
  frc = {model='Arrhenius', A=7.000000000000e+07, n=0.000000, C=0.000000000000e+00, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 3, 6,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 7,},
  prodCoeffs = { 2.000000e+00,},
  label = "r8",
}

reaction[9] = {
  equation = "H + OH + M <=> H2O + M",
  type = "anonymous_collider",
  frc = {model='Arrhenius', A=7.500000000000e+11, n=-2.600000, C=0.000000000000e+00, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 3, 7,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 5,},
  prodCoeffs = { 1.000000e+00,},
  efficiencies = {
    [0]=1.000000e+00,
    [1]=1.600000e+00,
    [2]=1.600000e+00,
    [3]=1.000000e+00,
    [4]=4.000000e+00,
    [5]=2.000000e+01,
    [6]=1.000000e+00,
    [7]=1.000000e+00,
    [8]=1.000000e+00,
  },
  label = "r9",
}

reaction[10] = {
  equation = "O + O + M <=> O2 + M",
  type = "anonymous_collider",
  frc = {model='Arrhenius', A=1.380000000000e+06, n=-1.000000, C=1.711122294917e+02, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 0,},
  reacCoeffs = { 2.000000e+00,},
  prodIdx = { 1,},
  prodCoeffs = { 1.000000e+00,},
  efficiencies = {
    [0]=1.000000e+00,
    [1]=1.000000e+00,
    [2]=1.000000e+00,
    [3]=1.000000e+00,
    [4]=1.000000e+00,
    [5]=1.000000e+00,
    [6]=1.000000e+00,
    [7]=1.000000e+00,
    [8]=1.000000e+00,
  },
  label = "r10",
}

reaction[11] = {
  equation = "O + H2O <=> OH + OH",
  type = "elementary",
  frc = {model='Arrhenius', A=5.750000000000e+07, n=0.000000, C=9.058882737796e+03, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 0, 5,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 7,},
  prodCoeffs = { 2.000000e+00,},
  label = "r11",
}

reaction[12] = {
  equation = "H2 + O2 <=> OH + OH",
  type = "elementary",
  frc = {model='Arrhenius', A=1.000000000000e+07, n=0.000000, C=2.164066431807e+04, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 1, 4,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 7,},
  prodCoeffs = { 2.000000e+00,},
  label = "r12",
}

reaction[13] = {
  equation = "HO2 + OH <=> H2O + O2",
  type = "elementary",
  frc = {model='Arrhenius', A=6.300000000000e+06, n=0.000000, C=0.000000000000e+00, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 6, 7,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 1, 5,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
  label = "r13",
}

reaction[14] = {
  equation = "HO2 + O <=> O2 + OH",
  type = "elementary",
  frc = {model='Arrhenius', A=6.000000000000e+06, n=0.000000, C=0.000000000000e+00, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 0, 6,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 1, 7,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
  label = "r14",
}

reaction[15] = {
  equation = "HO2 + HO2 <=> H2O2 + O2",
  type = "elementary",
  frc = {model='Arrhenius', A=1.800000000000e+06, n=0.000000, C=0.000000000000e+00, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 6,},
  reacCoeffs = { 2.000000e+00,},
  prodIdx = { 1, 8,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
  label = "r15",
}

reaction[16] = {
  equation = "OH + H2O2 <=> H2O + HO2",
  type = "elementary",
  frc = {model='Arrhenius', A=1.000000000000e+07, n=0.000000, C=9.058882737796e+02, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 7, 8,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 5, 6,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
  label = "r16",
}

reaction[17] = {
  equation = "O + H2O2 <=> OH + HO2",
  type = "elementary",
  frc = {model='Arrhenius', A=8.000000000000e+07, n=0.000000, C=5.032712632109e+02, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 0, 8,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 6, 7,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
  label = "r17",
}

reaction[18] = {
  equation = "H + H2O2 <=> H2O + OH",
  type = "elementary",
  frc = {model='Arrhenius', A=3.180000000000e+08, n=0.000000, C=4.529441368898e+03, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 3, 8,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 5, 7,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
  label = "r18",
}

