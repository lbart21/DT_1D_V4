species = {[0]='H', [1]='H2', [2]='O', [3]='O2', [4]='H2O2', [5]='H2O', [6]='HO2', [7]='OH', [8]='Ar', [9]='He', [10]='N2', }
config = {
  tempLimits = {lower=300.000000, upper=50000.000000},
  odeStep = {method='rkf', errTol=1.000000e-03},
  tightTempCoupling = false,
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
  reacIdx = { 0, 3,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 2, 7,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[2] = {
  equation = "O + H2 <=> H + OH",
  type = "elementary",
  frc = {model='Arrhenius', A=5.080000000000e-02, n=2.670000, C=3.165236634795e+03, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 1, 2,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 0, 7,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[3] = {
  equation = "H2 + OH <=> H2O + H",
  type = "elementary",
  frc = {model='Arrhenius', A=2.160000000000e+02, n=1.510000, C=1.726035239642e+03, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 1, 7,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 0, 5,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[4] = {
  equation = "O + H2O <=> OH + OH",
  type = "elementary",
  frc = {model='Arrhenius', A=2.970000000000e+00, n=2.020000, C=6.743111431836e+03, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 2, 5,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 7,},
  prodCoeffs = { 2.000000e+00,},
}

reaction[5] = {
  equation = "H2 + M <=> H + H + M",
  type = "anonymous_collider",
  frc = {model='Arrhenius', A=4.580000000000e+13, n=-1.400000, C=5.252581875038e+04, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 1,},
  reacCoeffs = { 1.000000e+00,},
  prodIdx = { 0,},
  prodCoeffs = { 2.000000e+00,},
  efficiencies = {
    [0]=1.000000e+00,
    [1]=2.500000e+00,
    [2]=1.000000e+00,
    [3]=1.000000e+00,
    [4]=1.000000e+00,
    [5]=1.200000e+01,
    [6]=1.000000e+00,
    [7]=1.000000e+00,
    [10]=1.000000e+00,
  },
}

reaction[6] = {
  equation = "H2 + Ar <=> H + H + Ar",
  type = "elementary",
  frc = {model='Arrhenius', A=5.840000000000e+12, n=-1.100000, C=5.252581875038e+04, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 1, 8,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 0, 8,},
  prodCoeffs = { 2.000000e+00, 1.000000e+00,},
}

reaction[7] = {
  equation = "H2 + He <=> H + H + He",
  type = "elementary",
  frc = {model='Arrhenius', A=5.840000000000e+12, n=-1.100000, C=5.252581875038e+04, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 1, 9,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 0, 9,},
  prodCoeffs = { 2.000000e+00, 1.000000e+00,},
}

reaction[8] = {
  equation = "O + O + M <=> O2 + M",
  type = "anonymous_collider",
  frc = {model='Arrhenius', A=6.160000000000e+03, n=-0.500000, C=0.000000000000e+00, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 2,},
  reacCoeffs = { 2.000000e+00,},
  prodIdx = { 3,},
  prodCoeffs = { 1.000000e+00,},
  efficiencies = {
    [0]=1.000000e+00,
    [1]=2.500000e+00,
    [2]=1.000000e+00,
    [3]=1.000000e+00,
    [4]=1.000000e+00,
    [5]=1.200000e+01,
    [6]=1.000000e+00,
    [7]=1.000000e+00,
    [10]=1.000000e+00,
  },
}

reaction[9] = {
  equation = "O + O + Ar <=> O2 + Ar",
  type = "elementary",
  frc = {model='Arrhenius', A=1.890000000000e+01, n=0.000000, C=-9.007589151482e+02, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 2, 8,},
  reacCoeffs = { 2.000000e+00, 1.000000e+00,},
  prodIdx = { 3, 8,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[10] = {
  equation = "O + O + He <=> O2 + He",
  type = "elementary",
  frc = {model='Arrhenius', A=1.890000000000e+01, n=0.000000, C=-9.007589151482e+02, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 2, 9,},
  reacCoeffs = { 2.000000e+00, 1.000000e+00,},
  prodIdx = { 3, 9,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[11] = {
  equation = "O + H + M <=> OH + M",
  type = "anonymous_collider",
  frc = {model='Arrhenius', A=4.710000000000e+06, n=-1.000000, C=0.000000000000e+00, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 0, 2,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 7,},
  prodCoeffs = { 1.000000e+00,},
  efficiencies = {
    [0]=1.000000e+00,
    [1]=2.500000e+00,
    [2]=1.000000e+00,
    [3]=1.000000e+00,
    [4]=1.000000e+00,
    [5]=1.200000e+01,
    [6]=1.000000e+00,
    [7]=1.000000e+00,
    [8]=7.500000e-01,
    [9]=7.500000e-01,
    [10]=1.000000e+00,
  },
}

reaction[12] = {
  equation = "H + OH + M <=> H2O + M",
  type = "anonymous_collider",
  frc = {model='Arrhenius', A=3.800000000000e+10, n=-2.000000, C=0.000000000000e+00, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 0, 7,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 5,},
  prodCoeffs = { 1.000000e+00,},
  efficiencies = {
    [0]=1.000000e+00,
    [1]=1.000000e+00,
    [2]=1.000000e+00,
    [3]=1.000000e+00,
    [4]=1.000000e+00,
    [5]=1.200000e+01,
    [6]=1.000000e+00,
    [7]=1.000000e+00,
    [8]=3.800000e-01,
    [9]=3.800000e-01,
    [10]=1.000000e+00,
  },
}

reaction[13] = {
  equation = "H + O2 + M <=> HO2 + M",
  type = "anonymous_collider",
  frc = {model='Troe',
   kInf={A=1.480000000000e+00, n=0.600000, C=0.000000000000e+00, rctIndex=-1},
   k0={A=9.040000000000e+01, n=-1.500000, C=2.465764628059e+02, rctIndex=-1},
 T3=1.000000000000e-32,  T1=1.000000000000e+32,  a=5.000000000000e-01,  T2=1.000000000000e+32, 
},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 0, 3,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 6,},
  prodCoeffs = { 1.000000e+00,},
  efficiencies = {
    [0]=1.000000e+00,
    [1]=3.000000e+00,
    [2]=1.000000e+00,
    [3]=1.100000e+00,
    [4]=1.000000e+00,
    [5]=1.600000e+01,
    [6]=1.000000e+00,
    [7]=1.000000e+00,
    [8]=1.000000e+00,
    [9]=1.200000e+00,
    [10]=1.000000e+00,
  },
}

reaction[14] = {
  equation = "HO2 + H <=> H2 + O2",
  type = "elementary",
  frc = {model='Arrhenius', A=1.660000000000e+07, n=0.000000, C=4.126381622467e+02, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 0, 6,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 1, 3,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[15] = {
  equation = "HO2 + H <=> OH + OH",
  type = "elementary",
  frc = {model='Arrhenius', A=7.080000000000e+07, n=0.000000, C=1.509651813098e+02, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 0, 6,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 7,},
  prodCoeffs = { 2.000000e+00,},
}

reaction[16] = {
  equation = "HO2 + O <=> OH + O2",
  type = "elementary",
  frc = {model='Arrhenius', A=3.250000000000e+07, n=0.000000, C=0.000000000000e+00, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 2, 6,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 3, 7,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[17] = {
  equation = "HO2 + OH <=> H2O + O2",
  type = "elementary",
  frc = {model='Arrhenius', A=2.890000000000e+07, n=0.000000, C=-2.516086355163e+02, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 6, 7,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 3, 5,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[18] = {
  equation = "HO2 + HO2 <=> H2O2 + O2",
  type = "elementary",
  frc = {model='Arrhenius', A=4.200000000000e+08, n=0.000000, C=6.028542906970e+03, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 6,},
  reacCoeffs = { 2.000000e+00,},
  prodIdx = { 3, 4,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[19] = {
  equation = "HO2 + HO2 <=> H2O2 + O2",
  type = "elementary",
  frc = {model='Arrhenius', A=1.300000000000e+05, n=0.000000, C=-8.202441517830e+02, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 6,},
  reacCoeffs = { 2.000000e+00,},
  prodIdx = { 3, 4,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[20] = {
  equation = "H2O2 + M <=> OH + OH + M",
  type = "anonymous_collider",
  frc = {model='Troe',
   kInf={A=2.950000000000e+08, n=0.000000, C=2.435571591797e+04, rctIndex=-1},
   k0={A=1.200000000000e+05, n=0.000000, C=2.289638583198e+04, rctIndex=-1},
 T3=1.000000000000e-32,  T1=1.000000000000e+32,  a=5.000000000000e-01,  T2=1.000000000000e+32, 
},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 4,},
  reacCoeffs = { 1.000000e+00,},
  prodIdx = { 7,},
  prodCoeffs = { 2.000000e+00,},
  efficiencies = {
    [0]=1.000000e+00,
    [1]=2.500000e+00,
    [2]=1.000000e+00,
    [3]=1.000000e+00,
    [4]=1.000000e+00,
    [5]=1.200000e+01,
    [6]=1.000000e+00,
    [7]=1.000000e+00,
    [8]=6.400000e-01,
    [9]=6.400000e-01,
    [10]=1.000000e+00,
  },
}

reaction[21] = {
  equation = "H2O2 + H <=> H2O + OH",
  type = "elementary",
  frc = {model='Arrhenius', A=2.410000000000e+07, n=0.000000, C=1.997772565999e+03, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 0, 4,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 5, 7,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[22] = {
  equation = "H2O2 + H <=> H2 + HO2",
  type = "elementary",
  frc = {model='Arrhenius', A=4.820000000000e+07, n=0.000000, C=4.000577304709e+03, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 0, 4,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 1, 6,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[23] = {
  equation = "H2O2 + O <=> OH + HO2",
  type = "elementary",
  frc = {model='Arrhenius', A=9.550000000000e+00, n=2.000000, C=1.997772565999e+03, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 2, 4,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 6, 7,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[24] = {
  equation = "H2O2 + OH <=> H2O + HO2",
  type = "elementary",
  frc = {model='Arrhenius', A=1.000000000000e+06, n=0.000000, C=0.000000000000e+00, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 4, 7,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 5, 6,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[25] = {
  equation = "H2O2 + OH <=> H2O + HO2",
  type = "elementary",
  frc = {model='Arrhenius', A=5.800000000000e+08, n=0.000000, C=4.810757111071e+03, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 4, 7,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 5, 6,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

