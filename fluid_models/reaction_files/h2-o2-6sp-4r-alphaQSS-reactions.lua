-- chain reaction mechanism for H2 + O2 reactions
-- Author: Luke Bartholomew
-- Date: 31/3/23
-- Place: UQ, Saint Lucia
-- Reference: An Updated Comprehensive Model of Hydrogen Combustion, 2004
-- url: https://onlinelibrary.wiley.com/doi/epdf/10.1002/kin.20026
--
-- History:
--   
C_factor = 4184.0 / 8.3145

Config{
   odeStep = {method='alpha-qss', eps1=0.001}, --default eps1=0.001
   tightTempCoupling = true
}

Reaction{
    'H + O2 <=> O + OH',
    fr={'Arrhenius', A=3.55e15,  n=-0.41, C=16.6*C_factor}
 }

Reaction{
    'O + H2 <=> H + OH',
    fr={'Arrhenius', A=5.08e4,  n=2.67, C=6.29*C_factor}
}

Reaction{
    'H2 + OH <=> H2O + H',
    fr={'Arrhenius', A=2.16e8,  n=1.51, C=3.43*C_factor}
}

Reaction{
    'O + H2O <=> OH + OH',
    fr={'Arrhenius', A=2.97e6,  n=2.02, C=13.4*C_factor}
}
