-- chain, dissociation, recombination and HO2/H2O2 formation and consumption 
-- reaction mechanism for H2 + O2 reactions.
-- Full 19 reaction model from reference.
-- Author: Luke Bartholomew
-- Date: 31/3/23
-- Place: UQ, Saint Lucia
-- Reference: An Updated Comprehensive Model of Hydrogen Combustion, 2004
-- url: https://onlinelibrary.wiley.com/doi/epdf/10.1002/kin.20026
--
-- History:
--   

C_factor = 4184.0 / 8.3145

-- Chain reactions

Reaction{
    'H + O2 <=> O + OH',
    fr={'Arrhenius', A=3.55e15, n=-0.41, C=16.6*C_factor}
}
-- In the right order as https://pubs.acs.org/doi/pdf/10.1021/jp971416a

Reaction{
    'O + H2 <=> H + OH',
    fr={'Arrhenius', A=5.08e4, n=2.67, C=6.29*C_factor}
}
-- In the right order as https://reader.elsevier.com/reader/sd/pii/S0082078488803254?token=F7BDB3DC5F42958BEA966258322DC84F1A0A364EA7A74969470BDE3C6D21BB904935375816DFD5A9F54B33C869093A98&originRegion=us-east-1&originCreation=20230404053410

Reaction{
    'H2 + OH <=> H2O + H',
    fr={'Arrhenius', A=2.16e8, n=1.51, C=3.43*C_factor}
}
-- In the right order as https://pubs.acs.org/doi/pdf/10.1021/j100324a035

Reaction{
    'O + H2O <=> OH + OH',
    fr={'Arrhenius', A=2.97e6, n=2.02, C=13.4*C_factor}
}
-- In the right order as https://reader.elsevier.com/reader/sd/pii/S0082078406802419?token=BD57E8D46921B93153D527790E2A5FD94515C4565A67EC5C237C1EE1F772B819CDD004AFC18D4369BAFB3ABDCED60C94&originRegion=us-east-1&originCreation=20230404060452

-- H2/O2 dissociation/recombination reactions

Reaction{
    'H2 + M <=> H + H + M',
    fr={'Arrhenius', A=4.58e19, n=-1.4, C=104.38*C_factor},
    efficiencies = {Ar=0.0, He=0.0, H2O=12.0, H2=2.5}
}
-- In the right order as https://aip.scitation.org/doi/pdf/10.1063/1.555759
-- Might need to change efficiencies so that only N2 has a non-zero efficiency. 
-- Paper has this reaction for M = N2

Reaction{
    'H2 + Ar <=> H + H + Ar',
    fr={'Arrhenius', A=5.84e18, n=-1.1, C=104.38*C_factor}
}
-- In the right order as https://aip.scitation.org/doi/pdf/10.1063/1.555759

Reaction{
    'H2 + He <=> H + H + He',
    fr={'Arrhenius', A=5.84e18, n=-1.1, C=104.38*C_factor}
}
-- In the right order as https://aip.scitation.org/doi/pdf/10.1063/1.555759

Reaction{
    'O + O + M <=> O2 + M',
    fr={'Arrhenius', A=6.16e15, n=-0.5, C=0.0*C_factor},
    efficiencies = {Ar=0.0, He=0.0, H2O=12.0, H2=2.5}
}
-- In the right order as https://aip.scitation.org/doi/pdf/10.1063/1.555759
-- Model is chosen for 2000-10000K, but 200-4000K model is available for M = Ar

Reaction{
    'O + O + Ar <=> O2 + Ar',
    fr={'Arrhenius', A=1.89e13, n=0.0, C=-1.79*C_factor}
}
-- In the right order as https://aip.scitation.org/doi/pdf/10.1063/1.555759

Reaction{
    'O + O + He <=> O2 + He',
    fr={'Arrhenius', A=1.89e13, n=0.0, C=-1.79*C_factor}
}
-- In the right order as https://aip.scitation.org/doi/pdf/10.1063/1.555759

Reaction{
    'O + H + M <=> OH + M',
    fr={'Arrhenius', A=4.71e18, n=-1.0, C=0.0*C_factor},
    efficiencies = {H2O=12.0, H2=2.5, Ar=0.75, He=0.75}
}
-- In the right order as https://aip.scitation.org/doi/pdf/10.1063/1.555759

Reaction{
    'H + OH + M <=> H2O + M',
    fr={'Arrhenius', A=3.8e22, n=-2.0, C=0.0*C_factor},
    efficiencies = {H2O=12.0, Ar=2.5, Ar=0.38, He=0.38}
}
-- In the right order as https://aip.scitation.org/doi/pdf/10.1063/1.555759

-- HO2 formation and consumption reactions

Reaction{
    'H + O2 + M <=> HO2 + M',
    fr={'pressure dependent', 
        kInf={A=1.48e12, n=0.6, C=0.0*C_factor},
        k0={A=6.37e20, n=-1.72, C=0.52*C_factor},
        Troe={a=0.8, T1=1e32, T2=1e32, T3=1e-32}},
    efficiencies = {H2O=11.0, H2=2.0, O2=0.78}
}
-- When main bath is N2
-- In the right order as https://pubs.acs.org/doi/pdf/10.1021/j100248a033 for  


--Reaction{
    --'H + O2 + M <=> HO2 + M',
    --fr={'pressure dependent', 
        --kInf={A=1.48e12, n=0.6, C=0.0*C_factor},
        --k0={A=9.04e19, n=-1.5, C=0.49*C_factor},
        --Troe={a=0.5, T1=1e32, T2=1e32, T3=1e-32}},
    --efficiencies = {H2O=16.0, H2=3.0, O2=1.1, He=1.2}
--}
-- When main bath is Ar or He

Reaction{
    'HO2 + H <=> H2 + O2',
    fr={'Arrhenius', A=1.66e13, n=0.0, C=0.82*C_factor}
}

Reaction{
    'HO2 + H <=> OH + OH',
    fr={'Arrhenius', A=7.08e13, n=0.0, C=0.3*C_factor}
}

Reaction{
    'HO2 + O <=> OH + O2',
    fr={'Arrhenius', A=3.25e13, n=0.0, C=0.0*C_factor}
}

Reaction{
    'HO2 + OH <=> H2O + O2',
    fr={'Arrhenius', A=2.89e13, n=0.0, C=-0.5*C_factor}
}

-- H2O2 formation and consumption reactions

Reaction{
    'HO2 + HO2 <=> H2O2 + O2',
    fr={'Arrhenius', A=4.2e14, n=0.0, C=11.98*C_factor}
}
Reaction{
    'HO2 + HO2 <=> H2O2 + O2',
    fr={'Arrhenius', A=1.3e11, n=0.0, C=-1.63*C_factor}
}

Reaction{
    'H2O2 + M <=> OH + OH + M',
    fr={'pressure dependent', 
        kInf={A=2.95e14, n=0.0, C=48.4*C_factor},
        k0={A=1.2e17, n=0.0, C=45.5*C_factor},
        Troe={a=0.5, T1=1e32, T2=1e32, T3=1e-32}},
    efficiencies = {H2O=12.0, H2=2.5, Ar=0.64, He=0.64}
}

Reaction{
    'H2O2 + H <=> H2O + OH',
    fr={'Arrhenius', A=2.41e13, n=0.0, C=3.97*C_factor}
}

Reaction{
    'H2O2 + H <=> H2 + HO2',
    fr={'Arrhenius', A=4.82e13, n=0.0, C=7.95*C_factor}
}

Reaction{
    'H2O2 + O <=> OH + HO2',
    fr={'Arrhenius', A=9.55e6, n=2.0, C=3.97*C_factor}
}

Reaction{
    'H2O2 + OH <=> H2O + HO2',
    fr={'Arrhenius', A=1.0e12, n=0.0, C=0.0*C_factor}
}

Reaction{
    'H2O2 + OH <=> H2O + HO2',
    fr={'Arrhenius', A=5.8e14, n=0.0, C=9.56*C_factor}
}