-- H2-O2-6sp-8r.lua
--
-- This chemical kinetic system provides
-- a simple hydrogen-oxygen reaction mechanism.
--
-- Author: Luke E. Bartholomew
-- Date: 2 March, 2023
-- Place: The University of Queensland

conv_factor = 4184.0 / 8.3145

Reaction{
    'H2 + M <=> H + H + M',
    fr = {'Arrhenius', A = 4.58e19,  n = -1.4, C = 104.38 * conv_factor},
    efficiencies = {H2O = 12.0, H2 = 2.5}
}

Reaction{
    'H + OH + M <=> H2O + M',
    fr = {'Arrhenius', A = 3.8e22,  n = -2.0, C = 0.0 * conv_factor},
    efficiencies = {H2O = 12.0, H2 = 2.5}
}

Reaction{
    'O + O + M <=> O2 + M',
    fr = {'Arrhenius', A = 6.16e15,  n = -0.5, C = 0.0 * conv_factor},
    efficiencies = {H2O = 12.0, H2 = 2.5}
}

Reaction{
    'O + H + M <=> OH + M',
    fr = {'Arrhenius', A = 4.71e18,  n = -1.0, C = 0.0 * conv_factor},
    efficiencies = {H2O = 12.0, H2 = 2.5}
}

Reaction{
    'O2 + H <=> O + OH',
    fr = {'Arrhenius', A=3.55e15,  n = -0.41, C = 16.6 * conv_factor}
}

Reaction{
    'H2 + O <=> H + OH',
    fr = {'Arrhenius', A = 5.08e4,  n = 2.67, C = 6.29 * conv_factor}
}

Reaction{
    'H2 + OH <=> H2O + H',
    fr = {'Arrhenius', A = 2.16e8,  n = 1.51, C = 3.43 * conv_factor}
}

Reaction{
    'O + H2O <=> OH + OH',
    fr = {'Arrhenius', A = 2.97e6,  n = 2.02, C = 13.4 * conv_factor}
}