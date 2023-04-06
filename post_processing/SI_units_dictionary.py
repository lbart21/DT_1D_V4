"""
Function:
Author: Luke Bartholomew
Edits:
"""
SI_UNITS = {
    "T"                         : r'$K$', #K
    "T_G"                       : r'$K$',
    "T_L"                       : r'$K$',
    "T_g"                       : r'$K$',
    "T_l"                       : r'$K$',
    "T_1"                       : r'$K$',
    "T_2"                       : r'$K$',

    "T_ave"                     : r'$K$',
    "T_t"                       : r'$K$',

    "p"                         : r'$Pa$', #Pa
    "p_G"                       : r'$Pa$',
    "p_L"                       : r'$Pa$',
    "p_g"                       : r'$Pa$',
    "p_l"                       : r'$Pa$',
    "p_1"                       : r'$Pa$',
    "p_2"                       : r'$Pa$',

    "p_t"                       : r'$Pa$',


    "vel_x"                     : r'$ms^{-1}$', #ms^-1
    "vel_x_G"                   : r'$ms^{-1}$',
    "vel_x_L"                   : r'$ms^{-1}$',
    "vel_x_g"                   : r'$ms^{-1}$',
    "vel_x_l"                   : r'$ms^{-1}$',
    "vel_x_1"                   : r'$ms^{-1}$',
    "vel_x_2"                   : r'$ms^{-1}$',

    "vel_x_ave"                 : r'$ms^{-1}$',


    "vel_y"                     : r'$ms^{-1}$', #ms^-1
    "vel_y_G"                   : r'$ms^{-1}$',
    "vel_y_L"                   : r'$ms^{-1}$',
    "vel_y_g"                   : r'$ms^{-1}$',
    "vel_y_l"                   : r'$ms^{-1}$',
    "vel_y_1"                   : r'$ms^{-1}$',
    "vel_y_2"                   : r'$ms^{-1}$',
    
    "vel_y_ave"                 : r'$ms^{-1}$',

    "vel_z"                     : r'$ms^{-1}$', #ms^-1
    "vel_z_G"                   : r'$ms^{-1}$',
    "vel_z_L"                   : r'$ms^{-1}$',
    "vel_z_g"                   : r'$ms^{-1}$',
    "vel_z_l"                   : r'$ms^{-1}$',
    "vel_z_1"                   : r'$ms^{-1}$',
    "vel_z_2"                   : r'$ms^{-1}$',

    "vel_z_ave"                 : r'$ms^{-1}$',

    "vel_mag"                   : r'$ms^{-1}$', #ms^-1
    
    "u"                         : r'$Jkg^{-1}$', #Jkg^-1
    "u_G"                       : r'$Jkg^{-1}$',
    "u_L"                       : r'$Jkg^{-1}$',
    "u_g"                       : r'$Jkg^{-1}$',
    "u_l"                       : r'$Jkg^{-1}$',
    "u_1"                       : r'$Jkg^{-1}$',
    "u_2"                       : r'$Jkg^{-1}$',
    

    "h"                         : r'$Jkg^{-1}$', #Jkg^-1
    "h_G"                       : r'$Jkg^{-1}$',
    "h_L"                       : r'$Jkg^{-1}$',
    "h_g"                       : r'$Jkg^{-1}$',
    "h_l"                       : r'$Jkg^{-1}$',
    "h_1"                       : r'$Jkg^{-1}$',
    "h_2"                       : r'$Jkg^{-1}$',


    "s"                         : r'$J(kgK)^{-1}$', #J(kgK)^-1
    "s_G"                       : r'$J(kgK)^{-1}$',
    "s_L"                       : r'$J(kgK)^{-1}$',
    "s_g"                       : r'$J(kgK)^{-1}$',
    "s_l"                       : r'$J(kgK)^{-1}$',
    "s_1"                       : r'$J(kgK)^{-1}$',
    "s_2"                       : r'$J(kgK)^{-1}$',


    "rho"                       : r'$kgm^{-3}$', #kgm^-3
    "rho_G"                     : r'$kgm^{-3}$',
    "rho_L"                     : r'$kgm^{-3}$',
    "rho_g"                     : r'$kgm^{-3}$',
    "rho_l"                     : r'$kgm^{-3}$',
    "rho_1"                     : r'$kgm^{-3}$',
    "rho_2"                     : r'$kgm^{-3}$',

    "rho_ave"                   : r'$kgm^{-3}$',

    "k"                         : r'$W(mK)^{-1}$', #W(mK)^-1
    "k_G"                       : r'$W(mK)^{-1}$',
    "k_L"                       : r'$W(mK)^{-1}$',
    "k_g"                       : r'$W(mK)^{-1}$',
    "k_l"                       : r'$W(mK)^{-1}$',
    "k_1"                       : r'$W(mK)^{-1}$',
    "k_2"                       : r'$W(mK)^{-1}$',


    "dynamic viscosity"         : r'$Pas$', #Pas
    "dynamic viscosity_G"       : r'$Pas$',
    "dynamic viscosity_L"       : r'$Pas$',
    "dynamic viscosity_g"       : r'$Pas$',
    "dynamic viscosity_l"       : r'$Pas$',
    "dynamic viscosity_1"       : r'$Pas$',
    "dynamic viscosity_2"       : r'$Pas$',


    "kinematic viscosity"       : r'$m^2s^{-1}$', #m^2s^-1
    "kinematic viscosity_G"     : r'$m^2s^{-1}$',
    "kinematic viscosity_L"     : r'$m^2s^{-1}$',
    "kinematic viscosity_g"     : r'$m^2s^{-1}$',
    "kinematic viscosity_l"     : r'$m^2s^{-1}$',
    "kinematic viscosity_1"     : r'$m^2s^{-1}$',
    "kinematic viscosity_2"     : r'$m^2s^{-1}$',


    "htc"                       : r'$Wm^{-2}K^{-1}$', #Wm^-2K^-1
    "htc_G"                     : r'$Wm^{-2}K^{-1}$',
    "htc_L"                     : r'$Wm^{-2}K^{-1}$',
    "htc_g"                     : r'$Wm^{-2}K^{-1}$',
    "htc_l"                     : r'$Wm^{-2}K^{-1}$',
    "htc_1"                     : r'$Wm^{-2}K^{-1}$',
    "htc_2"                     : r'$Wm^{-2}K^{-1}$',


    "a"                         : r'$ms^{-1}$', #ms^-1
    "a_G"                       : r'$ms^{-1}$',
    "a_L"                       : r'$ms^{-1}$',
    "a_g"                       : r'$ms^{-1}$',
    "a_l"                       : r'$ms^{-1}$',
    "a_1"                       : r'$ms^{-1}$',
    "a_2"                       : r'$ms^{-1}$',


    "Ma"                        : r'$-$', #-
    "Ma_G"                      : r'$-$',
    "Ma_L"                      : r'$-$',
    "Ma_g"                      : r'$-$',
    "Ma_l"                      : r'$-$',
    "Ma_1"                      : r'$-$',
    "Ma_2"                      : r'$-$',


    "psi"                       : r'$-$', #-
    "phi"                       : r'$m$', #m
 
    "alpha_G"                   : r'$-$', #-
    "alpha_L"                   : r'$-$',
    "alpha_g"                   : r'$-$',
    "alpha_l"                   : r'$-$',
    "alpha_1"                   : r'$-$',
    "alpha_2"                   : r'$-$',

    "mass_flux"                 : r'$kgs^{-1}$',
    "energy_flux"               : r'$Js^{-1}$',

    "A"                         : r'$m^{2}$',

    "massf"                     : r'$-$',

    "molef"                     : r'$-$',

    "conc"                      : r'$moleL^{-1}$',

    "D"                         : r'$m$',
    "R"                         : r'$m$'
}