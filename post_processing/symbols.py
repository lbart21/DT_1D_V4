"""
Function:
Author: Luke Bartholomew
Edits:
"""
SYMBOLS = {
    "T"                     : r'$T$',
    "T_G"                   : r'$T_{G}$',
    "T_L"                   : r'$T_{L}$',
    "T_g"                   : r'$T_{g}$',
    "T_l"                   : r'$T_{l}$',
    "T_1"                   : r'$T_{1}$',
    "T_2"                   : r'$T_{2}$',

    "T_ave"                 : r'$T_{ave}$',
    "T_t"                   : r'$T_{t}$',

    "p"                     : r'$p$',
    "p_G"                   : r'$p_{G}$',
    "p_L"                   : r'$p_{L}$',
    "p_g"                   : r'$p_{g}$',
    "p_l"                   : r'$p_{l}$',
    "p_1"                   : r'$p_{1}$',
    "p_2"                   : r'$p_{2}$',

    "p_t"                   : r'$p_{t}$',

    "vel_x"                 : r'$v_{x}$',
    "vel_x_G"               : r'$v_{x(G)}$',
    "vel_x_L"               : r'$v_{x(L)}$',
    "vel_x_g"               : r'$v_{x(g)}$',
    "vel_x_l"               : r'$v_{x(l)}$',
    "vel_x_1"               : r'$v_{x(1)}$',
    "vel_x_2"               : r'$v_{x(1)}$',

    "vel_x_ave"             : r'$v_{x(ave)}$',

    "vel_y"                 : r'$v_{z}$',
    "vel_y_G"               : r'$v_{z(G)}$',
    "vel_y_L"               : r'$v_{z(L)}$',
    "vel_y_g"               : r'$v_{z(g)}$',
    "vel_y_l"               : r'$v_{z(l)}$',
    "vel_y_1"               : r'$v_{z(1)}$',
    "vel_y_2"               : r'$v_{z(2)}$',

    "vel_y_ave"             : r'$v_{y(ave)}$',

    "vel_z"                 : r'$v_{z}$',
    "vel_z_G"               : r'$v_{z(G)}$',
    "vel_z_L"               : r'$v_{z(L)}$',
    "vel_z_g"               : r'$v_{z(g)}$',
    "vel_z_l"               : r'$v_{z(l)}$',
    "vel_z_1"               : r'$v_{z(1)}$',
    "vel_z_2"               : r'$v_{z(2)}$',

    "vel_z_ave"             : r'$v_{z(ave)}$',

    "vel_mag"               : r'$v$',

    "u"                     : r'$u$',
    "u_G"                   : r'$u_{G}$',
    "u_L"                   : r'$u_{L}$',
    "u_g"                   : r'$u_{g}$',
    "u_l"                   : r'$u_{l}$',
    "u_1"                   : r'$u_{1}$',
    "u_2"                   : r'$u_{2}$',
    
    "h"                     : r'$u$',
    "h_G"                   : r'$h_{G}$',
    "h_L"                   : r'$h_{L}$',
    "h_g"                   : r'$h_{g}$',
    "h_l"                   : r'$h_{l}$',
    "h_1"                   : r'$h_{1}$',
    "h_2"                   : r'$h_{2}$',

    "s"                     : r'$s$',
    "s_G"                   : r'$s_{G}$',
    "s_L"                   : r'$s_{L}$',
    "s_g"                   : r'$s_{g}$',
    "s_l"                   : r'$s_{l}$',
    "s_1"                   : r'$s_{1}$',
    "s_2"                   : r'$s_{2}$',

    "rho"                   : r'$\rho$',
    "rho_G"                 : r'$\rho_{G}$',
    "rho_L"                 : r'$\rho_{L}$',
    "rho_g"                 : r'$\rho_{g}$',
    "rho_l"                 : r'$\rho_{l}$',
    "rho_1"                 : r'$\rho_{1}$',
    "rho_2"                 : r'$\rho_{2}$',

    "rho_ave"                 : r'$\rho_{ave}$',

    "k"                     : r'$k$',
    "k_G"                   : r'$k_{G}$',
    "k_L"                   : r'$k_{L}$',
    "k_g"                   : r'$k_{g}$',
    "k_l"                   : r'$k_{l}$',
    "k_1"                   : r'$k_{1}$',
    "k_2"                   : r'$k_{2}$',

    "dynamic viscosity"     : r'$\mu$',
    "dynamic viscosity_G"   : r'$\mu_{G}$',
    "dynamic viscosity_L"   : r'$\mu_{L}$',
    "dynamic viscosity_g"   : r'$\mu_{g}$',
    "dynamic viscosity_l"   : r'$\mu_{l}$',
    "dynamic viscosity_1"   : r'$\mu_{1}$',
    "dynamic viscosity_2"   : r'$\mu_{2}$',

    "kinematic viscosity"   : r'$\nu$',
    "kinematic viscosity_G" : r'$\nu_{G}$',
    "kinematic viscosity_L" : r'$\nu_{L}$',
    "kinematic viscosity_g" : r'$\nu_{g}$',
    "kinematic viscosity_l" : r'$\nu_{l}$',
    "kinematic viscosity_1" : r'$\nu_{1}$',
    "kinematic viscosity_2" : r'$\nu_{2}$',

    "htc"                   : r'$htc$',
    "htc_G"                 : r'$htc_{G}$',
    "htc_L"                 : r'$htc_{L}$',
    "htc_g"                 : r'$htc_{g}$',
    "htc_l"                 : r'$htc_{l}$',
    "htc_1"                 : r'$htc_{1}$',
    "htc_2"                 : r'$htc_{2}$',

    "a"                     : r'$a$',
    "a_G"                   : r'$a_{G}$',
    "a_L"                   : r'$a_{L}$',
    "a_g"                   : r'$a_{g}$',
    "a_l"                   : r'$a_{l$',
    "a_1"                   : r'$a_{1}$',
    "a_2"                   : r'$a_{2}$',

    "Ma"                    : r'$Ma$',
    "Ma_G"                  : r'$Ma_{G}$',
    "Ma_L"                  : r'$Ma_{L}$',
    "Ma_g"                  : r'$Ma_{g}$',
    "Ma_l"                  : r'$Ma_{l}$',
    "Ma_1"                  : r'$Ma_{1}$',
    "Ma_2"                  : r'$Ma_{2}$',

    "psi"                   : r'$\psi$',
    "phi"                   : r'$\phi$',

    "alpha_G"               : r'$\alpha_{G}$',
    "alpha_L"               : r'$\alpha_{L}$',
    "alpha_g"               : r'$\alpha_{g}$',
    "alpha_l"               : r'$\alpha_{l}$',
    "alpha_1"               : r'$\alpha_{1}$',
    "alpha_2"               : r'$\alpha_{2}$',

    "mass_flux"             : r'$\dot{m}$',
    "energy_flux"           : r'$\dot{m}h_{t}$',

    "A"                     : r'$A$',

    "massf"                 : r'$f$',

    "molef"                 : r'$\chi$',

    "conc"                  : r'$C$',

    "D"                     : r'$D$',
    "R"                     : r'$R$'
}
