from fluidFluxes import AUSMPlusupORIGINAL
from gdtk.gas import GasModel, GasState
gm_L = GasModel("ideal-air-gas-model.lua")
gs_L = GasState(gm_L)
gs_L.p = 1e6
gs_L.T = 300.0
gs_L.update_thermo_from_pT()

print("a = ", gs_L.a, ", p = ", gs_L.p, ", rho = ", gs_L.rho, ", u = ", gs_L.u, ", vel_x_L = ", 2.0 * gs_L.a, ", vel_x_R = ", 0.0)
f = AUSMPlusupORIGINAL(a_L_forMa = gs_L.a, a_R_forMa = gs_L.a, \
                        p_L = gs_L.p, rho_L = gs_L.rho, \
                        u_L = gs_L.u, a_L = gs_L.a, vel_x_L = 2.0 * gs_L.a, \
                        p_R = gs_L.p/10.0, rho_R = gs_L.rho / 10.0, u_R = gs_L.u, \
                        a_R = gs_L.a, vel_x_R = 0.0)

print(f.fluxes)