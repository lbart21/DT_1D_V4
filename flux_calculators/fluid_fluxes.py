"""
Function:
Author: Luke Bartholomew
Edits:
"""
import math as m
import numpy as np

class AUSMPlusupORIGINAL():
    def __init__(self, lft_flow_state, rght_flow_state, multi_species_flux, multi_phase) -> None:
        a_L_for_Ma = lft_flow_state.fluid_state.a
        p_L = lft_flow_state.fluid_state.p
        rho_L = lft_flow_state.fluid_state.rho
        u_L = lft_flow_state.fluid_state.u
        a_L = lft_flow_state.fluid_state.a
        vel_x_L = lft_flow_state.vel_x

        a_R_for_Ma = rght_flow_state.fluid_state.a
        p_R = rght_flow_state.fluid_state.p
        rho_R = rght_flow_state.fluid_state.rho
        u_R = rght_flow_state.fluid_state.u
        a_R = rght_flow_state.fluid_state.a
        vel_x_R = rght_flow_state.vel_x

        h_L = u_L + p_L / rho_L
        h_R = u_R + p_R / rho_R
        #print(a_L_forMa, a_R_forMa, p_L, rho_L, h_L, a_L, vel_x_L, p_R, rho_R, h_R, a_R, vel_x_R)
        H_L = h_L + 0.5 * vel_x_L ** 2.0
        H_R = h_R + 0.5 * vel_x_R ** 2.0
        
        self.sigma = 1.0
        self.K_p = 0.25
        self.K_u = 0.75
        self.M_inf = 1e-2
        self.beta = 0.125

        self.a_half_for_Ma = 0.5 * (a_L_for_Ma + a_R_for_Ma)
        self.a_half = 0.5 * (a_L + a_R)
        self.rho_half = 0.5 * (rho_L + rho_R)
        self.M_L = vel_x_L / self.a_half_for_Ma
        self.M_R = vel_x_R / self.a_half_for_Ma

        self.Mbar_sq = (vel_x_L ** 2.0 + vel_x_R ** 2.0) / (2.0 * self.a_half ** 2.0)
        #print(self.Mbar_sq)
        self.M0_sq = min(1.0, max(self.Mbar_sq, self.M_inf ** 2.0))
        #print(self.M0_sq)
        self.fa = self.M0_sq ** 0.5 * (2.0 - self.M0_sq ** 0.5)
        #print(self.fa)
        self.alpha = 0.1875 * (-4.0 + 5.0 * self.fa ** 2.0)
        #print(self.alpha)
        M4plus_ML = self.M_4_p(M = self.M_L, beta = self.beta)

        P5plus_ML = self.P_5_p(M = self.M_L, alpha = self.alpha)
        #print(P5plus_ML)
    
        M4minus_MR = self.M_4_m(M = self.M_R, beta = self.beta)
        P5minus_MR = self.P_5_m(M = self.M_R, alpha = self.alpha)
        #print(P5minus_MR)

        Mp = -1.0 * self.K_p * max(1.0 - self.sigma * self.Mbar_sq, 0.0) * (p_R - p_L) / (self.fa * self.rho_half * self.a_half ** 2.0)
        #print("max function", max(1.0) - self.sigma * self.Mbar_sq, 0.0)))
        Pu = -1.0 * self.K_u * P5plus_ML * P5minus_MR * (rho_L + rho_R) * self.fa * self.a_half * (vel_x_R - vel_x_L)
        #print(Pu)
        Ma_half = M4plus_ML + M4minus_MR + Mp
        #print(M4plus_ML, M4minus_MR, Mp)
        vel_x_half = self.a_half_for_Ma * Ma_half
        #print(self.a_half_forMa, Ma_half, self.a_half)
        self.fluxes = {}
        self.fluxes["Ma"] = Ma_half
        self.fluxes["vel_x"] = vel_x_half
        p_half = P5plus_ML * p_L + P5minus_MR * p_R + Pu
        if vel_x_half >= 0.0:
            self.fluxes["mass"] = vel_x_half * rho_L
            self.fluxes["xMom"] = vel_x_half * rho_L * vel_x_L
            self.fluxes["energy"] = vel_x_half * rho_L * H_L
            if multi_species_flux:
                massf_lft = lft_flow_state.fluid_state.massf
                self.fluxes["massf"] = (vel_x_half * rho_L * np.array(massf_lft)).tolist()
            
        else:
            self.fluxes["mass"] = vel_x_half * rho_R
            self.fluxes["xMom"] = vel_x_half * rho_R * vel_x_R
            self.fluxes["energy"] = vel_x_half * rho_R * H_R
            if multi_species_flux:
                massf_rght = rght_flow_state.fluid_state.massf
                self.fluxes["massf"] = (vel_x_half * rho_R * np.array(massf_rght)).tolist()
        self.fluxes["p"] = p_half
        
    def M_1_p(self, M):
        return 0.5 * (M + abs(M))
    
    def M_1_m(self, M): 
        return 0.5 * (M - abs(M))
    
    def M_2_p(self, M):
        return 0.25 * (M + 1.0) ** 2.0
    
    def M_2_m(self, M):
        return -0.25 * (M - 1.0) ** 2.0

    def M_4_p(self, M, beta):
        if abs(M) >= 1.0:
            return self.M_1_p(M = M)
        else:
            return self.M_2_p(M = M) * (1.0 - 16.0 * beta * self.M_2_m(M = M))
    
    def M_4_m(self, M, beta):
        if abs(M) >= 1.0:
            return self.M_1_m(M = M)
        else:
            return self.M_2_m(M = M) * (1.0 + 16.0 * beta * self.M_2_p(M = M))

    def P_5_p(self, M, alpha):
        if abs(M) >= 1.0:
            return self.M_1_p(M = M) / M
        else:
            return self.M_2_p(M = M) * (2.0 - M - 16.0 * alpha * M * self.M_2_m(M = M))
    
    def P_5_m(self, M, alpha):
        if abs(M) >= 1.0:
            return self.M_1_m(M = M) / M
        else:
            return self.M_2_m(M = M) * (-2.0 - M + 16.0 * alpha * M * self.M_2_p(M = M))
    

class AUSMPlusupPAPER():
    def __init__(self, lft_flow_state, rght_flow_state, multi_species_flux) -> None:
        a_L = lft_flow_state.fluid_state.a
        a_L_forMa = lft_flow_state.fluid_state.a
        p_L = lft_flow_state.fluid_state.p
        rho_L = lft_flow_state.fluid_state.rho
        vel_x_L = lft_flow_state.vel_x
        u_L = lft_flow_state.fluid_state.u
        h_L = u_L + p_L / rho_L
        H_L = h_L + 0.5 * vel_x_L ** 2.0

        a_R = rght_flow_state.fluid_state.a
        a_R_forMa = rght_flow_state.fluid_state.a
        p_R = rght_flow_state.fluid_state.p
        rho_R = rght_flow_state.fluid_state.rho
        vel_x_R = rght_flow_state.vel_x
        u_R = rght_flow_state.fluid_state.u
        h_R = u_R + p_R / rho_R
        H_R = h_R + 0.5 * vel_x_R ** 2.0

        self.K_p = 0.25
        self.K_u = 0.75

        h_L = u_L + p_L / rho_L
        h_R = u_R + p_R / rho_R
        
        H_L = h_L + 0.5 * vel_x_L ** 2.0
        H_R = h_R + 0.5 * vel_x_R ** 2.0

        self.a_half_forMa = 0.5 * (a_L_forMa + a_R_forMa)
        self.a_bar = 0.5 * (a_L + a_R)
        self.rho_bar = 0.5 * (rho_L + rho_R)
        self.M_L = vel_x_L / self.a_half_forMa
        self.M_R = vel_x_R / self.a_half_forMa
        self.Mbar = 0.5 * (self.M_L + self.M_R)
        self.Mbar_sq = self.Mbar ** 2.0
        #self.Mbar_sq = (vel_x_L ** 2.0) + vel_x_R ** 2.0)) / (2.0) * self.a_bar ** 2.0))
        
        M_4_p_ML = self.M_4_p(M = self.M_L)
        M_1_p_ML = self.M_1_p(M = self.M_L)
        P_5_p_Mbar = self.P_5_p(M = self.Mbar)
        P_5_p_ML = self.P_5_p(M = self.M_L)
        
        M_4_m_MR = self.M_4_m(M = self.M_R)
        M_1_m_MR = self.M_1_m(M = self.M_R)
        P_5_m_Mbar = self.P_5_m(M = self.Mbar)
        P_5_m_MR = self.P_5_m(M = self.M_R)
        self.Delta_M = M_4_p_ML - M_1_p_ML - M_4_m_MR + M_1_m_MR
        D_p = self.K_p * self.Delta_M * max(1.0 - self.Mbar_sq, 0.0) * (p_L - p_R) / self.a_bar
        D_vel_x = self.K_u * P_5_p_Mbar * P_5_m_Mbar * self.rho_bar * self.a_bar * (vel_x_L - vel_x_R)
        self.fluxes = {}
        self.M_half = M_4_p_ML + M_4_m_MR
        vel_x_half = self.a_half_forMa * self.M_half
        if vel_x_half > 0.0:
            massFlux = vel_x_half * rho_L
        else:
            massFlux = vel_x_half * rho_R
        self.fluxes["mass"] = massFlux + D_p
        p_half = P_5_p_ML * p_L + P_5_m_MR * p_R + D_vel_x
        if self.fluxes["mass"] > 0.0:
            self.fluxes["xMom"] = self.fluxes["mass"] * vel_x_L
            self.fluxes["energy"] = self.fluxes["mass"] * H_L
        else:
            self.fluxes["xMom"] = self.fluxes["mass"] * vel_x_R
            self.fluxes["energy"] = self.fluxes["mass"] * H_R
        self.fluxes["p"] = p_half

    def M_1_p(self, M):
        return 0.5 * (M + abs(M))
    
    def M_1_m(self, M):
        return 0.5 * (M - abs(M))

    def M_2_p(self, M):
        return 0.25 * (M + 1.0) ** 2.0

    def M_2_m(self, M):
        return -0.25 * (M - 1.0) ** 2.0
    
    def M_4_p(self, M):
        if abs(M) > 1.0:
            return self.M_1_p(M = M)
        else:
            return self.M_2_p(M = M) * (1.0 - 2.0 * self.M_2_m(M = M))
    
    def M_4_m(self, M):
        if abs(M) > 1.0:
            return self.M_1_m(M = M)
        else:
            return self.M_2_m(M = M) * (1.0 + 2.0 * self.M_2_p(M = M))
        
    def P_5_p(self, M):
        if abs(M) >= 1.0:
            return self.M_1_p(M = M) / M
        else:
            return self.M_2_p(M = M) * (2.0 - M - 3.0 * M * self.M_2_m(M = M))
    
    def P_5_m(self, M):
        if abs(M) >= 1.0:
            return self.M_1_m(M = M) / M
        else:
            return self.M_2_m(M = M) * (-2.0 - M + 3.0 * M * self.M_2_p(M = M))
    
class AUSMPlusM():
    def __init__(self, lft_flow_state, rght_flow_state, multi_species_flux):
        pass

class AUSM_U_splitting():
    def __init__(self, lft_flow_state, rght_flow_state, multi_species_flux) -> None:
        a_L = lft_flow_state.fluid_state.a
        p_L = lft_flow_state.fluid_state.p
        rho_L = lft_flow_state.fluid_state.rho
        vel_x_L = lft_flow_state.vel_x
        u_L = lft_flow_state.fluid_state.u
        h_L = u_L + p_L / rho_L
        H_L = h_L + 0.5 * vel_x_L ** 2.0

        a_R = rght_flow_state.fluid_state.a
        p_R = rght_flow_state.fluid_state.p
        rho_R = rght_flow_state.fluid_state.rho
        vel_x_R = rght_flow_state.vel_x
        u_R = rght_flow_state.fluid_state.u
        h_R = u_R + p_R / rho_R
        H_R = h_R + 0.5 * vel_x_R ** 2.0

        if abs(vel_x_L) <= a_L:
            vel_x_L_plus = 0.25 / a_L * (vel_x_L + a_L) ** 2.0
            p_L_plus_factor = (2.0 - vel_x_L / a_L) / a_L
            
        else:
            vel_x_L_plus = 0.5 * (vel_x_L + abs(vel_x_L))
            p_L_plus_factor = 1.0 / vel_x_L
            

        if abs(vel_x_R) <= a_R:
            vel_x_R_minus = - 0.25 / a_R * (vel_x_R - a_R) ** 2.0
            p_R_minus_factor = (-2.0 - vel_x_R / a_R) / a_R
            
        else:
            vel_x_R_minus = 0.5 * (vel_x_R - abs(vel_x_R))
            p_R_minus_factor = 1.0 / vel_x_R

        p_L_plus = p_L * vel_x_L_plus * p_L_plus_factor
        p_R_minus = p_R * vel_x_R_minus * p_R_minus_factor

        vel_x_half = vel_x_L_plus + vel_x_R
        
        
        self.fluxes = {}
        rhoU_half = 0.5 * (vel_x_half * (rho_L + rho_R) - abs(vel_x_half) * (rho_R - rho_L)) 
        xMom_flux = 0.5 * (vel_x_half * (rho_L * vel_x_L + rho_R * vel_x_R) - abs(vel_x_half) * (rho_R * vel_x_R - rho_L * vel_x_L)) 
        enthalpy_flux = 0.5 * (vel_x_half * (rho_L * H_L + rho_R * H_R) - abs(vel_x_half) * (rho_R * H_R - rho_L * H_L))
        p_half = p_L_plus + p_R_minus
        
        self.fluxes["p"] = p_half
        self.fluxes["mass"] = rhoU_half
        self.fluxes["xMom"] = xMom_flux
        self.fluxes["energy"] = enthalpy_flux

        if multi_species_flux:
            if rhoU_half >= 0.0:
                rhoU_species_half = (rhoU_half * np.array(lft_flow_state.fluid_state.massf)).tolist()
            else:
                rhoU_species_half = (rhoU_half * np.array(rght_flow_state.fluid_state.massf)).tolist()
            self.fluxes["massf"] = rhoU_species_half

class AUSM_M_splitting():
    def __init__(self, lft_flow_state, rght_flow_state, multi_species_flux) -> None:
        a_L = lft_flow_state.fluid_state.a
        p_L = lft_flow_state.fluid_state.p
        rho_L = lft_flow_state.fluid_state.rho
        vel_x_L = lft_flow_state.vel_x
        u_L = lft_flow_state.fluid_state.u
        h_L = u_L + p_L / rho_L
        H_L = h_L + 0.5 * vel_x_L ** 2.0

        a_R = rght_flow_state.fluid_state.a
        p_R = rght_flow_state.fluid_state.p
        rho_R = rght_flow_state.fluid_state.rho
        vel_x_R = rght_flow_state.vel_x
        u_R = rght_flow_state.fluid_state.u
        h_R = u_R + p_R / rho_R
        H_R = h_R + 0.5 * vel_x_R ** 2.0

        if abs(vel_x_L) <= a_L:
            vel_x_L_plus = 0.25 / a_L * (vel_x_L + a_L) ** 2.0
            p_L_plus_factor = (2.0 - vel_x_L / a_L) / a_L
            
        else:
            vel_x_L_plus = 0.5 * (vel_x_L + abs(vel_x_L))
            p_L_plus_factor = 1.0 / vel_x_L
            

        if abs(vel_x_R) <= a_R:
            vel_x_R_minus = - 0.25 / a_R * (vel_x_R - a_R) ** 2.0
            p_R_minus_factor = (-2.0 - vel_x_R / a_R) / a_R
            
        else:
            vel_x_R_minus = 0.5 * (vel_x_R - abs(vel_x_R))
            p_R_minus_factor = 1.0 / vel_x_R
        
        p_L_plus = p_L * vel_x_L_plus * p_L_plus_factor
        p_R_minus = p_R * vel_x_R_minus * p_R_minus_factor

        self.fluxes = {}
        p_half = p_L_plus + p_R_minus
        M_half = vel_x_L_plus / a_L + vel_x_R_minus / a_R
        rhoU_half = 0.5 * (M_half * (rho_L * a_L + rho_R * a_R) - abs(M_half) * (rho_R * a_R - rho_L * a_L))
        xMom_flux = 0.5 * (M_half * (rho_L * a_L * vel_x_L + rho_R * a_R * vel_x_R) - abs(M_half) * (rho_R * a_R * vel_x_R - rho_L * a_L * vel_x_L))
        enthalpy_flux = 0.5 * (M_half * (rho_L * a_L * H_L + rho_R * a_R * H_R) - abs(M_half) * (rho_R * a_R * H_R - rho_L * a_L * H_L))

        self.fluxes["mass"] = rhoU_half
        self.fluxes["p"] = p_half
        self.fluxes["xMom"] = xMom_flux
        self.fluxes["energy"] = enthalpy_flux

        if multi_species_flux:
            if rhoU_half >= 0.0:
                rhoU_species_half = (rhoU_half * np.array(lft_flow_state.fluid_state.massf)).tolist()
            else:
                rhoU_species_half = (rhoU_half * np.array(rght_flow_state.fluid_state.massf)).tolist()
            self.fluxes["massf"] = rhoU_species_half

class AUSMD():
    def __init__(self, lft_flow_state, rght_flow_state, multi_species_flux) -> None:
        a_L = lft_flow_state.fluid_state.a
        p_L = lft_flow_state.fluid_state.p
        rho_L = lft_flow_state.fluid_state.rho
        vel_x_L = lft_flow_state.vel_x
        u_L = lft_flow_state.fluid_state.u
        h_L = u_L + p_L / rho_L
        H_L = h_L + 0.5 * vel_x_L ** 2.0

        a_R = rght_flow_state.fluid_state.a
        p_R = rght_flow_state.fluid_state.p
        rho_R = rght_flow_state.fluid_state.rho
        vel_x_R = rght_flow_state.vel_x
        u_R = rght_flow_state.fluid_state.u
        h_R = u_R + p_R / rho_R
        H_R = h_R + 0.5 * vel_x_R ** 2.0
        
        a_m = max(a_L, a_R)
        alpha_L = 2.0 * (p_L / rho_L) / ((p_L / rho_L) + (p_R / rho_R))
        alpha_R = 2.0 * (p_R / rho_R) / ((p_L / rho_L) + (p_R / rho_R))

        if abs(vel_x_L) <= a_m:
            vel_x_L_plus = alpha_L * (vel_x_L + a_m) ** 2.0 / (4.0 * a_m) \
                            + 0.5 * (1.0 - alpha_L) * (vel_x_L + abs(vel_x_L))
            p_L_plus = p_L * (vel_x_L + a_m) ** 2.0 * (2.0 - vel_x_L / a_m) / (4.0 * a_m ** 2.0)
        else:
            vel_x_L_plus = 0.5 * (vel_x_L + abs(vel_x_L))
            p_L_plus = 0.5 * p_L * (vel_x_L + abs(vel_x_L)) / vel_x_L

        if abs(vel_x_R) <= a_m:
            vel_x_R_minus = - alpha_R * (vel_x_R - a_m) ** 2.0 / (4.0 * a_m) \
                            + 0.5 * (1.0 - alpha_R) * (vel_x_R - abs(vel_x_R)) 
            p_R_minus = p_R * (vel_x_R - a_m) ** 2.0 * (2.0 + vel_x_R / a_m) / (4.0 * a_m ** 2.0)
        else:
            vel_x_R_minus = 0.5 * (vel_x_R - abs(vel_x_R))
            p_R_minus = 0.5 * p_R * (vel_x_R - abs(vel_x_R)) / vel_x_R
        
        self.fluxes = {}
        p_half = p_L_plus + p_R_minus
        rhoU_half = vel_x_L_plus * rho_L + vel_x_R_minus * rho_R
        xMom_flux = 0.5 * (rhoU_half * (vel_x_L + vel_x_R) - abs(rhoU_half) * (vel_x_R - vel_x_L))
        enthalpy_flux = 0.5 * (rhoU_half * (H_L + H_R) - abs(rhoU_half) * (H_R - H_L))

        self.fluxes["mass"] = rhoU_half
        self.fluxes["xMom"] = xMom_flux
        self.fluxes["energy"] = enthalpy_flux
        self.fluxes["p"] = p_half

        if multi_species_flux:
            if rhoU_half >= 0.0:
                rhoU_species_half = (rhoU_half * np.array(lft_flow_state.fluid_state.massf)).tolist()
            else:
                rhoU_species_half = (rhoU_half * np.array(rght_flow_state.fluid_state.massf)).tolist()
            self.fluxes["massf"] = rhoU_species_half

class AUSMV():
    def __init__(self, lft_flow_state, rght_flow_state, multi_species_flux) -> None:
        a_L = lft_flow_state.fluid_state.a
        p_L = lft_flow_state.fluid_state.p
        rho_L = lft_flow_state.fluid_state.rho
        vel_x_L = lft_flow_state.vel_x
        u_L = lft_flow_state.fluid_state.u
        h_L = u_L + p_L / rho_L
        H_L = h_L + 0.5 * vel_x_L ** 2.0

        a_R = rght_flow_state.fluid_state.a
        p_R = rght_flow_state.fluid_state.p
        rho_R = rght_flow_state.fluid_state.rho
        vel_x_R = rght_flow_state.vel_x
        u_R = rght_flow_state.fluid_state.u
        h_R = u_R + p_R / rho_R
        H_R = h_R + 0.5 * vel_x_R ** 2.0

        a_m = max(a_L, a_R)
        alpha_L = 2.0 * (p_L / rho_L) / ((p_L / rho_L) + (p_R / rho_R))
        alpha_R = 2.0 * (p_R / rho_R) / ((p_L / rho_L) + (p_R / rho_R))

        if abs(vel_x_L) <= a_m:
            vel_x_L_plus = alpha_L * (vel_x_L + a_m) ** 2.0 / (4.0 * a_m) \
                            + 0.5 * (1.0 - alpha_L) * (vel_x_L + abs(vel_x_L))
            p_L_plus = p_L * (vel_x_L + a_m) ** 2.0 * (2.0 - vel_x_L / a_m) / (4.0 * a_m ** 2.0)
        else:
            vel_x_L_plus = 0.5 * (vel_x_L + abs(vel_x_L))
            p_L_plus = 0.5 * p_L * (vel_x_L + abs(vel_x_L)) / vel_x_L

        if abs(vel_x_R) <= a_m:
            vel_x_R_minus = - alpha_R * (vel_x_R - a_m) ** 2.0 / (4.0 * a_m) \
                            + 0.5 * (1.0 - alpha_R) * (vel_x_R - abs(vel_x_R)) 
            p_R_minus = p_R * (vel_x_R - a_m) ** 2.0 * (2.0 + vel_x_R / a_m) / (4.0 * a_m ** 2.0)
        else:
            vel_x_R_minus = 0.5 * (vel_x_R - abs(vel_x_R))
            p_R_minus = 0.5 * p_R * (vel_x_R - abs(vel_x_R)) / vel_x_R

        self.fluxes = {}
        p_half = p_L_plus + p_R_minus
        rhoU_half = vel_x_L_plus * rho_L + vel_x_R_minus * rho_R
        xMom_flux = vel_x_L_plus * rho_L * vel_x_L + vel_x_R_minus * rho_R * vel_x_R
        enthalpy_flux = 0.5 * (rhoU_half * (H_L + H_R) - abs(rhoU_half) * (H_R - H_L))
        
        if multi_species_flux:
            if rhoU_half >= 0.0:
                rhoU_species_half = (rhoU_half * np.array(lft_flow_state.fluid_state.massf)).tolist()
            else:
                rhoU_species_half = (rhoU_half * np.array(rght_flow_state.fluid_state.massf)).tolist()
            self.fluxes["massf"] = rhoU_species_half
        
        self.fluxes["mass"] = rhoU_half
        self.fluxes["xMom"] = xMom_flux
        self.fluxes["energy"] = enthalpy_flux
        self.fluxes["p"] = p_half


class AUSMDV():
    def __init__(self, lft_flow_state, rght_flow_state, multi_species_flux) -> None:
        K = 10.0
        a_L = lft_flow_state.fluid_state.a
        p_L = lft_flow_state.fluid_state.p
        rho_L = lft_flow_state.fluid_state.rho
        vel_x_L = lft_flow_state.vel_x
        u_L = lft_flow_state.fluid_state.u
        h_L = u_L + p_L / rho_L
        H_L = h_L + 0.5 * vel_x_L ** 2.0

        a_R = rght_flow_state.fluid_state.a
        p_R = rght_flow_state.fluid_state.p
        rho_R = rght_flow_state.fluid_state.rho
        vel_x_R = rght_flow_state.vel_x
        u_R = rght_flow_state.fluid_state.u
        h_R = u_R + p_R / rho_R
        H_R = h_R + 0.5 * vel_x_R ** 2.0

        s = min(1.0, K * abs(p_R - p_L) / min(p_L, p_R))

        a_m = max(a_L, a_R)
        alpha_L = 2.0 * (p_L / rho_L) / ((p_L / rho_L) + (p_R / rho_R))
        alpha_R = 2.0 * (p_R / rho_R) / ((p_L / rho_L) + (p_R / rho_R))

        if abs(vel_x_L) <= a_m:
            vel_x_L_plus = alpha_L * (vel_x_L + a_m) ** 2.0 / (4.0 * a_m) \
                            + 0.5 * (1.0 - alpha_L) * (vel_x_L + abs(vel_x_L))
            p_L_plus = p_L * (vel_x_L + a_m) ** 2.0 * (2.0 - vel_x_L / a_m) / (4.0 * a_m ** 2.0)
        else:
            vel_x_L_plus = 0.5 * (vel_x_L + abs(vel_x_L))
            p_L_plus = 0.5 * p_L * (vel_x_L + abs(vel_x_L)) / vel_x_L

        if abs(vel_x_R) <= a_m:
            vel_x_R_minus = - alpha_R * (vel_x_R - a_m) ** 2.0 / (4.0 * a_m) \
                            + 0.5 * (1.0 - alpha_R) * (vel_x_R - abs(vel_x_R)) 
            p_R_minus = p_R * (vel_x_R - a_m) ** 2.0 * (2.0 + vel_x_R / a_m) / (4.0 * a_m ** 2.0)
        else:
            vel_x_R_minus = 0.5 * (vel_x_R - abs(vel_x_R))
            p_R_minus = 0.5 * p_R * (vel_x_R - abs(vel_x_R)) / vel_x_R
        
        self.fluxes = {}
        p_half = p_L_plus + p_R_minus
        rhoU_half = vel_x_L_plus * rho_L + vel_x_R_minus * rho_R
        xMom_flux_AUSMV = vel_x_L_plus * rho_L * vel_x_L + vel_x_R_minus * rho_R * vel_x_R
        xMom_flux_AUSMD = 0.5 * (rhoU_half * (vel_x_L + vel_x_R) - abs(rhoU_half) * (vel_x_R - vel_x_L))
        xMom_flux = 0.5 * (1.0 + s) * xMom_flux_AUSMV + 0.5 * (1.0 - s) * xMom_flux_AUSMD
        enthalpy_flux = 0.5 * (rhoU_half * (H_L + H_R) - abs(rhoU_half) * (H_R - H_L))
        
        self.fluxes["p"] = p_half
        self.fluxes["mass"] = rhoU_half
        self.fluxes["xMom"] = xMom_flux
        self.fluxes["energy"] = enthalpy_flux

        if multi_species_flux:
            if rhoU_half >= 0.0:
                rhoU_species_half = (rhoU_half * np.array(lft_flow_state.fluid_state.massf)).tolist()
            else:
                rhoU_species_half = (rhoU_half * np.array(rght_flow_state.fluid_state.massf)).tolist()
            self.fluxes["massf"] = rhoU_species_half

        
class HLLC():
    def __init__(self) -> None:
        pass
    
class EFMflx():
    def __init__(self) -> None:
        pass

class LDFSS():
    def __init__(self) -> None:
        pass

class Hanel():
    def __init__(self) -> None:
        pass

