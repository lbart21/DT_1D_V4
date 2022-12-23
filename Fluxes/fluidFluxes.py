import math as m
import numpy as np

class AUSMPlusupORIGINAL():
    def __init__(self, a_L_forMa, a_R_forMa, p_L, rho_L, u_L, a_L, vel_x_L, p_R, rho_R, u_R, a_R, vel_x_R) -> None:
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

        self.a_half_forMa = 0.5 * (a_L_forMa + a_R_forMa)
        self.a_half = 0.5 * (a_L + a_R)
        self.rho_half = 0.5 * (rho_L + rho_R)
        self.M_L = vel_x_L / self.a_half_forMa
        self.M_R = vel_x_R / self.a_half_forMa

        self.Mbar_sq = (vel_x_L ** 2.0 + vel_x_R ** 2.0) / (2.0 * self.a_half ** 2.0)
        self.M0_sq = min(1.0, max(self.Mbar_sq, self.M_inf))
        self.fa = self.M0_sq ** 0.5 * (2.0 - self.M0_sq ** 0.5)
        self.alpha = 0.1875 * (-4.0 + 5.0 * self.fa ** 2.0)
        M4plus_ML = self.M_4_p(M = self.M_L, beta = self.beta)
        P5plus_ML = self.P_5_p(M = self.M_L, alpha = self.alpha)
    
        M4minus_MR = self.M_4_m(M = self.M_R, beta = self.beta)
        P5minus_MR = self.P_5_m(M = self.M_R, alpha = self.alpha)

        Mp = -1.0 * self.K_p * max(1.0 - self.sigma * self.Mbar_sq, 0.0) * (p_R - p_L) / (self.fa * self.rho_half * self.a_half ** 2.0)
        #print("max function", max(1.0) - self.sigma * self.Mbar_sq, 0.0)))
        Pu = -1.0 * self.K_u * P5plus_ML * P5minus_MR * (rho_L + rho_R) * self.fa * self.a_half * (vel_x_R - vel_x_L)

        Ma_half = M4plus_ML + M4minus_MR + Mp
        vel_x_half = self.a_half_forMa * Ma_half
        self.fluxes = {}
        self.fluxes["Ma"] = Ma_half
        self.fluxes["vel_x"] = vel_x_half
        p_half = P5plus_ML * p_L + P5minus_MR * p_R + Pu
        if vel_x_half >= 0.0:
            self.fluxes["mass"] = vel_x_half * rho_L
            self.fluxes["xMom"] = vel_x_half * rho_L * vel_x_L
            self.fluxes["energy"] = vel_x_half * rho_L * H_L
            
        else:
            self.fluxes["mass"] = vel_x_half * rho_R
            self.fluxes["xMom"] = vel_x_half * rho_R * vel_x_R
            self.fluxes["energy"] = vel_x_half * rho_R * H_R        
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
    def __init__(self, a_L_forMa, a_L, p_L, h_L, rho_L, vel_x_L, a_R_forMa, a_R, p_R, h_R, rho_R, vel_x_R) -> None:
        self.K_p = 0.25
        self.K_u = 0.75

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
    def __init__(self, a_L_forMa, p_L, h_L, a_L, vel_x_L, gamma_L, a_R_forMa, p_R, h_R, a_R, vel_x_R, gamma_R, p_L_left, p_R_left, p_L_right, p_R_right):
        self.a_half_forMa = 0.5 * (a_L_forMa + a_R_forMa)
        self.M_L = vel_x_L / self.a_half_forMa
        self.M_R = vel_x_R / self.a_half_forMa
        self.M_inf = 1e-2
        gamma = 0.5 * (gamma_L + gamma_R)
        H_L = h_L + 0.5 * vel_x_L ** 2.0 #total enthalpy
        H_R = h_R + 0.5 * vel_x_R ** 2.0 #total enthalpy
        h_normal = 0.5 * (H_L + H_R - 0.5 * vel_x_L ** 2.0 - 0.5 * vel_x_R ** 2.0) #equation 29
        self.a_s = (2.0 * (gamma - 1.0) / (gamma + 1.0) * h_normal) ** 0.5 #equation 29
        if 0.5 * (vel_x_L + vel_x_R) >= 0.0:                      # equation 28
            self.a_half = self.a_s ** 2.0 / max(abs(vel_x_L), self.a_s)    # equation 28
        else:                                                                       # equatoin 28
            self.a_half = self.a_s ** 2.0 / max(abs(vel_x_R), self.a_s)    # equation 28
        M_for_f = min(1.0, max(abs(self.M_L), abs(self.M_R))) # equation 15
        f = 0.5 * (1.0 - m.cos(m.pi * M_for_f)) # equation 15
        f_0 = min(1.0, max(f, self.M_inf ** 2.0))
        h_left = min(p_L_left / p_R_left, p_R_left / p_L_left)
        h_mid = min(p_L / p_R, p_R / p_L)
        h_right = min(p_L_right / p_R_right, p_R_right / p_L_right)
        h = min(h_left, h_mid, h_right)
        g = 0.5 * (1.0 + m.cos(m.pi * h)) # equaiton 25
        rho_half = gamma * (p_L + p_R) / (2.0 * self.a_half ** 2.0)
        Mp = -0.5 * (1.0 - f) * (p_R - p_L) * (1.0 - g) / (self.rho_half * self.a_half ** 2.0) # equation 14
        P_L_p_ML = self.P_p(M = self.M_L)
        P_R_m_MR = self.P_m(M = self.M_R)
        p_s = 0.5 * (p_L + p_R) + 0.5 * (P_L_p_ML - P_R_m_MR) * (p_L - p_R) + 0.5 * f_0 * (P_L_p_ML + P_R_m_MR - 1.0) * (p_L + p_R) #equation 19
        p_vel_x = -1.0 * g * gamma * (p_L + p_R) * P_L_p_ML * P_R_m_MR * (vel_x_R - vel_x_L) / (2.0 * self.a_half)# equation 26
        self.p_half = p_s + p_vel_x #equation 27

        pass

    def P_p(self, M):
        if abs(M) >= 1.0:
            return 0.5 * (1.0 + np.sign(M))
        else:
            return 0.25 * (M + 1.0) ** 2.0 * (2.0 - M) + self.alpha * M * (M ** 2.0 - 1.0) ** 2.0

    def P_m(self, M):
        if abs(M) >= 1.0:
            return 0.5 * (1.0 - np.sign(M))
        else:
            return 0.25 * (M - 1.0) ** 2.0 * (2.0 + M) - self.alpha * M * (M ** 2.0 - 1.0) ** 2.0

    def M_p(self, M):
        if abs(M) >= 1.0:
            return 0.5 * (M + abs(M))
        else:
            return 0.25 * (M + 1.0) ** 2.0 + self.beta * (M ** 2.0 - 1.0) ** 2.0
        

    def M_m(self, M):
        if abs(M) >= 1.0:
            return 0.5 * (M - abs(M))
        else:
            return -0.25 * (M - 1.0) ** 2.0 - self.beta * (M ** 2.0 - 1.0) ** 2.0

class AUSM():
    def __init__(self) -> None:
        pass

class AUSMD():
    def __init__(self) -> None:
        pass

class AUSMV():
    def __init__(self) -> None:
        pass

class AUSMDV():
    def __init__(self, LftState, RghtState) -> None:
        self.fluxes = {}
        rL, rR = LftState.fs["rho"], RghtState.fs["rho"]
        pL, pR = LftState.fs["p"], RghtState.fs["p"]
        vel_xL, vel_xR = LftState.fs["vel_x"], RghtState.fs["vel_x"]
        aL, aR = LftState.fs["a"], RghtState.fs["a"]
        hL, hR = LftState.fs["h"], RghtState.fs["h"]
        HL, HR = hL + 0.5 * vel_xL ** 2.0, hR + 0.5 * vel_xR ** 2.0
        pLrL = pL / rL
        pRrR = pR / rR
        alphaL = 2.0 * pLrL / (pLrL + pRrR)
        alphaR = 2.0 * pRrR / (pLrL + pRrR)
        a_m = max(aL, aR)
        M_L = vel_xL / a_m
        M_R = vel_xR / a_m

        dvel_xL = 0.5 * (vel_xL + abs(vel_xL))
        
        if abs(M_L) <= 1.0:
            pLplus = 0.25 * pL * (M_L + 1.0) ** 2.0 * (2.0 - M_L)
            vel_xLplus = alphaL * ((vel_xL + a_m) * (vel_xL + a_m) / (4.0 * a_m) - dvel_xL) + dvel_xL
        else:
            pLplus = pL * dvel_xL / vel_xL
            vel_xLplus = dvel_xL

        dvel_xR = 0.5 * (vel_xR - abs(vel_xR))
        if abs(M_R) <= 1.0:
            pRminus = 0.25 * pR * (M_R - 1.0) ** 2.0 * (2.0 + M_R) 
            vel_xRminus = alphaR * (-(vel_xR - a_m) * (vel_xR - a_m) / (4.0 * a_m) - dvel_xR) + dvel_xR;
        else:
            pRminus = pR * dvel_xR / vel_xR
            vel_xRminus = dvel_xR
        
        ru_half = vel_xLplus * rL + vel_xRminus * rR
        p_half = pLplus + pRminus
        K = 10.0
        dp = pL - pR
        dp = K * abs(dp) / min(pL, pR)
        s = 0.5 * min(1.0, dp)
        ru2_AUSMV = vel_xLplus * rL * vel_xL + vel_xRminus * rR * vel_xR
        ru2_AUSMD = 0.5 * (ru_half * (vel_xL + vel_xR) - abs(ru_half) * (vel_xR - vel_xL))
        ru2_half = (0.5 + s) * ru2_AUSMV + (0.5 - s) * ru2_AUSMD

        self.fluxes["mass"] = ru_half

        if ru_half >= 0.0:
            self.fluxes["xMom"] = ru2_half + p_half
            self.fluxes["energy"] = ru_half * HL
        
        else:
            self.fluxes["xMom"] = ru2_half + p_half
            self.fluxes["energy"] = ru_half * HR
        
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

