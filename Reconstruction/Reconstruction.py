import math as m
from Algorithms.DesignToolAlgorithmV3_1D.Reconstruction.Limiters import getLimiters


class Copy_Recon():
    def __init__(self) -> None:
        pass
    
    def Recon(self, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0] = qL_stencil
        [q_R0] = qR_stencil
        
        Lft = q_L0
        Rght = q_R0
        return Lft, Rght

class eilmer_Linear_Recon():
    def __init__(self) -> None:
        pass
    
    def Recon(self, limiter, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil
        [len_L0, len_L1] = dxL_stencil
        [q_R0, q_R1] = qR_stencil
        [len_R0, len_R1] = dxR_stencil

        w_0 = len_L0
        w_1 = len_R0
        w_2 = 0.5 * len_L0 / (len_L1 + 2.0 * len_L0 + len_R0)
        w_3 = 0.5 * len_R0 / (len_L0 + 2.0 * len_R0 + len_R1)
        w_4 = (2.0 * len_L0 + len_L1)
        w_5 = (2.0 * len_R0 + len_R1)
        w_6 = 2.0 / (len_L0 + len_L1)
        w_7 = 2.0 / (len_R0 + len_L0)
        w_8 = 2.0 / (len_R1 + len_R0)
        DelLminus = (q_L0 - q_L1) * w_6
        Del = (q_R0 - q_L0) * w_7
        DelRplus = (q_R1 - q_R0) * w_8
        ### Apply limiter
        sL, sR = 1.0, 1.0
        if limiter != None:
            sL, sR = getLimiters(limiter = limiter, qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                                    qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
        
        Lft = q_L0 + sL * w_2 * (Del * w_4 + DelLminus * w_1)
        Rght = q_R0 - sR * w_3 * (DelRplus * w_0 + Del * w_5)
        ### Apply clipping
        Lft, Rght = clipping().clipToBounds(q_L0 = q_L0, q_R0 = q_R0, Lft = Lft, Rght = Rght)
        return Lft, Rght

class eilmer_Quadratic():
    def __init__(self) -> None:
        pass
    
    def Recon(self, limiter, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil
        [len_L0, len_L1] = dxL_stencil
        [q_R0, q_R1] = qR_stencil
        [len_R0, len_R1] = dxR_stencil
       
        x_L1 = - (len_L0 + 0.5 * len_L1)
        x_L0 = - 0.5 * len_L0
        x_R0 = 0.5 * len_R0
        x_R1 = len_R0 + 0.5 * len_R1

        w_0 = (x_L0 * x_R0) / ((x_L1 - x_L0) * (x_L1 - x_R0))
        w_1 = (x_L1 * x_R0) / ((x_L0 - x_L1) * (x_L0 - x_R0))
        w_2 = (x_R0 * x_R1) / ((x_L0 - x_R0) * (x_L0 - x_R1))
        w_3 = (x_L1 * x_L0) / ((x_R0 - x_L1) * (x_R0 - x_L0))
        w_4 = (x_L0 * x_R1) / ((x_R0 - x_L0) * (x_R0 - x_R1))
        w_5 = (x_L0 * x_R0) / ((x_R1 - x_L0) * (x_R1 - x_R0))

        ### Apply limiter
        sL, sR = 1.0, 1.0
        if limiter != None:
            sL, sR = getLimiters(limiter = limiter, qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                                    qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
          
        Lft = q_L0 + sL * (q_L1 * w_0 + q_L0 * (w_1 - 1) + q_R0 * w_3)
        Rght = q_R0 + sR * (q_L0 * w_2 + q_R0 * (w_4 - 1) + q_R1 * w_5)
        ### Apply clipping
        Lft, Rght = clipping().clipToBounds(q_L0 = q_L0, q_R0 = q_R0, Lft = Lft, Rght = Rght)
        return Lft, Rght

class eilmer_Cubic():
    def __init__(self) -> None:
        pass
    
    def Recon(self, limiter, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1, q_L2] = qL_stencil
        [len_L0, len_L1, len_L2] = dxL_stencil
        [q_R0, q_R1, q_R2] = qR_stencil
        [len_R0, len_R1, len_R2] = dxR_stencil
        
        ### Define node locations
        x_L2 = -1.0 * (len_L0 + len_L1 + 0.5 * len_L2)
        x_L1 = -1.0 * (len_L0 + 0.5 * len_L1)
        x_L0 = -0.5 * len_L0
        x_R0 = 0.5 * len_R0
        x_R1 = len_R0 + 0.5 * len_R1
        x_R2 = len_R0 + len_R1 + 0.5 * len_R2
        ### Form weights
        w_0 = x_L1 * x_L0 * x_R0 * x_R1 / ((x_L2 - x_L1) * (x_L2 - x_L0) * (x_L2 - x_R0) * (x_L2 - x_R1))
        w_1 = x_L2 * x_L0 * x_R0 * x_R1 / ((x_L1 - x_L2) * (x_L1 - x_L0) * (x_L1 - x_R0) * (x_L1 - x_R1))
        w_2 = x_L2 * x_L1 * x_R0 * x_R1 / ((x_L0 - x_L2) * (x_L0 - x_L1) * (x_L0 - x_R0) * (x_L0 - x_R1))
        w_3 = x_L2 * x_L1 * x_L0 * x_R1 / ((x_R0 - x_L2) * (x_R0 - x_L1) * (x_R0 - x_L0) * (x_R0 - x_R1))
        w_4 = x_L2 * x_L1 * x_L0 * x_R0 / ((x_R1 - x_L2) * (x_R1 - x_L1) * (x_R1 - x_L0) * (x_R1 - x_R0))
        w_5 = x_L0 * x_R0 * x_R1 * x_R2 / ((x_L1 - x_L0) * (x_L1 - x_R0) * (x_L1 - x_R1) * (x_L1 - x_R2))
        w_6 = x_L1 * x_R0 * x_R1 * x_R2 / ((x_L0 - x_L1) * (x_L0 - x_R0) * (x_L0 - x_R1) * (x_L0 - x_R2))
        w_7 = x_L1 * x_L0 * x_R1 * x_R2 / ((x_R0 - x_L1) * (x_R0 - x_L0) * (x_R0 - x_R1) * (x_R0 - x_R2))
        w_8 = x_L1 * x_L0 * x_R0 * x_R2 / ((x_R1 - x_L1) * (x_R1 - x_L0) * (x_R1 - x_R0) * (x_R1 - x_R2))
        w_9 = x_L1 * x_L0 * x_R0 * x_R1 / ((x_R2 - x_L1) * (x_R2 - x_L0) * (x_R2 - x_R0) * (x_R2 - x_R1))
        
        ### Apply limiter
        sL, sR = 1.0, 1.0
        if limiter != None:
            sL, sR = getLimiters(limiter = limiter, qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                                    qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

        Lft = q_L0 + sL * (w_0 * q_L2 + w_1 * q_L1 + (w_2 - 1.0) * q_L0 + w_3 * q_R0 + w_4 * q_R1)
        Rght = q_R0 + sR * (w_5 * q_L1 + w_6 * q_L0 + (w_7 - 1.0) * q_R0 + w_8 * q_R1 + w_9 * q_R2)

        ### Apply clipping
        Lft, Rght = clipping().clipToBounds(q_L0 = q_L0, q_R0 = q_R0, Lft = Lft, Rght = Rght)
        return Lft, Rght

class MUSCL():
    def __init__(self) -> None:
        self.k = 1/3
        pass
    
    def Recon(self, limiter, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        k = self.k
        [q_L0, q_L1] = qL_stencil
        [len_L0, len_L1] = dxL_stencil
        [q_R0, q_R1] = qR_stencil
        [len_R0, len_R1] = dxR_stencil
    
        ### Apply limiter
        sL, sR = 1.0, 1.0
        if limiter != None:
            sL, sR = getLimiters(limiter = limiter, qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                                    qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

        Lft = q_L0 + sL / 4.0 * ((1.0 - k) * (q_L0 - q_L1) + (1.0 + k) * (q_R0 - q_L0))
        Rght = q_R0 - sR / 4.0 * ((1.0 + k) * (q_R0 - q_L0) + (1.0 - k) * (q_R1 - q_R0))
        
        return Lft, Rght

class THINC():
    def __init__(self) -> None:
        self.epsilon = 1e-20
        self.beta = 1.6
        pass
    
    def Recon(self, limiter, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        epsilon = self.epsilon
        beta = self.beta
        [q_L0, q_L1] = qL_stencil
        [len_L0, len_L1] = dxL_stencil
        [q_R0, q_R1] = qR_stencil
        [len_R0, len_R1] = dxR_stencil
        
        if q_L1 == q_R0:
            Lft = q_L1
        else:
            q_min_L = min(q_L1, q_R0)
            q_max_L = max(q_L1, q_R0) - q_min_L
            if abs(q_max_L) < 1e-10:
                Lft = q_min_L
            elif abs(q_L1 - q_L0) < 1e-10 or abs(q_L0 - q_R0) < 1e-10:
                Lft = q_L0
            elif q_L0 > max(q_L1, q_R0) or q_L0 < min(q_L1, q_R0):
                Lft = q_L0
            else:
                THINC_Theta_L = (q_R0 - q_L1) / (abs(q_R0 - q_L1) + epsilon)

                THINC_C_L = (q_L0 - q_min_L + epsilon) / (q_max_L + epsilon)
                THINC_B_L = m.exp(THINC_Theta_L * beta * (2.0 * THINC_C_L - 1.0))
                THINC_A_L = (THINC_B_L / m.cosh(beta) - 1.0) / m.tanh(beta)
                Lft = q_min_L + 0.5 * q_max_L * (1.0 + THINC_Theta_L * (THINC_A_L + m.tanh(beta)) / (THINC_A_L * m.tanh(beta) + 1.0)) 
        if q_L0 == q_R1:
            Rght = q_L0
        else:
            q_min_R = min(q_L0, q_R1)
            q_max_R = max(q_L0, q_R1) - q_min_R
            if abs(q_max_R) < 1e-10:
                Rght = q_min_R
            elif abs(q_L0 - q_R0) < 1e-10 or abs(q_R0 - q_R1) < 1e-10:
                Rght = q_R0
            elif q_R0 > max(q_L0, q_R1) or q_R0 < min(q_L0, q_R1):
                Rght = q_R0
            else:
                THINC_Theta_R = (q_R1 - q_L0) / (abs(q_R1 - q_L0) + epsilon)

                THINC_C_R = (q_R0 - q_min_R + epsilon) / (q_max_R + epsilon)
                THINC_B_R = m.exp(THINC_Theta_R * beta * (2.0 * THINC_C_R - 1.0))
                THINC_A_R = (THINC_B_R / m.cosh(beta) - 1.0) / m.tanh(beta)
                Rght = q_min_R + 0.5 * q_max_R * (1.0 + THINC_Theta_R * THINC_A_R)
        ### Apply limiter
        sL, sR = 1.0, 1.0
        if limiter != None:
            sL, sR = getLimiters(limiter = limiter, qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                                    qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

        Lft = q_L0 + sL * (Lft - q_L0)
        Rght = q_R0 + sR * (Rght - q_R0)
        ### Apply clipping
        Lft, Rght = clipping().clipToBounds(q_L0 = q_L0, q_R0 = q_R0, Lft = Lft, Rght = Rght)
        return Lft, Rght

class MUSCLTHINCBVD():
    def __init__(self) -> None:
        pass
    
    def Recon(self, limiter, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        Lft, Rght = None, None
        return Lft, Rght

class clipping():
    def __init__(self) -> None:
        pass

    def clipToBounds(self, q_L0, q_R0, Lft, Rght):
        q_min = min(q_L0, q_R0)
        q_max = max(q_L0, q_R0)
        if Lft < q_min:
            Lft = q_min
        elif Lft > q_max:
            Lft = q_max
        if Rght < q_min:
            Rght = q_min
        elif Rght > q_max:
            Rght = q_max
        return Lft, Rght

def getReconstruction(reconstruction, limiter, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    if reconstruction == "Copy":
        Lft, Rght = Copy_Recon().Recon(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif reconstruction == "eilmer_linear":
        Lft, Rght = eilmer_Linear_Recon().Recon(limiter = limiter, qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif reconstruction == "eilmer_quadratic":
        Lft, Rght = eilmer_Quadratic().Recon(limiter = limiter, qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
    
    elif reconstruction == "eilmer_cubic":
        Lft, Rght = eilmer_Cubic().Recon(limiter = limiter, qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
               
    elif reconstruction == "MUSCL":
        Lft, Rght = MUSCL().Recon(limiter = limiter, qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
               
    elif reconstruction == "THINC":
        Lft, Rght = THINC().Recon(limiter = limiter, qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
                      
    elif reconstruction == "MUSCLTHINCBVD":
        Lft, Rght = MUSCLTHINCBVD().Recon(limiter = limiter, qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
                                                 
    return Lft, Rght