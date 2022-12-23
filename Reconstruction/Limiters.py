import math as m

class eilmer_vanAlbada():
    def __init__(self) -> None:
        self.epsilon = 1e-6

    def limiter(self, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil[:2]
        [len_L0, len_L1] = dxL_stencil[:2]
        [q_R0, q_R1] = qR_stencil[:2]
        [len_R0, len_R1] = dxR_stencil[:2]
        
        w_1 = 2.0 / (len_L0 + len_L1)
        w_2 = 2.0 / (len_R0 + len_L0)
        w_3 = 2.0 / (len_R1 + len_R0)
        van_DelLminus = (q_L0 - q_L1) * w_1
        van_Del = (q_R0 - q_L0) * w_2
        van_DelRplus = (q_R1 - q_R0) * w_3
        sL = (van_DelLminus * van_Del + abs(van_DelLminus * van_Del) + self.epsilon) / (van_DelLminus * van_DelLminus + van_Del * van_Del + self.epsilon)
        sR = (van_Del * van_DelRplus + abs(van_Del * van_DelRplus) + self.epsilon) / (van_Del * van_Del + van_DelRplus * van_DelRplus + self.epsilon)
        return sL, sR

class nonuniform_vanAlbada():
    def __init__(self) -> None:
        self.epsilon = 1e-6
        self.k = 2.0
        pass

    def limiter(self, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil[:2]
        [len_L0, len_L1] = dxL_stencil[:2]
        [q_R0, q_R1] = qR_stencil[:2]
        [len_R0, len_R1] = dxR_stencil[:2]
        
        A_L = (len_L1 + len_L0) / (len_L0 + len_R0)
        B_L = (2.0 * len_L0) / (len_L0 + len_R0)
        if q_L0 == q_R0:
            sL, sR = 0.0, 0.0
            return sL, sR
        else:
            theta_L = (q_L0 - q_L1) / (q_R0 - q_L0)
            A_R = (len_R0 + len_R1) / (len_L0 + len_R0)
            B_R = (2.0 * len_R0) / (len_L0 + len_R0)
            theta_R = (q_R1 - q_R0) / (q_R0 - q_L0)
            if theta_L <  0.0:
                sL = 0.0
            else:
                sL = (B_L * (theta_L ** self.k + theta_L)) / (theta_L ** self.k + A_L)
            if theta_R < 0.0:
                sR = 0.0
            else:
                sR = (B_R * (theta_R ** self.k + theta_R)) / (theta_R ** self.k + A_R)
            return sL, sR

class nonuniform_vanLeer():
    def __init__(self) -> None:
        self.k = 3.0
        pass

    def limiter(self, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil[:2]
        [len_L0, len_L1] = dxL_stencil[:2]
        [q_R0, q_R1] = qR_stencil[:2]
        [len_R0, len_R1] = dxR_stencil[:2]
        
        k = self.k
        if q_L0 == q_R0:
            sL, sR = 0.0, 0.0
            return sL, sR
        else:
            A_L = (len_L1 + len_L0) / (len_L0 + len_R0)
            B_L = (2.0 * len_L0) / (len_L0 + len_R0)
            theta_L = (q_L0 - q_L1) / (q_R0 - q_L0)
            A_R = (len_R0 + len_R1) / (len_L0 + len_R0)
            B_R = (2.0 * len_R0) / (len_L0 + len_R0)
            theta_R = (q_R1 - q_R0) / (q_R0 - q_L0)
              
            if theta_L < 0.0:
                sL = 0.0
            else:
                theta_L_num, A_L_denom = 0.0, 0.0
                for i in range(1, k+1):
                    theta_L_num += theta_L ** i
                    A_L_denom += A_L ** i
                sL = B_L * theta_L_num * (A_L_denom + 1.0)/ ((theta_L_num + 1.0) * A_L_denom)
            
            if theta_R < 0.0:
                sR = 0.0
            else:
                theta_R_num, A_R_denom = 0.0, 0.0
                for i in range(1, k+1):
                    theta_R_num += theta_R ** i
                    A_R_denom += A_R ** i
                sR = B_R * theta_R_num * (A_R_denom + 1.0) / ((theta_R_num + 1.0) * A_R_denom)
            
            return sL, sR

class nonuniform_minmod():
    def __init__(self) -> None:
        pass

    def limiter(self, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil[:2]
        [len_L0, len_L1] = dxL_stencil[:2]
        [q_R0, q_R1] = qR_stencil[:2]
        [len_R0, len_R1] = dxR_stencil[:2]
        
        
        if q_L0 == q_R0:
            sL, sR = 0.0, 0.0
            return sL, sR
        else:
            A_L = (len_L1 + len_L0) / (len_L0 + len_R0)
            B_L = (2.0 * len_L0) / (len_L0 + len_R0)
            theta_L = (q_L0 - q_L1) / (q_R0 - q_L0)
            A_R = (len_R0 + len_R1) / (len_L0 + len_R0)
            B_R = (2.0 * len_R0) / (len_L0 + len_R0)
            theta_R = (q_R1 - q_R0) / (q_R0 - q_L0)
            sL = max(0.0, B_L / A_L * min(theta_L, A_L))
            sR = max(0.0, B_R / A_R * min(theta_R, A_R))
            return sL, sR

class nonuniform_superbee():
    def __init__(self) -> None:
        pass

    def limiter(self, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil[:2]
        [len_L0, len_L1] = dxL_stencil[:2]
        [q_R0, q_R1] = qR_stencil[:2]
        [len_R0, len_R1] = dxR_stencil[:2]
        
        if q_L0 == q_R0:
            sL, sR = 0.0, 0.0
            return sL, sR
        else:
            A_L = (len_L1 + len_L0) / (len_L0 + len_R0)
            B_L = (2.0 * len_L0) / (len_L0 + len_R0)
            theta_L = (q_L0 - q_L1) / (q_R0 - q_L0)
            A_R = (len_R0 + len_R1) / (len_L0 + len_R0)
            B_R = (2.0 * len_R0) / (len_L0 + len_R0)
            theta_R = (q_R1 - q_R0) / (q_R0 - q_L0)
            sL = max(0.0, min(2.0 * theta_L, B_L), min(B_L * theta_L / A_L, 2.0))
            sR = max(0.0, min(2.0 * theta_R, B_R), min(B_R * theta_R / A_R, 2.0))
            return sL, sR

class nonuniform_MC():
    def __init__(self) -> None:
        pass

    def limiter(self, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil[:2]
        [len_L0, len_L1] = dxL_stencil[:2]
        [q_R0, q_R1] = qR_stencil[:2]
        [len_R0, len_R1] = dxR_stencil[:2]
        
        if q_L0 == q_R0:
            sL, sR = 0.0, 0.0
            return sL, sR
        else:
            A_L = (len_L1 + len_L0) / (len_L0 + len_R0)
            B_L = (2.0 * len_L0) / (len_L0 + len_R0)
            theta_L = (q_L0 - q_L1) / (q_R0 - q_L0)
            A_R = (len_R0 + len_R1) / (len_L0 + len_R0)
            B_R = (2.0 * len_R0) / (len_L0 + len_R0)
            theta_R = (q_R1 - q_R0) / (q_R0 - q_L0)
            sL = max(0.0, min(2.0 * theta_L, B_L * (theta_L + 1.0) / (A_L + 1.0), 2.0))
            sR = max(0.0, min(2.0 * theta_R, B_R * (theta_R + 1.0) / (A_R + 1.0), 2.0))
            return sL, sR

class uniform_Koren():
    def __init__(self) -> None:
        pass

    def limiter(self, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil[:2]
        [len_L0, len_L1] = dxL_stencil[:2]
        [q_R0, q_R1] = qR_stencil[:2]
        [len_R0, len_R1] = dxR_stencil[:2]
        
        if abs(q_R0 - q_L0) < 1e-9:
            sL, sR = 0.0, 0.0
        else:
            theta_L = (q_L0 - q_L1) / (q_R0 - q_L0)
            theta_R = (q_R1 - q_R0) / (q_R0 - q_L0)
            if theta_L < 0.0:
                sL = 0.0
            else:
                sL = max(0.0, min(2.0 * theta_L, min((1.0 + 2.0 * theta_L) / 3.0, 2.0)))
            if theta_R < 0.0:
                sR = 0.0
            else:
                sR = max(0.0, min(2.0 * theta_R, min((1.0 + 2.0 * theta_R) / 3.0, 2.0)))
        return sL, sR
        

class uniform_minmod():
    def __init__(self) -> None:
        pass

    def limiter(self, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil[:2]
        [len_L0, len_L1] = dxL_stencil[:2]
        [q_R0, q_R1] = qR_stencil[:2]
        [len_R0, len_R1] = dxR_stencil[:2]
        
        if abs(q_R0 - q_L0) < 1e-9:
            sL, sR = 0.0, 0.0
        else:
            theta_L = (q_L0 - q_L1) / (q_R0 - q_L0)
            theta_R = (q_R1 - q_R0) / (q_R0 - q_L0)
            if theta_L < 0.0:
                sL = 0.0
            else:
                sL = max(0.0, min(1.0, theta_L))
            if theta_R < 0.0:
                sR = 0.0
            else:
                sR = max(0.0, min(1.0, theta_R))
        return sL, sR

class uniform_MC():
    def __init__(self) -> None:
        pass

    def limiter(self, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil[:2]
        [len_L0, len_L1] = dxL_stencil[:2]
        [q_R0, q_R1] = qR_stencil[:2]
        [len_R0, len_R1] = dxR_stencil[:2]
        
        if abs(q_R0 - q_L0) < 1e-9:
            sL, sR = 0.0, 0.0
        else:
            theta_L = (q_L0 - q_L1) / (q_R0 - q_L0)
            theta_R = (q_R1 - q_R0) / (q_R0 - q_L0)
            if theta_L < 0.0:
                sL = 0.0
            else:
                sL = max(0.0, min(2.0 * theta_L, 0.5 * (1.0 + theta_L), 2.0))
            if theta_R < 0.0:
                sR = 0.0
            else:
                sR = max(0.0, min(2.0 * theta_R, 0.5 * (1.0 + theta_R), 2.0))
        return sL, sR

class uniform_Osher():
    def __init__(self) -> None:
        self.beta = 1.5 # any number between 1 and 2
        pass

    def limiter(self, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil[:2]
        [len_L0, len_L1] = dxL_stencil[:2]
        [q_R0, q_R1] = qR_stencil[:2]
        [len_R0, len_R1] = dxR_stencil[:2]
        beta = self.beta
        
        if abs(q_R0 - q_L0) < 1e-9:
            sL, sR = 0.0, 0.0
        else:
            theta_L = (q_L0 - q_L1) / (q_R0 - q_L0)
            theta_R = (q_R1 - q_R0) / (q_R0 - q_L0)
            if theta_L < 0.0:
                sL = 0.0
            else:
                sL = max(0.0, min(theta_L, beta))
            if theta_R < 0.0:
                sR = 0.0
            else:
                sR = max(0.0, min(theta_R, beta))
        return sL, sR

class uniform_ospre():
    def __init__(self) -> None:
        pass

    def limiter(self, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil[:2]
        [len_L0, len_L1] = dxL_stencil[:2]
        [q_R0, q_R1] = qR_stencil[:2]
        [len_R0, len_R1] = dxR_stencil[:2]
        
        if abs(q_R0 - q_L0) < 1e-9:
            sL, sR = 0.0, 0.0
        else:
            theta_L = (q_L0 - q_L1) / (q_R0 - q_L0)
            theta_R = (q_R1 - q_R0) / (q_R0 - q_L0)
            if theta_L < 0.0:
                sL = 0.0
            else:
                sL = (1.5 * (theta_L ** 2.0 + theta_L)) / (theta_L ** 2.0 + theta_L + 1.0)
            if theta_R < 0.0:
                sR = 0.0
            else:
                sR = (1.5 * (theta_R ** 2.0 + theta_R)) / (theta_R ** 2.0 + theta_R + 1.0)
        return sL, sR

class uniform_superbee():
    def __init__(self) -> None:
        pass

    def limiter(self, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil[:2]
        [len_L0, len_L1] = dxL_stencil[:2]
        [q_R0, q_R1] = qR_stencil[:2]
        [len_R0, len_R1] = dxR_stencil[:2]
        
        
        if abs(q_R0 - q_L0) < 1e-9:
            sL, sR = 0.0, 0.0
        else:
            theta_L = (q_L0 - q_L1) / (q_R0 - q_L0)
            theta_R = (q_R1 - q_R0) / (q_R0 - q_L0)
            if theta_L < 0.0:
                sL = 0.0
            else:
                sL = max(0.0, min(2.0 * theta_L, 1.0), min(theta_L, 2.0))
            if theta_R < 0.0:
                sR = 0.0
            else:
                sR = max(0.0, min(2.0 * theta_R, 1.0), min(theta_R, 2.0))
        return sL, sR

class uniform_Sweby():
    def __init__(self) -> None:
        self.beta = 1.5 # any number between 1 and 2
        pass

    def limiter(self, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil[:2]
        [len_L0, len_L1] = dxL_stencil[:2]
        [q_R0, q_R1] = qR_stencil[:2]
        [len_R0, len_R1] = dxR_stencil[:2]
        beta = self.beta
        
        
        if abs(q_R0 - q_L0) < 1e-9:
            sL, sR = 0.0, 0.0
        else:
            theta_L = (q_L0 - q_L1) / (q_R0 - q_L0)
            theta_R = (q_R1 - q_R0) / (q_R0 - q_L0)
            if theta_L < 0.0:
                sL = 0.0
            else:
                sL = max(0.0, min(beta * theta_L, 1.0), min(theta_L, beta))
            if theta_R < 0.0:
                sR = 0.0
            else:
                sR = max(0.0, min(beta * theta_R, 1.0), min(theta_R, beta))
        return sL, sR

class uniform_UMIST():
    def __init__(self) -> None:
        pass

    def limiter(self, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil[:2]
        [len_L0, len_L1] = dxL_stencil[:2]
        [q_R0, q_R1] = qR_stencil[:2]
        [len_R0, len_R1] = dxR_stencil[:2]
        
        
        if abs(q_R0 - q_L0) < 1e-9:
            sL, sR = 0.0, 0.0
        else:
            theta_L = (q_L0 - q_L1) / (q_R0 - q_L0)
            theta_R = (q_R1 - q_R0) / (q_R0 - q_L0)
            if theta_L < 0.0:
                sL = 0.0
            else:
                sL = max(0.0, min(2.0 * theta_L, 0.25 + 0.75 * theta_L, 0.75 + 0.25 * theta_L, 2.0))
            if theta_R < 0.0:
                sR = 0.0
            else:
                sR = max(0.0, min(2.0 * theta_R, 0.25 + 0.75 * theta_R, 0.75 + 0.25 * theta_R, 2.0))
        return sL, sR

class uniform_vanAlbada1():
    def __init__(self) -> None:
        pass

    def limiter(self, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil[:2]
        [len_L0, len_L1] = dxL_stencil[:2]
        [q_R0, q_R1] = qR_stencil[:2]
        [len_R0, len_R1] = dxR_stencil[:2]
        
        if abs(q_R0 - q_L0) < 1e-9:
            sL, sR = 0.0, 0.0
        else:
            theta_L = (q_L0 - q_L1) / (q_R0 - q_L0)
            theta_R = (q_R1 - q_R0) / (q_R0 - q_L0)
            if theta_L < 0.0:
                sL = 0.0
            else:
                sL = (theta_L ** 2.0 + theta_L) / (theta_L ** 2.0 + 1.0)
            if theta_R < 0.0:
                sR = 0.0
            else:
                sR = (theta_R ** 2.0 + theta_R) / (theta_R ** 2.0 + 1.0)
        return sL, sR

class uniform_vanAlbada2():
    def __init__(self) -> None:
        pass

    def limiter(self, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil[:2]
        [len_L0, len_L1] = dxL_stencil[:2]
        [q_R0, q_R1] = qR_stencil[:2]
        [len_R0, len_R1] = dxR_stencil[:2]
        
        
        if abs(q_R0 - q_L0) < 1e-9:
            sL, sR = 0.0, 0.0
        else:
            theta_L = (q_L0 - q_L1) / (q_R0 - q_L0)
            theta_R = (q_R1 - q_R0) / (q_R0 - q_L0)
            if theta_L < 0.0:
                sL = 0.0
            else:
                sL = (2.0 * theta_L) / (theta_L ** 2.0 + 1.0)
            if theta_R < 0.0:
                sR = 0.0
            else:
                sR = (2.0 * theta_R) / (theta_R ** 2.0 + 1.0)
        return sL, sR

class uniform_vanLeer():
    def __init__(self) -> None:
        pass

    def limiter(self, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil[:2]
        [len_L0, len_L1] = dxL_stencil[:2]
        [q_R0, q_R1] = qR_stencil[:2]
        [len_R0, len_R1] = dxR_stencil[:2]
        
        if abs(q_R0 - q_L0) < 1e-9:
            sL, sR = 0.0, 0.0
        else:
            theta_L = (q_L0 - q_L1) / (q_R0 - q_L0)
            theta_R = (q_R1 - q_R0) / (q_R0 - q_L0)
            if theta_L < 0.0:
                sL = 0.0
            else:
                sL = (theta_L + abs(theta_L)) / (1.0 + abs(theta_L))
            if theta_R < 0.0:
                sR = 0.0
            else:
                sR = (theta_R + abs(theta_R)) / (1.0 + abs(theta_R))
        return sL, sR

class uniform_Gregori():
    def __init__(self) -> None:
        pass

    def limiter(self, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
        [q_L0, q_L1] = qL_stencil[:2]
        [len_L0, len_L1] = dxL_stencil[:2]
        [q_R0, q_R1] = qR_stencil[:2]
        [len_R0, len_R1] = dxR_stencil[:2]

        if abs(q_R0 - q_L0) < 1e-9:
            sL, sR = 0.0, 0.0
        else:
            theta_L = (q_L0 - q_L1) / (q_R0 - q_L0)
            theta_R = (q_R1 - q_R0) / (q_R0 - q_L0)
            if theta_L < 0.0:
                sL = 0.0
            else:
                sL = max(0.0, min(2.0 * theta_L, theta_L ** 0.5, 2.0))
            if theta_R < 0.0:
                sR = 0.0
            else:
                sR = max(0.0, min(2.0 * theta_R, theta_R ** 0.5, 2.0))
        return sL, sR


class example():
    def __init__(self) -> None:
        pass

    def limiter(self):
        pass

def getLimiters(limiter, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    if limiter == "eilmer_vanAlbada":
        sL, sR = eilmer_vanAlbada().limiter(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "nonuniform_vanAlbada":
        sL, sR = nonuniform_vanAlbada().limiter(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "nonuniform_vanLeer":
        sL, sR = nonuniform_vanLeer().limiter(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "nonuniform_minmod":
        sL, sR = nonuniform_minmod().limiter(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "nonuniform_superbee":
        sL, sR = nonuniform_superbee().limiter(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "nonuniform_MC":
        sL, sR = nonuniform_MC().limiter(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_Koren":
        sL, sR = uniform_Koren().limiter(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_minmod":
        sL, sR = uniform_minmod().limiter(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_MC":
        sL, sR = uniform_MC().limiter(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_Osher":
        sL, sR = uniform_Osher().limiter(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_ospre":
        sL, sR = uniform_ospre().limiter(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_superbee":
        sL, sR = uniform_superbee().limiter(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_Sweby":
        sL, sR = uniform_Sweby().limiter(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
    
    elif limiter == "uniform_UMIST":
        sL, sR = uniform_UMIST().limiter(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_vanAlbada1":
        sL, sR = uniform_vanAlbada1().limiter(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_vanAlbada2":
        sL, sR = uniform_vanAlbada2().limiter(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_vanLeer":
        sL, sR = uniform_vanLeer().limiter(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_Gregori":
        sL, sR = uniform_Gregori().limiter(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    return sL, sR