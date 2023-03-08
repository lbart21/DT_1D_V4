"""
Function:
Author: Luke Bartholomew
Edits:
"""
import math as m

def eilmer_vanAlbada_L3R3_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    epsilon = 1e-6
    [q_L0, q_L1] = qL_stencil[:2]
    [q_R0, q_R1] = qR_stencil[:2]

    delta = q_R0 - q_L0

    if side == "left":
        delta_L_minus = q_L0 - q_L1
        s = (delta_L_minus * delta + abs(delta_L_minus * delta) + epsilon) \
                / (delta_L_minus * delta_L_minus + delta * delta + epsilon)
    
    elif side == "right":
        delta_R_plus = q_R1 - q_R0
        s = (delta * delta_R_plus + abs(delta * delta_R_plus) + epsilon) \
                / (delta * delta + delta_R_plus * delta_R_plus + epsilon)
    
    return s

def eilmer_vanAlbada_L2R2_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    epsilon = 1e-6
    
    if side == "left":
        [q_L0, q_L1] = qL_stencil[:2]
        [len_L0, len_L1] = dxL_stencil[:2]

        [q_R0] = qR_stencil[:1]
        [len_R0] = dxR_stencil[:1]

        w_1 = 2.0 / (len_L0 + len_L1)
        w_2 = 2.0 / (len_R0 + len_L0)
        delta = (q_R0 - q_L0) * w_2
        delta_L_minus = (q_L0 - q_L1) * w_1
        s = (delta_L_minus * delta + abs(delta_L_minus * delta) + epsilon) \
                / (delta_L_minus * delta_L_minus + delta * delta + epsilon)
    
    elif side == "right":
        [q_R0, q_R1] = qR_stencil[:2]
        [len_R0, len_R1] = dxR_stencil[:2]

        [q_L0] = qL_stencil[:1]
        [len_L0] = dxL_stencil[:1]

        w_2 = 2.0 / (len_R0 + len_L0)
        w_3 = 2.0 / (len_R1 + len_R0)
        delta = (q_R0 - q_L0) * w_2
        delta_R_plus = (q_R1 - q_R0) * w_3
        s = (delta * delta_R_plus + abs(delta * delta_R_plus) + epsilon) \
                / (delta * delta + delta_R_plus * delta_R_plus + epsilon)
    
    return s

def nonuniform_vanAlbada_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    print("This limit is due to be redone.")
    # TODO Redo this function
    epsilon = 1e-6
    k = 2.0
    [q_L0, q_L1] = qL_stencil[:2]
    [len_L0, len_L1] = dxL_stencil[:2]
    [q_R0, q_R1] = qR_stencil[:2]
    [len_R0, len_R1] = dxR_stencil[:2]
        
    A_L = (len_L1 + len_L0) / (len_L0 + len_R0)
    B_L = (2.0 * len_L0) / (len_L0 + len_R0)
    return 1.0

def nonuniform_vanLeer_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    print("This limit is due to be redone.")
    # TODO Redo this function
    [q_L0, q_L1] = qL_stencil[:2]
    [len_L0, len_L1] = dxL_stencil[:2]
    [q_R0, q_R1] = qR_stencil[:2]
    [len_R0, len_R1] = dxR_stencil[:2]
        
    k = 3.0
    return 1.0 

def nonuniform_minmod_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    print("This limit is due to be redone.")
    # TODO Redo this function
    [q_L0, q_L1] = qL_stencil[:2]
    [len_L0, len_L1] = dxL_stencil[:2]
    [q_R0, q_R1] = qR_stencil[:2]
    [len_R0, len_R1] = dxR_stencil[:2]
        
    return 1.0

def nonuniform_superbee_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    print("This limit is due to be redone.")
    # TODO Redo this function
    [q_L0, q_L1] = qL_stencil[:2]
    [len_L0, len_L1] = dxL_stencil[:2]
    [q_R0, q_R1] = qR_stencil[:2]
    [len_R0, len_R1] = dxR_stencil[:2]
        
    return 1.0

def nonuniform_MC_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    print("This limit is due to be redone.")
    # TODO Redo this function
    [q_L0, q_L1] = qL_stencil[:2]
    [len_L0, len_L1] = dxL_stencil[:2]
    [q_R0, q_R1] = qR_stencil[:2]
    [len_R0, len_R1] = dxR_stencil[:2]

    return 1.0

def uniform_Koren_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    print("This limit is due to be redone.")
    # TODO Redo this function
    [q_L0, q_L1] = qL_stencil[:2]
    [len_L0, len_L1] = dxL_stencil[:2]
    [q_R0, q_R1] = qR_stencil[:2]
    [len_R0, len_R1] = dxR_stencil[:2]
        
    return 1.0

def uniform_minmod_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    print("This limit is due to be redone.")
    # TODO Redo this function
    [q_L0, q_L1] = qL_stencil[:2]
    [len_L0, len_L1] = dxL_stencil[:2]
    [q_R0, q_R1] = qR_stencil[:2]
    [len_R0, len_R1] = dxR_stencil[:2]
        
    return 1.0
        
def uniform_MC_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    print("This limit is due to be redone.")
    # TODO Redo this function
    [q_L0, q_L1] = qL_stencil[:2]
    [len_L0, len_L1] = dxL_stencil[:2]
    [q_R0, q_R1] = qR_stencil[:2]
    [len_R0, len_R1] = dxR_stencil[:2]
        
    return 1.0
        
def uniform_Osher_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    print("This limit is due to be redone.")
    # TODO Redo this function
    [q_L0, q_L1] = qL_stencil[:2]
    [len_L0, len_L1] = dxL_stencil[:2]
    [q_R0, q_R1] = qR_stencil[:2]
    [len_R0, len_R1] = dxR_stencil[:2]
    beta = 1.5
        
    return 1.0

def uniform_ospre_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    print("This limit is due to be redone.")
    # TODO Redo this function
    [q_L0, q_L1] = qL_stencil[:2]
    [len_L0, len_L1] = dxL_stencil[:2]
    [q_R0, q_R1] = qR_stencil[:2]
    [len_R0, len_R1] = dxR_stencil[:2]
        
    return 1.0

def uniform_superbee_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    print("This limit is due to be redone.")
    # TODO Redo this function
    [q_L0, q_L1] = qL_stencil[:2]
    [len_L0, len_L1] = dxL_stencil[:2]
    [q_R0, q_R1] = qR_stencil[:2]
    [len_R0, len_R1] = dxR_stencil[:2]
        
    return 1.0

def uniform_Sweby_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    print("This limit is due to be redone.")
    # TODO Redo this function
    [q_L0, q_L1] = qL_stencil[:2]
    [len_L0, len_L1] = dxL_stencil[:2]
    [q_R0, q_R1] = qR_stencil[:2]
    [len_R0, len_R1] = dxR_stencil[:2]
    beta = 1.5
        
    return 1.0

def uniform_UMIST_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    print("This limit is due to be redone.")
    # TODO Redo this function
    [q_L0, q_L1] = qL_stencil[:2]
    [len_L0, len_L1] = dxL_stencil[:2]
    [q_R0, q_R1] = qR_stencil[:2]
    [len_R0, len_R1] = dxR_stencil[:2]
        
    return 1.0

def uniform_vanAlbada1_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    print("This limit is due to be redone.")
    # TODO Redo this function
    [q_L0, q_L1] = qL_stencil[:2]
    [len_L0, len_L1] = dxL_stencil[:2]
    [q_R0, q_R1] = qR_stencil[:2]
    [len_R0, len_R1] = dxR_stencil[:2]
        
    return 1.0

def uniform_vanAlbada2_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    [q_L0, q_L1] = qL_stencil[:2]
    [len_L0, len_L1] = dxL_stencil[:2]
    [q_R0, q_R1] = qR_stencil[:2]
    [len_R0, len_R1] = dxR_stencil[:2]
    
    return 1.0
        
def uniform_vanLeer_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    print("This limit is due to be redone.")
    # TODO Redo this function
    [q_L0, q_L1] = qL_stencil[:2]
    [len_L0, len_L1] = dxL_stencil[:2]
    [q_R0, q_R1] = qR_stencil[:2]
    [len_R0, len_R1] = dxR_stencil[:2]
        
    return 1.0

def uniform_Gregori_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    print("This limit is due to be redone.")
    # TODO Redo this function
    [q_L0, q_L1] = qL_stencil[:2]
    [len_L0, len_L1] = dxL_stencil[:2]
    [q_R0, q_R1] = qR_stencil[:2]
    [len_R0, len_R1] = dxR_stencil[:2]

    return 1.0

def example_limiter(side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    pass

def get_limiter(limiter, side, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    if limiter == "eilmer_vanAlbada_L3R3":
        s = eilmer_vanAlbada_L3R3_limiter(side = side, \
                                            qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
    
    if limiter == "eilmer_vanAlbada_L2R2":
        s = eilmer_vanAlbada_L2R2_limiter(side = side, \
                                            qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "nonuniform_vanAlbada":
        s = nonuniform_vanAlbada_limiter(side = side, \
                                            qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "nonuniform_vanLeer":
        s = nonuniform_vanLeer_limiter(side = side, \
                                        qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "nonuniform_minmod":
        s = nonuniform_minmod_limiter(side = side, \
                                        qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "nonuniform_superbee":
        s = nonuniform_superbee_limiter(side = side, \
                                        qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "nonuniform_MC":
        s = nonuniform_MC_limiter(side = side, \
                                    qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                    qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_Koren":
        s = uniform_Koren_limiter(side = side, \
                                    qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                    qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_minmod":
        s = uniform_minmod_limiter(side = side, \
                                        qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_MC":
        s = uniform_MC_limiter(side = side, \
                                    qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                    qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_Osher":
        s = uniform_Osher_limiter(side = side, \
                                        qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_ospre":
        s = uniform_ospre_limiter(side = side, \
                                        qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_superbee":
        s = uniform_superbee_limiter(side = side, \
                                            qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_Sweby":
        s = uniform_Sweby_limiter(side = side, \
                                        qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
    
    elif limiter == "uniform_UMIST":
        s = uniform_UMIST_limiter(side = side, \
                                        qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_vanAlbada1":
        s = uniform_vanAlbada1_limiter(side = side, \
                                            qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_vanAlbada2":
        s = uniform_vanAlbada2_limiter(side = side, \
                                            qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_vanLeer":
        s = uniform_vanLeer_limiter(side = side, \
                                        qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif limiter == "uniform_Gregori":
        s = uniform_Gregori_limiter(side = side, \
                                        qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    return s