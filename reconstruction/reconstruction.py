"""
Function:
Author: Luke Bartholomew
Edits:
"""
import math as m
from Algorithms.DT_1D_V4.reconstruction.limiters import get_limiter

def recon_eilmer_L3R3(limiter, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
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
    if limiter is not None:

        sL = get_limiter(limiter = limiter, side = "left", \
                            qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                            qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
        sR = get_limiter(limiter = limiter, side = "right", \
                                qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    lft = q_L0 + sL * (w_0 * q_L2 + w_1 * q_L1 + (w_2 - 1.0) * q_L0 + w_3 * q_R0 + w_4 * q_R1)
    rght = q_R0 + sR * (w_5 * q_L1 + w_6 * q_L0 + (w_7 - 1.0) * q_R0 + w_8 * q_R1 + w_9 * q_R2)

    ### Apply clipping
    lft, rght = clip_to_bounds(q_L0 = q_L0, q_R0 = q_R0, lft = lft, rght = rght)
    return lft, rght

def recon_eilmer_L2R2(limiter, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
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
    if limiter is not None:
        sL = get_limiter(limiter = limiter, side = "left", \
                                qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
        sR = get_limiter(limiter = limiter, side = "right", \
                                qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
        
    lft = q_L0 + sL * w_2 * (Del * w_4 + DelLminus * w_1)
    rght = q_R0 - sR * w_3 * (DelRplus * w_0 + Del * w_5)
    ### Apply clipping
    lft, rght = clip_to_bounds(q_L0 = q_L0, q_R0 = q_R0, lft = lft, rght = rght)
    return lft, rght

def recon_eilmer_L2R1(limiter, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    [q_L0, q_L1] = qL_stencil
    [len_L0, len_L1] = dxL_stencil
    [q_R0] = qR_stencil
    [len_R0] = dxR_stencil

    a_L0 = 0.5 * len_L0 / (len_L1 + 2.0 * len_L0 + len_R0)

    delta = (q_R0 - q_L0) * 2.0 / (len_L0 + len_R0)
    delta_L_minus = (q_L0 - q_L1) * 2.0 / (len_L1 + len_L0)

    ### Apply limiter
    sL = 1.0
    if limiter is not None:
        sL = get_limiter(limiter = limiter, side = "left", \
                                qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    lft = q_L0 + sL * a_L0 * (delta * (2.0 * len_L0 + len_L1) + delta_L_minus * len_R0)

    if delta * delta_L_minus >= 0.0:
        rght = q_L0 * len_R0 / (len_L0 + len_R0)
    
    else:
        rght = q_R0
    ### Apply clipping
    lft, rght = clip_to_bounds(q_L0 = q_L0, q_R0 = q_R0, lft = lft, rght = rght)
    
    return lft, rght

def recon_eilmer_L1R2(limiter, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    [q_L0] = qL_stencil
    [len_L0] = dxL_stencil
    [q_R0, q_R1] = qR_stencil
    [len_R0, len_R1] = dxR_stencil

    delta = (q_R0 - q_L0) * 2.0 / (len_L0 + len_R0)
    delta_R_plus = (q_R1 - q_R0) * 2.0 / (len_R0 + len_R1)

    a_R0 = 0.5 * len_R0 / (len_L0 + 2.0 * len_R0 + len_R1)

    ### Apply limiter
    sR = 1.0
    if limiter is not None:
        sR = get_limiter(limiter = limiter, side = "right", \
                                qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
    rght = q_R0 - sR * a_R0 * (delta_R_plus * len_L0 + delta * (2.0 * len_R0 + len_R1))
    if delta * delta_R_plus >= 0.0:
        lft = q_L0 * len_R0 / (len_L0 + len_R0) + q_R0 * len_L0 / (len_L0 + len_R0)

    else:
        lft = q_L0

    ### Apply clipping
    lft, rght = clip_to_bounds(q_L0 = q_L0, q_R0 = q_R0, lft = lft, rght = rght)
    
    return lft, rght

def recon_eilmer_L1R1(qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    [q_L0] = qL_stencil
    [q_R0] = qR_stencil
        
    lft = q_L0
    rght = q_R0
    return lft, rght

def clip_to_bounds(q_L0, q_R0, lft, rght):
        q_min = min(q_L0, q_R0)
        q_max = max(q_L0, q_R0)
        if lft < q_min:
            lft = q_min
        elif lft > q_max:
            lft = q_max
        if rght < q_min:
            rght = q_min
        elif rght > q_max:
            rght = q_max
        return lft, rght

def get_reconstruction(reconstruction, limiter, qL_stencil, dxL_stencil, qR_stencil, dxR_stencil):
    if reconstruction == "eilmer_L1R1":
        lft, rght = recon_eilmer_L1R1(qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif reconstruction == "eilmer_L2R2":
        lft, rght = recon_eilmer_L2R2(limiter = limiter, qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
    
    elif reconstruction == "eilmer_L3R3":
        lft, rght = recon_eilmer_L3R3(limiter = limiter, qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)

    elif reconstruction == "eilmer_L2R1":
        lft, rght = recon_eilmer_L2R1(limiter = limiter, qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
        
    elif reconstruction == "eilmer_L1R2":
        lft, rght = recon_eilmer_L1R2(limiter = limiter, qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
                                         
    return lft, rght