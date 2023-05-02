"""
Function:
Author: Luke Bartholomew
Edits:
"""
import numpy as np
def multi_species_decode_to_primative_properties(cqs):
    rho_tol = 0.1

    rho = cqs["mass"]
    d_inv = 1.0 / rho

    vel_x = cqs["xMom"] * d_inv

    rU = cqs["energy"]
    u = rU * d_inv
    u -= 0.5 * vel_x ** 2.0

    for ind, val in enumerate(cqs["spcs_mass"]):
        if val < 0.0:
            cqs["spcs_mass"][ind] = 0.0
    
    rho_sum = sum(cqs["spcs_mass"])

    if abs(rho - rho_sum) > rho_tol:
        print("Conserved quantity: rho too different from sum of species mass.")
    
    scale = 1.0
    if abs(rho - rho_sum) > 0.0:
        scale = rho / rho_sum
        for ind, val in enumerate(cqs["spcs_mass"]):
            cqs["spcs_mass"][ind] *= scale

    massf = (np.array(cqs["spcs_mass"]) * d_inv).tolist()
    
    return rho, u, massf, vel_x