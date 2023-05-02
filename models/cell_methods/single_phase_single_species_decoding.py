"""
Function:
Author: Luke Bartholomew
Edits:
"""

def single_species_decode_to_primative_properties(cqs):
    rho = cqs["mass"]
    xMom = cqs["xMom"]
    rhoU = cqs["energy"]
    vel_x = xMom / rho
    u = rhoU / rho - 0.5 * vel_x ** 2.0

    return rho, u, vel_x