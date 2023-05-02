import numpy as np
def encode_multi_species_cqs(flow_state):
    rho = flow_state.fluid_state.rho
    u = flow_state.fluid_state.u
    massf = flow_state.fluid_state.massf
    vel_x = flow_state.vel_x
    cqs = {
        "mass"      : rho,
        "xMom"      : rho * vel_x,
        "energy"    : rho * (u + 0.5 * vel_x ** 2.0),
        "spcs_mass" : (rho * np.array(massf)).tolist()
    }
    return cqs
