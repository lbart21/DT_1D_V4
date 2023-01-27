"""
Function:
Author: Luke Bartholomew
Edits:
"""
from copy import deepcopy

def FromStagnationWithMassFlowRateInFlow_BC(mesh, BC, onWestBoundaryBool):
    """
    [p_stag, T_stag] = BC[1][2]
    Take p_stag and T_stag to create a stagnation state. Use velocity on inside of 
    boundary and stagnation enthalpy to find static enthalpy. Update properties from 
    stagnation entropy and static enthalpy. All cells have the interior cell's velocity.
    Velocity is set to 0 if flow is out of boundary.
    Then update the flow properties to complete the state.
    """

    [stag_gas_state, mass_flux] = BC[1][2]
    p_stag = stag_gas_state.p
    T_stag = stag_gas_state.T
    relaxation_factor = 0.1

    p0_min = 0.1 * p_stag
    p0_max = 10.0 * p_stag
    vel_x_boundary = mesh.cellArray[0].fs.vel_x
    p_boundary = mesh.cellArray[0].fs.fluid_state.p
    rho_boundary = mesh.cellArray[0].fs.fluid_state.rho
    rhoU_boundary = rho_boundary * vel_x_boundary
    A_boundary = mesh.interfaceArray[0].GEO["A"]
    dp_over_p = 0.5 * relaxation_factor / rho_boundary * ( (mass_flux / A_boundary) ** 2.0 \
                                                            - rhoU_boundary * abs(rhoU_boundary) ) / p_boundary
    new_p0 = (1.0 + dp_over_p) * p_stag
    new_p0 = min(max(new_p0, p0_min), p0_max)
    flowState = deepcopy(mesh.cellArray[0].fs)
    flowState.fluid_state.p = new_p0
    flowState.fluid_state.T = T_stag
    flowState.fluid_state.update_thermo_from_pT()
    stagnation_enthalpy = flowState.fluid_state.enthalpy
    
    if onWestBoundaryBool:
        bulk_speed = max(0.0, vel_x_boundary) # 0.0 if vel_x_boundary negative, vel_x_boundary if positive, 
        
    else:
        bulk_speed = min(0.0, vel_x_boundary)

    static_enthalpy = stagnation_enthalpy - 0.5 * bulk_speed ** 2.0
    
    flowState.vel_x = bulk_speed
    flowState.fluid_state.update_thermo_from_hs(h = static_enthalpy, s = flowState.fluid_state.entropy)

    for cell in range(len(mesh.cell_array)):
        mesh.cellArray[cell].fs.vel_x = flowState.vel_x
        mesh.cellArray[cell].fs.fluid_state.copy_values(flowState.fluid_state)
        
    return mesh