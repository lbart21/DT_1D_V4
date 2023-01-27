"""
Function:
Author: Luke Bartholomew
Edits:
"""
from copy import deepcopy

def FromStagnationInFlow_BC(mesh, BC, onWestBoundaryBool):
    """
    GasState = BC[1][2]
    Take p_stag and T_stag to create a stagnation state. Use velocity on inside of 
    boundary and stagnation enthalpy to find static enthalpy. Update properties from 
    stagnation entropy and static enthalpy. All cells have the interior cell's velocity.
    Velocity is set to 0 if flow is out of boundary.
    Then update the flow properties to complete the state.
    """
    stag_gas_state = BC[1][2]
    p_stag = stag_gas_state.p
    T_stag = stag_gas_state.T
    vel_x_boundary = mesh.cellArray[0].fs.vel_x
    
    if onWestBoundaryBool:
        bulk_speed = max(0.0, vel_x_boundary) # 0.0 if vel_x_boundary negative, vel_x_boundary if positive, 
        
    else:
        bulk_speed = min(0.0, vel_x_boundary)

    flowState = deepcopy(mesh.cellArray[0].fs)
    flowState.fluid_state.p = p_stag
    flowState.fluid_state.T = T_stag
    flowState.vel_x = bulk_speed
    flowState.fluid_state.update_thermo_from_pT()
    stagnation_enthalpy = flowState.fluid_state.enthalpy
    static_enthalpy = stagnation_enthalpy - 0.5 * bulk_speed ** 2.0

    flowState.fluid_state.update_thermo_from_hs(h = static_enthalpy, s = flowState.fluid_state.entropy)

    for cell in range(len(mesh.cell_array)):
        mesh.cellArray[cell].fs.vel_x = flowState.vel_x
        mesh.cellArray[cell].fs.fluid_state.copy_values(flowState.fluid_state)
        
    return mesh