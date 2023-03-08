"""
Function:
Author: Luke Bartholomew
Edits:
"""
from copy import deepcopy

def from_stagnation_inflow_bc(mesh, BC, on_west_boundary_bool):
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
    vel_x_boundary = mesh.cell_array[0].fs.vel_x
    
    if on_west_boundary_bool:
        bulk_speed = max(0.0, vel_x_boundary) # 0.0 if vel_x_boundary negative, vel_x_boundary if positive, 
        
    else:
        bulk_speed = min(0.0, vel_x_boundary)

    flow_state_interior = deepcopy(mesh.cellArray[0].flow_state)
    flow_state_interior.fluid_state.p = p_stag
    flow_state_interior.fluid_state.T = T_stag
    flow_state_interior.vel_x = bulk_speed
    flow_state_interior.fluid_state.update_thermo_from_pT()
    stagnation_enthalpy = flow_state_interior.fluid_state.enthalpy
    static_enthalpy = stagnation_enthalpy - 0.5 * bulk_speed ** 2.0

    flow_state_interior.fluid_state.update_thermo_from_hs(h = static_enthalpy, s = flow_state_interior.fluid_state.entropy)

    for cell in mesh.cell_array:
        cell.flow_state.vel_x = flow_state_interior.vel_x
        cell.flow_state.fluid_state.copy_values(flow_state_interior.fluid_state)
        
    return mesh