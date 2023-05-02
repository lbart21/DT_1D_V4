"""
Function:
Author: Luke Bartholomew
Edits:
"""

def mdot_from_stagnation_inflow_bc(mesh, BC, on_west_boundary_bool):
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
    #print("BC massf:", stag_gas_state.massf)
    #print("BC pressure:", stag_gas_state.p)

    p0_min = 0.1 * p_stag
    p0_max = 10.0 * p_stag
    vel_x_boundary = mesh.cell_array[0].flow_state.vel_x
    p_boundary = mesh.cell_array[0].flow_state.fluid_state.p
    rho_boundary = mesh.cell_array[0].flow_state.fluid_state.rho
    #print("rho_boundary", rho_boundary)
    rho_v_boundary = rho_boundary * vel_x_boundary
    A_boundary = mesh.interface_array[0].geo["A"]
    #print("A_boundary", A_boundary)
    dp_over_p = 0.5 * relaxation_factor / rho_boundary * ( (mass_flux / A_boundary) ** 2.0 \
                                                            - rho_v_boundary * abs(rho_v_boundary) ) / p_boundary
    #print("dp_over_p", dp_over_p)
    new_p0 = (1.0 + dp_over_p) * p_stag
    #print("new_p0", new_p0)
    new_p0 = min(max(new_p0, p0_min), p0_max)
    #print("new_p0", new_p0)
    #flow_state_interior = deepcopy(mesh.cell_array[0].flow_state)

    gm_object = mesh.cell_array[0].flow_state.fluid_state.gmodel.__class__
    gm_filename = mesh.cell_array[0].flow_state.fluid_state.gmodel.file_name
    gs_object = mesh.cell_array[0].flow_state.fluid_state.__class__
    fs_object = mesh.cell_array[0].flow_state.__class__

    gm = gm_object(gm_filename)
    gs = gs_object(gm)
    fs = fs_object(state = gs)

    fs.fluid_state.copy_values(stag_gas_state)
    fs.fluid_state.p = new_p0
    fs.fluid_state.T = T_stag
    fs.fluid_state.update_thermo_from_pT()

    stagnation_enthalpy = fs.fluid_state.enthalpy
    stagnation_entropy = fs.fluid_state.entropy


    #flow_state_interior.fluid_state.p = new_p0
    #flow_state_interior.fluid_state.T = T_stag
    #print("Stagnation conditions:", flow_state_interior.fluid_state.massf, new_p0, T_stag)
    #flow_state_interior.fluid_state.update_thermo_from_pT()
    #stagnation_enthalpy = flow_state_interior.fluid_state.enthalpy
    
    if on_west_boundary_bool:
        bulk_speed = max(0.0, vel_x_boundary) # 0.0 if vel_x_boundary negative, vel_x_boundary if positive, 
        
    else:
        bulk_speed = min(0.0, vel_x_boundary)

    static_enthalpy = stagnation_enthalpy - 0.5 * bulk_speed ** 2.0

    fs.vel_x = bulk_speed
    fs.fluid_state.update_thermo_from_hs(h = static_enthalpy, s = stagnation_entropy)
    
    #flow_state_interior.vel_x = bulk_speed
    #flow_state_interior.fluid_state.update_thermo_from_hs(h = static_enthalpy, s = flow_state_interior.fluid_state.entropy)

    for ind, cell in enumerate(mesh.cell_array):
        mesh.cell_array[ind].flow_state.vel_x = fs.vel_x
        mesh.cell_array[ind].flow_state.fluid_state.copy_values(fs.fluid_state)
        
    return mesh