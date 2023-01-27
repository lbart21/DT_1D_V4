"""
Function:
Author: Luke Bartholomew
Edits:
"""
def supersonic_inflow_bc(mesh, b_c):
    """
    Simply put specified flow state which is in index [1][2] of BC list in order
    [vel_x, p, T] = BC[1][2]
    Then update the flow properties to complete the state.
    All cells have the same properties
    """
    flow_state = b_c[1][2]
    vel_x = flow_state.vel_x
    for cell in mesh.cell_array:
        cell.flow_state.vel_x = vel_x
        cell.flow_state.fluid_state.copy_values(flow_state.fluid_state)

    return mesh
