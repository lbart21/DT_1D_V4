"""
Function:
Author: Luke Bartholomew
Edits:
"""
def fixed_p_outflow_bc(mesh, BC):
    """
    Ghost cells are mirrored values of the interior cells. Specify outlet pressure
    and update from PT.
    Cells will have same pressure but not necessarily the same velocity, temperature etc.
    """
    gs = BC[1][2]
    p = gs.p

    for cell in mesh.cell_array:
        cell.flow_state.fluid_state.p = p
        cell.flow_state.fluid_state.update_thermo_from_pT()

    return mesh
    