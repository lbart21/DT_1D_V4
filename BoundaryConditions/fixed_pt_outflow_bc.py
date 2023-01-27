"""
Function:
Author: Luke Bartholomew
Edits:
"""
def fixed_pt_outflow_bc(mesh, BC):
    """
    Ghost cells are mirrored values of the interior cells. Specify outlet pressure
    and temperature and update from PT.
    Cells will have same pressure and temperature (and thus all thermodynamic properties), 
    but not necessarily the same velocity.
    """
    gs = BC[1][2]
    p = gs.p
    T = gs.T

    for cell in mesh.cell_array:
        cell.flow_state.fluid_state.p = p
        cell.flow_state.fluid_state.T = T
        cell.flow_stiate.fluid_state.update_thermo_from_pT()

    return mesh