"""
Function:
Author: Luke Bartholomew
Edits:
"""
def wall_with_slip_bc(mesh):
    """
    Ghost cells are mirrored values of the interior cells. Specify outlet pressure
    and update from PT.
    Cells will have same pressure but not necessarily the same velocity, temperature etc.
    """
    for cell in mesh.cell_array:
        cell.flow_state.vel_x *= -1.0
   
    return mesh