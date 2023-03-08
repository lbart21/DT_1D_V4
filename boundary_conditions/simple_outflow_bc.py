"""
Function:
Author: Luke Bartholomew
Edits:
"""
def simple_outflow_bc(mesh):
    """
    Zero gradient BC, ghost cells are already generated in a zero-gradient sense, 
    so no need to change anything
    """
    return mesh