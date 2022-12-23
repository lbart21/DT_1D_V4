def SimpleOutFlow_BC(mesh):
    """
    Zero gradient BC, ghost cells are already generated in a zero-gradient sense, 
    so no need to change anything
    """
    return mesh