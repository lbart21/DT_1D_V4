def WallWithSlip_BC(mesh):
    """
    Ghost cells are mirrored values of the interior cells. Specify outlet pressure
    and update from PT.
    Cells will have same pressure but not necessarily the same velocity, temperature etc.
    """
    for cell in range(len(mesh.cellArray)):
        mesh.cellArray[cell].fs.vel_x *= -1
   
    return mesh