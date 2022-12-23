def FixedPOutFlow_BC(mesh, BC):
    """
    Ghost cells are mirrored values of the interior cells. Specify outlet pressure
    and update from PT.
    Cells will have same pressure but not necessarily the same velocity, temperature etc.
    """
    gs = BC[1][2]
    p = gs.p

    for cell in range(len(mesh.cellArray)):
        mesh.cellArray[cell].fs.fluid_state.p = p
        mesh.cellArray[cell].fs.fluid_state.update_thermo_from_pT()

    return mesh