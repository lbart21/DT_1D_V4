def FixedPTOutFlow_BC(mesh, BC):
    """
    Ghost cells are mirrored values of the interior cells. Specify outlet pressure
    and temperature and update from PT.
    Cells will have same pressure and temperature (and thus all thermodynamic properties), 
    but not necessarily the same velocity.
    """
    gs = BC[1][2]
    p = gs.p
    T = gs.T

    for cell in range(len(mesh.cellArray)):
        mesh.cellArray[cell].fs.fluid_state.p = p
        mesh.cellArray[cell].fs.fluid_state.T = T
        mesh.cellArray[cell].flowState.update_thermo_from_pT()

    return mesh