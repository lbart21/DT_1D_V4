def SupersonicInFlow_BC(mesh, BC):
    """
    Simply put specified flow state which is in index [1][2] of BC list in order
    [vel_x, p, T] = BC[1][2]
    Then update the flow properties to complete the state.
    All cells have the same properties
    """
    fs = BC[1][2]
    vel_x = fs.vel_x
    for cell in range(len(mesh.cellArray)):
        mesh.cellArray[cell].fs.vel_x = vel_x
        mesh.cellArray[cell].fs.fluid_state.copy_values(fs.fluid_state)

    return mesh