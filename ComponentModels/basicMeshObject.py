class basic1DMeshObject():
    """
    Basically the same as the meshObject class, but an assumed cell/interface 
    order allows an easier definition of blocks.
    reversed decides if cells are expected to be generated in a normal fashion
    of west to east or in reverse from east to west. This is only really used in 
    boundary cell generation where west boundary ghost cell layers are generated 
    in reverse order"""
    def __init__(self, nCells, reversed = False) -> None:
        self.cellArray = [None] * nCells
        self.interfaceArray = [None] * (nCells + 1)
        if reversed:
            self.mapCellIDToWestInterfaceIdx = [i + 1 for i in range(nCells)]
            self.mapCellIDToEastInterfaceIdx = [i for i in range(nCells)]
            self.mapInterfaceIDToWestCellIdx = [i for i in range(nCells)] + [None] 
            self.mapInterfaceIDToEastCellIdx = [None] + [i for i in range(nCells)] 
            
        else:
            self.mapCellIDToWestInterfaceIdx = [i for i in range(nCells)]
            self.mapCellIDToEastInterfaceIdx = [i + 1 for i in range(nCells)]
            self.mapInterfaceIDToWestCellIdx = [None] + [i for i in range(nCells)]
            self.mapInterfaceIDToEastCellIdx = [i for i in range(nCells)] + [None]
        self.componentLabels = []
        self.boundaryConditions = []
        self.boundaryInterfaceIDs = [0, nCells]
        self.cellIdxToTrack = []
        self.interfaceIdxToTrack = []
        