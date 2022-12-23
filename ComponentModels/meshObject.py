class meshObject():
    def __init__(self, nCells):
        self.cellArray = [None] * nCells
        self.interfaceArray = [None] * (nCells + 1)
        self.mapCellIDToWestInterfaceIdx = [None] * nCells
        self.mapCellIDToEastInterfaceIdx = [None] * nCells
        self.mapInterfaceIDToWestCellIdx = [None] * (nCells + 1)
        self.mapInterfaceIDToEastCellIdx = [None] * (nCells + 1)
        self.boundaryConditions = []
        self.boundaryInterfaceIDs = []
        self.componentLabels = []
        self.cellIdxToTrack = []
        self.interfaceIdxToTrack = []
        
