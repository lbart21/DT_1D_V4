class jointBlock():
    def __init__(self, meshObject1, meshObject2, block1InterfaceIDBeingReplaced, block2InterfaceIDBeingReplaced, newInterface, AddingGhostCellsBool = False, newComponentLabel = None) -> None: 
        """ 
        Given two meshes, makes a joined mesh, including cell array, interface array 
        and mapping arrays from cellIDs to interface indicies and vice versa.
        Assume block 1 is always the larger one, so we'll be adding 2 onto 1 rather than 1 onto 2 to make the resizing process faster.
        joiningBlocksBool = True if joining separate components, in which case the boundary interface indices have to be updated.
                          = False if just adding boundary cells, hence boundary interface indices don't have to be updated.
        """
                #--------------------Preliminary Calcs----------------------------------
        ### Find if joining to east or west edge of meshObject1
        joiningToWestEdge = False
        joiningToEastEdge = False
        if meshObject1.mapInterfaceIDToWestCellIdx[block1InterfaceIDBeingReplaced] == None:
            joiningToWestEdge = True
        elif meshObject1.mapInterfaceIDToEastCellIdx[block1InterfaceIDBeingReplaced] == None:
            joiningToEastEdge = True
        else:
            print("Index given for meshObject1 interface is not a mesh boundary interface.")

        ### Will need to know how many cells and interfaces are in each mesh
        nCellsMesh1 = len(meshObject1.cellArray)
        nInterfacesMesh1 = len(meshObject1.interfaceArray)
        nCellsMesh2 = len(meshObject2.cellArray)
        nInterfacesMesh2 = len(meshObject2.interfaceArray)

        ### When joining blocks, cell locations need to be updated. Only want to do this when joining blocks, not when adding interface
        # Will do this by always updating the east block, so the positions are always positive.
        # Restricting it to always updating mesh2 locations runs into issues when joining to west face of mesh1 because 
        # locations become negative and the calculation of these locations is inefficient 
        # Process is find location of joining cell of west block. Add half of that cell's width to find location of interface.
        # For every cell in east block, add this interface location to the pos_x values
        if not AddingGhostCellsBool:
            if joiningToWestEdge: #Need to find location of joining cell of mesh2. Need to look west of joining interface of mesh2
                mesh2JoiningBoundaryCellIdx = meshObject2.mapInterfaceIDToEastCellIdx[block2InterfaceIDBeingReplaced]
                mesh2JoiningBoundaryCell = meshObject2.cellArray[mesh2JoiningBoundaryCellIdx]
                mesh2JoiningBoundaryCellLocation = mesh2JoiningBoundaryCell.GEO["pos_x"]
                mesh2JoiningBoundaryCellWidth = mesh2JoiningBoundaryCell.GEO["dx"]
                for mesh1Cell in meshObject1.cellArray:
                    mesh1Cell.GEO["pos_x"] += mesh2JoiningBoundaryCellLocation + 0.5 * mesh2JoiningBoundaryCellWidth
                
            elif joiningToEastEdge: #Need to find location of joining block of block1. Need to look west of interface of mesh2.
                mesh1JoiningBoundaryCellIdx = meshObject1.mapInterfaceIDToWestCellIdx[block1InterfaceIDBeingReplaced]
                mesh1JoiningBoundaryCell = meshObject1.cellArray[mesh1JoiningBoundaryCellIdx]
                mesh1JoiningBoundaryCellLocation = mesh1JoiningBoundaryCell.GEO["pos_x"]
                mesh1JoiningBoundaryCellWidth = mesh1JoiningBoundaryCell.GEO["dx"]
                for mesh2Cell in meshObject2.cellArray:
                    mesh2Cell.GEO["pos_x"] += mesh1JoiningBoundaryCellLocation + 0.5 * mesh1JoiningBoundaryCellWidth      
            if newComponentLabel != None:
                self.componentLabels = [newComponentLabel]
            else:    
                self.componentLabels = meshObject1.componentLabels + meshObject2.componentLabels
            print("Individual meshs' boundary interface indices were: ", meshObject1.boundaryInterfaceIDs, " and ", meshObject2.boundaryInterfaceIDs)
        elif AddingGhostCellsBool:
            self.componentLabels = meshObject1.componentLabels
        #-----------------------------------------------------------------------

        # Join cell Arrays
        self.cellArray = meshObject1.cellArray + meshObject2.cellArray
        if newComponentLabel != None: #Check if something has been put in newComponentLabel
            for cell in self.cellArray:
                cell.label = newComponentLabel
                            
        # Join interfaceArrays
        meshObject1.interfaceArray[block1InterfaceIDBeingReplaced] = newInterface #Put new interface object into mesh1.interfaceArray
        del meshObject2.interfaceArray[block2InterfaceIDBeingReplaced]
        if joiningToWestEdge:   # On west edge of mesh1, so need to increase values of nL in mesh1.interfaceArray by number of cells in mesh2. 
                                # Then need to increase values of nR in mesh2.interfaceArray by number of cells in mesh1.
            for interface in meshObject1.interfaceArray:
                interface.nL += nCellsMesh2
            for interface in meshObject2.interfaceArray:
                interface.nR += nCellsMesh1
        if joiningToEastEdge:   # On east edge of mesh1, so need to increase values of nR in mesh1.interfaceArray by number of cells in mesh2.
                                # Then need to increase values of nL in mesh2.interfaceArray by number of cells in mesh1.
            for interface in meshObject1.interfaceArray:
                interface.nR += nCellsMesh2
            for interface in meshObject2.interfaceArray:
                interface.nL += nCellsMesh1

        self.interfaceArray = meshObject1.interfaceArray + meshObject2.interfaceArray

        # Join mapCellIDToWestInterfaceIdx
        for ind, value in enumerate(meshObject2.mapCellIDToWestInterfaceIdx):
            if value <= block2InterfaceIDBeingReplaced:
                meshObject2.mapCellIDToWestInterfaceIdx[ind] += nInterfacesMesh1
            else:
                meshObject2.mapCellIDToWestInterfaceIdx[ind] -= 1
                meshObject2.mapCellIDToWestInterfaceIdx[ind] += nInterfacesMesh1

        if joiningToEastEdge:
            mesh2BoundaryCellIdx = meshObject2.mapInterfaceIDToEastCellIdx[block2InterfaceIDBeingReplaced]
            meshObject2.mapCellIDToWestInterfaceIdx[mesh2BoundaryCellIdx] = block1InterfaceIDBeingReplaced

        self.mapCellIDToWestInterfaceIdx = meshObject1.mapCellIDToWestInterfaceIdx + meshObject2.mapCellIDToWestInterfaceIdx

        # Join mapCellIDToEastInterfaceIdx
        for ind, value in enumerate(meshObject2.mapCellIDToEastInterfaceIdx):
            if value <= block2InterfaceIDBeingReplaced:
                meshObject2.mapCellIDToEastInterfaceIdx[ind] += nInterfacesMesh1
            else:
                meshObject2.mapCellIDToEastInterfaceIdx[ind] -= 1
                meshObject2.mapCellIDToEastInterfaceIdx[ind] += nInterfacesMesh1

        if joiningToWestEdge:
            mesh2BoundaryCellIdx = meshObject2.mapInterfaceIDToWestCellIdx[block2InterfaceIDBeingReplaced]
            meshObject2.mapCellIDToEastInterfaceIdx[mesh2BoundaryCellIdx] = block1InterfaceIDBeingReplaced

        self.mapCellIDToEastInterfaceIdx = meshObject1.mapCellIDToEastInterfaceIdx + meshObject2.mapCellIDToEastInterfaceIdx

        # Join mapInterfaceIDToWestCellIdx
        if joiningToWestEdge:
            mesh2BoundaryCellIdx = meshObject2.mapInterfaceIDToWestCellIdx[block2InterfaceIDBeingReplaced]
            mesh2BoundaryCellIdx += nCellsMesh1
            meshObject1.mapInterfaceIDToWestCellIdx[block1InterfaceIDBeingReplaced] = mesh2BoundaryCellIdx
        
        del meshObject2.mapInterfaceIDToWestCellIdx[block2InterfaceIDBeingReplaced]

        for ind, value in enumerate(meshObject2.mapInterfaceIDToWestCellIdx):
            if value == None:
                pass
            else:
                meshObject2.mapInterfaceIDToWestCellIdx[ind] += nCellsMesh1

        self.mapInterfaceIDToWestCellIdx = meshObject1.mapInterfaceIDToWestCellIdx + meshObject2.mapInterfaceIDToWestCellIdx

        # Join mapInterfaceIDToEastCellIdx
        if joiningToEastEdge:
            mesh2BoundaryCellIdx = meshObject2.mapInterfaceIDToEastCellIdx[block2InterfaceIDBeingReplaced]
            mesh2BoundaryCellIdx += nCellsMesh1
            meshObject1.mapInterfaceIDToEastCellIdx[block1InterfaceIDBeingReplaced] = mesh2BoundaryCellIdx
        
        del meshObject2.mapInterfaceIDToEastCellIdx[block2InterfaceIDBeingReplaced]

        for ind, value in enumerate(meshObject2.mapInterfaceIDToEastCellIdx):
            if value == None:
                pass
            else:
                meshObject2.mapInterfaceIDToEastCellIdx[ind] += nCellsMesh1

        self.mapInterfaceIDToEastCellIdx = meshObject1.mapInterfaceIDToEastCellIdx + meshObject2.mapInterfaceIDToEastCellIdx
        
        for cell_ID in range(len(self.cellArray)):
            self.cellArray[cell_ID].cell_ID = cell_ID

        for interface_ID in range(len(self.interfaceArray)):
            self.interfaceArray[interface_ID].interface_ID = interface_ID
        
        meshObject1.boundaryInterfaceIDs.remove(block1InterfaceIDBeingReplaced)
        meshObject2.boundaryInterfaceIDs.remove(block2InterfaceIDBeingReplaced)

        meshObject2.boundaryInterfaceIDs[0] += nInterfacesMesh1 - 1 #meshObject1.boundaryInterfaceIDs doesn't change 

        if meshObject2.boundaryConditions: #Check if not empty
            meshObject2.boundaryConditions[0][0] += nInterfacesMesh1 - 1    # Update the index of the boundary condition of mesh2.
                                                                            # There can only be 1 BC, so we can safely index at 0 
                                                                            # in meshObject2.boundaryCondition list
      
        self.boundaryInterfaceIDs = meshObject1.boundaryInterfaceIDs + meshObject2.boundaryInterfaceIDs 
        if not AddingGhostCellsBool:
            print("New mesh boundary interface indices are: ", self.boundaryInterfaceIDs)
        self.boundaryConditions = meshObject1.boundaryConditions + meshObject2.boundaryConditions

        if not AddingGhostCellsBool:
            self.cellIdxToTrack = meshObject1.cellIdxToTrack + [i + nCellsMesh1 for i in meshObject2.cellIdxToTrack]
            self.interfaceIdxToTrack = meshObject1.interfaceIdxToTrack + [j + nInterfacesMesh1 - 1 for j in meshObject2.interfaceIdxToTrack]

        else:
            self.interfaceIdxToTrack = meshObject1.interfaceIdxToTrack
            self.cellIdxToTrack = meshObject1.cellIdxToTrack