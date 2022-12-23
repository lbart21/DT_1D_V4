def findIdxOfCellRecursively(interfaceID, recursionDepth, direction, \
                                mapInterfaceIDToEastCellIdx, cellArray, \
                                mapCellIDToEastInterfaceIdx, mapInterfaceIDToWestCellIdx, \
                                mapCellIDToWestInterfaceIdx, interfaceArray):
        """ Given interface ID, direction you want to look in and how far you want to look, return the index of the cell that you end at
        Eg. recursionDepth = 1, direction = "East" will give you the index of the East cell of the interface [ | | | | | !0| | | | | | ]
            recursionDepth = 3, direction = "West" will give you the index of the West cell of the West cell of the West cell of the interface [ | | | |0| | ! | | | | | ]
        """
        currentInterfaceID = interfaceID
        currentRecursionDepth = 0
        while currentRecursionDepth < recursionDepth:
            if direction == "East":
                newCellIdx = mapInterfaceIDToEastCellIdx[currentInterfaceID]
                newCellID = cellArray[newCellIdx].cell_ID
                currentRecursionDepth += 1
                if currentRecursionDepth != recursionDepth:
                    currentInterfaceID = interfaceArray[mapCellIDToEastInterfaceIdx[newCellID]].interface_ID
            
            if direction == "West":
                newCellIdx = mapInterfaceIDToWestCellIdx[currentInterfaceID]
                newCellID = cellArray[newCellIdx].cell_ID
                currentRecursionDepth += 1
                if currentRecursionDepth != recursionDepth:
                    currentInterfaceID = interfaceArray[mapCellIDToWestInterfaceIdx[newCellID]].interface_ID
        return newCellID