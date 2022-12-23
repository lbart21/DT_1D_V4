def findIdxOfInterfaceRecursively(interfaceID, recursionDepth, direction, \
                                mapInterfaceIDToEastCellIdx, cellArray, \
                                mapCellIDToEastInterfaceIdx, mapInterfaceIDToWestCellIdx, \
                                mapCellIDToWestInterfaceIdx, interfaceArray):
    currentInterfaceID = interfaceID

   
    currentRecursionDepth = 0
    while currentRecursionDepth < recursionDepth:
        if direction == "East":
            neighbouringCellIdx = mapInterfaceIDToEastCellIdx[currentInterfaceID]
            neighbouringCellID = cellArray[neighbouringCellIdx].cell_ID
            currentInterfaceID = mapCellIDToEastInterfaceIdx[neighbouringCellID]
            currentRecursionDepth += 1
            
        elif direction == "West":
            neighbouringCellIdx = mapInterfaceIDToWestCellIdx[currentInterfaceID]
            neighbouringCellID = cellArray[neighbouringCellIdx].cell_ID
            currentInterfaceID = mapCellIDToWestInterfaceIdx[neighbouringCellID]
            currentRecursionDepth += 1
            
    return currentInterfaceID