"""
Function:
Author: Luke Bartholomew
Edits:
"""
from Algorithms.DesignToolAlgorithmV4_1D.ComponentModels.SinglePhaseStraightPipeCell import SinglePhaseStraightPipeCell
from Algorithms.DesignToolAlgorithmV4_1D.ComponentModels.SinglePhaseInterface import SinglePhaseInterface
from Algorithms.DesignToolAlgorithmV4_1D.ComponentModels.meshObject import meshObject

import numpy as np

class SinglePhaseStraightPipe(meshObject):
    def __init__(self, nCells, Geometry, flowState, reconstructionProperties, \
                        reconstructionScheme, limiter, updateFrom, \
                        componentLabel, fluxScheme) -> None:
        """
        Geometry = [D, L]
        """
        super().__init__(nCells)

        self.componentLabels = [componentLabel]
        ### Define cells
        self.initialiseCells(Geometry = Geometry, flowState = flowState, \
                            componentLabel = componentLabel, nCells = nCells)

        ### Define interfaces
        self.initialiseInterfaces(Geometry = Geometry, limiter = limiter, \
                                reconstructionScheme = reconstructionScheme, \
                                reconstructionProperties = reconstructionProperties, \
                                updateFrom = updateFrom, fluxScheme = fluxScheme, \
                                nCells = nCells, flowState = flowState)

        ### Form connections
        self.connectCellsToInterfaces()

    def initialiseCells(self, Geometry, componentLabel, nCells, flowState):
        [D, L] = Geometry
    
        for cell in range(nCells):
            flowState_object = flowState.__class__
            fluidState_object = flowState.fluid_state.__class__
            fluidModel_filename = flowState.fluid_state.gmodel.file_name
            fluidModel_object = flowState.fluid_state.gmodel.__class__
            gm = fluidModel_object(fluidModel_filename)
            gs = fluidState_object(gm)
            fs = flowState_object(gs)

            cellObject = SinglePhaseStraightPipeCell(cell_ID = cell, label = componentLabel)
            GEO = {
                "dx"    :   L / nCells,
                "dV"    :   0.25 * np.pi * D ** 2 * L / nCells,
                "A_c"   :   0.25 * np.pi * D ** 2,
                "A_s"   :   np.pi * D * L / nCells,
                "pos_x" :   (0.5 + cell) * L / nCells
            }
            cellObject.fillGeometry(Geometry = GEO)
            cellObject.fs = fs
            cellObject.fs.fluid_state.copy_values(flowState.fluid_state)
            #print(empty_flowState.fluid_state.p)
            cellObject.fs.vel_x = flowState.vel_x
            cellObject.initialiseConservedProperties()            
            self.cellArray[cell] = cellObject

    def initialiseInterfaces(self, Geometry, reconstructionScheme, limiter, reconstructionProperties, updateFrom, fluxScheme, nCells, flowState):
        [D, L] = Geometry
        for interface in range(nCells + 1):
            flowState_object = flowState.__class__
            fluidState_object = flowState.fluid_state.__class__
            fluidModel_filename = flowState.fluid_state.gmodel.file_name
            fluidModel_object = flowState.fluid_state.gmodel.__class__
            gm_Lft = fluidModel_object(fluidModel_filename)
            gm_Rght = fluidModel_object(fluidModel_filename)
            gs_Lft = fluidState_object(gm_Lft)
            gs_Rght = fluidState_object(gm_Rght)
            fs_Lft = flowState_object(gs_Lft)
            fs_Rght = flowState_object(gs_Rght)

            interfaceObject = SinglePhaseInterface(interface_ID = interface, nL = interface, \
                                                    nR = nCells - interface, \
                                                    fluxScheme = fluxScheme, \
                                                    reconstructionScheme = reconstructionScheme, \
                                                    limiter = limiter, \
                                                    reconstructionProperties = reconstructionProperties, \
                                                    updateFrom = updateFrom)
            GEO = {"A"  : 0.25 * np.pi * D ** 2}
            interfaceObject.fillGeometry(Geometry = GEO)
            interfaceObject.LftState = fs_Lft
            interfaceObject.RghtState = fs_Rght
            self.interfaceArray[interface] = interfaceObject
        
        self.boundaryInterfaceIDs = [0, nCells]

    def connectCellsToInterfaces(self):
        nCells = len(self.cellArray)
        self.mapCellIDToWestInterfaceIdx = [None] * nCells
        self.mapCellIDToEastInterfaceIdx = [None] * nCells
        self.mapInterfaceIDToWestCellIdx = [None] * (nCells + 1)
        self.mapInterfaceIDToEastCellIdx = [None] * (nCells + 1)

        for array_idx in range(nCells):
            self.mapCellIDToWestInterfaceIdx[array_idx] = array_idx
            self.mapCellIDToEastInterfaceIdx[array_idx] = array_idx + 1
            self.mapInterfaceIDToWestCellIdx[array_idx + 1] = array_idx
            self.mapInterfaceIDToEastCellIdx[array_idx] = array_idx
          
    def addBoundaryConditions(self, BC):
        self.boundaryConditions.append(BC)






 
            


        


