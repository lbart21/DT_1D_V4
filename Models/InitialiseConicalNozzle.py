"""
Function:
Author: Luke Bartholomew
Edits:
"""
from Algorithms.DesignToolAlgorithmV4_1D.ComponentModels.SinglePhaseInterface import SinglePhaseInterface
from Algorithms.DesignToolAlgorithmV4_1D.ComponentModels.SinglePhaseConicalNozzleCell import SinglePhaseConicalNozzleCell
from Algorithms.DesignToolAlgorithmV4_1D.ComponentModels.meshObject import meshObject
import numpy as np
from copy import deepcopy

class SinglePhaseConicalNozzle(meshObject):
    def __init__(self, nCells, Geometry, flowState, reconstructionScheme, limiter, \
                reconstructionProperties, updateFrom, componentLabel, fluxScheme) -> None:
        """
        Geometry = [D1, D2, L]
        """
        super().__init__(nCells)

        self.componentLabels = [componentLabel]

        #print(empty_flowState_to_copy_in.fluid_state.p)

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
        pass

    def initialiseCells(self, Geometry, nCells, flowState, componentLabel):
        [D1, D2, L] = Geometry 
        #print(empty_flowState.fluid_state.p)
        
        for cell in range(nCells):
            flowState_object = flowState.__class__
            fluidState_object = flowState.fluid_state.__class__
            fluidModel_filename = flowState.fluid_state.gmodel.file_name
            fluidModel_object = flowState.fluid_state.gmodel.__class__
            gm = fluidModel_object(fluidModel_filename)
            gs = fluidState_object(gm)
            fs = flowState_object(gs)

            cellObject = SinglePhaseConicalNozzleCell(cell_ID = cell, label = componentLabel)
            dx = L / nCells
            D_L = self.D_at_x(D1 = D1, D2 = D2, x = cell * dx, L = L)
            D_c = self.D_at_x(D1 = D1, D2 = D2, x = (0.5 + cell) * dx, L = L)
            D_R = self.D_at_x(D1 = D1, D2 = D2, x = (1.0 + cell) * dx, L = L)

            GEO = {
                "dx"    :   dx,
                "dV"    :   np.pi * (D_L ** 2.0 + D_L * D_R + D_R ** 2.0) * dx / 12.0,
                "A_c"   :   0.25 * np.pi * D_c ** 2.0,
                "A_s"   :   0.5 * np.pi * (D_L + D_R) * dx * (1.0 + (0.5 * (D_R - D_L) / dx) ** 2.0) ** 0.5,
                "pos_x" :   (0.5 + cell) * dx
            }
            cellObject.fillGeometry(Geometry = GEO)
            #print(empty_flowState.fluid_state.p)
            cellObject.fs = fs
            cellObject.fs.fluid_state.copy_values(flowState.fluid_state)
            #print(empty_flowState.fluid_state.p)
            cellObject.fs.vel_x = flowState.vel_x
            cellObject.initialiseConservedProperties()
            self.cellArray[cell] = cellObject

    def initialiseInterfaces(self, Geometry, nCells, fluxScheme, reconstructionScheme, limiter, reconstructionProperties, updateFrom, flowState):
        [D1, D2, L] = Geometry
        #print(empty_flowState.fluid_state.p)
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

            interfaceObject = SinglePhaseInterface(interface_ID = interface, nL = interface, updateFrom = updateFrom, \
                                                    nR = nCells - interface, fluxScheme = fluxScheme, \
                                                    reconstructionScheme = reconstructionScheme, limiter = limiter, \
                                                    reconstructionProperties = reconstructionProperties)
            D = self.D_at_x(D1 = D1, D2 = D2, x = interface * L / nCells, L = L)
            GEO = {"A"  : 0.25 * np.pi * D ** 2}
            interfaceObject.fillGeometry(Geometry = GEO)
            #print(empty_flowState.fluid_state.p)
            interfaceObject.LftState = fs_Lft
            #interfaceObject.LftState.fluid_state.copy_values(empty_flowState.fluid_state)
            interfaceObject.RghtState = fs_Rght
            #interfaceObject.RghtState.fluid_state.copy_values(empty_flowState.fluid_state)
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
    
    def D_at_x(self, D1, D2, x, L):
        return D1 + (D2 - D1) * x / L
