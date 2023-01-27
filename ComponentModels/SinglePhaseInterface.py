"""
Function:
Author: Luke Bartholomew
Edits:
"""
import math as m
import importlib
import sys
from Algorithms.DesignToolAlgorithmV4_1D.Reconstruction.Reconstruction import getReconstruction
from Algorithms.DesignToolAlgorithmV4_1D.Fluxes.fluidFluxes import *
from Algorithms.DesignToolAlgorithmV4_1D.Reconstruction.LocateNeighbouringCellIndices import findIdxOfCellRecursively
from Algorithms.DesignToolAlgorithmV4_1D.FluidModel.FlowState import FlowState
from gdtk.gas import GasModel, GasState

class SinglePhaseInterface():
    def __init__(self, interface_ID, nL, nR, reconstructionScheme, limiter, reconstructionProperties, updateFrom, fluxScheme) -> None:
        self.interface_ID = interface_ID
        self.nL = nL
        self.nR = nR
        
        self.fluxScheme = fluxScheme
        self.FluxFlag = True
        self.reconstructionScheme = reconstructionScheme
        self.limiter = limiter
        self.reconstructionProperties = reconstructionProperties
        self.updateFrom = updateFrom
        self.stencilBeenFormed = False
        
        self.LftState = None
        self.RghtState = None
        

    def completeInterfaceMethods(self, cellArray, interfaceArray, mapCellIDToWestInterfaceIdx, \
                                        mapCellIDToEastInterfaceIdx, mapInterfaceIDToWestCellIdx, \
                                        mapInterfaceIDToEastCellIdx, dt_inv):
        if self.FluxFlag:
            if not self.stencilBeenFormed: 
                self.FormListsOfStencilIndices(mapCellIDToWestInterfaceIdx = mapCellIDToWestInterfaceIdx, \
                                                mapCellIDToEastInterfaceIdx = mapCellIDToEastInterfaceIdx, \
                                                mapInterfaceIDToWestCellIdx = mapInterfaceIDToWestCellIdx, \
                                                mapInterfaceIDToEastCellIdx = mapInterfaceIDToEastCellIdx, \
                                                cellArray = cellArray, interfaceArray = interfaceArray)
                self.stencilBeenFormed = True

            self.reconstructStates(cellArray = cellArray)
            self.calculateFluxes()
            self.updateNeighbouringCellsCPs(dt_inv = dt_inv, cellArray = cellArray)
    
    def fillGeometry(self, Geometry):
        self.GEO = Geometry

    def reconstructStates(self, cellArray):
        """
        reconstructionScheme = ["Copy", [rnL, rnR]] etc etc where rnL and rnR are the REQUIRED stencil lengths in the left and right directions for the chosen scheme
        reconstructionProperties = list of properties to reconstruct eg ["p", "T", "vel_x", "alpha_g"]
        updateFrom = str(fluidPropertyUpdatePair) eg "pT", "rhoP", "rhoU"
        """
        #print(self.interface_ID, self.LftStencilIdxs, self.RghtStencilIdxs)
        if self.FluxFlag == True: #Do reconstruction
            for prop in self.reconstructionProperties:
                qL_stencil, dxL_stencil, \
                qR_stencil, dxR_stencil = self.getStencils(property = prop, cellArray = cellArray)
                LftProp, RghtProp = getReconstruction(reconstruction = self.reconstructionScheme[0], limiter = self.limiter, \
                                                        qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
                if prop == "vel_x":
                    self.LftState.vel_x = LftProp
                    self.RghtState.vel_x = RghtProp
                else:
                    setattr(self.LftState.fluid_state, prop, LftProp)
                    setattr(self.RghtState.fluid_state, prop, RghtProp)

            #print(self.LftState, self.RghtState)
            if self.updateFrom == "pT":
                self.LftState.fluid_state.update_thermo_from_pT()
                self.RghtState.fluid_state.update_thermo_from_pT()
            #print(self.LftState.fs["p"], self.RghtState.fs["p"])
        else: #Not an interface to calculate fluxes over, so skip
            pass

    def calculateFluxes(self):
        if self.FluxFlag == True:
            #print(self.LftStencilIdxs, self.RghtStencilIdxs)
            if self.fluxScheme == "AUSMPlusUPOriginal":
                self.boundaryFluxes = AUSMPlusupORIGINAL(   a_L_forMa = self.LftState.fluid_state.a,        a_R_forMa = self.RghtState.fluid_state.a, \
                                                            p_L = self.LftState.fluid_state.p,              rho_L = self.LftState.fluid_state.rho, \
                                                            u_L = self.LftState.fluid_state.u,              a_L = self.LftState.fluid_state.a, \
                                                            vel_x_L = self.LftState.vel_x,                  p_R = self.RghtState.fluid_state.p, \
                                                            rho_R = self.RghtState.fluid_state.rho,         u_R = self.RghtState.fluid_state.u, \
                                                            a_R = self.RghtState.fluid_state.a,             vel_x_R = self.RghtState.vel_x).fluxes

            elif self.fluxScheme == "AUSMPlusUPPaper":
                self.boundaryFluxes = AUSMPlusupPAPER(  a_L_forMa = self.LftState.fs["a"],  a_L = self.LftState.fs["a"], \
                                                        p_L = self.LftState.fs["p"],        h_L = self.LftState.fs["h"], \
                                                        rho_L = self.LftState.fs["rho"],    vel_x_L = self.LftState.fs["vel_x"], \
                                                        a_R_forMa = self.RghtState.fs["a"], a_R = self.RghtState.fs["a"], \
                                                        p_R = self.RghtState.fs["p"],       h_R = self.RghtState.fs["h"], \
                                                        rho_R = self.RghtState.fs["rho"],   vel_x_R = self.RghtState.fs["vel_x"]).fluxes
        
        else:
            pass

        #print("Fluxes: ", self.interface_ID, self.boundaryFluxes)
        #print("Lft and Rght interface values:" , self.LftState.fluid_state.p, self.RghtState.fluid_state.p, self.LftState.fluid_state.rho, self.RghtState.fluid_state.rho)

    def updateNeighbouringCellsCPs(self, dt_inv, cellArray):
        if self.FluxFlag == True:
            westCell = cellArray[self.LftStencilIdxs[0]]
            if westCell.InteriorCellFlag: # Check if cell is a ghost cell as we don't want to update the conserved properties of ghost cells
                #print(dt_inv, self.GEO["A"], self.boundaryFluxes["mass"], westCell.GEO["dV"])
                westCell.conservedProperties["mass"] -= dt_inv * self.GEO["A"] * self.boundaryFluxes["mass"] / westCell.GEO["dV"]
                westCell.conservedProperties["xMom"] -= dt_inv * self.GEO["A"] * self.boundaryFluxes["xMom"] / westCell.GEO["dV"]
                westCell.conservedProperties["xMom"] -= dt_inv * westCell.GEO["A_c"] * self.boundaryFluxes["p"] / westCell.GEO["dV"]
                westCell.conservedProperties["energy"] -= dt_inv * self.GEO["A"] * self.boundaryFluxes["energy"] / westCell.GEO["dV"]

            eastCell = cellArray[self.RghtStencilIdxs[0]]
            if eastCell.InteriorCellFlag: # Check if cell is a ghost cell as we don't want to update the conserved properties of ghost cells
                eastCell.conservedProperties["mass"] += dt_inv * self.GEO["A"] * self.boundaryFluxes["mass"] / eastCell.GEO["dV"]
                eastCell.conservedProperties["xMom"] += dt_inv * self.GEO["A"] * self.boundaryFluxes["xMom"] / eastCell.GEO["dV"]
                eastCell.conservedProperties["xMom"] += dt_inv * eastCell.GEO["A_c"] * self.boundaryFluxes["p"] / eastCell.GEO["dV"]
                eastCell.conservedProperties["energy"] += dt_inv * self.GEO["A"] * self.boundaryFluxes["energy"] / eastCell.GEO["dV"]
                
        else:
            pass

    def getStencils(self, property, cellArray):
        ### Fill left stencil
        qL_stencil = [None] * len(self.LftStencilIdxs)
        dxL_stencil = [None] * len(self.LftStencilIdxs)
        for ind, LftCellIdx in enumerate(self.LftStencilIdxs):
            cell_L = cellArray[LftCellIdx]
            if property == "vel_x":
                qL_stencil[ind] = cell_L.fs.vel_x
            else:
                qL_stencil[ind] = getattr(cell_L.fs.fluid_state, property)
            dxL_stencil[ind] = cell_L.GEO["dx"]

        ### Fill right stencil
        qR_stencil = [None] * len(self.RghtStencilIdxs)
        dxR_stencil = [None] * len(self.RghtStencilIdxs)
        for ind, RghtCellIdx in enumerate(self.RghtStencilIdxs):
            cell_R = cellArray[RghtCellIdx]
            if property == "vel_x":
                qR_stencil[ind] = cell_R.fs.vel_x
            else:
                qR_stencil[ind] = getattr(cell_R.fs.fluid_state, property)
            dxR_stencil[ind] = cell_R.GEO["dx"]
            
        return qL_stencil, dxL_stencil, qR_stencil, dxR_stencil
        

    def FormListsOfStencilIndices(self, mapCellIDToWestInterfaceIdx, mapCellIDToEastInterfaceIdx, \
                                        mapInterfaceIDToWestCellIdx, mapInterfaceIDToEastCellIdx, cellArray, interfaceArray):
        self.LftStencilIdxs = [None] * self.reconstructionScheme[1][0]
        self.RghtStencilIdxs = [None] * self.reconstructionScheme[1][1]
        for LftStencilIdx in range(self.reconstructionScheme[1][0]):
            cellIdx = findIdxOfCellRecursively(interfaceID = self.interface_ID, recursionDepth = LftStencilIdx + 1, \
                                                direction = "West", mapInterfaceIDToEastCellIdx = mapInterfaceIDToEastCellIdx, \
                                                cellArray = cellArray, mapCellIDToEastInterfaceIdx = mapCellIDToEastInterfaceIdx, \
                                                mapInterfaceIDToWestCellIdx = mapInterfaceIDToWestCellIdx, \
                                                mapCellIDToWestInterfaceIdx = mapCellIDToWestInterfaceIdx, \
                                                interfaceArray = interfaceArray)
            self.LftStencilIdxs[LftStencilIdx] = cellIdx #In order of [L_Idx0, L_Idx1, ...]
        
        for RghtStencilIdx in range(self.reconstructionScheme[1][1]):
            cellIdx = findIdxOfCellRecursively(interfaceID = self.interface_ID, recursionDepth = RghtStencilIdx + 1, \
                                                direction = "East", mapInterfaceIDToEastCellIdx = mapInterfaceIDToEastCellIdx, \
                                                cellArray = cellArray, mapCellIDToEastInterfaceIdx = mapCellIDToEastInterfaceIdx, \
                                                mapInterfaceIDToWestCellIdx = mapInterfaceIDToWestCellIdx, \
                                                mapCellIDToWestInterfaceIdx = mapCellIDToWestInterfaceIdx, \
                                                interfaceArray = interfaceArray)
            self.RghtStencilIdxs[RghtStencilIdx] = cellIdx #In order of [R_Idx0, R_Idx1, ...]
   