from Algorithms.DesignToolAlgorithmV3_1D.ComponentModels.basicMeshObject import basic1DMeshObject
from Algorithms.DesignToolAlgorithmV3_1D.Reconstruction.LocateNeighbouringCellIndices import findIdxOfCellRecursively
from Algorithms.DesignToolAlgorithmV3_1D.Reconstruction.LocateNeighbouringInterfaceIndices import findIdxOfInterfaceRecursively
from Algorithms.DesignToolAlgorithmV3_1D.BoundaryConditions.SimpleOutFlow_BC import SimpleOutFlow_BC
from Algorithms.DesignToolAlgorithmV3_1D.BoundaryConditions.SupersonicInflow_BC import SupersonicInFlow_BC
from Algorithms.DesignToolAlgorithmV3_1D.BoundaryConditions.FromStagnationInFlow_BC import FromStagnationInFlow_BC
from Algorithms.DesignToolAlgorithmV3_1D.BoundaryConditions.FromStagnationWithMassFlowRateInFlow_BC import FromStagnationWithMassFlowRateInFlow_BC
from Algorithms.DesignToolAlgorithmV3_1D.BoundaryConditions.FixedPOutFlow_BC import FixedPOutFlow_BC
from Algorithms.DesignToolAlgorithmV3_1D.BoundaryConditions.FixedPTOutFlow_BC import FixedPTOutFlow_BC
from Algorithms.DesignToolAlgorithmV3_1D.BoundaryConditions.WallWithSlip_BC import WallWithSlip_BC
from Algorithms.DesignToolAlgorithmV3_1D.BoundaryConditions.WallNoSlip_BC import WallNoSlip_BC
from copy import deepcopy

class FormBoundaryConditionInformation():
    def __init__(self, mesh, BC) -> None:
        ### Find out if BC is on east or west boundary
        onWestBoundaryBool = False
        onEastBoundaryBool = False
        if mesh.mapInterfaceIDToWestCellIdx[BC[0]] == None:
            onWestBoundaryBool = True

        elif mesh.mapInterfaceIDToEastCellIdx[BC[0]] == None:
            onEastBoundaryBool = True
        
        else:
            print("Boundary is not on an edge of the mesh")
        ### Wall BCs
        if BC[1][0] == "WallNoSlip_BC":
            if onEastBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = False)
            elif onWestBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = True)
            # Copy cells from inner layer to get correct geometry props 
            self.MirrorCopyCells(BC = BC, onEastBoundaryBool = onEastBoundaryBool, \
                                    onWestBoundaryBool = onWestBoundaryBool, mesh = mesh)
            self.ghostCellLayer = WallNoSlip_BC(mesh = self.ghostCellLayer)
        
        elif BC[1][0] == "WallWithSlip_BC":
            if onEastBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = False)
            elif onWestBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = True)
            # Copy cells from inner layer to get correct geometry props 
            self.MirrorCopyCells(BC = BC, onEastBoundaryBool = onEastBoundaryBool, \
                                    onWestBoundaryBool = onWestBoundaryBool, mesh = mesh)
            self.ghostCellLayer = WallWithSlip_BC(mesh = self.ghostCellLayer)

        ### InFlow BCs        
        elif BC[1][0] == "SupersonicInFlow_BC":
            if onEastBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = False)
            elif onWestBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = True)
            # Copy cells from inner layer to get correct geometry props 
            self.MirrorCopyCells(BC = BC, onEastBoundaryBool = onEastBoundaryBool, \
                                    onWestBoundaryBool = onWestBoundaryBool, mesh = mesh)
            self.ghostCellLayer = SupersonicInFlow_BC(mesh = self.ghostCellLayer, BC = BC)
        
        elif BC[1][0] == "FromStagnationInFlow_BC":
            if onEastBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = False)
            elif onWestBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = True)
            # Copy cells from inner layer to get correct geometry props 
            self.MirrorCopyCells(BC = BC, onEastBoundaryBool = onEastBoundaryBool, \
                                    onWestBoundaryBool = onWestBoundaryBool, mesh = mesh)
            self.ghostCellLayer = FromStagnationInFlow_BC(mesh = self.ghostCellLayer, BC = BC, onWestBoundaryBool = onWestBoundaryBool)
    
        elif BC[1][0] == "FromStagnationWithMassFlowRateInFlow_BC":
            if onEastBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = False)
            elif onWestBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = True)
            # Copy cells from inner layer to get correct geometry props 
            self.MirrorCopyCells(BC = BC, onEastBoundaryBool = onEastBoundaryBool, \
                                    onWestBoundaryBool = onWestBoundaryBool, mesh = mesh)
            self.ghostCellLayer = FromStagnationWithMassFlowRateInFlow_BC(mesh = self.ghostCellLayer, BC = BC, onWestBoundaryBool = onWestBoundaryBool)
            
        ### OutFlow BCs
        elif BC[1][0] == "SimpleOutFlow_BC":
            if onEastBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = False)
            elif onWestBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = True)

            # Copy cells and interfaces from inner layer to get correct geometry props 
            self.MirrorCopyCells(BC = BC, onEastBoundaryBool = onEastBoundaryBool, \
                                    onWestBoundaryBool = onWestBoundaryBool, mesh = mesh)
            self.ghostCellLayer = SimpleOutFlow_BC(mesh = self.ghostCellLayer)

        elif BC[1][0] == "SimpleExtrapolateOutFlow_BC":
            pass

        elif BC[1][0] == "FixedPOutFlow_BC":
            if onEastBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = False)
            elif onWestBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = True)

            # Copy cells and interfaces from inner layer to get correct geometry props 
            self.MirrorCopyCells(BC = BC, onEastBoundaryBool = onEastBoundaryBool, \
                                    onWestBoundaryBool = onWestBoundaryBool, mesh = mesh)
            self.ghostCellLayer = FixedPOutFlow_BC(mesh = self.ghostCellLayer, BC = BC)

        elif BC[1][0] == "FixedPTOutFlow_BC":
            if onEastBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = False)
            elif onWestBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = True)

            # Copy cells and interfaces from inner layer to get correct geometry props 
            self.MirrorCopyCells(BC = BC, onEastBoundaryBool = onEastBoundaryBool, \
                                    onWestBoundaryBool = onWestBoundaryBool, mesh = mesh)
            self.ghostCellLayer = FixedPTOutFlow_BC(mesh = self.ghostCellLayer, BC = BC)
        


    def MirrorCopyCells(self, BC, onEastBoundaryBool, onWestBoundaryBool, mesh):
        for cell in range(BC[1][1]): #Index [1][1] of BC gives how many ghost cells get generated
            if onEastBoundaryBool:
                cellIdx = findIdxOfCellRecursively(interfaceID = BC[0], recursionDepth = cell + 1, \
                                                        direction = "West", mapInterfaceIDToEastCellIdx = mesh.mapInterfaceIDToEastCellIdx, \
                                                        mapCellIDToEastInterfaceIdx = mesh.mapCellIDToEastInterfaceIdx, \
                                                        mapInterfaceIDToWestCellIdx = mesh.mapInterfaceIDToWestCellIdx, \
                                                        mapCellIDToWestInterfaceIdx = mesh.mapCellIDToWestInterfaceIdx, \
                                                        cellArray = mesh.cellArray, interfaceArray = mesh.interfaceArray)
            elif onWestBoundaryBool:
                cellIdx = findIdxOfCellRecursively(interfaceID = BC[0], recursionDepth = cell + 1, \
                                                        direction = "East", mapInterfaceIDToEastCellIdx = mesh.mapInterfaceIDToEastCellIdx, \
                                                        mapCellIDToEastInterfaceIdx = mesh.mapCellIDToEastInterfaceIdx, \
                                                        mapInterfaceIDToWestCellIdx = mesh.mapInterfaceIDToWestCellIdx, \
                                                        mapCellIDToWestInterfaceIdx = mesh.mapCellIDToWestInterfaceIdx, \
                                                        cellArray = mesh.cellArray, interfaceArray = mesh.interfaceArray)
            
            interior_cell = mesh.cellArray[cellIdx]
            interior_cell_gs_class = interior_cell.fs.fluid_state.__class__
            ghost_cell_gs = interior_cell_gs_class(interior_cell.fs.fluid_state.gmodel)
            ghost_cell_gs.copy_values(interior_cell.fs.fluid_state)
            ghostCell = deepcopy(interior_cell)
            ghostCell.fs.fluid_state = ghost_cell_gs
            
            ghostCell.InteriorFlag = False
            self.ghostCellLayer.cellArray[cell] = ghostCell
            
        for interface in range(BC[1][1] + 1):
            if onEastBoundaryBool:
                interfaceIdx = findIdxOfInterfaceRecursively(interfaceID = BC[0], recursionDepth = interface, \
                                                                direction = "West", mapInterfaceIDToEastCellIdx = mesh.mapInterfaceIDToEastCellIdx, \
                                                                mapCellIDToEastInterfaceIdx = mesh.mapCellIDToEastInterfaceIdx, \
                                                                mapInterfaceIDToWestCellIdx = mesh.mapInterfaceIDToWestCellIdx, \
                                                                mapCellIDToWestInterfaceIdx = mesh.mapCellIDToWestInterfaceIdx, \
                                                                cellArray = mesh.cellArray, interfaceArray = mesh.interfaceArray)
            elif onWestBoundaryBool:
                interfaceIdx = findIdxOfInterfaceRecursively(interfaceID = BC[0], recursionDepth = interface, \
                                                                direction = "East", mapInterfaceIDToEastCellIdx = mesh.mapInterfaceIDToEastCellIdx, \
                                                                mapCellIDToEastInterfaceIdx = mesh.mapCellIDToEastInterfaceIdx, \
                                                                mapInterfaceIDToWestCellIdx = mesh.mapInterfaceIDToWestCellIdx, \
                                                                mapCellIDToWestInterfaceIdx = mesh.mapCellIDToWestInterfaceIdx, \
                                                                cellArray = mesh.cellArray, interfaceArray = mesh.interfaceArray)
            interior_interface = mesh.interfaceArray[interfaceIdx]
            interior_interface_gsL_class = interior_interface.LftState.fluid_state.__class__
            interior_interface_gsR_class = interior_interface.RghtState.fluid_state.__class__
            LftState_gs = interior_interface_gsL_class(interior_interface.LftState.fluid_state.gmodel)
            RghtState_gs = interior_interface_gsR_class(interior_interface.RghtState.fluid_state.gmodel)
            LftState_gs.copy_values(interior_interface.LftState.fluid_state)
            RghtState_gs.copy_values(interior_interface.RghtState.fluid_state)
            ghostInterface = deepcopy(interior_interface)
            ghostInterface.LftState.fluid_state = LftState_gs
            ghostInterface.RghtState.fluid_state = RghtState_gs
            ghostInterface.FluxFlag = False
            self.ghostCellLayer.interfaceArray[interface] = ghostInterface
        
