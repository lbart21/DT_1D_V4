"""
Function:
Author: Luke Bartholomew
Edits:
"""
from copy import deepcopy

from Algorithms.DesignToolAlgorithmV4_1D.ComponentModels.basicMeshObject import basic1DMeshObject
from Algorithms.DesignToolAlgorithmV4_1D.Reconstruction.LocateNeighbouringCellIndices import find_idx_of_cell_recursively
from Algorithms.DesignToolAlgorithmV4_1D.Reconstruction.LocateNeighbouringInterfaceIndices import find_idx_of_interface_recursively
from simple_outflow_bc import SimpleOutFlow_BC
from Algorithms.DesignToolAlgorithmV4_1D.BoundaryConditions.supersonic_inflow_bc import supersonic_inflow_bc
from from_stagnation_inflow_bc import FromStagnationInFlow_BC
from from_stagnation_with_mass_flow_rate_inflow_bc import FromStagnationWithMassFlowRateInFlow_BC
from Algorithms.DesignToolAlgorithmV4_1D.BoundaryConditions.fixed_p_outflow_bc import fixed_p_outflow_bc
from Algorithms.DesignToolAlgorithmV4_1D.BoundaryConditions.fixed_pt_outflow_bc import fixed_pt_outflow_bc
from wall_with_slip_bc import WallWithSlip_BC
from wall_no_slip_bc import WallNoSlip_BC


class FormBoundaryConditionInformation():
    def __init__(self, mesh, BC) -> None:
        ### Find out if BC is on east or west boundary
        onWestBoundaryBool = False
        onEastBoundaryBool = False
        if mesh.map_interface_id_to_west_cell_idx[BC[0]] is None:
            onWestBoundaryBool = True

        elif mesh.map_interface_id_to_east_cell_idx[BC[0]] is None:
            onEastBoundaryBool = True
        
        else:
            print("Boundary is not on an edge of the mesh")
            print(BC)
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
            self.ghostCellLayer = supersonic_inflow_bc(mesh = self.ghostCellLayer, b_c = BC)
            #for cell in self.ghostCellLayer.cell_array:
                #print(cell.flow_state.fluid_state.p, cell.flow_state.vel_x)
        
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
            self.ghostCellLayer = fixed_p_outflow_bc(mesh = self.ghostCellLayer, BC = BC)

        elif BC[1][0] == "FixedPTOutFlow_BC":
            if onEastBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = False)
            elif onWestBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = True)

            # Copy cells and interfaces from inner layer to get correct geometry props 
            self.MirrorCopyCells(BC = BC, onEastBoundaryBool = onEastBoundaryBool, \
                                    onWestBoundaryBool = onWestBoundaryBool, mesh = mesh)
            self.ghostCellLayer = fixed_pt_outflow_bc(mesh = self.ghostCellLayer, BC = BC)
        


    def MirrorCopyCells(self, BC, onEastBoundaryBool, onWestBoundaryBool, mesh):
        for cell in range(BC[1][1]): #Index [1][1] of BC gives how many ghost cells get generated
            if onEastBoundaryBool:
                cellIdx = find_idx_of_cell_recursively(interface_id = BC[0], recursion_depth = cell + 1, \
                                                        direction = "West", map_interface_id_to_east_cell_idx = mesh.map_interface_id_to_east_cell_idx, \
                                                        map_cell_id_to_east_interface_idx = mesh.map_cell_id_to_east_interface_idx, \
                                                        map_interface_id_to_west_cell_idx = mesh.map_interface_id_to_west_cell_idx, \
                                                        map_cell_id_to_west_interface_idx = mesh.map_cell_id_to_west_interface_idx, \
                                                        cell_array = mesh.cell_array, interface_array = mesh.interface_array)
            elif onWestBoundaryBool:
                cellIdx = find_idx_of_cell_recursively(interface_id = BC[0], recursion_depth = cell + 1, \
                                                        direction = "East", map_interface_id_to_east_cell_idx = mesh.map_interface_id_to_east_cell_idx, \
                                                        map_cell_id_to_east_interface_idx = mesh.map_cell_id_to_east_interface_idx, \
                                                        map_interface_id_to_west_cell_idx = mesh.map_interface_id_to_west_cell_idx, \
                                                        map_cell_id_to_west_interface_idx = mesh.map_cell_id_to_west_interface_idx, \
                                                        cell_array = mesh.cell_array, interface_array = mesh.interface_array)
            
            interior_cell = mesh.cell_array[cellIdx]
            interior_cell_gs_class = interior_cell.flow_state.fluid_state.__class__
            ghost_cell_gs = interior_cell_gs_class(interior_cell.flow_state.fluid_state.gmodel)
            ghost_cell_gs.copy_values(interior_cell.flow_state.fluid_state)
            ghostCell = deepcopy(interior_cell)
            ghostCell.flow_state.fluid_state = ghost_cell_gs
            
            ghostCell.InteriorFlag = False
            self.ghostCellLayer.cell_array[cell] = ghostCell
            
        for interface in range(BC[1][1] + 1):
            if onEastBoundaryBool:
                interfaceIdx = find_idx_of_interface_recursively(interface_id = BC[0], recursion_depth = interface, \
                                                                direction = "West", map_interface_id_to_east_cell_idx = mesh.map_interface_id_to_east_cell_idx, \
                                                                map_cell_id_to_east_interface_idx = mesh.map_cell_id_to_east_interface_idx, \
                                                                map_interface_id_to_west_cell_idx = mesh.map_interface_id_to_west_cell_idx, \
                                                                map_cell_id_to_west_interface_idx = mesh.map_cell_id_to_west_interface_idx, \
                                                                cell_array = mesh.cell_array)
            elif onWestBoundaryBool:
                interfaceIdx = find_idx_of_interface_recursively(interface_id = BC[0], recursion_depth = interface, \
                                                                direction = "East", map_interface_id_to_east_cell_idx = mesh.map_interface_id_to_east_cell_idx, \
                                                                map_cell_id_to_east_interface_idx = mesh.map_cell_id_to_east_interface_idx, \
                                                                map_interface_id_to_west_cell_idx = mesh.map_interface_id_to_west_cell_idx, \
                                                                map_cell_id_to_west_interface_idx = mesh.map_cell_id_to_west_interface_idx, \
                                                                cell_array = mesh.cell_array)
            interior_interface = mesh.interface_array[interfaceIdx]
            interior_interface_gsL_class = interior_interface.lft_state.fluid_state.__class__
            interior_interface_gsR_class = interior_interface.rght_state.fluid_state.__class__
            lft_state_gs = interior_interface_gsL_class(interior_interface.lft_state.fluid_state.gmodel)
            rght_state_gs = interior_interface_gsR_class(interior_interface.rght_state.fluid_state.gmodel)
            lft_state_gs.copy_values(interior_interface.lft_state.fluid_state)
            rght_state_gs.copy_values(interior_interface.rght_state.fluid_state)
            ghostInterface = deepcopy(interior_interface)
            ghostInterface.lft_state.fluid_state = lft_state_gs
            ghostInterface.rght_state.fluid_state = rght_state_gs
            ghostInterface.flux_flag = False
            self.ghostCellLayer.interface_array[interface] = ghostInterface
        
