"""
Function:
Author: Luke Bartholomew
Edits:
"""
from copy import deepcopy

from Algorithms.DT_1D_V4.models.prefilled_single_inlet_mesh_object import SingleInlet1DMeshObject
from Algorithms.DT_1D_V4.reconstruction.locate_neighbouring_cell_indices import find_idx_of_cell_recursively
from Algorithms.DT_1D_V4.reconstruction.locate_neighbouring_interface_indices import find_idx_of_interface_recursively
from Algorithms.DT_1D_V4.boundary_conditions.simple_outflow_bc import simple_outflow_bc
from Algorithms.DT_1D_V4.boundary_conditions.supersonic_inflow_bc import supersonic_inflow_bc
from Algorithms.DT_1D_V4.boundary_conditions.from_stagnation_inflow_bc import from_stagnation_inflow_bc
from Algorithms.DT_1D_V4.boundary_conditions.mdot_from_stagnation_inflow_bc import mdot_from_stagnation_inflow_bc
from Algorithms.DT_1D_V4.boundary_conditions.fixed_p_outflow_bc import fixed_p_outflow_bc
from Algorithms.DT_1D_V4.boundary_conditions.fixed_pt_outflow_bc import fixed_pt_outflow_bc
from Algorithms.DT_1D_V4.boundary_conditions.wall_with_slip_bc import wall_with_slip_bc
from Algorithms.DT_1D_V4.boundary_conditions.wall_no_slip_bc import wall_no_slip_bc


class FormBoundaryConditionInformation():
    def __init__(self, mesh, BC) -> None:
        ### Find out if BC is on east or west boundary
        on_west_boundary_bool = False
        on_east_boundary_bool = False
        if mesh.map_interface_id_to_west_cell_idx[BC[0]] is None:
            on_west_boundary_bool = True

        elif mesh.map_interface_id_to_east_cell_idx[BC[0]] is None:
            on_east_boundary_bool = True
        
        else:
            print("Boundary is not on an edge of the mesh")
            print(BC)
            print("This is what the west cell ID is:", mesh.map_interface_id_to_west_cell_idx[BC[0]])
            print("This is what the east cell ID is:", mesh.map_interface_id_to_east_cell_idx[BC[0]])
        ### Wall BCs
        if BC[1][0] == "WallNoSlip_BC":
            if on_east_boundary_bool:
                self.ghost_cell_layer = SingleInlet1DMeshObject(n_cells = BC[1][1], reverse_direction_for_ghost_cells = False)
            elif on_west_boundary_bool:
                self.ghost_cell_layer = SingleInlet1DMeshObject(n_cells = BC[1][1], reverse_direction_for_ghost_cells = True)
            # Copy cells from inner layer to get correct geometry props 
            self.mirror_copy_cells(BC = BC, on_east_boundary_bool = on_east_boundary_bool, \
                                    on_west_boundary_bool = on_west_boundary_bool, mesh = mesh)
            self.ghost_cell_layer = wall_no_slip_bc(mesh = self.ghost_cell_layer)
        
        elif BC[1][0] == "WallWithSlip_BC":
            if on_east_boundary_bool:
                self.ghost_cell_layer = SingleInlet1DMeshObject(n_cells = BC[1][1], reverse_direction_for_ghost_cells = False)
            elif on_west_boundary_bool:
                self.ghost_cell_layer = SingleInlet1DMeshObject(n_cells = BC[1][1], reverse_direction_for_ghost_cells = True)
            # Copy cells from inner layer to get correct geometry props 
            self.mirror_copy_cells(BC = BC, on_east_boundary_bool = on_east_boundary_bool, \
                                    on_west_boundary_bool = on_west_boundary_bool, mesh = mesh)
            self.ghost_cell_layer = wall_with_slip_bc(mesh = self.ghost_cell_layer)

        ### InFlow BCs        
        elif BC[1][0] == "SupersonicInFlow_BC":
            if on_east_boundary_bool:
                self.ghost_cell_layer = SingleInlet1DMeshObject(n_cells = BC[1][1], reverse_direction_for_ghost_cells = False)
            elif on_west_boundary_bool:
                self.ghost_cell_layer = SingleInlet1DMeshObject(n_cells = BC[1][1], reverse_direction_for_ghost_cells = True)
            # Copy cells from inner layer to get correct geometry props 
            self.mirror_copy_cells(BC = BC, on_east_boundary_bool = on_east_boundary_bool, \
                                    on_west_boundary_bool = on_west_boundary_bool, mesh = mesh)
            self.ghost_cell_layer = supersonic_inflow_bc(mesh = self.ghost_cell_layer, b_c = BC)
            #for cell in self.ghostCellLayer.cell_array:
                #print(cell.flow_state.fluid_state.p, cell.flow_state.vel_x)
        
        elif BC[1][0] == "FromStagnationInFlow_BC":
            if on_east_boundary_bool:
                self.ghost_cell_layer = SingleInlet1DMeshObject(n_cells = BC[1][1], reverse_direction_for_ghost_cells = False)
            elif on_west_boundary_bool:
                self.ghost_cell_layer = SingleInlet1DMeshObject(n_cells = BC[1][1], reverse_direction_for_ghost_cells = True)
            # Copy cells from inner layer to get correct geometry props 
            self.mirror_copy_cells(BC = BC, on_east_boundary_bool = on_east_boundary_bool, \
                                    on_west_boundary_bool = on_west_boundary_bool, mesh = mesh)
            self.ghost_cell_layer = from_stagnation_inflow_bc(mesh = self.ghost_cell_layer, BC = BC, on_west_boundary_bool = on_west_boundary_bool)
    
        elif BC[1][0] == "FromStagnationWithMassFlowRateInFlow_BC":
            if on_east_boundary_bool:
                self.ghost_cell_layer = SingleInlet1DMeshObject(n_cells = BC[1][1], reverse_direction_for_ghost_cells = False)
            elif on_west_boundary_bool:
                self.ghost_cell_layer = SingleInlet1DMeshObject(n_cells = BC[1][1], reverse_direction_for_ghost_cells = True)
            # Copy cells from inner layer to get correct geometry props 
            self.mirror_copy_cells(BC = BC, on_east_boundary_bool = on_east_boundary_bool, \
                                    on_west_boundary_bool = on_west_boundary_bool, mesh = mesh)
            self.ghost_cell_layer = mdot_from_stagnation_inflow_bc(mesh = self.ghost_cell_layer, BC = BC, on_west_boundary_bool = on_west_boundary_bool)
            
        ### OutFlow BCs
        elif BC[1][0] == "SimpleOutFlow_BC":
            if on_east_boundary_bool:
                self.ghost_cell_layer = SingleInlet1DMeshObject(n_cells = BC[1][1], reverse_direction_for_ghost_cells = False)
            elif on_west_boundary_bool:
                self.ghost_cell_layer = SingleInlet1DMeshObject(n_cells = BC[1][1], reverse_direction_for_ghost_cells = True)

            # Copy cells and interfaces from inner layer to get correct geometry props 
            self.mirror_copy_cells(BC = BC, on_east_boundary_bool = on_east_boundary_bool, \
                                    on_west_boundary_bool = on_west_boundary_bool, mesh = mesh)
            self.ghost_cell_layer = simple_outflow_bc(mesh = self.ghost_cell_layer)

        elif BC[1][0] == "SimpleExtrapolateOutFlow_BC":
            print("SimpleExtrapolateOutFlow_BC has not been implemented yet.")

        elif BC[1][0] == "FixedPOutFlow_BC":
            if on_east_boundary_bool:
                self.ghost_cell_layer = SingleInlet1DMeshObject(n_cells = BC[1][1], reverse_direction_for_ghost_cells = False)
            elif on_west_boundary_bool:
                self.ghost_cell_layer = SingleInlet1DMeshObject(n_cells = BC[1][1], reverse_direction_for_ghost_cells = True)

            # Copy cells and interfaces from inner layer to get correct geometry props 
            self.mirror_copy_cells(BC = BC, on_east_boundary_bool = on_east_boundary_bool, \
                                    on_west_boundary_bool = on_west_boundary_bool, mesh = mesh)
            self.ghost_cell_layer = fixed_p_outflow_bc(mesh = self.ghost_cell_layer, BC = BC)

        elif BC[1][0] == "FixedPTOutFlow_BC":
            if on_east_boundary_bool:
                self.ghost_cell_layer = SingleInlet1DMeshObject(n_cells = BC[1][1], reverse_direction_for_ghost_cells = False)
            elif on_west_boundary_bool:
                self.ghost_cell_layer = SingleInlet1DMeshObject(n_cells = BC[1][1], reverse_direction_for_ghost_cells = True)

            # Copy cells and interfaces from inner layer to get correct geometry props 
            self.mirror_copy_cells(BC = BC, on_east_boundary_bool = on_east_boundary_bool, \
                                    on_west_boundary_bool = on_west_boundary_bool, mesh = mesh)
            self.ghost_cell_layer = fixed_pt_outflow_bc(mesh = self.ghost_cell_layer, BC = BC)
        


    def mirror_copy_cells(self, BC, on_east_boundary_bool, on_west_boundary_bool, mesh):
        for cell in range(BC[1][1]): #Index [1][1] of BC gives how many ghost cells get generated
            if on_east_boundary_bool:
                cell_idx, hit_bad_interface = find_idx_of_cell_recursively(interface_id = BC[0], recursion_depth = cell + 1, \
                                                        direction = "West", map_interface_id_to_east_cell_idx = mesh.map_interface_id_to_east_cell_idx, \
                                                        map_cell_id_to_east_interface_idx = mesh.map_cell_id_to_east_interface_idx, \
                                                        map_interface_id_to_west_cell_idx = mesh.map_interface_id_to_west_cell_idx, \
                                                        map_cell_id_to_west_interface_idx = mesh.map_cell_id_to_west_interface_idx, \
                                                        cell_array = mesh.cell_array, interface_array = mesh.interface_array)
            elif on_west_boundary_bool:
                cell_idx, hit_bad_interface = find_idx_of_cell_recursively(interface_id = BC[0], recursion_depth = cell + 1, \
                                                        direction = "East", map_interface_id_to_east_cell_idx = mesh.map_interface_id_to_east_cell_idx, \
                                                        map_cell_id_to_east_interface_idx = mesh.map_cell_id_to_east_interface_idx, \
                                                        map_interface_id_to_west_cell_idx = mesh.map_interface_id_to_west_cell_idx, \
                                                        map_cell_id_to_west_interface_idx = mesh.map_cell_id_to_west_interface_idx, \
                                                        cell_array = mesh.cell_array, interface_array = mesh.interface_array)
            
            interior_cell = mesh.cell_array[cell_idx]
            interior_cell_gs_class = interior_cell.flow_state.fluid_state.__class__
            ghost_cell_gs = interior_cell_gs_class(interior_cell.flow_state.fluid_state.gmodel)
            ghost_cell_gs.copy_values(interior_cell.flow_state.fluid_state)
            ghost_cell = deepcopy(interior_cell)
            ghost_cell.flow_state.fluid_state = ghost_cell_gs
            
            ghost_cell.InteriorFlag = False
            self.ghost_cell_layer.cell_array[cell] = ghost_cell
            
        for interface in range(BC[1][1] + 1):
            if on_east_boundary_bool:
                interface_idx = find_idx_of_interface_recursively(interface_id = BC[0], recursion_depth = interface, \
                                                                direction = "West", map_interface_id_to_east_cell_idx = mesh.map_interface_id_to_east_cell_idx, \
                                                                map_cell_id_to_east_interface_idx = mesh.map_cell_id_to_east_interface_idx, \
                                                                map_interface_id_to_west_cell_idx = mesh.map_interface_id_to_west_cell_idx, \
                                                                map_cell_id_to_west_interface_idx = mesh.map_cell_id_to_west_interface_idx, \
                                                                cell_array = mesh.cell_array)
            elif on_west_boundary_bool:
                interface_idx = find_idx_of_interface_recursively(interface_id = BC[0], recursion_depth = interface, \
                                                                direction = "East", map_interface_id_to_east_cell_idx = mesh.map_interface_id_to_east_cell_idx, \
                                                                map_cell_id_to_east_interface_idx = mesh.map_cell_id_to_east_interface_idx, \
                                                                map_interface_id_to_west_cell_idx = mesh.map_interface_id_to_west_cell_idx, \
                                                                map_cell_id_to_west_interface_idx = mesh.map_cell_id_to_west_interface_idx, \
                                                                cell_array = mesh.cell_array)
            interior_interface = mesh.interface_array[interface_idx]
            interior_interface_lft_gs_class = interior_interface.lft_state.fluid_state.__class__
            interior_interface_rght_gs_class = interior_interface.rght_state.fluid_state.__class__
            lft_state_gs = interior_interface_lft_gs_class(interior_interface.lft_state.fluid_state.gmodel)
            rght_state_gs = interior_interface_rght_gs_class(interior_interface.rght_state.fluid_state.gmodel)
            lft_state_gs.copy_values(interior_interface.lft_state.fluid_state)
            rght_state_gs.copy_values(interior_interface.rght_state.fluid_state)
            ghost_interface = deepcopy(interior_interface)
            ghost_interface.lft_state.fluid_state = lft_state_gs
            ghost_interface.rght_state.fluid_state = rght_state_gs
            ghost_interface.flux_flag = False
            self.ghost_cell_layer.interface_array[interface] = ghost_interface
        
