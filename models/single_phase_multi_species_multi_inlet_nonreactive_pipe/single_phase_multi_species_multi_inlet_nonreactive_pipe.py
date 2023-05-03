"""
Function:
Author: Luke Bartholomew
Edits:
"""
import numpy as np

from Algorithms.DT_1D_V4.models.prefilled_multi_inlet_mesh_object \
    import MultiInlet1DMeshObject
from Algorithms.DT_1D_V4.models.single_phase_multi_species_multi_inlet_nonreactive_pipe.\
        single_phase_multi_species_nonreactive_pipe_cell \
            import SinglePhaseMultiSpeciesNonReactivePipeCell
from Algorithms.DT_1D_V4.models.interface_models.\
        single_phase_uniform_massf_constant_area_interface \
            import SinglePhaseUniformMassfConstantAreaInterface
class SinglePhaseMultiSpeciesMultiInletNonReactivePipe(MultiInlet1DMeshObject):
    def __init__(self, n_cells, n_inlets, bulk_geometry, inlet_areas, init_flow_state, \
                    comp_label, limiter, recon_scheme, recon_props, update_from, \
                    flux_scheme,  reverse_direction_for_mirrored_flow = False) -> None:
        super().__init__(n_cells = n_cells, n_inlets = n_inlets, reverse_direction_for_mirrored_flow = reverse_direction_for_mirrored_flow)
        self.component_labels = [comp_label]
        area_tol = 1e-6
        if abs(sum(inlet_areas) - 0.25 * np.pi * bulk_geometry[0] ** 2.0) > area_tol:
            print("Inlet areas don't sum to the bulk pipe area")
        """
        bulk_geometry = [D, L]
        inlet_areas = [A_1, A_2, ..., A_n]
        """
        self.initialise_cells(geometry = bulk_geometry, \
                                init_flow_state = init_flow_state, \
                                comp_label = comp_label, \
                                n_cells = n_cells)
        
        self.initialise_interfaces(geometry = bulk_geometry, \
                                    limiter = limiter, \
                                    recon_scheme = recon_scheme, \
                                    n_cells = n_cells, \
                                    recon_props = recon_props, \
                                    update_from = update_from, \
                                    flux_scheme = flux_scheme, \
                                    init_flow_state = init_flow_state, \
                                    inlet_areas = inlet_areas, \
                                    reverse_direction = reverse_direction_for_mirrored_flow)
    
    def initialise_cells(self, geometry, init_flow_state, comp_label, n_cells):
        [D, L] = geometry
        for cell in range(n_cells):
            flow_state_object = init_flow_state.__class__
            fluid_state_object = init_flow_state.fluid_state.__class__
            fluid_model_filename = init_flow_state.fluid_state.gmodel.file_name
            fluid_model_object = init_flow_state.fluid_state.gmodel.__class__
            gm = fluid_model_object(fluid_model_filename)
            gs = fluid_state_object(gm)
            fs = flow_state_object(gs)

            cell_object = SinglePhaseMultiSpeciesNonReactivePipeCell(cell_id = cell, \
                                                                    label = comp_label)
            geo = {
                "dx"    :   L / n_cells,
                "dV"    :   0.25 * np.pi * D ** 2 * L / n_cells,
                "A_c"   :   0.25 * np.pi * D ** 2,
                "A_s"   :   np.pi * D * L / n_cells,
                "pos_x" :   (0.5 + cell) * L / n_cells
            }
            cell_object.fill_geometry(geometry = geo)
            cell_object.flow_state = fs
            cell_object.flow_state.fluid_state.copy_values(init_flow_state.fluid_state)
            cell_object.flow_state.vel_x = init_flow_state.vel_x
            cell_object.initialise_conserved_quantities()            
            self.cell_array[cell] = cell_object

    def initialise_interfaces(self, geometry, limiter, recon_scheme, n_cells, recon_props, \
                                update_from, flux_scheme, init_flow_state, inlet_areas, reverse_direction):
        [D, _] = geometry
        inlet_interface_array = [None] * len(inlet_areas)
        bulk_interface_array = [None] * n_cells
        for inlet_ind, area in enumerate(inlet_areas): #Doing the inlets that are all on one side of the block
            flow_state_object = init_flow_state.__class__
            fluid_state_object = init_flow_state.fluid_state.__class__
            fluid_model_filename = init_flow_state.fluid_state.gmodel.file_name
            fluid_model_object = init_flow_state.fluid_state.gmodel.__class__
            gm_lft = fluid_model_object(fluid_model_filename)
            gm__rght = fluid_model_object(fluid_model_filename)
            gs_lft = fluid_state_object(gm_lft)
            gs_rght = fluid_state_object(gm__rght)
            fs_lft = flow_state_object(gs_lft)
            fs_rght = flow_state_object(gs_rght)
            if reverse_direction: #Inlets are on the east boundary
                interface_id = n_cells + inlet_ind
                
            else:
                interface_id = inlet_ind

            interface_object = SinglePhaseUniformMassfConstantAreaInterface(\
                                    interface_id = interface_id, \
                                    flux_scheme = flux_scheme, \
                                    recon_scheme = recon_scheme, \
                                    limiter = limiter, \
                                    recon_props = recon_props, \
                                    update_from = update_from)
            geo = {"A"  : area}
            interface_object.fill_geometry(geometry = geo)
            interface_object.lft_state = fs_lft
            interface_object.rght_state = fs_rght
            inlet_interface_array[inlet_ind] = interface_object

        for interface in range(n_cells):
            flow_state_object = init_flow_state.__class__
            fluid_state_object = init_flow_state.fluid_state.__class__
            fluid_model_filename = init_flow_state.fluid_state.gmodel.file_name
            fluid_model_object = init_flow_state.fluid_state.gmodel.__class__
            gm_lft = fluid_model_object(fluid_model_filename)
            gm__rght = fluid_model_object(fluid_model_filename)
            gs_lft = fluid_state_object(gm_lft)
            gs_rght = fluid_state_object(gm__rght)
            fs_lft = flow_state_object(gs_lft)
            fs_rght = flow_state_object(gs_rght)
            if reverse_direction: #Inlets are on the east boundary
                interface_id = interface
            else:
                interface_id = len(inlet_areas) + interface

            interface_object = SinglePhaseUniformMassfConstantAreaInterface(\
                                    interface_id = interface_id, \
                                    flux_scheme = flux_scheme, \
                                    recon_scheme = recon_scheme, \
                                    limiter = limiter, \
                                    recon_props = recon_props, \
                                    update_from = update_from)
            geo = {"A"  : 0.25 * np.pi * D ** 2}
            interface_object.fill_geometry(geometry = geo)
            interface_object.lft_state = fs_lft
            interface_object.rght_state = fs_rght
            bulk_interface_array[interface] = interface_object

        if reverse_direction:
            self.interface_array = bulk_interface_array + inlet_interface_array
        else:
            self.interface_array = inlet_interface_array + bulk_interface_array

    def add_boundary_conditions(self, BC):
        self.boundary_conditions.append(BC)
