"""
Function:
Author: Luke Bartholomew
Edits:
"""
from Algorithms.DT_1D_V4.models.prefilled_single_inlet_mesh_object import SingleInlet1DMeshObject
from Algorithms.DT_1D_V4.models.single_phase_multi_species_nonreactive_pipe.\
        single_phase_multi_species_nonreactive_pipe_cell \
            import SinglePhaseMultiSpeciesNonReactivePipeCell
from Algorithms.DT_1D_V4.models.interface_models.single_phase_uniform_massf_interface \
            import SinglePhaseUniformMassfInterface

import numpy as np
class SinglePhaseMultiSpeciesNonReactivePipe(SingleInlet1DMeshObject):
    def __init__(self, n_cells, geometry, init_flow_state, \
                recon_props, recon_scheme, limiter, update_from, comp_label,\
                flux_scheme, reverse_direction_for_ghost_cells = False):
        super().__init__(n_cells = n_cells, reverse_direction_for_ghost_cells = reverse_direction_for_ghost_cells)
        """
        geometry = [D, L]
        """
        self.component_labels = [comp_label]
        self.initialise_cells(geometry = geometry, \
                                init_flow_state = init_flow_state, \
                                comp_label = comp_label, \
                                n_cells = n_cells)

        self.initialise_interfaces(geometry = geometry, \
                                    limiter = limiter, \
                                    recon_scheme = recon_scheme, \
                                    n_cells = n_cells, \
                                    recon_props = recon_props, \
                                    update_from = update_from, \
                                    flux_scheme = flux_scheme, \
                                    init_flow_state = init_flow_state)

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
        

    def initialise_interfaces(self, geometry, limiter, recon_scheme, n_cells, \
                                    recon_props, update_from, flux_scheme, \
                                    init_flow_state):
        [D, _] = geometry
        for interface in range(n_cells + 1):
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

            interface_object = SinglePhaseUniformMassfInterface(\
                                    interface_id = interface, \
                                    flux_scheme = flux_scheme, \
                                    recon_scheme = recon_scheme, \
                                    limiter = limiter, \
                                    recon_props = recon_props, \
                                    update_from = update_from)
            geo = {"A"  : 0.25 * np.pi * D ** 2}
            interface_object.fill_geometry(geometry = geo)
            interface_object.lft_state = fs_lft
            interface_object.rght_state = fs_rght
            self.interface_array[interface] = interface_object

    def add_boundary_conditions(self, BC):
        self.boundary_conditions.append(BC)
