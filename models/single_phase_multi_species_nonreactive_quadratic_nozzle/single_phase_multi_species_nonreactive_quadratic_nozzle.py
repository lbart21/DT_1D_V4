"""
Function: 
Author: Luke Bartholomew
Edits: 
"""
import numpy as np

from Algorithms.DT_1D_V4.models.interface_models.single_phase_uniform_massf_interface \
            import SinglePhaseUniformMassfInterface
from Algorithms.DT_1D_V4.models.single_phase_multi_species_nonreactive_quadratic_nozzle.\
        single_phase_multi_species_nonreactive_quadratic_nozzle_cell \
            import SinglePhaseMultiSpeciesNonReactiveQuadraticNozzleCell
from Algorithms.DT_1D_V4.models.prefilled_single_inlet_mesh_object import SingleInlet1DMeshObject

class SinglePhaseMultiSpeciesNonReactiveQuadraticNozzle(SingleInlet1DMeshObject):
    def __init__(self, n_cells, geometry, init_flow_state, \
                recon_props, recon_scheme, limiter, update_from, comp_label,\
                flux_scheme, reverse_direction_for_ghost_cells = False) -> None:
        super().__init__(n_cells = n_cells, reverse_direction_for_ghost_cells = reverse_direction_for_ghost_cells)
        self.component_labels = [comp_label]
        """
        geometry = [D1, D2, L]
        """
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
        [_, _,  L] = geometry
        for cell in range(n_cells):
            flow_state_object = init_flow_state.__class__
            fluid_state_object = init_flow_state.fluid_state.__class__
            fluid_model_filename = init_flow_state.fluid_state.gmodel.file_name
            fluid_model_object = init_flow_state.fluid_state.gmodel.__class__
            gm = fluid_model_object(fluid_model_filename)
            gs = fluid_state_object(gm)
            fs = flow_state_object(gs)

            cell_object = SinglePhaseMultiSpeciesNonReactiveQuadraticNozzleCell(cell_id = cell, \
                                                                            label = comp_label)
            pos_x_w = cell * L / n_cells
            pos_x_c = (cell + 0.5) * L / n_cells
            pos_x_e = (cell + 1.0) * L / n_cells

            D_w = self.find_diameter_at_x(geometry = geometry, x = pos_x_w)
            D_c = self.find_diameter_at_x(geometry = geometry, x = pos_x_c)
            D_e = self.find_diameter_at_x(geometry = geometry, x = pos_x_e)

            dx = L / n_cells
            geo = {
                "dx"    :   dx,
                "dV"    :   np.pi * (D_w ** 2.0 + D_w * D_e + D_e ** 2.0) * dx / 12.0,
                "A_c"   :   0.25 * np.pi * D_c ** 2.0,
                "A_s"   :   0.25 * np.pi * (D_w + D_e) * (4.0 * dx ** 2.0 + (D_w - D_e) ** 2.0) ** 0.5,
                "pos_x" :   pos_x_c
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
        [_, _, L] = geometry
        for interface in range(n_cells + 1):
            D = self.find_diameter_at_x(geometry = geometry, x = interface * L / n_cells)
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

    def find_diameter_at_x(self, geometry, x):
        [D_1, D_2, L] = geometry
        return D_1 + (D_2 - D_1) * x / L + (D_2 - D_1) * x * (L - x) / L ** 2.0
    
    def add_boundary_conditions(self, BC):
        self.boundary_conditions.append(BC)