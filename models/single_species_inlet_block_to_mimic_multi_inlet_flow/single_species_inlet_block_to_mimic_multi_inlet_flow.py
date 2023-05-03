"""
Function:
Author: Luke Bartholomew
Edits:
"""
from Algorithms.DT_1D_V4.models.prefilled_single_inlet_mesh_object import SingleInlet1DMeshObject
from Algorithms.DT_1D_V4.models.single_species_inlet_block_to_mimic_multi_inlet_flow.single_species_inlet_cell_to_mimic_multi_inlet_flow \
                import SingleSpeciesCellToMimicMultiInletFlow
from Algorithms.DT_1D_V4.models.interface_models.single_phase_uniform_massf_interface \
                import SinglePhaseUniformMassfInterface

import numpy as np
class SingleSpeciesInletBlockToMimicMultiInletFlow(SingleInlet1DMeshObject):
    def __init__(self, inlet_flow_state, init_flow_state, geometry, inlet_area, comp_label, \
                        recon_scheme, limiter, recon_props, update_from, flux_scheme) -> None:
        super().__init__(n_cells = 1, reverse_direction_for_ghost_cells = False)
        self.component_labels = [comp_label]
        
        self.initialise_cell(comp_label = comp_label, inlet_flow_state = inlet_flow_state, \
                            init_flow_state = init_flow_state, geometry = geometry, inlet_area = inlet_area)

        self.initialise_interfaces(recon_scheme = recon_scheme, limiter = limiter, \
                            recon_props = recon_props, update_from = update_from, \
                            flux_scheme = flux_scheme, geometry = geometry, init_flow_state = init_flow_state)

    def initialise_cell(self, comp_label, inlet_flow_state, init_flow_state, geometry, inlet_area):
        cell = SingleSpeciesCellToMimicMultiInletFlow(cell_id = 0, label = comp_label, \
                                            source_fs = inlet_flow_state, source_area = inlet_area)
        [D, L] = geometry
        geo = {
                "dx"    :   L,
                "dV"    :   0.25 * np.pi * D ** 2 * L,
                "A_c"   :   0.25 * np.pi * D ** 2,
                "A_s"   :   np.pi * D * L,
                "pos_x" :   0.5 * L
        }
        cell.fill_geometry(geometry = geo)
        cell.flow_state = init_flow_state
        cell.initialise_conserved_quantities()
        self.cell_array[0] = cell

    def initialise_interfaces(self, recon_scheme, limiter, recon_props, update_from, flux_scheme, geometry, init_flow_state):
        [D, _] = geometry
        for i in range(2):
            interface = SinglePhaseUniformMassfInterface(interface_id = i, \
                                        recon_scheme = recon_scheme, limiter = limiter, \
                                        recon_props = recon_props, \
                                        update_from = update_from, flux_scheme = flux_scheme)
            
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

            geo = {"A"  : 0.25 * np.pi * D ** 2}
            interface.fill_geometry(geometry = geo)
            interface.lft_state = fs_lft
            interface.rght_state = fs_rght
            self.interface_array[i] = interface

    def add_boundary_conditions(self, BC):
        self.boundary_conditions.append(BC)
