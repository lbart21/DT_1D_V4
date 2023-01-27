"""
Function:
Author: Luke Bartholomew
Edits:
"""
from prefilled_single_inlet_mesh_object import basic1DMeshObject
from Algorithms.DesignToolAlgorithmV4_1D.ComponentModels.inlet_cell_to_mimic_multi_inlet_flow import CellToMimicMultiInletFlow
from Algorithms.DesignToolAlgorithmV4_1D.ComponentModels.single_phase_multi_species_nonreactive_interface import SinglePhaseMultiSpeciesNonReactiveInterface

import numpy as np
class InletBlockToMimicMultiInletFlow(basic1DMeshObject):
    def __init__(self, inlet_flow_state, init_flow_state, geometry, inlet_area, comp_label, \
                        reconstruction_scheme, limiter, reconstruction_properties, update_from, flux_scheme) -> None:
        super().__init__(nCells = 1, reversed = False)
        self.component_labels = [comp_label]
        
        self.initialise_cell(comp_label = comp_label, inlet_flow_state = inlet_flow_state, \
                            init_flow_state = init_flow_state, Geometry = geometry, inlet_area = inlet_area)

        self.initialise_interfaces(reconstruction_scheme = reconstruction_scheme, limiter = limiter, \
                            reconstruction_properties = reconstruction_properties, update_from = update_from, \
                            flux_scheme = flux_scheme, Geometry = geometry, init_flow_state = init_flow_state)

    def initialise_cell(self, comp_label, inlet_flow_state, init_flow_state, Geometry, inlet_area):
        cell = CellToMimicMultiInletFlow(cell_id = 0, label = comp_label, source_fs = inlet_flow_state, source_area = inlet_area)
        [D, L] = Geometry
        GEO = {
                "dx"    :   L,
                "dV"    :   0.25 * np.pi * D ** 2 * L,
                "A_c"   :   0.25 * np.pi * D ** 2,
                "A_s"   :   np.pi * D * L,
                "pos_x" :   0.5 * L
        }
        cell.fill_geometry(Geometry = GEO)
        cell.flow_state = init_flow_state
        cell.initialise_conserved_properties()
        self.cell_array[0] = cell

    def initialise_interfaces(self, reconstruction_scheme, limiter, reconstruction_properties, update_from, flux_scheme, Geometry, init_flow_state):
        [D, L] = Geometry
        for i in range(2):
            interface = SinglePhaseMultiSpeciesNonReactiveInterface(interface_id = i, nL = i, nR = 1 - i, \
                                        reconstruction_scheme = reconstruction_scheme, limiter = limiter, \
                                        reconstruction_properties = reconstruction_properties, \
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

            GEO = {"A"  : 0.25 * np.pi * D ** 2}
            interface.fill_geometry(geometry = GEO)
            interface.lft_state = fs_lft
            interface.rght_state = fs_rght
            self.interface_array[i] = interface

    def add_boundary_conditions(self, BC):
        self.boundary_conditions.append(BC)
