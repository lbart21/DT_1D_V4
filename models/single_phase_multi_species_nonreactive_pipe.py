"""
Function:
Author: Luke Bartholomew
Edits:
"""
from prefilled_single_inlet_mesh_object \
                                                import basic1DMeshObject
from Algorithms.DesignToolAlgorithmV4_1D.ComponentModels.single_phase_multi_species_nonreactive_pipe_cell \
                                        import SinglePhaseMultiSpeciesNonReactivePipeCell
from Algorithms.DesignToolAlgorithmV4_1D.ComponentModels.single_phase_multi_species_nonreactive_interface \
                import SinglePhaseMultiSpeciesNonReactiveInterface




import numpy as np
class SinglePhaseMultiSpeciesNonReactivePipe(basic1DMeshObject):
    def __init__(self, nCells, Geometry, init_flow_state, \
                RECON_PROPS, RECON_SCHEME, LIMITER, UPDATE_FROM, COMP_LABEL,\
                FLUX_SCHEME, reversed = False):
        super().__init__(nCells = nCells, reversed = reversed)
        """
        Geometry = [D, L]
        """
        self.component_labels = [COMP_LABEL]
        self.initialise_cells(Geometry = Geometry, \
                                init_flow_state = init_flow_state, \
                                COMP_LABEL = COMP_LABEL, \
                                nCells = nCells)

        self.initialise_interfaces(Geometry = Geometry, \
                                    LIMITER = LIMITER, \
                                    RECON_SCHEME = RECON_SCHEME, \
                                    nCells = nCells, \
                                    RECON_PROPS = RECON_PROPS, \
                                    UPDATE_FROM = UPDATE_FROM, \
                                    FLUX_SCHEME = FLUX_SCHEME, \
                                    init_flow_state = init_flow_state)

    def initialise_cells(self, Geometry, init_flow_state, COMP_LABEL, nCells):
        [D, L] = Geometry

        for cell in range(nCells):
            flowState_object = init_flow_state.__class__
            fluidState_object = init_flow_state.fluid_state.__class__
            fluidModel_filename = init_flow_state.fluid_state.gmodel.file_name
            fluidModel_object = init_flow_state.fluid_state.gmodel.__class__
            gm = fluidModel_object(fluidModel_filename)
            gs = fluidState_object(gm)
            fs = flowState_object(gs)

            cellObject = SinglePhaseMultiSpeciesNonReactivePipeCell(cell_id = cell, \
                                                            label = COMP_LABEL)
            GEO = {
                "dx"    :   L / nCells,
                "dV"    :   0.25 * np.pi * D ** 2 * L / nCells,
                "A_c"   :   0.25 * np.pi * D ** 2,
                "A_s"   :   np.pi * D * L / nCells,
                "pos_x" :   (0.5 + cell) * L / nCells
            }
            cellObject.fill_geometry(Geometry = GEO)
            cellObject.flow_state = fs
            cellObject.flow_state.fluid_state.copy_values(init_flow_state.fluid_state)
            cellObject.flow_state.vel_x = init_flow_state.vel_x
            cellObject.initialise_conserved_properties()            
            self.cell_array[cell] = cellObject
        

    def initialise_interfaces(self, Geometry, LIMITER, RECON_SCHEME, nCells, \
                                    RECON_PROPS, UPDATE_FROM, FLUX_SCHEME, \
                                    init_flow_state):
        [D, L] = Geometry
        for interface in range(nCells + 1):
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

            interface_object = SinglePhaseMultiSpeciesNonReactiveInterface(\
                                    interface_id = interface, nL = interface, \
                                    nR = nCells - interface, \
                                    flux_scheme = FLUX_SCHEME, \
                                    reconstruction_scheme = RECON_SCHEME, \
                                    limiter = LIMITER, \
                                    reconstruction_properties = RECON_PROPS, \
                                    update_from = UPDATE_FROM)
            GEO = {"A"  : 0.25 * np.pi * D ** 2}
            interface_object.fill_geometry(geometry = GEO)
            interface_object.lft_state = fs_lft
            interface_object.rght_state = fs_rght
            self.interface_array[interface] = interface_object

    def add_boundary_conditions(self, BC):
        self.boundary_conditions.append(BC)
