"""
Function:
Author: Luke Bartholomew
Edits:
"""

import numpy as np

from Algorithms.DT_1D_V4.models.cell_methods.goal_massf_solver \
            import form_goal_massf_solver_matrices, goal_massf_solver
from Algorithms.DT_1D_V4.models.cell_methods.bulk_reaction_fast_chemistry_method \
            import bulk_reaction_fast_chemsitry
from Algorithms.DT_1D_V4.models.cell_methods.single_phase_multi_species_decoding \
            import multi_species_decode_to_primative_properties
from Algorithms.DT_1D_V4.models.cell_methods.single_phase_multi_species_encode_cqs \
            import encode_multi_species_cqs

class SinglePhaseMultiSpeciesBulkChemistryForGoalMassfUniformCell():
    def __init__(self, cell_id, label) -> None:
        self.geo = {}
        self.flow_state = None
        self.reactor_model = None
        self.cell_id = cell_id
        self.label = label
        self.phase = "Single"
        self.interior_cell_flag = True

        self.formed_goal_massf_matrices = False
        self.LHS_A = None
        self.RHS_A = None
        self.goal_massf_source = None

        self.bulk_reaction_parameters = []


        self.dt_suggest = 1e-11

        self.cqs = {
            "mass"      : 0.0,
            "xMom"      : 0.0,
            "energy"    : 0.0,
            "spcs_mass" : []
        }

    def complete_cell_methods(self, **kwargs):
        if self.interior_cell_flag:
            
            if not self.formed_goal_massf_matrices:
                self.LHS_A, self.RHS_A, self.goal_massf_source = \
                    form_goal_massf_solver_matrices(massf_current = self.flow_state.fluid_state.massf, \
                                                    bulk_reaction_parameters = self.bulk_reaction_parameters, \
                                                    molar_masses = self.flow_state.fluid_state.gmodel.mol_masses, \
                                                    species_names = self.flow_state.fluid_state.gmodel.species_names)

                self.formed_goal_massf_matrices = True
            
            goal_massf = goal_massf_solver(massf_current = self.flow_state.fluid_state.massf, \
                                            bulk_reaction_parameters = self.bulk_reaction_parameters, \
                                            species_names = self.flow_state.fluid_state.gmodel.species_names, \
                                            LHS_A = self.LHS_A, RHS_A = self.RHS_A, source = self.goal_massf_source)
            
            new_massf = bulk_reaction_fast_chemsitry(massf_current = self.flow_state.fluid_state.massf, \
                                                    x = self.bulk_reaction_parameters[0], \
                                                    massf_goal = goal_massf)
            
            self.flow_state.fluid_state.massf = new_massf
            self.reinitialise_massf_conserved_quantity()

    def fill_geometry(self, geometry):
        self.geo = geometry

    def decode_to_primative_properties(self):
        if self.interior_cell_flag:
            rho, u, massf, vel_x = multi_species_decode_to_primative_properties(cqs = self.cqs)

            self.flow_state.fluid_state.rho = rho
            self.flow_state.fluid_state.u = u
            self.flow_state.fluid_state.massf = massf
            self.flow_state.vel_x = vel_x

            self.flow_state.fluid_state.update_thermo_from_rhou()
    
    def initialise_conserved_quantities(self):
        self.cqs = encode_multi_species_cqs(flow_state = self.flow_state)

    def max_allowable_dt(self, cfl):
        return cfl * self.geo["dx"] / (abs(self.flow_state.vel_x) + self.flow_state.fluid_state.a)

    def decode_to_primative_properties_after_interface_methods(self):
        self.decode_to_primative_properties()

    def decode_to_primative_properties_after_cell_methods(self):
        self.decode_to_primative_properties()
    
    def reinitialise_massf_conserved_quantity(self):
        self.cqs["spcs_mass"] = (self.flow_state.fluid_state.rho * np.array(self.flow_state.fluid_state.massf)).tolist()