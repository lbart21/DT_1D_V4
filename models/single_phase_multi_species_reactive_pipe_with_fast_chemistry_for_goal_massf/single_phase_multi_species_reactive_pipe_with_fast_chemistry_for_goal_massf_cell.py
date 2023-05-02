"""
Function:
Author: Luke Bartholomew
Edits:
"""
import numpy as np

from Algorithms.DT_1D_V4.models.cell_methods.single_phase_multi_species_decoding \
        import multi_species_decode_to_primative_properties

from Algorithms.DT_1D_V4.models.cell_methods.single_phase_multi_species_encode_cqs \
        import encode_multi_species_cqs
class SinglePhaseMultiSpeciesReactivePipeWithFastChemistryForGoalMassfCell():
    def __init__(self, cell_id, label) -> None:
        self.geo = {}
        self.flow_state = None
        self.cell_id = cell_id
        self.label = label
        self.phase = "Single"
        self.interior_cell_flag = True
        self.combustion_parameters = [] # Will be [sp_name, f_goal, x, reaction]

        self.cqs = {
            "mass"      : 0.0,
            "xMom"      : 0.0,
            "energy"    : 0.0,
            "spcs_mass" : []
        }

    def fill_geometry(self, geometry):
        self.geo = geometry
    
    def max_allowable_dt(self, cfl):
        return cfl * self.geo["dx"] / (abs(self.flow_state.vel_x) + self.flow_state.fluid_state.a)

    def initialise_conserved_quantities(self):
        self.cqs = encode_multi_species_cqs(flow_state = self.flow_state)

    def decode_to_primative_properties(self):
        if self.interior_cell_flag:
            rho, u, massf, vel_x = multi_species_decode_to_primative_properties(cqs = self.cqs)

            self.flow_state.fluid_state.rho = rho
            self.flow_state.fluid_state.u = u
            self.flow_state.fluid_state.massf = massf

            self.flow_state.vel_x = vel_x

            self.flow_state.fluid_state.update_thermo_from_rhou()
            
    
    def complete_cell_methods(self, **kwargs):
        if self.interior_cell_flag:
            [sp_name, f_goal, x, reaction] = self.combustion_parameters
            reactant_stoichimetric_coefficients = reaction.reactants_stoichiometric_coefficients
            product_stoichimetric_coefficients = reaction.products_stoichiometric_coefficients
            # Find index of sp_name in gmodel.species_names
            species_names = self.flow_state.fluid_state.gmodel.species_names
            species_molar_masses = self.flow_state.fluid_state.gmodel.mol_masses
            current_massf = self.flow_state.fluid_state.massf
            new_massf = self.flow_state.fluid_state.massf
            sp_ind = species_names.index(sp_name)
            delta_f = (f_goal - current_massf[sp_ind])
            for ind, species in enumerate(species_names):
                if species != sp_name:
                    alpha_i = reactant_stoichimetric_coefficients[species]# For current species in loop
                    M_i = species_molar_masses[ind]
                    beta_j = product_stoichimetric_coefficients[sp_name]# For specified species
                    M_j = species_molar_masses[sp_ind]
                    x_trimmed = min(beta_j * M_j * current_massf[ind] / (alpha_i * M_i * (delta_f + 1e-10)), x)
                    if x_trimmed != x:
                        print("Insufficient amount of species:", species)
                        x = x_trimmed
            for ind, species in enumerate(species_names):
                if species == sp_name:
                    new_massf[ind] = current_massf[ind] + delta_f * x
                else:
                    alpha_i = reactant_stoichimetric_coefficients[species]# For current species in loop
                    M_i = species_molar_masses[ind]
                    beta_j = product_stoichimetric_coefficients[sp_name]# For specified species
                    M_j = species_molar_masses[sp_ind]
                    new_massf[ind] = current_massf[ind] - x * delta_f * (alpha_i * M_i) / (beta_j * M_j)
            self.flow_state.fluid_state.massf = new_massf
            self.reinitialise_massf_conserved_quantity()

    def decode_to_primative_properties_after_interface_methods(self):
        self.decode_to_primative_properties()

    def decode_to_primative_properties_after_cell_methods(self):
        self.decode_to_primative_properties()

    def reinitialise_massf_conserved_quantity(self):
        self.cqs["spcs_mass"] = (self.flow_state.fluid_state.rho * np.array(self.flow_state.fluid_state.massf)).tolist()