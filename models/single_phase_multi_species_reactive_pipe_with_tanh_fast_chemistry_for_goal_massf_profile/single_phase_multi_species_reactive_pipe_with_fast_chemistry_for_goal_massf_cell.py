"""
Function:
Author: Luke Bartholomew
Edits:
"""
import numpy as np

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
        self.cqs["mass"] = self.flow_state.fluid_state.rho
        self.cqs["xMom"] = self.flow_state.fluid_state.rho * self.flow_state.vel_x
        self.cqs["energy"] = self.flow_state.fluid_state.rho * \
                                                    (self.flow_state.fluid_state.u + 0.5 * self.flow_state.vel_x ** 2)
        self.cqs["spcs_mass"] = (self.flow_state.fluid_state.rho * np.array(self.flow_state.fluid_state.massf)).tolist()

    def decode_to_primative_properties(self):
        if self.interior_cell_flag:
            rho = sum(self.cqs["spcs_mass"])
            rho_tol = 1e-6
            if abs(rho - self.cqs["mass"]) > rho_tol:
                print("Too large of an error between conserved quantity mass and sum of species mass")
            massf = (np.array(self.cqs["spcs_mass"]) / rho).tolist()
            self.flow_state.fluid_state.massf = massf
            self.flow_state.fluid_state.rho = rho
            vel_x = self.cqs["xMom"] / rho
            self.flow_state.vel_x = vel_x
            self.flow_state.fluid_state.u = self.cqs["energy"] / rho - 0.5 * vel_x ** 2.0
            #print("Printing fluid state props before prop update, after decoding")
            #print("rho: ", self.flow_state.fluid_state.rho, "p: ", self.flow_state.fluid_state.p, \
                        #"T: ", self.flow_state.fluid_state.T, "massf: ", self.flow_state.fluid_state.massf)
            self.flow_state.fluid_state.update_thermo_from_rhou()
            #print("Printing fluid state props after prop update, after decoding")
            #print("rho: ", self.flow_state.fluid_state.rho, "p: ", self.flow_state.fluid_state.p, \
                        #"T: ", self.flow_state.fluid_state.T, "massf: ", self.flow_state.fluid_state.massf)
            self.flow_state.fluid_state.update_sound_speed()
    
    def complete_cell_methods(self, **kwargs):
        if self.interior_cell_flag:
            [sp_name, f_goal, x, reaction] = self.combustion_parameters
            #print(self.cell_id, f_goal)
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
                    if delta_f == 0.0:
                        x = 0.0
                    elif delta_f > 0.0:
                        x_trimmed = min(beta_j * M_j * current_massf[ind] / (alpha_i * M_i * delta_f), x)
                        if x_trimmed != x:
                            print("Cell ID:", self.cell_id, "Insufficient amount of species:", \
                                    sp_name, " for reaction. Full step would result in mass fraction > 1.0")
                            x = x_trimmed

                    elif delta_f < 0.0:
                        x_trimmed = min((current_massf[ind] - 1.0) * beta_j * M_j / (alpha_i * M_i * delta_f), x)
                        if x_trimmed != x:
                            print("Cell ID:", self.cell_id, "Insufficient amount of species:", \
                                    species, " for reaction. Full step would result in negative mass fraction")
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