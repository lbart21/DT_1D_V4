"""
Function:
Author: Luke Bartholomew
Edits:
"""
import numpy as np

class SinglePhaseMultiSpeciesReactivePipeWithCombustionFractionCell():
    def __init__(self, cell_id, label) -> None:
        self.geo = {}
        self.flow_state = None
        self.cell_id = cell_id
        self.label = label
        self.phase = "Single"
        self.interior_cell_flag = True
        self.combustion_fraction = [] #When filled, will be [[x, "element", reaction], [y, "element", reaction], ...]

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
            species_names = self.flow_state.fluid_state.gmodel.species_names
            species_molar_masses = self.flow_state.fluid_state.gmodel.mol_masses

            for reaction in self.combustion_fraction:
                [x, chosen_molecule, reaction_object] = reaction
                reactant_stoichimetric_coefficients = reaction_object.reactants_stoichiometric_coefficients
                product_stoichimetric_coefficients = reaction_object.products_stoichiometric_coefficients
                old_massf = self.flow_state.fluid_state.massf
                #print("Original mass fractions:", old_massf)
                new_massf = self.flow_state.fluid_state.massf
                chosen_molecule_ind = species_names.index(chosen_molecule)
                C_0 = reactant_stoichimetric_coefficients[chosen_molecule]
                M_0 = species_molar_masses[chosen_molecule_ind]
                f_0 = old_massf[chosen_molecule_ind]
                for molecule in reaction_object.reactants_molecule_names:
                    cqs_ind = species_names.index(molecule)
                    C_i = reactant_stoichimetric_coefficients[molecule]
                    M_i = species_molar_masses[cqs_ind]
                    new_massf[cqs_ind] = new_massf[cqs_ind] - x * C_i / C_0 * M_i / M_0 * f_0
                #print(new_massf)
                for molecule in reaction_object.products_molecule_names:
                    cqs_ind = species_names.index(molecule)
                    C_i = product_stoichimetric_coefficients[molecule]
                    M_i = species_molar_masses[cqs_ind]
                    new_massf[cqs_ind] = new_massf[cqs_ind] + x * C_i / C_0 * M_i / M_0 * f_0
                #print(new_massf)
                #print(sum(new_massf))
                self.flow_state.fluid_state.massf = new_massf
        #print(self.flow_state.fluid_state.massf)
        self.reinitialise_massf_conserved_quantity()

    def decode_to_primative_properties_after_interface_methods(self):
        self.decode_to_primative_properties()

    def decode_to_primative_properties_after_cell_methods(self):
        self.decode_to_primative_properties()

    def reinitialise_massf_conserved_quantity(self):
        self.cqs["spcs_mass"] = (self.flow_state.fluid_state.rho * np.array(self.flow_state.fluid_state.massf)).tolist()