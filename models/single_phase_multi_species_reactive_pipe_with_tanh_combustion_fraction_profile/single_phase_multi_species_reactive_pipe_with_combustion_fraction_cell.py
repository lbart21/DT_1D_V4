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