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

class SinglePhaseMultiSpeciesReactivePipeWithFastChemistryCell():
    def __init__(self, cell_id, label) -> None:
        self.geo = {}
        self.flow_state = None
        self.cell_id = cell_id
        self.label = label
        self.phase = "Single"
        self.interior_cell_flag = True
        self.fast_chemistry_fraction = 0.0 
        self.fast_chemistry_gs = None

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
    
    def decode_to_primative_properties_after_interface_methods(self):
        self.decode_to_primative_properties()

    def decode_to_primative_properties_after_cell_methods(self):
        self.decode_to_primative_properties()

    def reinitialise_massf_conserved_quantity(self):
        self.cqs["spcs_mass"] = (self.flow_state.fluid_state.rho * np.array(self.flow_state.fluid_state.massf)).tolist()

    def complete_cell_methods(self, **kwargs):
        if self.interior_cell_flag:
            self.fast_chemistry_gs.rho = self.flow_state.fluid_state.rho
            self.fast_chemistry_gs.u = self.flow_state.fluid_state.u
            self.fast_chemistry_gs.update_thermo_from_rhou()
            #print(self.flow_state.fluid_state.gmodel.species_names)
            #print(list(self.fast_chemistry_gs.ceaSavedData.values()))
            equil_massf = np.array([list(self.fast_chemistry_gs.ceaSavedData["massf"].values())[ind] for ind in \
                           [list(self.fast_chemistry_gs.ceaSavedData["massf"].keys()).index(name) \
                                for name in self.flow_state.fluid_state.gmodel.species_names]])
            #print(equil_massf)
            current_massf = np.array(self.flow_state.fluid_state.massf)
            #print(current_massf)
            delta_massf = equil_massf - current_massf
            new_massf = (current_massf + self.fast_chemistry_fraction * delta_massf).tolist()
            self.flow_state.fluid_state.massf = new_massf
            
        self.reinitialise_massf_conserved_quantity()
