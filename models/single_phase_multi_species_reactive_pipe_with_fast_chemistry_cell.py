"""
Function:
Author: Luke Bartholomew
Edits:
"""
import numpy as np

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
