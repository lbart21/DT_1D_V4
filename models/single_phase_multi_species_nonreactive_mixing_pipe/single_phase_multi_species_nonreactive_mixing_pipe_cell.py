"""
Function:
Author: Luke Bartholomew
Edits:
"""
import numpy as np
class SinglePhaseMultiSpeciesNonReactiveMixingPipeCell():
    def __init__(self, cell_id, label) -> None:
        self.geo = {}
        self.flow_state = None
        self.reactor_model = None
        self.cell_id = cell_id
        self.label = label
        self.phase = "Single"
        self.interior_cell_flag = True

        self.dt_suggest = 1e-11

        self.cqs = {
            "mass"      : 0.0,
            "xMom"      : 0.0,
            "energy"    : 0.0,
            "spcs_mass" : []
        }
    
    def fill_geometry(self, geometry):
        self.geo = geometry
    
    def decode_to_primative_properties(self):
        if self.interior_cell_flag:
            rho = sum(self.cqs["spcs_mass"])
            rho_tol = 1e-6
            if abs(rho - self.cqs["mass"]) > rho_tol:
                #print("Too large of an error between conserved quantity mass and sum of species mass")
                #print("Cell: ", self.cell_id, "Error: ", abs(rho - self.cqs["mass"]))
                self.cqs["spcs_mass"] = (np.array(self.cqs["spcs_mass"]) * self.cqs["mass"] / rho).tolist()

            massf = (np.array(self.cqs["spcs_mass"]) / self.cqs["mass"]).tolist()
            massf = (np.array(massf) / sum(massf)).tolist()
            self.flow_state.fluid_state.massf = massf
            self.flow_state.fluid_state.rho = self.cqs["mass"]
            vel_x = self.cqs["xMom"] / self.cqs["mass"]
            self.flow_state.vel_x = vel_x
            self.flow_state.fluid_state.u = self.cqs["energy"] / self.cqs["mass"] - 0.5 * vel_x ** 2.0
            #print("Printing fluid state props before prop update, after decoding")
            #print("rho: ", self.flow_state.fluid_state.rho, "p: ", self.flow_state.fluid_state.p, \
                        #"T: ", self.flow_state.fluid_state.T, "massf: ", self.flow_state.fluid_state.massf)
            self.flow_state.fluid_state.update_thermo_from_rhou()
            #print("Printing fluid state props after prop update, after decoding")
            #print("rho: ", self.flow_state.fluid_state.rho, "p: ", self.flow_state.fluid_state.p, \
                        #"T: ", self.flow_state.fluid_state.T, "massf: ", self.flow_state.fluid_state.massf)
            self.flow_state.fluid_state.update_sound_speed()
        
    def initialise_conserved_quantities(self):
        self.cqs["mass"] = self.flow_state.fluid_state.rho
        self.cqs["xMom"] = self.flow_state.fluid_state.rho * self.flow_state.vel_x
        self.cqs["energy"] = self.flow_state.fluid_state.rho * \
                                                    (self.flow_state.fluid_state.u + 0.5 * self.flow_state.vel_x ** 2)
        self.cqs["spcs_mass"] = (self.flow_state.fluid_state.rho * np.array(self.flow_state.fluid_state.massf)).tolist()

    def max_allowable_dt(self, cfl):
        return cfl * self.geo["dx"] / (abs(self.flow_state.vel_x) + self.flow_state.fluid_state.a)

    def complete_cell_methods(self, **kwargs):
        pass
        
    def decode_to_primative_properties_after_interface_methods(self):
        self.decode_to_primative_properties()

    def decode_to_primative_properties_after_cell_methods(self):
        pass
        