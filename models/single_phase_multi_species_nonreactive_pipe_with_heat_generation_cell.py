"""
Function:
Author: Luke Bartholomew
Edits:
"""

class SinglePhaseMultiSpeciesNonReactivePipeWithHeatGenCell():
    def __init__(self, cell_id, label, heat_gen) -> None:
        self.geo = {}
        self.flow_state = None
        self.cell_id = cell_id
        self.label = label
        self.phase = "Single"
        self.interior_cell_flag = True

        self.heat_gen = heat_gen

        self.cqs = {
            "mass"      : 0.0,
            "xMom"      : 0.0,
            "energy"    : 0.0
        }

        self.flow_state = None
    
    def fill_geometry(self, geometry):
        self.geo = geometry
    
    def decode_to_primative_properties(self):
        new_density = self.cqs["mass"] 
        new_vel_x = self.cqs["xMom"] / self.cqs["mass"]
        new_u = self.cqs["energy"] / self.cqs["mass"] \
                                                        - 0.5 * new_vel_x ** 2.0
        self.flow_state.vel_x = new_vel_x
        self.flow_state.fluid_state.rho = new_density
        self.flow_state.fluid_state.u = new_u
        self.flow_state.fluid_state.update_thermo_from_rhou()
        self.flow_state.fluid_state.update_sound_speed()
    
    def initialise_conserved_quantities(self):
        self.cqs["mass"] = self.flow_state.fluid_state.rho
        self.cqs["xMom"] = self.flow_state.fluid_state.rho * self.flow_state.vel_x
        self.cqs["energy"] = self.flow_state.fluid_state.rho * \
                                                    (self.flow_state.fluid_state.u + 0.5 * self.flow_state.vel_x ** 2)
    
    def max_allowable_dt(self, cfl):
        return cfl * self.geo["dx"] / (abs(self.flow_state.vel_x) + self.flow_state.fluid_state.a)
    
    def complete_cell_methods(self, **kwargs):
        dt_inv = kwargs["dt_inv"]
        self.cqs["energy"] += dt_inv * self.heat_gen
    
    def decode_to_primative_properties_after_interface_methods(self):
        self.decode_to_primative_properties()

    def decode_to_primative_properties_after_cell_methods(self):
        self.decode_to_primative_properties()