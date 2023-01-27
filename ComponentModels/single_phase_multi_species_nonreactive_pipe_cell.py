"""
Function:
Author: Luke Bartholomew
Edits:
"""
class SinglePhaseMultiSpeciesNonReactivePipeCell():
    def __init__(self, cell_id, label) -> None:
        self.GEO = {}
        self.flow_state = None
        self.cell_id = cell_id
        self.label = label
        self.phase = "Single"
        self.interior_cell_flag = True

        self.cqs = {
            "mass"      : 0.0,
            "xMom"      : 0.0,
            "energy"    : 0.0
        }

        self.flow_state = None
    
    def fill_geometry(self, Geometry):
        self.GEO = Geometry
    
    def update_primative_properties(self):
        new_density = self.cqs["mass"] 
        new_vel_x = self.cqs["xMom"] / self.cqs["mass"]
        new_u = self.cqs["energy"] / self.cqs["mass"] \
                                                        - 0.5 * new_vel_x ** 2.0
        self.flow_state.vel_x = new_vel_x
        self.flow_state.fluid_state.rho = new_density
        self.flow_state.fluid_state.u = new_u
        self.flow_state.fluid_state.update_thermo_from_rhou()
        self.flow_state.fluid_state.update_sound_speed()
    
    def initialise_conserved_properties(self):
        self.cqs["mass"] = self.flow_state.fluid_state.rho
        self.cqs["xMom"] = self.flow_state.fluid_state.rho * self.flow_state.vel_x
        self.cqs["energy"] = self.flow_state.fluid_state.rho * \
                                                    (self.flow_state.fluid_state.u + 0.5 * self.flow_state.vel_x ** 2)

    def max_allowable_dt(self, cfl):
        return cfl * self.GEO["dx"] / (abs(self.flow_state.vel_x) + self.flow_state.fluid_state.a)
    
    def update_properties(self):
        self.update_primative_properties()

    def complete_cell_methods(self, cell_array, interface_array, map_cell_id_to_west_interface_idx, \
                                map_cell_id_to_east_interface_idx, map_interface_id_to_west_cell_idx, \
                                map_interface_id_to_east_cell_idx, dt_inv):
        pass