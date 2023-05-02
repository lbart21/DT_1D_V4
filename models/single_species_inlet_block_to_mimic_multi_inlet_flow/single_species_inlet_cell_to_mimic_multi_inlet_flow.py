"""
Function:
Author: Luke Bartholomew
Edits:
"""

from Algorithms.DT_1D_V4.models.cell_methods.single_phase_single_species_decoding \
        import single_species_decode_to_primative_properties

from Algorithms.DT_1D_V4.models.cell_methods.single_phase_single_species_encode_cqs \
        import encode_single_species_cqs

class SingleSpeciesCellToMimicMultiInletFlow():
    def __init__(self, cell_id, label, source_fs, source_area) -> None:
        self.geo = {}
        self.flow_state = None
        self.cell_id = cell_id
        self.label = label
        self.phase = "Single"
        self.interior_cell_flag = True

        self.source_fs = source_fs
        self.source_area = source_area

        self.cqs = {
            "mass"      : 0.0,
            "xMom"      : 0.0,
            "energy"    : 0.0
        }

        self.flow_state = None
    
    def fill_geometry(self, geometry):
        self.geo = geometry
    
    def decode_to_primative_properties(self):
        if self.interior_cell_flag:
            rho, u, vel_x = single_species_decode_to_primative_properties(cqs = self.cqs)
            
            self.flow_state.vel_x = vel_x
            self.flow_state.fluid_state.rho = rho
            self.flow_state.fluid_state.u = u
            self.flow_state.fluid_state.update_thermo_from_rhou()
    
    def initialise_conserved_quantities(self):
        self.cqs = encode_single_species_cqs(flow_state = self.flow_state)
        
    def max_allowable_dt(self, cfl):
        return cfl * self.geo["dx"] / (abs(self.flow_state.vel_x) + self.flow_state.fluid_state.a)
    
    def complete_cell_methods(self, **kwargs):
        interface_array = kwargs["interface_array"]
        map_cell_id_to_west_interface_idx = kwargs["map_cell_id_to_west_interface_idx"]
        dt_inv = kwargs["dt_inv"]

        if self.interior_cell_flag:
            #print("Start of cell methods within cell object")
            dV = self.geo["dV"]
            A_inlet = self.source_area
            rho_in = self.source_fs.fluid_state.rho
            vel_x_in = self.source_fs.vel_x
            u_in = self.source_fs.fluid_state.u
            p_in = self.source_fs.fluid_state.p
            west_interface = interface_array[map_cell_id_to_west_interface_idx[self.cell_id][0]]

            self.cqs["mass"] += dt_inv * A_inlet * rho_in * vel_x_in / dV
            self.cqs["energy"] += dt_inv * A_inlet * rho_in * vel_x_in * (u_in + p_in / rho_in + 0.5 * vel_x_in ** 2) / dV
            self.cqs["xMom"] -= dt_inv * A_inlet * west_interface.boundary_fluxes["p"] / dV
            self.cqs["xMom"] += dt_inv * A_inlet * (p_in + rho_in * vel_x_in ** 2) / dV

        else:
            pass
        
    def decode_to_primative_properties_after_interface_methods(self):
        self.decode_to_primative_properties()

    def decode_to_primative_properties_after_cell_methods(self):
        self.decode_to_primative_properties()