"""
Function:
Author: Luke Bartholomew
Edits:
"""
import numpy as np
class SinglePhaseMultiSpeciesFixedFluxInterface():
    def __init__(self, interface_id, flow_state, on_west_boundary, multi_species) -> None:
        self.interface_id = interface_id
        self.flux_flag = True
        self.on_west_boundary = on_west_boundary
        self.flow_state = flow_state
        self.fluxes_been_calculated = False 
        self.multi_species = multi_species

        self.boundary_fluxes = {}
        
    def fill_geometry(self, geometry):
        self.geo = geometry

    def complete_interface_methods(self, **kwargs):
        cell_array = kwargs["cell_array"]
        
        map_interface_id_to_west_cell_idx = kwargs["map_interface_id_to_west_cell_idx"]
        map_interface_id_to_east_cell_idx = kwargs["map_interface_id_to_east_cell_idx"]

        dt_inv = kwargs["dt_inv"]
        if self.flux_flag:
            if not self.fluxes_been_calculated:
                self.calculate_fluxes()
            self.update_neighbouring_cells_cqs(dt_inv = dt_inv, cell_array = cell_array, \
                            map_interface_id_to_west_cell_idx = map_interface_id_to_west_cell_idx, \
                            map_interface_id_to_east_cell_idx = map_interface_id_to_east_cell_idx)
            
    def update_neighbouring_cells_cqs(self, dt_inv, cell_array, map_interface_id_to_west_cell_idx, \
                                                    map_interface_id_to_east_cell_idx):
        if self.on_west_boundary: # On west boundary, so only want to update cqs in cell to the east
            east_cell = cell_array[map_interface_id_to_east_cell_idx[self.interface_id]]
            east_cell.cqs["mass"] += dt_inv * self.geo["A"] * self.boundary_fluxes["mass"] / east_cell.geo["dV"]
            east_cell.cqs["xMom"] += dt_inv * self.geo["A"] * self.boundary_fluxes["xMom"] / east_cell.geo["dV"]
            east_cell.cqs["xMom"] += dt_inv * east_cell.geo["A_c"] * self.boundary_fluxes["p"] / east_cell.geo["dV"]
            east_cell.cqs["energy"] += dt_inv * self.geo["A"] * self.boundary_fluxes["energy"] / east_cell.geo["dV"]
            if self.multi_species:
                for ind in range(len(self.boundary_fluxes["massf"])):
                    east_cell.cqs["spcs_mass"][ind] += dt_inv * self.geo["A"] * self.boundary_fluxes["massf"][ind] / east_cell.geo["dV"]
            
        else: # On east boundary, so only want to update cqs in cell to the west
            west_cell = cell_array[map_interface_id_to_west_cell_idx[self.interface_id]]
            west_cell.cqs["mass"] -= dt_inv * self.geo["A"] * self.boundary_fluxes["mass"] / west_cell.geo["dV"]
            west_cell.cqs["xMom"] -= dt_inv * self.geo["A"] * self.boundary_fluxes["xMom"] / west_cell.geo["dV"]
            west_cell.cqs["xMom"] -= dt_inv * west_cell.geo["A_c"] * self.boundary_fluxes["p"] / west_cell.geo["dV"]
            west_cell.cqs["energy"] -= dt_inv * self.geo["A"] * self.boundary_fluxes["energy"] / west_cell.geo["dV"]
            if self.multi_species:
                for ind in range(len(self.boundary_fluxes["massf"])):
                    west_cell.cqs["spcs_mass"][ind] -= dt_inv * self.geo["A"] * self.boundary_fluxes["massf"][ind] / west_cell.geo["dV"]
            
    def calculate_fluxes(self):
        rho = self.flow_state.fluid_state.rho
        vel_x = self.flow_state.vel_x
        u = self.flow_state.fluid_state.u
        p = self.flow_state.fluid_state.p
        a = self.flow_state.fluid_state.a
        h_t = u + p / rho + 0.5 * vel_x ** 2.0

        self.boundary_fluxes["mass"] = vel_x * rho
        self.boundary_fluxes["xMom"] = vel_x * vel_x * rho
        self.boundary_fluxes["energy"] = vel_x * rho * h_t
        self.boundary_fluxes["p"] = p
        self.boundary_fluxes["vel_x"] = vel_x
        self.boundary_fluxes["Ma"] = vel_x / a

        if self.multi_species:
            massf = self.flow_state.fluid_state.massf
            self.boundary_fluxes["massf"] = (vel_x * np.array(massf) * rho).tolist()

        self.fluxes_been_calculated = True
        