"""
Function:
Author: Luke Bartholomew
Edits:
"""
from Algorithms.DesignToolAlgorithmV4_1D.Fluxes.fluidFluxes \
                        import AUSMPlusupORIGINAL, AUSMPlusupPAPER
from Algorithms.DesignToolAlgorithmV4_1D.Reconstruction.Reconstruction import getReconstruction
from Algorithms.DesignToolAlgorithmV4_1D.Reconstruction.LocateNeighbouringCellIndices import find_idx_of_cell_recursively



class SinglePhaseMultiSpeciesNonReactiveInterface():
    def __init__(self, interface_id, nL, nR, reconstruction_scheme, limiter, reconstruction_properties, update_from, flux_scheme) -> None:
        self.interface_id = interface_id
        self.nL = nL
        self.nR = nR
        
        self.flux_scheme = flux_scheme
        self.flux_flag = True
        self.reconstruction_scheme = reconstruction_scheme
        self.limiter = limiter
        self.reconstruction_properties = reconstruction_properties
        self.update_from = update_from
        self.stencil_been_formed = False
        
        self.lft_state = None
        self.rght_state = None

        self.boundary_fluxes = {}
    
    def fill_geometry(self, geometry):
        self.GEO = geometry
    
    def complete_interface_methods(self, cell_array, interface_array, map_cell_id_to_west_interface_idx, \
                                        map_cell_id_to_east_interface_idx, map_interface_id_to_west_cell_idx, \
                                        map_interface_id_to_east_cell_idx, dt_inv):
        if self.flux_flag:
            if not self.stencil_been_formed:
                self.form_lists_of_stencil_indices(map_cell_id_to_west_interface_idx = map_cell_id_to_west_interface_idx, \
                                                map_cell_id_to_east_interface_idx = map_cell_id_to_east_interface_idx, \
                                                map_interface_id_to_west_cell_idx = map_interface_id_to_west_cell_idx, \
                                                map_interface_id_to_east_cell_idx = map_interface_id_to_east_cell_idx, \
                                                cell_array = cell_array, interface_array = interface_array)
                self.stencil_been_formed = True

            self.reconstruct_states(cell_array = cell_array)
            self.calculate_fluxes()
            self.update_neighbouring_cells_cqs(dt_inv = dt_inv, cell_array = cell_array)
    
    def reconstruct_states(self, cell_array):
        if self.flux_flag is True: #Do reconstruction
            for prop in self.reconstruction_properties:
                qL_stencil, dxL_stencil, \
                qR_stencil, dxR_stencil = self.get_stencils(property = prop, cell_array = cell_array)
                if prop == "massf":
                    qL_stencil = list(map(list, zip(*qL_stencil)))
                    qR_stencil = list(map(list, zip(*qR_stencil)))
                    n_species = len(qL_stencil)
                    lft_reconstructed_massf = [None] * n_species
                    rght_reconstructed_massf = [None] * n_species
                    
                    for species_idx in range(len(qL_stencil)): 
                        lft_prop, rght_prop = getReconstruction(reconstruction = self.reconstruction_scheme[0], limiter = self.limiter, \
                                                        qL_stencil = qL_stencil[species_idx], dxL_stencil = dxL_stencil, \
                                                        qR_stencil = qR_stencil[species_idx], dxR_stencil = dxR_stencil)
                        lft_reconstructed_massf[species_idx] = lft_prop
                        rght_reconstructed_massf[species_idx] = rght_prop

                    setattr(self.lft_state.fluid_state, prop, lft_reconstructed_massf)
                    setattr(self.rght_state.fluid_state, prop, rght_reconstructed_massf)

                elif prop == "vel_x":
                    lft_prop, rght_prop = getReconstruction(reconstruction = self.reconstruction_scheme[0], limiter = self.limiter, \
                                                        qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
                    self.lft_state.vel_x = lft_prop
                    self.rght_state.vel_x = rght_prop
                else:
                    lft_prop, rght_prop = getReconstruction(reconstruction = self.reconstruction_scheme[0], limiter = self.limiter, \
                                                        qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
                    setattr(self.lft_state.fluid_state, prop, lft_prop)
                    setattr(self.rght_state.fluid_state, prop, rght_prop)

            if self.update_from == "pT":
                self.lft_state.fluid_state.update_thermo_from_pT()
                self.rght_state.fluid_state.update_thermo_from_pT()

    def calculate_fluxes(self):
        if self.flux_flag is True:
            #print(self.LftStencilIdxs, self.RghtStencilIdxs)
            if self.flux_scheme == "AUSMPlusUPOriginal":
                self.boundary_fluxes = AUSMPlusupORIGINAL(  a_L_forMa = self.lft_state.fluid_state.a,        a_R_forMa = self.rght_state.fluid_state.a, \
                                                            p_L = self.lft_state.fluid_state.p,              rho_L = self.lft_state.fluid_state.rho, \
                                                            u_L = self.lft_state.fluid_state.u,              a_L = self.lft_state.fluid_state.a, \
                                                            vel_x_L = self.lft_state.vel_x,                  p_R = self.rght_state.fluid_state.p, \
                                                            rho_R = self.rght_state.fluid_state.rho,         u_R = self.rght_state.fluid_state.u, \
                                                            a_R = self.rght_state.fluid_state.a,             vel_x_R = self.rght_state.vel_x).fluxes

            elif self.flux_scheme == "AUSMPlusUPPaper":
                self.boundary_fluxes = AUSMPlusupPAPER(  a_L_forMa = self.lft_state.fluid_state.a,  a_L = self.lft_state.fluid_state.a, \
                                                        p_L = self.lft_state.fluid_state.p,        u_L = self.lft_state.fluid_state.u, \
                                                        rho_L = self.lft_state.fluid_state.rho,    vel_x_L = self.lft_state.vel_x, \
                                                        a_R_forMa = self.rght_state.fluid_state.a, a_R = self.rght_state.fluid_state.a, \
                                                        p_R = self.rght_state.fluid_state.p,       u_R = self.rght_state.fluid_state.u, \
                                                        rho_R = self.rght_state.fluid_state.rho,   vel_x_R = self.rght_state.vel_x).fluxes
        
        else:
            pass

    def update_neighbouring_cells_cqs(self, dt_inv, cell_array):
        if self.flux_flag is True:
            westCell = cell_array[self.lft_stencil_idxs[0]]
            if westCell.interior_cell_flag: # Check if cell is a ghost cell as we don't want to update the conserved properties of ghost cells
                #print(dt_inv, self.GEO["A"], self.boundaryFluxes["mass"], westCell.GEO["dV"])
                westCell.cqs["mass"] -= dt_inv * self.GEO["A"] * self.boundary_fluxes["mass"] / westCell.GEO["dV"]
                westCell.cqs["xMom"] -= dt_inv * self.GEO["A"] * self.boundary_fluxes["xMom"] / westCell.GEO["dV"]
                westCell.cqs["xMom"] -= dt_inv * westCell.GEO["A_c"] * self.boundary_fluxes["p"] / westCell.GEO["dV"]
                westCell.cqs["energy"] -= dt_inv * self.GEO["A"] * self.boundary_fluxes["energy"] / westCell.GEO["dV"]

            eastCell = cell_array[self.rght_stencil_idxs[0]]
            if eastCell.interior_cell_flag: # Check if cell is a ghost cell as we don't want to update the conserved properties of ghost cells
                eastCell.cqs["mass"] += dt_inv * self.GEO["A"] * self.boundary_fluxes["mass"] / eastCell.GEO["dV"]
                eastCell.cqs["xMom"] += dt_inv * self.GEO["A"] * self.boundary_fluxes["xMom"] / eastCell.GEO["dV"]
                eastCell.cqs["xMom"] += dt_inv * eastCell.GEO["A_c"] * self.boundary_fluxes["p"] / eastCell.GEO["dV"]
                eastCell.cqs["energy"] += dt_inv * self.GEO["A"] * self.boundary_fluxes["energy"] / eastCell.GEO["dV"]
                
        else:
            pass
    
    def get_stencils(self, property, cell_array):
        ### Fill left stencil
        qL_stencil = [None] * len(self.lft_stencil_idxs)
        dxL_stencil = [None] * len(self.lft_stencil_idxs)
        for ind, LftCellIdx in enumerate(self.lft_stencil_idxs):
            cell_L = cell_array[LftCellIdx]
            if property == "vel_x":
                qL_stencil[ind] = cell_L.flow_state.vel_x
            else:
                qL_stencil[ind] = getattr(cell_L.flow_state.fluid_state, property)
            dxL_stencil[ind] = cell_L.GEO["dx"]

        ### Fill right stencil
        qR_stencil = [None] * len(self.rght_stencil_idxs)
        dxR_stencil = [None] * len(self.rght_stencil_idxs)
        for ind, RghtCellIdx in enumerate(self.rght_stencil_idxs):
            cell_R = cell_array[RghtCellIdx]
            if property == "vel_x":
                qR_stencil[ind] = cell_R.flow_state.vel_x
            else:
                qR_stencil[ind] = getattr(cell_R.flow_state.fluid_state, property)
            dxR_stencil[ind] = cell_R.GEO["dx"]
            
        return qL_stencil, dxL_stencil, qR_stencil, dxR_stencil
    
    def form_lists_of_stencil_indices(self, map_cell_id_to_west_interface_idx, map_cell_id_to_east_interface_idx, \
                                        map_interface_id_to_west_cell_idx, map_interface_id_to_east_cell_idx, cell_array, interface_array):
        self.lft_stencil_idxs = [None] * self.reconstruction_scheme[1][0]
        self.rght_stencil_idxs = [None] * self.reconstruction_scheme[1][1]
        for lft_stencil_idx in range(self.reconstruction_scheme[1][0]):
            cell_idx = find_idx_of_cell_recursively(interface_id = self.interface_id, recursion_depth = lft_stencil_idx + 1, \
                                                direction = "West", map_interface_id_to_east_cell_idx = map_interface_id_to_east_cell_idx, \
                                                cell_array = cell_array, map_cell_id_to_east_interface_idx = map_cell_id_to_east_interface_idx, \
                                                map_interface_id_to_west_cell_idx = map_interface_id_to_west_cell_idx, \
                                                map_cell_id_to_west_interface_idx = map_cell_id_to_west_interface_idx, \
                                                interface_array = interface_array)
            self.lft_stencil_idxs[lft_stencil_idx] = cell_idx #In order of [L_Idx0, L_Idx1, ...]
        
        for rght_stencil_idx in range(self.reconstruction_scheme[1][1]):
            cell_idx = find_idx_of_cell_recursively(interface_id = self.interface_id, recursion_depth = rght_stencil_idx + 1, \
                                                direction = "East", map_interface_id_to_east_cell_idx = map_interface_id_to_east_cell_idx, \
                                                cell_array = cell_array, map_cell_id_to_east_interface_idx = map_cell_id_to_east_interface_idx, \
                                                map_interface_id_to_west_cell_idx = map_interface_id_to_west_cell_idx, \
                                                map_cell_id_to_west_interface_idx = map_cell_id_to_west_interface_idx, \
                                                interface_array = interface_array)
            self.rght_stencil_idxs[rght_stencil_idx] = cell_idx #In order of [R_Idx0, R_Idx1, ...]