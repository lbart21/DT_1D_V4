"""
Function:
Author: Luke Bartholomew
Edits:
    Jan 17 2023: Reworked to generate cell_to_interface mapping as list of lists.
"""
class MultiInlet1DMeshObject():
    def __init__(self, n_cells, n_inlets, reverse_direction_for_mirrored_flow = False) -> None:
        """
        reverse_direction = False -> multiple inlets are on west edge of block
                            True  -> multiple inlets are on east edge of block
        """
        n_interfaces = n_cells + n_inlets
        self.cell_array = [None] * n_cells
        self.interface_array = [None] * n_interfaces
        if reverse_direction_for_mirrored_flow:
            self.map_cell_id_to_west_interface_idx = [[i] for i in range(n_cells)] # Done
            self.map_cell_id_to_east_interface_idx = [[i + 1] for i in range(n_cells - 1)] + [list(range(n_cells, n_cells + n_inlets))] # Done
            self.map_interface_id_to_west_cell_idx = [None] + list(range(n_cells - 1)) + [n_cells - 1] * n_inlets # Done
            self.map_interface_id_to_east_cell_idx = list(range(n_cells)) + [None] * n_inlets # Done
            self.boundary_interface_ids = [0] + [n_cells + i for i in range(n_inlets)] # Done
        else:
            self.map_cell_id_to_west_interface_idx = [list(range(n_inlets))] + [[i + n_inlets] for i in range(n_cells - 1)] # Done
            self.map_cell_id_to_east_interface_idx = [[i + n_inlets] for i in range(n_cells)] # Done
            self.map_interface_id_to_west_cell_idx = [None] * n_inlets + list(range(n_cells)) # Done
            self.map_interface_id_to_east_cell_idx = [0] * n_inlets + list(range(1, n_cells)) + [None] # Done
            self.boundary_interface_ids = list(range(n_inlets)) + [n_interfaces - 1] # Done

        self.component_labels = []
        self.boundary_conditions = []
        self.cell_idx_to_track = []
        self.interface_idx_to_track = []
        