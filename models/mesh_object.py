"""
Function:
Author: Luke Bartholomew
Edits:
"""
class MeshObject():
    def __init__(self, n_cells):
        self.cell_array = [None] * n_cells
        self.interface_array = [None] * (n_cells + 1)
        self.map_cell_id_to_west_interface_idx = [[None]] * n_cells
        self.map_cell_id_to_east_interface_idx = [[None]] * n_cells
        self.map_interface_id_to_west_cell_idx = [None] * (n_cells + 1)
        self.map_interface_id_to_east_cell_idx = [None] * (n_cells + 1)
        self.boundary_conditions = []
        self.boundary_interface_ids = []
        self.component_labels = []
        self.cell_idx_to_track = []
        self.interface_idx_to_track = []
        
