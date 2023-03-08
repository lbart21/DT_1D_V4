"""
Function:
Author: Luke Bartholomew
Edits:
    Jan 17 2023: Reworked to generate cell_to_interface mapping as list of lists.
"""
class SingleInlet1DMeshObject():
    """
    Basically the same as the meshObject class, but an assumed cell/interface
    order allows an easier definition of blocks.
    reversed decides if cells are expected to be generated in a normal fashion
    of west to east or in reverse from east to west. This is only really used in
    boundary cell generation where west boundary ghost cell layers are generated
    in reverse order.
    """
    def __init__(self, n_cells, reverse_direction_for_ghost_cells = False) -> None:
        self.cell_array = [None] * n_cells
        self.interface_array = [None] * (n_cells + 1)
        if reverse_direction_for_ghost_cells:
            self.map_cell_id_to_west_interface_idx = [[i + 1] for i in range(n_cells)]
            self.map_cell_id_to_east_interface_idx = [[i] for i in range(n_cells)]
            self.map_interface_id_to_west_cell_idx = list(range(n_cells)) + [None] 
            self.map_interface_id_to_east_cell_idx = [None] + list(range(n_cells))
            
        else:
            self.map_cell_id_to_west_interface_idx = [[i] for i in range(n_cells)]
            self.map_cell_id_to_east_interface_idx = [[i + 1] for i in range(n_cells)]
            self.map_interface_id_to_west_cell_idx = [None] + list(range(n_cells))
            self.map_interface_id_to_east_cell_idx = list(range(n_cells)) + [None]

        self.component_labels = []
        self.boundary_conditions = []
        self.boundary_interface_ids = [0, n_cells]
        self.cell_idx_to_track = []
        self.interface_idx_to_track = []
        