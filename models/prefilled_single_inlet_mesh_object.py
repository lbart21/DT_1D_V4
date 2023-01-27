"""
Function:
Author: Luke Bartholomew
Edits:
    Jan 17 2023: Reworked to generate cell_to_interface mapping as list of lists.
"""
class basic1DMeshObject():
    """
    Basically the same as the meshObject class, but an assumed cell/interface
    order allows an easier definition of blocks.
    reversed decides if cells are expected to be generated in a normal fashion
    of west to east or in reverse from east to west. This is only really used in
    boundary cell generation where west boundary ghost cell layers are generated
    in reverse order.
    """
    def __init__(self, nCells, reversed = False) -> None:
        self.cell_array = [None] * nCells
        self.interface_array = [None] * (nCells + 1)
        if reversed:
            self.map_cell_id_to_west_interface_idx = [[i + 1] for i in range(nCells)]
            self.map_cell_id_to_east_interface_idx = [[i] for i in range(nCells)]
            self.map_interface_id_to_west_cell_idx = [i for i in range(nCells)] + [None] 
            self.map_interface_id_to_east_cell_idx = [None] + [i for i in range(nCells)] 
            
        else:
            self.map_cell_id_to_west_interface_idx = [[i] for i in range(nCells)]
            self.map_cell_id_to_east_interface_idx = [[i + 1] for i in range(nCells)]
            self.map_interface_id_to_west_cell_idx = [None] + [i for i in range(nCells)]
            self.map_interface_id_to_east_cell_idx = [i for i in range(nCells)] + [None]
        self.component_labels = []
        self.boundary_conditions = []
        self.boundary_interface_ids = [0, nCells]
        self.cell_idx_to_track = []
        self.interface_idx_to_track = []
        