"""
Function:
Author: Luke Bartholomew
Edits:
"""
def find_idx_of_interface_recursively(interface_id, recursion_depth, direction, \
                                map_interface_id_to_east_cell_idx, cell_array, \
                                map_cell_id_to_east_interface_idx, \
                                map_interface_id_to_west_cell_idx, \
                                map_cell_id_to_west_interface_idx):

    current_interface_id = interface_id
    current_recursion_depth = 0

    while current_recursion_depth < recursion_depth:
        if direction == "East":
            neighbouring_cell_idx = map_interface_id_to_east_cell_idx[current_interface_id]
            neighbouring_cell_id = cell_array[neighbouring_cell_idx].cell_id
            if len(map_cell_id_to_east_interface_idx[neighbouring_cell_id]) > 1:
                print("Ghost cell generation that enters a multi-interface cell boundary is not implemented yet.")
            else:
                current_interface_id = map_cell_id_to_east_interface_idx[neighbouring_cell_id][0]
                current_recursion_depth += 1

        elif direction == "West":
            neighbouring_cell_idx = map_interface_id_to_west_cell_idx[current_interface_id]
            neighbouring_cell_id = cell_array[neighbouring_cell_idx].cell_id
            if len(map_cell_id_to_west_interface_idx[neighbouring_cell_id]) > 1:
                print("Ghost cell generation that enters a multi-interface cell boundary is not implemented yet.")
            else:
                current_interface_id = map_cell_id_to_west_interface_idx[neighbouring_cell_id][0]
                current_recursion_depth += 1

    return current_interface_id
