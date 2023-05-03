"""
Function:
Author: Luke Bartholomew
Edits:
    Jan 17 2023: Reworked from V3 model to work with list of lists for cell to interface mapping
             and optimised by reducing number of repeated loops over arrays
"""
# pylint: disable = too-many-locals
from Algorithms.DT_1D_V4.models.empty_mesh_object import EmptyMeshObject
class JointBlock(EmptyMeshObject):
    """
    Take 2 mesh objects and return a joint mesh.
    """
    def __init__(self, mesh_object_1, mesh_object_2, block_1_interface_id_being_replaced, \
                        block_2_interface_id_being_replaced, new_interface, \
                        adding_ghost_cells_bool = False, new_component_label = None):
        """
        Given two meshes, makes a joined mesh, including cell array, interface array
        and mapping arrays from cellIDs to interface indicies and vice versa.
        Assume block 1 is always the larger one, so we'll be adding 2 onto 1,
        rather than 1 onto 2 to make the resizing process faster.
        joiningBlocksBool = True if joining separate components,
                            in which case the boundary interface indices have to be updated.
                          = False if just adding boundary cells,
                            hence boundary interface indices don't have to be updated.
        """
        super().__init__()
        #--------------------Preliminary Calcs----------------------------------
        ### Find if joining to east or west edge of meshObject1
        joining_to_west_edge = False
        joining_to_east_edge = False
        if mesh_object_1.map_interface_id_to_west_cell_idx[block_1_interface_id_being_replaced] is None:
            joining_to_west_edge = True
            # TODO also check opposite direction of meshObject2, if not None -> throw error
        elif mesh_object_1.map_interface_id_to_east_cell_idx[block_1_interface_id_being_replaced] is None:
            joining_to_east_edge = True
            # TODO also check opposite direction of meshObject2, if not None -> throw error
        else:
            print("Index given for meshObject1 interface is not a mesh boundary interface.")

        ### Will need to know how many cells and interfaces are in each mesh
        ncells_mesh_1 = len(mesh_object_1.cell_array)
        ninterfaces_mesh_1 = len(mesh_object_1.interface_array)
        ncells_mesh_2 = len(mesh_object_2.cell_array)

        """
        Update cell locations
        Update cell and component labels
        Update cell IDs
        Join cell_array and component_labels

        Remove boundary interfaces in mesh_object_1 and mesh_object_2.boundary_interface_ids
        Update mesh_object_2.boundary_interface_ids
        Join boundary_interface_ids arrays or assign new boundary_interface_ids as mesh_object_1.boundary_interface_ids
            if adding ghost cells

        Update mesh_object_2 boundary_condition locations if not adding ghost cells then join boundary_condition arrays
        Assign new boundary_condition as mesh_object_1.boundary_condition if adding ghost cells

        Update cell and interface IDs to track if not adding ghost cells
        Assign new cell and interface IDs to track as mesh_object_1 values if adding ghost cells
        """
        if not adding_ghost_cells_bool:
            if joining_to_west_edge: #Need to find location of joining cell of mesh2. Need to look west of joining interface of mesh2
                mesh_2_boundary_cell_idx = mesh_object_2.map_interface_id_to_west_cell_idx[block_2_interface_id_being_replaced]
                mesh_2_boundary_cell = mesh_object_2.cell_array[mesh_2_boundary_cell_idx]
                mesh_2_boundary_cell_location = mesh_2_boundary_cell.geo["pos_x"]
                mesh_2_boundary_cell_width = mesh_2_boundary_cell.geo["dx"]
                for mesh_1_cell in mesh_object_1.cell_array:
                    if new_component_label is not None:
                        mesh_1_cell.label = new_component_label
                    mesh_1_cell.geo["pos_x"] += mesh_2_boundary_cell_location + 0.5 * mesh_2_boundary_cell_width
                if new_component_label is not None:
                    for mesh_2_cell in mesh_object_2.cell_array:
                        mesh_2_cell.label = new_component_label
                
            elif joining_to_east_edge: #Need to find location of joining block of block1. Need to look west of interface of mesh2.
                mesh_1_boundary_cell_idx = mesh_object_1.map_interface_id_to_west_cell_idx[block_1_interface_id_being_replaced]
                mesh_1_boundary_cell = mesh_object_1.cell_array[mesh_1_boundary_cell_idx]
                mesh_1_boundary_cell_location = mesh_1_boundary_cell.geo["pos_x"]
                mesh_1_boundary_cell_width = mesh_1_boundary_cell.geo["dx"]
                if new_component_label is not None:
                    for mesh_1_cell in mesh_object_1.cell_array:
                        mesh_1_cell.label = new_component_label
                for mesh_2_cell in mesh_object_2.cell_array:
                    mesh_2_cell.geo["pos_x"] += mesh_1_boundary_cell_location + 0.5 * mesh_1_boundary_cell_width
                    mesh_2_cell.label = new_component_label
                
            print("Individual meshs' boundary interface indices were: ", mesh_object_1.boundary_interface_ids, " and ", mesh_object_2.boundary_interface_ids)

            if new_component_label is not None:
                self.component_labels = [new_component_label]
            else:
                self.component_labels = mesh_object_1.component_labels + mesh_object_2.component_labels

            mesh_object_1.boundary_interface_ids.remove(block_1_interface_id_being_replaced)
            mesh_object_2.boundary_interface_ids.remove(block_2_interface_id_being_replaced)

            #print("Interfaces in mesh 1:", ninterfaces_mesh_1)

            for ind, val in enumerate(mesh_object_2.boundary_interface_ids):
                if val > block_2_interface_id_being_replaced:
                    mesh_object_2.boundary_interface_ids[ind] += ninterfaces_mesh_1 - 1
                else:
                    mesh_object_2.boundary_interface_ids[ind] += ninterfaces_mesh_1

            self.boundary_interface_ids = mesh_object_1.boundary_interface_ids + mesh_object_2.boundary_interface_ids

            print("New mesh boundary interface indices are: ", self.boundary_interface_ids)
            #print("Mesh 2 boundary condition indices:", [i[0] for i in mesh_object_2.boundary_conditions])

            for ind, bc in enumerate(mesh_object_2.boundary_conditions):
                if bc[0] > block_2_interface_id_being_replaced:
                    mesh_object_2.boundary_conditions[ind] = [bc[0] + ninterfaces_mesh_1 - 1] + bc[1::]
                else:
                    mesh_object_2.boundary_conditions[ind] = [bc[0] + ninterfaces_mesh_1] + bc[1::]
            

            self.boundary_conditions = mesh_object_1.boundary_conditions + mesh_object_2.boundary_conditions
            #print("New boundary condition indices:", [i[0] for i in self.boundary_conditions])

            self.cell_idx_to_track = mesh_object_1.cell_idx_to_track + [i + ncells_mesh_1 for i in mesh_object_2.cell_idx_to_track]
            self.interface_idx_to_track = mesh_object_1.interface_idx_to_track + [j + ninterfaces_mesh_1 - 1 for j in mesh_object_2.interface_idx_to_track]

        else:
            self.component_labels = mesh_object_1.component_labels

            self.boundary_interface_ids = mesh_object_1.boundary_interface_ids

            self.boundary_conditions = mesh_object_1.boundary_conditions

            self.cell_idx_to_track = mesh_object_1.cell_idx_to_track
            self.interface_idx_to_track = mesh_object_1.interface_idx_to_track
            for cell in mesh_object_2.cell_array:
                cell.interior_cell_flag = False

        self.cell_array = mesh_object_1.cell_array + mesh_object_2.cell_array

        for ind, cell in enumerate(self.cell_array):
            cell.cell_id = ind

        # Update mesh_object_1.interface_array with new_interface
        # Delete interface from mesh_object_2.interface_array
        # Update interface IDs
        # Join interface_arrays

        mesh_object_1.interface_array[block_1_interface_id_being_replaced] = new_interface #Put new interface object into mesh1.interfaceArray

        del mesh_object_2.interface_array[block_2_interface_id_being_replaced]

        if joining_to_west_edge:   
            for ind, interface in enumerate(mesh_object_1.interface_array):
                interface.interface_id = ind
            for ind, interface in enumerate(mesh_object_2.interface_array):
                interface.interface_id = ind + ninterfaces_mesh_1

        if joining_to_east_edge:   
            for ind, interface in enumerate(mesh_object_1.interface_array):
                interface.interface_id = ind
            for ind, interface in enumerate(mesh_object_2.interface_array):
                interface.interface_id = ind + ninterfaces_mesh_1

        self.interface_array = mesh_object_1.interface_array + mesh_object_2.interface_array

        for ind, interface in enumerate(self.interface_array):
            interface.interface_ID = ind

        # Update values in mesh_object_2.map_cell_id_to_west_interface_idx
        # Join map_cell_id_to_west_interface_idx arrays

        for cell_boundary in mesh_object_2.map_cell_id_to_west_interface_idx:
            for ind, interface_idx in enumerate(cell_boundary):
                if interface_idx > block_2_interface_id_being_replaced:
                    cell_boundary[ind] -= 1
                cell_boundary[ind] += ninterfaces_mesh_1

        if joining_to_east_edge:
            mesh_2_boundary_cell_idx = mesh_object_2.map_interface_id_to_east_cell_idx[block_2_interface_id_being_replaced]
            mesh_object_2.map_cell_id_to_west_interface_idx[mesh_2_boundary_cell_idx] = [block_1_interface_id_being_replaced]

        self.map_cell_id_to_west_interface_idx = mesh_object_1.map_cell_id_to_west_interface_idx \
                                                + mesh_object_2.map_cell_id_to_west_interface_idx

        """
        Update values in mesh_object_2.map_cell_id_to_east_interface_idx
        Join map_cell_id_to_east_interface_idx arrays
        """

        for cell_boundary in mesh_object_2.map_cell_id_to_east_interface_idx:
            for ind, interface_idx in enumerate(cell_boundary):
                if interface_idx > block_2_interface_id_being_replaced:
                    cell_boundary[ind] -= 1
                cell_boundary[ind] += ninterfaces_mesh_1

        if joining_to_west_edge:
            mesh_2_boundary_cell_idx = mesh_object_2.map_interface_id_to_west_cell_idx[block_2_interface_id_being_replaced]
            mesh_object_2.map_cell_id_to_east_interface_idx[mesh_2_boundary_cell_idx] = [block_1_interface_id_being_replaced]

        self.map_cell_id_to_east_interface_idx = mesh_object_1.map_cell_id_to_east_interface_idx + mesh_object_2.map_cell_id_to_east_interface_idx

        """
        Remap mesh_object_1 boundary interface's west cell to mesh_object_2's bounadry cell
            if joining to west edge of mesh_object_1.
        Delete entry in mesh_object_2.map_interface_id_to_west_cell_idx for interface that
            is getting removed.
        Update mesh_object_2.map_interface_id_to_west_cell_idx values
        Join map_interface_id_to_west_cell_idx arrays
        """

        if joining_to_west_edge:
            mesh_2_boundary_cell_idx = mesh_object_2.map_interface_id_to_west_cell_idx[block_2_interface_id_being_replaced]
            mesh_2_boundary_cell_idx += ncells_mesh_1
            mesh_object_1.map_interface_id_to_west_cell_idx[block_1_interface_id_being_replaced] = mesh_2_boundary_cell_idx

        del mesh_object_2.map_interface_id_to_west_cell_idx[block_2_interface_id_being_replaced]

        for ind, value in enumerate(mesh_object_2.map_interface_id_to_west_cell_idx):
            if value is not None:
                mesh_object_2.map_interface_id_to_west_cell_idx[ind] += ncells_mesh_1

        self.map_interface_id_to_west_cell_idx = mesh_object_1.map_interface_id_to_west_cell_idx \
                                            + mesh_object_2.map_interface_id_to_west_cell_idx

        """
        Remap mesh_object_1 boundary interface's east cell to mesh_object_2's bounadry cell
            if joining to east edge of mesh_object_1.
        Delete entry in mesh_object_2.map_interface_id_to_east_cell_idx for interface that
            is getting removed.
        """

        if joining_to_east_edge:
            mesh_2_boundary_cell_idx = mesh_object_2.map_interface_id_to_east_cell_idx[block_2_interface_id_being_replaced]
            mesh_2_boundary_cell_idx += ncells_mesh_1
            mesh_object_1.map_interface_id_to_east_cell_idx[block_1_interface_id_being_replaced] = mesh_2_boundary_cell_idx

        del mesh_object_2.map_interface_id_to_east_cell_idx[block_2_interface_id_being_replaced]

        for ind, value in enumerate(mesh_object_2.map_interface_id_to_east_cell_idx):
            if value is not None:
                mesh_object_2.map_interface_id_to_east_cell_idx[ind] += ncells_mesh_1

        self.map_interface_id_to_east_cell_idx = mesh_object_1.map_interface_id_to_east_cell_idx \
                                            + mesh_object_2.map_interface_id_to_east_cell_idx
