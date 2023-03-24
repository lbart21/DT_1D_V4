"""
Function:
Author: Luke Bartholomew
Edits:
"""
import math as m

from copy import deepcopy

from Algorithms.DT_1D_V4.extras.join_blocks import JointBlock
from Algorithms.DT_1D_V4.boundary_conditions.form_boundary_condition_information import FormBoundaryConditionInformation
class Integrate():
    def __init__(self, mesh, cfl_flag, t_current, current_step) -> None:
        """
        cfl_flag = [Bool, float] -> Bool = True if want cfl criteria to be used, false if fixed dt
                                 -> float = cfl_max if Bool = True, dt_max if Bool = False
        """
        self.total_ghost_cells = 0
        self.mesh = mesh

        ### Want to make copies of mapping arrays because it's inefficient to remake them after removing the ghost cells at the end of the time step
        old_map_cell_id_to_west_interface_idx = deepcopy(mesh.map_cell_id_to_west_interface_idx)
        old_map_cell_id_to_east_interface_idx = deepcopy(mesh.map_cell_id_to_east_interface_idx)
        old_map_interface_id_to_west_cell_idx = deepcopy(mesh.map_interface_id_to_west_cell_idx)
        old_map_interface_id_to_east_cell_idx = deepcopy(mesh.map_interface_id_to_east_cell_idx)
        boundary_conditions = deepcopy(mesh.boundary_conditions)
        boundary_interface_ids = deepcopy(mesh.boundary_interface_ids)
        
        ### Add boundary conditions
        self.add_boundary_conditions(BCs = self.mesh.boundary_conditions)
        
        ### Find out inviscid flux dt from CFL
        if not cfl_flag[0]: #Using fixed dt, cfl_flag has form [False, float dt_fixed]
            self.dt_total = cfl_flag[1]
        else:
            if cfl_flag[1] == "no_ramping":
                # Using fixed cfl
                # cfl_flag has form [True, "no_ramping", float cfl_fixed]
                self.dt_total = 1e6 #Large number to be initialised, won't be used
                for cell in self.mesh.cell_array:
                    self.dt_total = min(self.dt_total, cell.max_allowable_dt(cfl = cfl_flag[2]))
                
            elif cfl_flag[1] == "cfl_ramp_from_tCurrent": 
                # Ramping cfl based on tCurrent
                # cfl_flag has form [True, "cfl_ramp_from_tCurrent", [[cfl_1, t_1], [cfl_2, t_2], ..., [cfl_n, t_n]]]
                # If we are outside of t_n, use cfl_n. Else, find first instance tCurrent < t_i, then use cfl_i
                if t_current > cfl_flag[2][-1][1]:
                    cfl_max = cfl_flag[2][-1][0]
                    #print("cfl value used: ", cfl_max)
                    self.dt_total = 1e6 #Large number to be initialised, won't be used
                    for cell in self.mesh.cell_array:
                        self.dt_total = min(self.dt_total, cell.max_allowable_dt(cfl = cfl_max))
                else:
                    for cfl_pair in range(len(cfl_flag[2])):
                        if t_current > cfl_flag[2][cfl_pair][1]:
                            pass
                        else:
                            cfl_current = cfl_flag[2][cfl_pair][0]
                            #print("cfl value used: ", cfl_current)
                            self.dt_total = 1e6 #Large number to be initialised, won't be used
                            for cell in self.mesh.cell_array:
                                self.dt_total = min(self.dt_total, cell.max_allowable_dt(cfl = cfl_current))
                            break
                
            elif cfl_flag[1] == "cfl_ramp_from_step":
                # Ramping cfl based on currentStep
                # cfl_flag has form [True, "cfl_ramp_from_step", [[cfl_1, step_1], [cfl_2, step_2], ..., [cfl_n, step_n]]]
                # If we are outside of step_n, use cfl_n. Else, find first instance currentStep < step_i, then use cfl_i
                if current_step > cfl_flag[2][-1][1]: 
                    cfl_max = cfl_flag[2][-1][0]
                    #print("step: ", currentStep, " cfl used: ", cfl_max)
                    self.dt_total = 1e6 #Large number to be initialised, won't be used
                    for cell in self.mesh.cell_array:
                        self.dt_total = min(self.dt_total, cell.max_allowable_dt(cfl = cfl_max))

                else:
                    for cfl_pair in range(len(cfl_flag[2])):
                        if current_step > cfl_flag[2][cfl_pair][1]:
                            pass
                        
                        else:
                            cfl_current = cfl_flag[2][cfl_pair][0]
                            #print("step: ", currentStep, " cfl used: ", cfl_current)
                            self.dt_total = 1e6 #Large number to be initialised, won't be used
                            for cell in self.mesh.cell_array:
                                self.dt_total = min(self.dt_total, cell.max_allowable_dt(cfl = cfl_current))
                            break
                
            elif cfl_flag[1] == "dt_ramp_then_cfl_ramp_from_step":
                # Ramping dt initially, then use cfl based on current step number
                # cfl_flag has the form [True, "dt_ramp_then_cfl_ramp_from_step", [dt_init, scale], [[cfl_1, step_1], [cfl_2, step_2], ..., [cfl_n, step_n]]]
                # Starting at dt = dt_init, ramp this every time step by scale until the time step is greater than what comes from
                # cfl_1. Then ramp cfl similarly to "cfl_ramp_from_step" method.
                [dt_init, scale] = cfl_flag[2]
                min_dt_from_cfl = 1e6
                for cell in self.mesh.cell_array:
                    min_dt_from_cfl = min(min_dt_from_cfl, cell.max_allowable_dt(cfl = cfl_flag[3][0][0]))
                start_step = m.ceil( m.log(min_dt_from_cfl / dt_init, scale) )
                if current_step < start_step:
                    dt_ramped = dt_init * scale ** current_step
                    self.dt_total = dt_ramped
                
                else:
                    if current_step > cfl_flag[3][-1][1]:
                        self.dt_total = 1e6
                        for cell in self.mesh.cell_array:
                            self.dt_total = min(self.dt_total, cell.max_allowable_dt(cfl = cfl_flag[3][-1][0]))
                        print("step = ", current_step, " final cfl criteria used, cfl = ", cfl_flag[3][-1][0])
                    else:
                        for cfl_pair in range(len(cfl_flag[3])):
                            if current_step > cfl_flag[3][cfl_pair][1]:
                                pass
                            else:
                                cfl_current = cfl_flag[3][cfl_pair][0]
                                print("step = ", current_step, " Intermediate cfl value used = ", cfl_current)
                                self.dt_total = 1e6 #Large number to be initialised, won't be used
                                for cell in self.mesh.cell_array:
                                    self.dt_total = min(self.dt_total, cell.max_allowable_dt(cfl = cfl_current))
                                break

            elif cfl_flag[1] == "dt_ramp_then_cfl_ramp_from_tCurrent":
                # Ramping dt initially, then use cfl based on current time.
                # cfl_flag has the form [True, "dt_ramp_then_cfl_ramp_from_tCurrent", [dt_init, scale], [[cfl_1, t_1], [cfl_2, t_2], ..., [cfl_n, t_n]]]
                # Starting at dt = dt_init, ramp this every time step by scale until the time step is greater than what comes from
                # cfl_1. Then ramp cfl similarly to "cfl_ramp_from_tCurrent".
                [dt_init, scale] = cfl_flag[2]
                min_dt_from_cfl = 1e6
                for cell in self.mesh.cell_array:
                    min_dt_from_cfl = min(min_dt_from_cfl, cell.max_allowable_dt(cfl = cfl_flag[3][0][0]))
                start_step = m.ceil( m.log(min_dt_from_cfl / dt_init, scale) )
                if current_step < start_step:
                    dt_ramped = dt_init * scale ** current_step
                    self.dt_total = dt_ramped
                
                else:
                    if t_current > cfl_flag[3][-1][1]:
                        self.dt_total = 1e6
                        for cell in self.mesh.cell_array:
                            self.dt_total = min(self.dt_total, cell.max_allowable_dt(cfl = cfl_flag[3][-1][0]))
                    else:
                        for cfl_pair in range(len(cfl_flag[3])):
                            if t_current > cfl_flag[3][cfl_pair][1]:
                                pass
                            else:
                                cfl_current = cfl_flag[3][cfl_pair][0]
                                self.dt_total = 1e6 #Large number to be initialised, won't be used
                                for cell in self.mesh.cell_array:
                                    self.dt_total = min(self.dt_total, cell.max_allowable_dt(cfl = cfl_current))
                                break
                            
            else:
                print("Invalid CFL criteria given.")
        
        #print("This is the pre-convection properties")
        #for ind, cell in enumerate(self.mesh.cellArray):
            #print("index: ", ind, "p", cell.fs.fluid_state.p, ", T", cell.fs.fluid_state.T, ", rho", cell.fs.fluid_state.rho, ", vel_x", cell.fs.vel_x, "a", cell.fs.fluid_state.a)
            #print("index: ", ind, "Mass:", cell.conservedProperties["mass"], "xMom:", cell.conservedProperties["xMom"], "Energy:", cell.conservedProperties["energy"])
        ### Perform iteration over interfaces and update conserved properties
        #print("Pre-interface method conserved quantities")
        #print(self.mesh.cell_array[0].cqs, self.mesh.cell_array[0].geo["pos_x"])

        for interface in self.mesh.interface_array:
            ### Do all the interface specific methods that are defined within each interface (reconstruction, fluxes, perform conserved property updates etc etc)
            interface.complete_interface_methods(cell_array = self.mesh.cell_array, \
                                                interface_array = self.mesh.interface_array, \
                                                map_cell_id_to_west_interface_idx = self.mesh.map_cell_id_to_west_interface_idx, \
                                                map_cell_id_to_east_interface_idx = self.mesh.map_cell_id_to_east_interface_idx, \
                                                map_interface_id_to_west_cell_idx = self.mesh.map_interface_id_to_west_cell_idx, \
                                                map_interface_id_to_east_cell_idx = self.mesh.map_interface_id_to_east_cell_idx, \
                                                dt_inv = self.dt_total)
        #print("Post-interface method conserved quantities")
        #print(self.mesh.cell_array[0].cqs, self.mesh.cell_array[0].geo["pos_x"])
        ### Update properties after doing interface methods. Don't have to do this for boundary cells but we will for now
        #print("This is the post-convection properties")
        
        for ind, cell in enumerate(self.mesh.cell_array):
            cell.decode_to_primative_properties_after_interface_methods()
        
            #print("index: ", ind, "p", cell.fs.fluid_state.p, ", T", cell.fs.fluid_state.T, ", rho", cell.fs.fluid_state.rho, ", vel_x", cell.fs.vel_x, "a", cell.fs.fluid_state.a)
            #print("index: ", ind, "Mass:", cell.conservedProperties["mass"], "xMom:", cell.conservedProperties["xMom"], "Energy:", cell.conservedProperties["energy"])
        

        ### Perform iteration over cells and update conserved properties
        for cell in self.mesh.cell_array:
            cell.complete_cell_methods(cell_array = self.mesh.cell_array, \
                                        interface_array = self.mesh.interface_array, \
                                        map_cell_id_to_west_interface_idx = self.mesh.map_cell_id_to_west_interface_idx, \
                                        map_cell_id_to_east_interface_idx = self.mesh.map_cell_id_to_east_interface_idx, \
                                        map_interface_id_to_west_cell_idx = self.mesh.map_interface_id_to_west_cell_idx, \
                                        map_interface_id_to_east_cell_idx = self.mesh.map_interface_id_to_east_cell_idx, \
                                        dt_inv = self.dt_total)
        
        for cell in self.mesh.cell_array:
            cell.decode_to_primative_properties_after_cell_methods()
        #print("Post-cell method conserved quantities")
        #print(self.mesh.cell_array[100].cqs, self.mesh.cell_array[100].geo["pos_x"])
        
        
        #print("Printing cell properties at end of time step")
        #for cell in self.mesh.cellArray:
            #print(cell.cell_ID, cell.fs.fluid_state.p, cell.fs.fluid_state.T, cell.fs.fluid_state.rho, cell.fs.vel_x)
        
        ### Remove boundary cells and reset mapping arrays to what they were before adding boundary conditions
        self.mesh.map_cell_id_to_west_interface_idx = old_map_cell_id_to_west_interface_idx
        self.mesh.map_cell_id_to_east_interface_idx = old_map_cell_id_to_east_interface_idx
        self.mesh.map_interface_id_to_west_cell_idx = old_map_interface_id_to_west_cell_idx
        self.mesh.map_interface_id_to_east_cell_idx = old_map_interface_id_to_east_cell_idx
        self.mesh.boundary_conditions = boundary_conditions
        self.mesh.boundary_interface_ids = boundary_interface_ids
        self.mesh.cell_array = self.mesh.cell_array[:-1 * self.total_ghost_cells]
        self.mesh.interface_array = self.mesh.interface_array[:-1 * self.total_ghost_cells]

        #print(self.mesh.interface_array[0].boundary_fluxes)

        
        #print("This is the end of time step conserved quantities")
        #for ind, cell in enumerate(self.mesh.cellArray):
            #print("index: ", ind, "p", cell.fs.fluid_state.p, ", T", cell.fs.fluid_state.T, ", rho", cell.fs.fluid_state.rho, ", vel_x", cell.fs.vel_x, "a", cell.fs.fluid_state.a)
            #print("index: ", ind, "Mass:", cell.conservedProperties["mass"], "xMom:", cell.conservedProperties["xMom"], "Energy:", cell.conservedProperties["energy"])

    def add_boundary_conditions(self, BCs):

        #print(self.mesh.map_interface_id_to_west_cell_idx, self.mesh.map_interface_id_to_east_cell_idx, "\n")
        
        for boundary_condition in BCs:
            self.total_ghost_cells += boundary_condition[1][1]
            if boundary_condition[1][0] == "ConstantFlux_BC":
                pass
            else:
                ghost_cell_layer = FormBoundaryConditionInformation(mesh = self.mesh, BC = boundary_condition).ghost_cell_layer
                self.mesh = JointBlock(mesh_object_1 = self.mesh, mesh_object_2 = ghost_cell_layer, \
                    block_1_interface_id_being_replaced = boundary_condition[0], block_2_interface_id_being_replaced = 0, \
                    new_interface = deepcopy(self.mesh.interface_array[boundary_condition[0]]), adding_ghost_cells_bool = True)
                #print(self.mesh.map_interface_id_to_west_cell_idx, self.mesh.map_interface_id_to_east_cell_idx, "\n")