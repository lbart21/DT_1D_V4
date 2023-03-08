"""
Function:
Author: Luke Bartholomew
Edits:
"""
import os

def write_to_data_file(cell_array, time, labels, flow_property_variables, sim_number) -> None:
    cwd = os.getcwd()
    for comp_label in labels:
        file_name = "Sim" + str(sim_number) + "DataAt" + str(format(time, ".9f")) +"ForComponent" + comp_label + ".txt"
        file = open(cwd + "/data/" + file_name, "w")
        file.write("Sim: " + str(sim_number) + "\n")
        file.write("Label: " + comp_label + "\n")
        file.write("Time: " + str(time) + "\n")
        n_cells = len(cell_array)
        first_cell_in_component = True
        for cell_idx in range(n_cells):
            if cell_array[cell_idx].label == comp_label:
                if cell_array[cell_idx].phase == "Single":
                    cell_flow_data = {}
                    if "vel_x" in flow_property_variables:
                        cell_flow_data["vel_x"] = cell_array[cell_idx].flow_state.vel_x
                    if "Ma" in flow_property_variables:
                        cell_flow_data["Ma"] = cell_array[cell_idx].flow_state.vel_x / cell_array[cell_idx].flow_state.fluid_state.a
                    if "p" in flow_property_variables:
                        cell_flow_data["p"] = cell_array[cell_idx].flow_state.fluid_state.p
                    if "T" in flow_property_variables:
                        cell_flow_data["T"] = cell_array[cell_idx].flow_state.fluid_state.T
                    if "rho" in flow_property_variables:
                        cell_flow_data["rho"] = cell_array[cell_idx].flow_state.fluid_state.rho
                    if "gamma" in flow_property_variables:
                        cell_flow_data["gamma"] = cell_array[cell_idx].flow_state.fluid_state.gamma
                    if "a" in flow_property_variables:
                        cell_flow_data["a"] = cell_array[cell_idx].flow_state.fluid_state.a
                    if "u" in flow_property_variables:
                        cell_flow_data["u"] = cell_array[cell_idx].flow_state.fluid_state.u
                    if "h" in flow_property_variables:
                        cell_flow_data["h"] = cell_array[cell_idx].flow_state.fluid_state.enthalpy
                    if "p_t" in flow_property_variables:
                        p = cell_array[cell_idx].flow_state.fluid_state.p
                        gamma = cell_array[cell_idx].flow_state.fluid_state.gamma
                        Ma = cell_array[cell_idx].flow_state.vel_x / cell_array[cell_idx].flow_state.fluid_state.a
                        cell_flow_data["p_t"] = p * (1.0 + 0.5 * (gamma - 1.0) * Ma ** 2.0) ** (gamma / (gamma - 1.0))
                    if "T_t" in flow_property_variables:
                        T = cell_array[cell_idx].flow_state.fluid_state.T
                        gamma = cell_array[cell_idx].flow_state.fluid_state.gamma
                        Ma = cell_array[cell_idx].flow_state.vel_x / cell_array[cell_idx].flow_state.fluid_state.a
                        cell_flow_data["T_t"] = T * (1.0 + 0.5 * (gamma - 1.0) * Ma ** 2.0)
                    if "massf" in flow_property_variables:
                        species_names = cell_array[cell_idx].flow_state.fluid_state.gmodel.species_names
                        massf_names = ["massf_" + name for name in species_names]
                        for ind, name in enumerate(massf_names):
                            cell_flow_data[name] = cell_array[cell_idx].flow_state.fluid_state.massf[ind]

                    #cellFlowData = {var : getattr(cellArray[cell_idx].fs.fluid_state, var) for var in flow_property_variables if var != "vel_x"}
                    #if "vel_x" in flow_property_variables:
                        #cellFlowData["vel_x"] = getattr(cellArray[cell_idx].fs, "vel_x")
                    joint_data_for_cell = {**cell_flow_data, **cell_array[cell_idx].geo}
                    
                else: #Two phase cell, add "_g" to all keys in gas flowState and "_l" to keys in liquid flowState
                    print("Two-phase data extraction not implemented yet.")
                    
                variable_names = list(joint_data_for_cell.keys())
                if first_cell_in_component is True:
                    # Add variable names to file
                    file.write("Variables: " + str(len(variable_names)) + "\n")
                    file.write(" ".join(variable_names) + "\n")
                    first_cell_in_component = False
                for name in variable_names:
                    file.write(str(format(joint_data_for_cell[name], ".9f")) + " ")
                file.write("\n")
        file.close()