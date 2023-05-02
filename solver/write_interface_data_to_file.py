"""
Function:
Author: Luke Bartholomew
Edits:
"""
import os
def write_interface_data_to_file(interface, time, flow_property_variables, sim_number, simulation_description):
    cwd = os.getcwd()
    interface_id = interface.interface_id
    file_name = "Sim" + str(sim_number) + "DataAtInterfaceID" + str(interface_id) + ".txt"
    if not os.path.exists(cwd + "/data/" + file_name):
        with open(cwd + "/data/" + file_name, "w") as file:
            interface_data = {}
            if "A" in flow_property_variables:
                interface_data["A"] = interface.geo["A"]
            if "mass_flux" in flow_property_variables:
                interface_data["mass_flux"] = interface.boundary_fluxes["mass"]
            if "energy_flux" in flow_property_variables:
                interface_data["energy_flux"] = interface.boundary_fluxes["energy"]
            if "xMom_flux" in flow_property_variables:
                interface_data["xMom_flux"] = interface.boundary_fluxes["xMom"]
            if "p" in flow_property_variables:
                interface_data["p"] = interface.boundary_fluxes["p"]
            if "Ma" in flow_property_variables:
                interface_data["Ma"] = interface.boundary_fluxes["Ma"]
            if "vel_x" in flow_property_variables:
                interface_data["vel_x"] = interface.boundary_fluxes["vel_x"]
            if "gamma" in flow_property_variables:
                interface_data["gamma"] = interface.lft_state.fluid_state.gamma

            variable_names = list(interface_data.keys()) #Does not include time
            file.write("Sim: " + str(sim_number) + "\n")
            file.write("ID: " + str(interface_id) + "\n")
            file.write("Description: " + simulation_description + "\n")
            file.write("Variables: " + str(len(variable_names) + 1) + "\n")
            file.write("time " + " ".join(variable_names) + "\n")
            file.write(str(format(time, ".9f")) + " ")
            for name in variable_names:
                file.write(str(format(interface_data[name], ".9f")) + " ")
            file.write("\n")
            
    else:
        with open(cwd + "/data/" + file_name, "a") as file:
            interface_data = {}
            if "A" in flow_property_variables:
                interface_data["A"] = interface.geo["A"]
            if "mass_flux" in flow_property_variables:
                interface_data["mass_flux"] = interface.boundary_fluxes["mass"]
            if "energy_flux" in flow_property_variables:
                interface_data["energy_flux"] = interface.boundary_fluxes["energy"]
            if "xMom_flux" in flow_property_variables:
                interface_data["xMom_flux"] = interface.boundary_fluxes["xMom"]
            if "p" in flow_property_variables:
                interface_data["p"] = interface.boundary_fluxes["p"]
            if "Ma" in flow_property_variables:
                interface_data["Ma"] = interface.boundary_fluxes["Ma"]
            if "vel_x" in flow_property_variables:
                interface_data["vel_x"] = interface.boundary_fluxes["vel_x"]
            if "gamma" in flow_property_variables:
                interface_data["gamma"] = interface.lft_state.fluid_state.gamma

            variable_names = list(interface_data.keys()) #Does not include time

            file.write(str(format(time, ".9f")) + " ")
            for name in variable_names:
                file.write(str(format(interface_data[name], ".9f")) + " ")
            file.write("\n")
    file.close()

    
