"""
Function:
Author: Luke Bartholomew
Edits:
"""
import os
def write_interface_data_to_file(interface, time, flow_property_variables):
    cwd = os.getcwd()
    interface_id = interface.interface_id
    filename = "DataAtInterfaceID" + str(interface_id) + ".txt"
    if not os.path.exists(cwd + "/data/" + filename):
        with open(cwd + "/data/" + filename, "w") as file:
            interfaceData = {}
            if "A" in flow_property_variables:
                interfaceData["A"] = interface.GEO["A"]
            if "mass_flux" in flow_property_variables:
                interfaceData["mass_flux"] = interface.boundary_fluxes["mass"]
            if "energy_flux" in flow_property_variables:
                interfaceData["energy_flux"] = interface.boundary_fluxes["energy"]
            if "xMom_flux" in flow_property_variables:
                interfaceData["xMom_flux"] = interface.boundary_fluxes["xMom"]
            if "p" in flow_property_variables:
                interfaceData["p"] = interface.boundary_fluxes["p"]
            if "Ma" in flow_property_variables:
                interfaceData["Ma"] = interface.boundary_fluxes["Ma"]
            if "vel_x" in flow_property_variables:
                interfaceData["vel_x"] = interface.boundary_fluxes["vel_x"]
            if "gamma" in flow_property_variables:
                interfaceData["gamma"] = interface.lft_state.fluid_state.gamma

            variableNames = list(interfaceData.keys()) #Does not include time
            file.write("ID: " + str(interface_id) + "\n")
            file.write("Variables: " + str(len(variableNames) + 1) + "\n")
            file.write("time " + " ".join(variableNames) + "\n")
            file.write(str(format(time, ".9f")) + " ")
            for name in variableNames:
                file.write(str(format(interfaceData[name], ".9f")) + " ")
            file.write("\n")
            
    else:
        with open(cwd + "/data/" + filename, "a") as file:
            interfaceData = {}
            if "A" in flow_property_variables:
                interfaceData["A"] = interface.GEO["A"]
            if "mass_flux" in flow_property_variables:
                interfaceData["mass_flux"] = interface.boundary_fluxes["mass"]
            if "energy_flux" in flow_property_variables:
                interfaceData["energy_flux"] = interface.boundary_fluxes["energy"]
            if "xMom_flux" in flow_property_variables:
                interfaceData["xMom_flux"] = interface.boundary_fluxes["xMom"]
            if "p" in flow_property_variables:
                interfaceData["p"] = interface.boundary_fluxes["p"]
            if "Ma" in flow_property_variables:
                interfaceData["Ma"] = interface.boundary_fluxes["Ma"]
            if "vel_x" in flow_property_variables:
                interfaceData["vel_x"] = interface.boundary_fluxes["vel_x"]
            if "gamma" in flow_property_variables:
                interfaceData["gamma"] = interface.lft_state.fluid_state.gamma

            variableNames = list(interfaceData.keys()) #Does not include time

            file.write(str(format(time, ".9f")) + " ")
            for name in variableNames:
                file.write(str(format(interfaceData[name], ".9f")) + " ")
            file.write("\n")
    file.close()

    
