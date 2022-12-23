import os
def WriteInterfaceDataToFile(interface, time, flow_property_variables):
    cwd = os.getcwd()
    interface_ID = interface.interface_ID
    filename = "DataAtInterfaceID" + str(interface_ID) + ".txt"
    if not os.path.exists(cwd + "/data/" + filename):
        with open(cwd + "/data/" + filename, "w") as file:
            interfaceData = {}
            if "A" in flow_property_variables:
                interfaceData["A"] = interface.GEO["A"]
            if "mass_flux" in flow_property_variables:
                interfaceData["mass_flux"] = interface.boundaryFluxes["mass"]
            if "energy_flux" in flow_property_variables:
                interfaceData["energy_flux"] = interface.boundaryFluxes["energy"]
            if "xMom_flux" in flow_property_variables:
                interfaceData["xMom_flux"] = interface.boundaryFluxes["xMom"]
            if "p" in flow_property_variables:
                interfaceData["p"] = interface.boundaryFluxes["p"]
            if "Ma" in flow_property_variables:
                interfaceData["Ma"] = interface.boundaryFluxes["Ma"]
            if "vel_x" in flow_property_variables:
                interfaceData["vel_x"] = interface.boundaryFluxes["vel_x"]
            if "gamma" in flow_property_variables:
                interfaceData["gamma"] = interface.LftState.fluid_state.gamma

            variableNames = list(interfaceData.keys()) #Does not include time
            file.write("ID: " + str(interface_ID) + "\n")
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
                interfaceData["mass_flux"] = interface.boundaryFluxes["mass"]
            if "energy_flux" in flow_property_variables:
                interfaceData["energy_flux"] = interface.boundaryFluxes["energy"]
            if "xMom_flux" in flow_property_variables:
                interfaceData["xMom_flux"] = interface.boundaryFluxes["xMom"]
            if "p" in flow_property_variables:
                interfaceData["p"] = interface.boundaryFluxes["p"]
            if "Ma" in flow_property_variables:
                interfaceData["Ma"] = interface.boundaryFluxes["Ma"]
            if "vel_x" in flow_property_variables:
                interfaceData["vel_x"] = interface.boundaryFluxes["vel_x"]
            if "gamma" in flow_property_variables:
                interfaceData["gamma"] = interface.LftState.fluid_state.gamma

            variableNames = list(interfaceData.keys()) #Does not include time

            file.write(str(format(time, ".9f")) + " ")
            for name in variableNames:
                file.write(str(format(interfaceData[name], ".9f")) + " ")
            file.write("\n")
    file.close()

    
