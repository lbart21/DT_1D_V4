import os
def WriteCellDataToFile(cell, time, flow_property_variables):
    cwd = os.getcwd()
    cell_ID = cell.cell_ID
    filename = "DataAtCellID" + str(cell_ID) + ".txt"
    if not os.path.exists(cwd + "/data/" + filename):
        with open(cwd + "/data/" + filename, "w") as file:
            cellFlowData = {}
            if "vel_x" in flow_property_variables:
                cellFlowData["vel_x"] = cell.fs.vel_x
            if "Ma" in flow_property_variables:
                cellFlowData["Ma"] = cell.fs.vel_x / cell.fs.fluid_state.a
            if "p" in flow_property_variables:
                cellFlowData["p"] = cell.fs.fluid_state.p
            if "T" in flow_property_variables:
                cellFlowData["T"] = cell.fs.fluid_state.T
            if "rho" in flow_property_variables:
                cellFlowData["rho"] = cell.fs.fluid_state.rho
            if "u" in flow_property_variables:
                cellFlowData["u"] = cell.fs.fluid_state.u
            if "a" in flow_property_variables:
                cellFlowData["a"] = cell.fs.fluid_state.a
            if "h" in flow_property_variables:
                cellFlowData["h"] = cell.fs.fluid_state.enthalpy
            if "A_c" in flow_property_variables:
                cellFlowData["A_c"] = cell.GEO["A_c"]
            if "dV" in flow_property_variables:
                cellFlowData["dV"] = cell.GEO["dV"]
            if "p_t" in flow_property_variables:
                gamma = cell.fs.fluid_state.gamma
                Ma = cell.fs.vel_x / cell.fs.fluid_state.a
                p = cell.fs.fluid_state.p
                cellFlowData["p_t"] = p * (1.0 + 0.5 * (gamma - 1.0) * Ma ** 2.0) ** (gamma / (gamma - 1.0))
            if "T_t" in flow_property_variables:
                gamma = cell.fs.fluid_state.gamma
                Ma = cell.fs.vel_x / cell.fs.fluid_state.a
                T = cell.fs.fluid_state.T
                cellFlowData["T_t"] = T * (1.0 + 0.5 * (gamma - 1.0) * Ma ** 2.0)

            #cellFlowData = {var : getattr(cell.fs.fluid_state, var) for var in flow_property_variables if var != "vel_x"}
            #if "vel_x" in flow_property_variables:
                #cellFlowData["vel_x"] = getattr(cell.fs, "vel_x")

            variableNames = list(cellFlowData.keys()) #Does not include time
            file.write("ID: " + str(cell_ID) + "\n")
            file.write("Variables: " + str(len(variableNames) + 1) + "\n")
            file.write("time " + " ".join(variableNames) + "\n")
            file.write(str(format(time, ".9f")) + " ")
            for name in variableNames:
                file.write(str(format(cellFlowData[name], ".9f")) + " ")
            file.write("\n")
    else:
        with open(cwd + "/data/" + filename, "a") as file:

            cellFlowData = {}
            if "vel_x" in flow_property_variables:
                cellFlowData["vel_x"] = cell.fs.vel_x
            if "Ma" in flow_property_variables:
                cellFlowData["Ma"] = cell.fs.vel_x / cell.fs.fluid_state.a
            if "p" in flow_property_variables:
                cellFlowData["p"] = cell.fs.fluid_state.p
            if "T" in flow_property_variables:
                cellFlowData["T"] = cell.fs.fluid_state.T
            if "rho" in flow_property_variables:
                cellFlowData["rho"] = cell.fs.fluid_state.rho
            if "u" in flow_property_variables:
                cellFlowData["u"] = cell.fs.fluid_state.u
            if "a" in flow_property_variables:
                cellFlowData["a"] = cell.fs.fluid_state.a
            if "h" in flow_property_variables:
                cellFlowData["h"] = cell.fs.fluid_state.enthalpy
            if "A_c" in flow_property_variables:
                cellFlowData["A_c"] = cell.GEO["A_c"]
            if "dV" in flow_property_variables:
                cellFlowData["dV"] = cell.GEO["dV"]
            if "p_t" in flow_property_variables:
                gamma = cell.fs.fluid_state.gamma
                Ma = cell.fs.vel_x / cell.fs.fluid_state.a
                p = cell.fs.fluid_state.p
                cellFlowData["p_t"] = p * (1.0 + 0.5 * (gamma - 1.0) * Ma ** 2.0) ** (gamma / (gamma - 1.0))
            if "T_t" in flow_property_variables:
                gamma = cell.fs.fluid_state.gamma
                Ma = cell.fs.vel_x / cell.fs.fluid_state.a
                T = cell.fs.fluid_state.T
                cellFlowData["T_t"] = T * (1.0 + 0.5 * (gamma - 1.0) * Ma ** 2.0)

            #cellFlowData = {var : getattr(cell.fs.fluid_state, var) for var in flow_property_variables if var != "vel_x"}
            #if "vel_x" in flow_property_variables:
                #cellFlowData["vel_x"] = getattr(cell.fs, "vel_x")

            variableNames = list(cellFlowData.keys()) #Does not include time
            file.write(str(format(time, ".9f")) + " ")
            for name in variableNames:
                file.write(str(format(cellFlowData[name], ".9f")) + " ")
            file.write("\n")
    file.close()