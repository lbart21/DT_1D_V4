import os

def WriteToDataFile(cellArray, time, labels, flow_property_variables) -> None:
    cwd = os.getcwd()
    for compLabel in labels:
        filename = "dataAt" + str(format(time, ".9f")) +"forComponent" + compLabel + ".txt"
        file = open(cwd + "/data/" + filename, "w")
        file.write("Label: " + compLabel + "\n")
        file.write("Time: " + str(time) + "\n")
        nCells = len(cellArray)
        firstCellInComponent = True
        for cell_idx in range(nCells):
            if cellArray[cell_idx].label == compLabel:
                if cellArray[cell_idx].phase == "Single":
                    cellFlowData = {}
                    if "vel_x" in flow_property_variables:
                        cellFlowData["vel_x"] = cellArray[cell_idx].fs.vel_x
                    if "Ma" in flow_property_variables:
                        cellFlowData["Ma"] = cellArray[cell_idx].fs.vel_x / cellArray[cell_idx].fs.fluid_state.a
                    if "p" in flow_property_variables:
                        cellFlowData["p"] = cellArray[cell_idx].fs.fluid_state.p
                    if "T" in flow_property_variables:
                        cellFlowData["T"] = cellArray[cell_idx].fs.fluid_state.T
                    if "rho" in flow_property_variables:
                        cellFlowData["rho"] = cellArray[cell_idx].fs.fluid_state.rho
                    if "gamma" in flow_property_variables:
                        cellFlowData["gamma"] = cellArray[cell_idx].fs.fluid_state.gamma
                    if "a" in flow_property_variables:
                        cellFlowData["a"] = cellArray[cell_idx].fs.fluid_state.a
                    if "u" in flow_property_variables:
                        cellFlowData["u"] = cellArray[cell_idx].fs.fluid_state.u
                    if "h" in flow_property_variables:
                        cellFlowData["h"] = cellArray[cell_idx].fs.fluid_state.enthalpy
                    if "p_t" in flow_property_variables:
                        p = cellArray[cell_idx].fs.fluid_state.p
                        gamma = cellArray[cell_idx].fs.fluid_state.gamma
                        Ma = cellArray[cell_idx].fs.vel_x / cellArray[cell_idx].fs.fluid_state.a
                        cellFlowData["p_t"] = p * (1.0 + 0.5 * (gamma - 1.0) * Ma ** 2.0) ** (gamma / (gamma - 1.0))
                    if "T_t" in flow_property_variables:
                        T = cellArray[cell_idx].fs.fluid_state.T
                        gamma = cellArray[cell_idx].fs.fluid_state.gamma
                        Ma = cellArray[cell_idx].fs.vel_x / cellArray[cell_idx].fs.fluid_state.a
                        cellFlowData["T_t"] = T * (1.0 + 0.5 * (gamma - 1.0) * Ma ** 2.0)


                    #cellFlowData = {var : getattr(cellArray[cell_idx].fs.fluid_state, var) for var in flow_property_variables if var != "vel_x"}
                    #if "vel_x" in flow_property_variables:
                        #cellFlowData["vel_x"] = getattr(cellArray[cell_idx].fs, "vel_x")
                    jointDataForCell = {**cellFlowData, **cellArray[cell_idx].GEO}
                    pass
                else: #Two phase cell, add "_g" to all keys in gas flowState and "_l" to keys in liquid flowState
                    for key in cellArray[cell_idx].gasFlowState.fs.keys():
                        cellArray[cell_idx].gasFlowState.fs[key + "_g"] = cellArray[cell_idx].gasFlowState.fs.pop(key)
                    for key in cellArray[cell_idx].liquidFlowState.fs.keys():
                        cellArray[cell_idx].liquidFlowState.fs[key + "_l"] = cellArray[cell_idx].liquidFlowState.fs.pop(key)
                    jointDataForCell = {**cellArray[cell_idx].gasFlowState.fs, **cellArray[cell_idx].liquidFlowState.fs, **cellArray[cell_idx].GEO}
                    
                variableNames = list(jointDataForCell.keys())
                if firstCellInComponent == True:
                    # Add variable names to file
                    file.write("Variables: " + str(len(variableNames)) + "\n")
                    file.write(" ".join(variableNames) + "\n")
                    firstCellInComponent = False
                for name in variableNames:
                    file.write(str(format(jointDataForCell[name], ".9f")) + " ")
                file.write("\n")
        file.close()