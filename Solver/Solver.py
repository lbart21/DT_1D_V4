from Algorithms.DesignToolAlgorithmV3_1D.Solver.WriteToDataFile import WriteToDataFile
from Algorithms.DesignToolAlgorithmV3_1D.Solver.WriteCellDataToFile import WriteCellDataToFile
from Algorithms.DesignToolAlgorithmV3_1D.Solver.WriteInterfaceDataToFile import WriteInterfaceDataToFile
from Algorithms.DesignToolAlgorithmV3_1D.Integrate.Integrate import Integrate
import objsize
from copy import deepcopy
class Solver():
    def __init__(self, meshObject, cfl_flag, tFinal, dataSaveDt, transient_cell_flow_property_variables_to_write, \
                        transient_interface_flow_property_variables_to_write, \
                        spatial_cell_flow_property_variables_to_write, rapidDataSavenSteps) -> None:
        """
        meshObject = object with attributes: cellArray, interfaceArray, mapCellIDToWestInterfaceIdx, 
        cfl_flag = [Bool, float]
        labels = 
        """
        tCurrent = 0.0
        self.totalDataDict = {} #We'll store all data for performance calculations but will only write to text files at certain times
        #self.addToData(time = tCurrent, data = meshObject.cellArray)
        WriteToDataFile(cellArray = meshObject.cellArray, time = tCurrent, labels = meshObject.componentLabels, \
                        flow_property_variables = spatial_cell_flow_property_variables_to_write)
        for cellID in meshObject.cellIdxToTrack:
            #print("Writing cell data")
            WriteCellDataToFile(cell = meshObject.cellArray[cellID], time = tCurrent, \
                                flow_property_variables = transient_cell_flow_property_variables_to_write)
        
        time_tol = 1e-9
        tWriteTol = 1e-9
        tWrite = dataSaveDt
        writtenData = False
        rapidDataWritten = False
        currentMeshObject = meshObject
        #currentMeshObject = deepcopy(meshObject)
        currentStep = 0
        while tCurrent < tFinal and abs(tCurrent - tFinal) > time_tol:
            #print(currentStep)
            newData = Integrate(mesh = currentMeshObject, cfl_flag = cfl_flag, tCurrent = tCurrent, currentStep = currentStep)
            tCurrent += newData.dtTotal
            #print("t: ", tCurrent)
            writtenData = False
            rapidDataWritten = False
            if tCurrent > tWrite - tWriteTol:
                print("Writing data, t = ", tCurrent)
                WriteToDataFile(cellArray = newData.mesh.cellArray, time = tCurrent, labels = newData.mesh.componentLabels, \
                                flow_property_variables = spatial_cell_flow_property_variables_to_write)
                writtenData = True
                tWrite += dataSaveDt
            #self.addToData(time = tCurrent, data = newData.mesh.cellArray)
            
            currentMeshObject = newData.mesh
            #currentMeshObject = deepcopy(newData.mesh)
            currentStep += 1
            if currentStep % rapidDataSavenSteps == 0:
                #print(newData.mesh.cellIdxToTrack)
                #print(newData.mesh.interfaceIdxToTrack)
                for cellID in newData.mesh.cellIdxToTrack:
                    #print(cellID)
                    #print("Writing cell data")
                    WriteCellDataToFile(cell = newData.mesh.cellArray[cellID], time = tCurrent, \
                                        flow_property_variables = transient_cell_flow_property_variables_to_write)
                for interfaceID in newData.mesh.interfaceIdxToTrack:
                    #print(interfaceID)
                    WriteInterfaceDataToFile(interface = newData.mesh.interfaceArray[interfaceID], time = tCurrent - newData.dtTotal, \
                                                flow_property_variables = transient_interface_flow_property_variables_to_write)
                rapidDataWritten = True
            
        if not writtenData:
            WriteToDataFile(cellArray = currentMeshObject.cellArray, time = tCurrent, \
                            labels = currentMeshObject.componentLabels, \
                            flow_property_variables = spatial_cell_flow_property_variables_to_write)
        
        if not rapidDataWritten:
            #print(newData.mesh.cellIdxToTrack)
            #print(newData.mesh.interfaceIdxToTrack)
            for cellID in newData.mesh.cellIdxToTrack:
                #print(cellID)
                #print("Writing cell data")
                WriteCellDataToFile(cell = newData.mesh.cellArray[cellID], time = tCurrent, \
                                    flow_property_variables = transient_cell_flow_property_variables_to_write)
            for interfaceID in newData.mesh.interfaceIdxToTrack:
                #print(interfaceID)
                WriteInterfaceDataToFile(interface = newData.mesh.interfaceArray[interfaceID], \
                                    time = tCurrent - newData.dtTotal, \
                                    flow_property_variables = transient_interface_flow_property_variables_to_write)
        
    def addToData(self, time, data):
        self.totalDataDict[str(time)] = data
        print("totalDataDict size:", objsize.get_deep_size(self.totalDataDict) / 1e6, "Mb")