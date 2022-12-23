import os
import pandas

class ProcessEilmerData():
    def __init__(self, dataFiles) -> None:
        """
        dataFiles = list of raw data files
        """
        self.Data = {
            "nCells"    : 0
        }
        
        self.tFinal = None
        self.nVariables = None
        self.variableNames = []
        self.ReadEilmerFiles(dataFiles = dataFiles)

    def ReadEilmerFiles(self, dataFiles):
        firstFile = True
        for file in dataFiles:
            cwd = os.getcwd()
            f = open(cwd + "/data/" + file)
            dataStart = False
            currentFileData = {}
            if firstFile:
                variableFlag = -1
                for rowInd, row in enumerate(f):
                    splitRow = row.split(": ")
                    if splitRow[0] == "sim_time":
                        self.tFinal = float(splitRow[1])
                    elif splitRow[0] == "nicell":
                        self.Data["nCells"] += int(splitRow[1])
                        nCellsInRow = int(splitRow[1])
                        for name in self.variableNames:
                            self.Data[name] = [] 
                            currentFileData[name] = [None] * nCellsInRow
                    elif splitRow[0] == "variables":
                        self.nVariables = int(splitRow[1])
                        variableFlag = rowInd + 1
                    elif splitRow[0] == "nkcell":
                        dataStart = True
                        dataStartInd = rowInd + 1 
                        continue
                    if rowInd == variableFlag:
                        self.variableNames = row[1:-1].replace('"', "").split(" ")
                        self.variableNames = ["pos_x" if name == "pos.x" else name for name in self.variableNames]
                        self.variableNames = ["vel_x" if name == "vel.x" else name for name in self.variableNames]
                    if dataStart:
                        if rowInd in range(dataStartInd, dataStartInd+nCellsInRow):
                            for ind, value in enumerate(row[1:-1].split(" ")):
                                currentFileData[self.variableNames[ind]][rowInd - dataStartInd] = float(value)
                        else:
                            dataStart = False
                firstFile = False
                for name in self.variableNames:
                    self.Data[name] += currentFileData[name]
            else:
                for rowInd, row in enumerate(f):
                    splitRow = row.split(": ")
                    if splitRow[0] == "nicell":
                        self.Data["nCells"] += int(splitRow[1])
                        nCellsInRow = int(splitRow[1])
                        for name in self.variableNames:
                            currentFileData[name] = [None] * nCellsInRow
                    elif splitRow[0] == "nkcell":
                        dataStart = True
                        dataStartInd = rowInd + 1
                        continue
                    if dataStart:
                        if rowInd in range(dataStartInd, dataStartInd+nCellsInRow):
                            for ind, value in enumerate(row[1:-1].split(" ")):
                                currentFileData[self.variableNames[ind]][rowInd - dataStartInd] = float(value)
                        else:
                            dataStart = False
                for name in self.variableNames:
                    self.Data[name] += currentFileData[name]
            f.close()
            self.componentData = pandas.DataFrame(columns = self.variableNames, index = range(self.Data["nCells"]))
            for cell in range(self.Data["nCells"]):
                for var in self.variableNames:
                    self.componentData[var][cell] = self.Data[var][cell]
            self.componentData.sort_values(by = ["pos_x"])