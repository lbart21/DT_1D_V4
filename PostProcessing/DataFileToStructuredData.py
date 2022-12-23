import pandas
import os
class GenerateDataObject():
    def __init__(self, dataFileName) -> None:
        cwd = os.getcwd()
        fileLength = sum(1 for line in open(cwd + "/data/" + dataFileName))
        file = open(cwd + "/data/" + dataFileName)
        variableFlag = -1
        dataStart = False
        for rowInd, row in enumerate(file):
            splitRow = row.split(": ")
            if splitRow[0] == "Time":
                self.tFinal = float(splitRow[1])
            elif splitRow[0] == "Label":
                self.label = splitRow[1]
            elif splitRow[0] == "Variables":
                variableFlag = rowInd + 1
            if rowInd == variableFlag:
                variableNames = [name for name in row[:-1].split(" ")]
                dataStartInd = rowInd + 1
                self.nCells = fileLength - dataStartInd
                self.componentData = pandas.DataFrame(columns = variableNames, index = range(self.nCells))
                dataStart = True
                continue
            if dataStart:
                for ind, name in enumerate(variableNames):
                    self.componentData[name][rowInd - dataStartInd] = float(row.split(' ')[ind])
        file.close()
        self.componentData.sort_values(by = ["pos_x"])
        
        

