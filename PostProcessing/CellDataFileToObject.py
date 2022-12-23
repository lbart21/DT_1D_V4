import pandas
import os
class FormCellDataFromFile():
    def __init__(self, dataFileName) -> None:
        cwd = os.getcwd()
        fileLength = sum(1 for line in open(cwd + "/data/" + dataFileName))
        file = open(cwd + "/data/" + dataFileName)
        variableFlag = -1
        dataStart = False
        for rowInd, row in enumerate(file):
            splitRow = row.split(": ")
            if splitRow[0] == "Variables":
                variableFlag = rowInd + 1
            if splitRow[0] == "ID":
                self.cell_ID = int(splitRow[1])
            if rowInd == variableFlag:
                variableNames = [name for name in row[:-1].split(" ")]
                dataStartInd = rowInd + 1
                nTimes = fileLength - dataStartInd
                self.cellData = pandas.DataFrame(columns = variableNames, index = range(nTimes))
                dataStart = True
                continue
            if dataStart:
                for ind, name in enumerate(variableNames):
                    self.cellData[name][rowInd - dataStartInd] = float(row.split(' ')[ind])
        file.close()
        self.cellData.sort_values(by = ["time"])