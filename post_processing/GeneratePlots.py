"""
Function:
Author: Luke Bartholomew
Edits:
"""
from Algorithms.DesignToolAlgorithmV4_1D.PostProcessing.DataFileToStructuredData import GenerateDataObject
from Algorithms.DesignToolAlgorithmV4_1D.PostProcessing.CellDataFileToObject import FormCellDataFromFile
from Algorithms.DesignToolAlgorithmV4_1D.PostProcessing.InterfaceDataFileToObject import FormInterfaceDataFromFile
from Algorithms.DesignToolAlgorithmV4_1D.PostProcessing.ProcessEilmerData import ProcessEilmerData
from Algorithms.DesignToolAlgorithmV4_1D.PostProcessing.SIUnitsDictionary import SIUnits
from Algorithms.DesignToolAlgorithmV4_1D.PostProcessing.Symbols import symbols

import Algorithms.DesignToolAlgorithmV2_0D.PostProcessing.InterfaceDataFileToObject as ZeroDInterfaceData

import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import PillowWriter

import os


class GenerateSinglePlots():
    def __init__(self, dataFile, plotVars) -> None:
        self.dataObject = GenerateDataObject(dataFileName = dataFile)
        
        tFinal = self.dataObject.tFinal
 
        for var in plotVars:
            fig = plt.figure(figsize=(15, 5))
            formattedTitleTime = '{:.3f}'.format(tFinal / 1e-6)
            formattedFileNameTime = '{:.9f}'.format(tFinal)
            plt.title("Distribution of " + symbols[var] + " at t = " \
                                                    + formattedTitleTime + r'$\mu$' + "s")
            plt.ylabel(symbols[var] + " (" + SIUnits[var] +")", \
                        rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            plt.scatter(self.dataObject.componentData["pos_x"], self.dataObject.componentData[var], marker = '.')
            filename = var + " distribution at t = " + formattedFileNameTime + ".jpg"
            plt.grid()
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()
        pass

class GenerateWaterfallPlots():
    def __init__(self, dataFiles, plotVars) -> None:
        dataFromFiles = {}
        t_list = []

        for file in dataFiles:
            data = GenerateDataObject(dataFileName = file)
            dataFromFiles[str(data.tFinal)] = data.componentData
            t_list.append(data.tFinal)

        for var in plotVars:
            fig = plt.figure(figsize=(15, 5))
            plt.title("Distribution of " + symbols[var] + " at multiple time values")
            plt.ylabel(symbols[var] + " (" + SIUnits[var] +")", \
                        rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            for time in t_list:
                formattedTime = '{:.3f}'.format(time / 1e-6)
                plt.scatter(dataFromFiles[str(time)]["pos_x"], dataFromFiles[str(time)][var], \
                                label = "Distribution at t = " + formattedTime + r'$\mu$' + "s", \
                                marker = ".")
            plt.legend()
            plt.grid()
            filename = var + " distribution at multiple times.jpg"
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()

class generateSingleComponentAnimation():
    def __init__(self, dataFiles, slowDownFactor, plotVars) -> None:
        self.data = {}
        self.dataTimes = [0.0]
        for file in dataFiles:
            dataObject = GenerateDataObject(dataFileName = file)
            self.data[str(dataObject.tFinal)] = dataObject.componentData
            self.dataTimes.append(dataObject.tFinal)
        
        timeStepList = []
        for i in range(len(self.dataTimes)):
            timeStepList.append(self.dataTimes[i+1] - self.dataTimes[i])
        
        for var in plotVars:
            self.fig, self.ax = plt.subplots()

            anim = animation.FuncAnimation(self.fig, self.update, frames = self.dataTimes, \
                                            blit = True, repeat = False)
            writerVideo = animation.PillowWriter(fps = 30)
            fileName = "AnimationOf" + var + ".gif"
            currentDir = os.getcwd()
            anim.save(currentDir + "/plots/" + fileName, writer = writerVideo)

    def update(self, i, fargs):
        (var,) = fargs
        x = self.data[str(i)]["pos_x"]
        y = self.data[str(i)][var]
        self.scat = self.ax.scatter(x, y)
        return self.scat, 

        
class GenerateSinglePlotsFromEilmerData():
    def __init__(self, EilmerDataNames, plotVars) -> None:
        EilmerData = ProcessEilmerData(dataFiles = EilmerDataNames)
        t = EilmerData.tFinal
        formattedTitleTime = '{:.3f}'.format(t / 1e-6)
        formattedFileNameTime = '{:.9f}'.format(t)
        if "Ma" in plotVars:
            EilmerData.componentData["Ma"] = EilmerData.componentData["vel_x"] / EilmerData.componentData["a"]

        if "p_t" in plotVars:
            gamma = 1.4
            p = EilmerData.componentData["p"]
            Ma = EilmerData.componentData["Ma"]
            EilmerData.componentData["p_t"] = p * (1.0 + 0.5 * (gamma - 1.0) * Ma ** 2.0) ** (gamma / (gamma - 1.0))
            
        if "T_t" in plotVars:
            gamma = 1.4
            T = EilmerData.componentData["T"]
            Ma = EilmerData.componentData["Ma"]
            EilmerData.componentData["T_t"] = T * (1.0 + 0.5 * (gamma - 1.0) * Ma ** 2.0)
            
        for var in plotVars:
            fig = plt.figure(figsize=(15, 5))
            plt.title("Eilmer Simulation Distribution of " + symbols[var] + " at t = " \
                                                    + formattedTitleTime + r'$\mu$' + "s")
            plt.ylabel(symbols[var] + " (" + SIUnits[var] +")", \
                        rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            plt.scatter(EilmerData.componentData["pos_x"], EilmerData.componentData[var], marker = '.')
            filename = var + " distribution at t = " + formattedFileNameTime + "WithEilmerSimulation.jpg"
            plt.grid()
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()
        pass

class GenerateThrustPlot():
    def __init__(self, interfaceFileName) -> None:
        interfaceData = FormInterfaceDataFromFile(dataFileName = interfaceFileName).interfaceData
        m_dot = interfaceData["mass_flux"]
        p_exit = interfaceData["p"]

        vel_x_e = interfaceData["vel_x"]
        
        A_e = interfaceData["A"]

        interfaceData["Thrust"] = m_dot * vel_x_e * A_e + p_exit * A_e

        fig = plt.figure(figsize=(15, 5))
        plt.title("Transient Thrust Profile")
        plt.ylabel("Thrust (N)", rotation = "horizontal", ha = "right")
        plt.xlabel("Time (ms)")
        plt.scatter(interfaceData["time"] * 1e3, interfaceData["Thrust"], marker = '.')
        filename = "ThrustProfile.jpg"
        plt.grid()
        mng = plt.get_current_fig_manager()
        mng.full_screen_toggle()
        currentDir = os.getcwd()
        plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
        plt.close()
        pass

class GenerateTransientCellPropertyPlots():
    def __init__(self, cellFileName, plotVars) -> None:
        cellData = FormCellDataFromFile(dataFileName = cellFileName)

        for var in plotVars:
            fig = plt.figure(figsize=(15, 5))
            plt.title("Transient Development of " + symbols[var] + " at Cell " + str(cellData.cell_ID))
            plt.ylabel(symbols[var] + " (" + SIUnits[var] +")", \
                        rotation = "horizontal", ha = "right")
            plt.xlabel("Time (ms)")
            plt.scatter(cellData.cellData["time"] * 1e3, cellData.cellData[var], marker = '.')
            filename = "Transient Development of " + var + " at Cell " + str(cellData.cell_ID) + ".jpg"
            plt.grid()
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()
        pass

class GenerateTransientInterfacePropertyPlots():
    def __init__(self, interfaceFileName, plotVars) -> None:
        interfaceData = FormInterfaceDataFromFile(dataFileName = interfaceFileName)

        for var in plotVars:
            fig = plt.figure(figsize=(15, 5))
            if var == "mass_flux":
                interfaceData.interfaceData["mass_flux"] *= interfaceData.interfaceData["A"]
            if var == "energy_flux":
                interfaceData.interfaceData["energy_flux"] *= interfaceData.interfaceData["A"]
            plt.title("Transient Development of " + symbols[var]+ " at Interface " + str(interfaceData.interface_ID))
            plt.ylabel(symbols[var] + " (" + SIUnits[var] +")", \
                        rotation = "horizontal", ha = "right")
            plt.xlabel("Time (ms)")
            plt.scatter(interfaceData.interfaceData["time"] * 1e3, interfaceData.interfaceData[var], marker = '.')
            filename = "Transient Development of " + var + " at Interface " + str(interfaceData.interface_ID) + ".jpg"
            plt.grid()
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()
        pass

class Compare1DTo0DThurstProfiles():
    def __init__(self, OneDInterfaceFileName, OneDCellFileName, ZeroDInterfaceFileName) -> None:
        ZeroDinterfaceDataObject = ZeroDInterfaceData.FormInterfaceDataFromFile(dataFileName = ZeroDInterfaceFileName).interfaceData
        ### Zero dim calculations
        ZeroDm_dot = ZeroDinterfaceDataObject["mass_flux"]
        ZeroDp_exit = ZeroDinterfaceDataObject["p"]
        ZeroDvel_x_exit = ZeroDinterfaceDataObject["vel_x"]

        ZeroDA_e = ZeroDinterfaceDataObject["A"]

        ZeroDinterfaceDataObject["Thrust"] = ZeroDm_dot * ZeroDvel_x_exit * ZeroDA_e + ZeroDp_exit * ZeroDA_e

        ### One dim calculations
        OneDcellData = FormCellDataFromFile(dataFileName = OneDCellFileName).cellData
        OneDinterfaceData = FormInterfaceDataFromFile(dataFileName = OneDInterfaceFileName).interfaceData
        OneDm_dot = OneDinterfaceData["mass_flux"]
        OneDp_exit = OneDinterfaceData["p"]

        OneDvel_x_e = OneDcellData["vel_x"][:-1]
        
        OneDA_e = OneDinterfaceData["A"]

        OneDinterfaceData["Thrust"] = OneDm_dot * OneDvel_x_e * OneDA_e + OneDp_exit * OneDA_e

        fig = plt.figure(figsize=(15, 5))
        plt.title("Transient Thrust Profiles of 0- and 1-D simulations")
        plt.ylabel("Thrust (N)", rotation = "horizontal", ha = "right")
        plt.xlabel("Time (ms)")
        plt.scatter(OneDinterfaceData["time"] * 1e3, OneDinterfaceData["Thrust"], marker = '.', label = "1-D Thrust")
        plt.scatter(ZeroDinterfaceDataObject["time"] * 1e3, ZeroDinterfaceDataObject["Thrust"], marker = '.', label = "0-D Thrust")
        filename = "ZeroAndOneDThrustProfileComparison.jpg"
        plt.grid()
        plt.legend()
        mng = plt.get_current_fig_manager()
        mng.full_screen_toggle()
        currentDir = os.getcwd()
        plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
        plt.close()
        pass

        
