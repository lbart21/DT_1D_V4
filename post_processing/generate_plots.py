"""
Function:
Author: Luke Bartholomew
Edits:
"""
import os
import re

from Algorithms.DT_1D_V4.post_processing.data_file_to_structured_data import GenerateDataObject
from Algorithms.DT_1D_V4.post_processing.cell_data_file_to_object import FormCellDataFromFile
from Algorithms.DT_1D_V4.post_processing.interface_data_file_to_object \
            import FormInterfaceDataFromFile
from Algorithms.DT_1D_V4.post_processing.process_eilmer_data import ProcessEilmerData
from Algorithms.DT_1D_V4.post_processing.SI_units_dictionary import SI_UNITS
from Algorithms.DT_1D_V4.post_processing.symbols import SYMBOLS

from Algorithms.DT_0D_V2.post_processing.interface_data_file_to_object \
            import FormInterfaceDataFromFile as interface_data_object_0d

import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import PillowWriter

class GenerateSinglePlots():
    def __init__(self, data_file, plot_vars) -> None:
        self.data_object = GenerateDataObject(data_file_name = data_file)
        
        t_final = self.data_object.t_final
        sim_number = self.data_object.sim_number
 
        for var in plot_vars:
            fig = plt.figure(figsize=(15, 5))
            formatted_title_time = '{:.3f}'.format(t_final / 1e-6)
            formatted_file_name_time = '{:.9f}'.format(t_final)
            
            if var == "massf": # Plot all mass fractions
                column_names = list(self.data_object.component_data.columns)
                massf_species_names = [species for species in column_names \
                                                if len(species) > 5 and species[:5] == "massf"]
                species_names = [species.split("_")[1] for species in massf_species_names] 
                for index, species in enumerate(species_names):
                    split_name = re.findall('(\d+|[A-Za-z]+)', species)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    plt.scatter(self.data_object.component_data["pos_x"], \
                                self.data_object.component_data[massf_species_names[index]], \
                                marker = '.', label = r"$" + formatted_species_name + "$")
                plt.legend()

            elif len(var) > 5 and var[:5] == "massf": # Plot specific mass fractions
                species_name = var[6:] # Trim off starting massf_
                split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
                for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                SI_UNITS[var] = r'$-$'
                SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mass$' + ' ' + r'$Fraction$'
                plt.scatter(self.data_object.component_data["pos_x"], \
                            self.data_object.component_data[var], marker = '.')
            else: # All other properties
                plt.scatter(self.data_object.component_data["pos_x"], self.data_object.component_data[var], marker = '.')
            plt.xlabel("Position (m)")
            plt.ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
            plt.title("Distribution of " + SYMBOLS[var] + " at t = " \
                                                    + formatted_title_time + r'$\mu$' + "s")
            
            filename = "Sim" + str(sim_number) + ' ' + var + " distribution at t = " + formatted_file_name_time + ".jpg"
            plt.grid()
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            current_dir = os.getcwd()
            plt.savefig(current_dir + "/plots/" + filename, bbox_inches="tight")
            plt.close()

class GenerateWaterfallPlots():
    def __init__(self, data_files, plot_vars) -> None:
        data_from_files = {}
        t_list = []

        for file in data_files:
            data = GenerateDataObject(data_file_name = file)
            self.sim_number = data.sim_number
            data_from_files[str(data.t_final)] = data.component_data
            t_list.append(data.t_final)

        for var in plot_vars:
            fig = plt.figure(figsize=(15, 5))
            if len(var) > 5 and var[:5] == "massf":
                species_name = var[6:] # Trim off starting massf_
                split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
                for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                SI_UNITS[var] = r'$-$'
                SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mass$' + ' ' + r'$Fraction$'
            for time in t_list:
                formatted_time = '{:.3f}'.format(time / 1e-6)
                if var == "massf":
                    column_names = list(data_from_files[str(time)].columns)
                    massf_species_names = [species for species in column_names \
                                                if len(species) > 5 and species[:5] == "massf"]
                    species_names = [species.split("_")[1] for species in massf_species_names] 
                    for index, species in enumerate(species_names):
                        split_name = re.findall('(\d+|[A-Za-z]+)', species)
                        for ind, val in enumerate(split_name):
                            if val.isnumeric():
                                split_name[ind] = "_{" + split_name[ind] + "}"
                        formatted_species_name = ''.join(split_name)
                        plt.scatter(data_from_files[str(time)]["pos_x"], \
                                    data_from_files[str(time)][massf_species_names[index]], \
                                    marker = '.', label = r"$" + formatted_species_name + "$" \
                                                    + " at t = " + formatted_time + r'$\mu$' + "s")
                    plt.legend()
                else:
                    plt.scatter(data_from_files[str(time)]["pos_x"], data_from_files[str(time)][var], \
                                label = "Distribution at t = " + formatted_time + r'$\mu$' + "s", \
                                marker = ".")
            plt.title("Distribution of " + SYMBOLS[var] + " at Various Times")
            plt.ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            plt.legend()
            plt.grid()
            filename = "Sim " + str(self.sim_number) + ' ' + var + " distribution at multiple times.jpg"
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            current_dir = os.getcwd()
            plt.savefig(current_dir + "/plots/" + filename, bbox_inches="tight")
            plt.close()

class GenerateSingleComponentAnimation():
    def __init__(self, data_files, slow_down_factor, plot_vars) -> None:
        self.data = {}
        self.data_times = [0.0]
        for file in data_files:
            data_object = GenerateDataObject(data_file_name = file)
            self.data[str(data_object.t_final)] = data_object.component_data
            self.data_times.append(data_object.t_final)
        
        time_step_list = []
        for i in range(len(self.data_times)):
            time_step_list.append(self.data_times[i+1] - self.data_times[i])
        
        for var in plot_vars:
            self.fig, self.ax = plt.subplots()

            anim = animation.FuncAnimation(self.fig, self.update, frames = self.data_times, \
                                            blit = True, repeat = False)
            writer_video = animation.PillowWriter(fps = 30)
            filename = "AnimationOf" + var + ".gif"
            current_dir = os.getcwd()
            anim.save(current_dir + "/plots/" + filename, writer = writer_video)

    def update(self, i, fargs):
        (var,) = fargs
        x = self.data[str(i)]["pos_x"]
        y = self.data[str(i)][var]
        self.scat = self.ax.scatter(x, y)
        return self.scat, 

        
class GenerateSinglePlotsFromEilmerData():
    def __init__(self, eilmer_data_names, plot_vars) -> None:
        eilmer_data = ProcessEilmerData(data_files = eilmer_data_names)
        t = eilmer_data.t_final
        formatted_title_time = '{:.3f}'.format(t / 1e-6)
        formatted_file_name_time = '{:.9f}'.format(t)
        if "Ma" in plot_vars:
            eilmer_data.component_data["Ma"] = eilmer_data.component_data["vel_x"] / eilmer_data.component_data["a"]

        if "p_t" in plot_vars:
            gamma = 1.4
            p = eilmer_data.component_data["p"]
            Ma = eilmer_data.component_data["Ma"]
            eilmer_data.component_data["p_t"] = p * (1.0 + 0.5 * (gamma - 1.0) * Ma ** 2.0) ** (gamma / (gamma - 1.0))
            
        if "T_t" in plot_vars:
            gamma = 1.4
            T = eilmer_data.component_data["T"]
            Ma = eilmer_data.component_data["Ma"]
            eilmer_data.component_data["T_t"] = T * (1.0 + 0.5 * (gamma - 1.0) * Ma ** 2.0)
            
        for var in plot_vars:
            fig = plt.figure(figsize=(15, 5))
            plt.title("Eilmer Simulation Distribution of " + SYMBOLS[var] + " at t = " \
                                                    + formatted_title_time + r'$\mu$' + "s")
            plt.ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            plt.scatter(eilmer_data.component_data["pos_x"], eilmer_data.component_data[var], marker = '.')
            filename = var + " distribution at t = " + formatted_file_name_time + "WithEilmerSimulation.jpg"
            plt.grid()
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()

class GenerateThrustPlot():
    def __init__(self, interface_file_name) -> None:
        interface_data = FormInterfaceDataFromFile(data_file_name = interface_file_name).interface_data
        m_dot = interface_data["mass_flux"]
        p_exit = interface_data["p"]

        vel_x_e = interface_data["vel_x"]
        
        A_e = interface_data["A"]

        interface_data["Thrust"] = m_dot * vel_x_e * A_e + p_exit * A_e

        fig = plt.figure(figsize=(15, 5))
        plt.title("Transient Thrust Profile")
        plt.ylabel("Thrust (N)", rotation = "horizontal", ha = "right")
        plt.xlabel("Time (ms)")
        plt.scatter(interface_data["time"] * 1e3, interface_data["Thrust"], marker = '.')
        filename = "ThrustProfile.jpg"
        plt.grid()
        mng = plt.get_current_fig_manager()
        mng.full_screen_toggle()
        currentDir = os.getcwd()
        plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
        plt.close()

class GenerateTransientCellPropertyPlots():
    def __init__(self, cell_file_name, plot_vars) -> None:
        cell_data = FormCellDataFromFile(data_file_name = cell_file_name)
        sim_number = cell_data.sim_number
        for var in plot_vars:
            if len(var) > 5 and var[:5] == "massf":
                species_name = var[6:] # Trim off starting massf_
                split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
                for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                SI_UNITS[var] = r'$-$'
                SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mass$' + ' ' + r'$Fraction$'
            if var == "massf":
                massf_names = [name for name in cell_data.cell_data.columns if "massf" in name]

                fig = plt.figure(figsize=(15, 5))
                plt.title("Transient Development of " + SYMBOLS[var] + " at Cell " + str(cell_data.cell_id))
                plt.ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                            rotation = "horizontal", ha = "right")
                plt.xlabel("Time (ms)")
                for species in massf_names:
                    name = species[6:]
                    split_name = re.findall('(\d+|[A-Za-z]+)', name)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_name = ''.join(split_name)
                    plt.scatter(cell_data.cell_data["time"] * 1e3, cell_data.cell_data[species], marker = '.', label = r"$" + formatted_name + "$")
                filename = "Sim " + str(sim_number) + " Transient Development of " + var + " at Cell " + str(cell_data.cell_id) + ".jpg"
                plt.grid()
                plt.legend()
                mng = plt.get_current_fig_manager()
                mng.full_screen_toggle()
                current_dir = os.getcwd()
                plt.savefig(current_dir + "/plots/" + filename, bbox_inches="tight")
                plt.close()
            else:
                fig = plt.figure(figsize=(15, 5))
                plt.title("Transient Development of " + SYMBOLS[var] + " at Cell " + str(cell_data.cell_id))
                plt.ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                            rotation = "horizontal", ha = "right")
                plt.xlabel("Time (ms)")
                plt.scatter(cell_data.cell_data["time"] * 1e3, cell_data.cell_data[var], marker = '.')
                filename = "Sim " + str(sim_number) + " Transient Development of " + var + " at Cell " + str(cell_data.cell_id) + ".jpg"
                plt.grid()
                mng = plt.get_current_fig_manager()
                mng.full_screen_toggle()
                current_dir = os.getcwd()
                plt.savefig(current_dir + "/plots/" + filename, bbox_inches="tight")
                plt.close()

class GenerateTransientInterfacePropertyPlots():
    def __init__(self, interface_file_name, plot_vars) -> None:
        interface_data = FormInterfaceDataFromFile(data_file_name = interface_file_name)

        for var in plot_vars:
            fig = plt.figure(figsize=(15, 5))
            if var == "mass_flux":
                interface_data.interface_data["mass_flux"] *= interface_data.interface_data["A"]
            if var == "energy_flux":
                interface_data.interface_data["energy_flux"] *= interface_data.interface_data["A"]
            plt.title("Transient Development of " + SYMBOLS[var]+ " at Interface " + str(interface_data.interface_id))
            plt.ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
            plt.xlabel("Time (ms)")
            plt.scatter(interface_data.interface_data["time"] * 1e3, interface_data.interface_data[var], marker = '.')
            filename = "Transient Development of " + var + " at Interface " + str(interface_data.interface_id) + ".jpg"
            plt.grid()
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()

class Compare1DTo0DThurstProfiles():
    def __init__(self, interface_file_name_1d, cell_file_name_1d, interface_file_name_0d) -> None:
        interface_data_0d = interface_data_object_0d(data_file_name = interface_file_name_0d).interface_data
        ### Zero dim calculations
        m_dot_0d = interface_data_0d["mass_flux"]
        p_exit_0d = interface_data_0d["p"]
        vel_x_exit_0d = interface_data_0d["vel_x"]

        A_e_0d = interface_data_0d["A"]

        interface_data_0d["Thrust"] = m_dot_0d * vel_x_exit_0d * A_e_0d + p_exit_0d * A_e_0d

        ### One dim calculations
        cell_data_1d = FormCellDataFromFile(dataFileName = cell_file_name_1d).cellData
        interface_data_1d = FormInterfaceDataFromFile(dataFileName = interface_file_name_1d).interface_data
        m_dot_1d = interface_data_1d["mass_flux"]
        p_exit_1d = interface_data_1d["p"]

        vel_x_e_1d = cell_data_1d["vel_x"][:-1]
        
        A_e_1d = interface_data_1d["A"]

        interface_data_1d["Thrust"] = m_dot_1d * vel_x_e_1d * A_e_1d + p_exit_1d * A_e_1d

        fig = plt.figure(figsize=(15, 5))
        plt.title("Transient Thrust Profiles of 0- and 1-D simulations")
        plt.ylabel("Thrust (N)", rotation = "horizontal", ha = "right")
        plt.xlabel("Time (ms)")
        plt.scatter(interface_data_1d["time"] * 1e3, interface_data_1d["Thrust"], marker = '.', label = "1-D Thrust")
        plt.scatter(interface_data_0d["time"] * 1e3, interface_data_0d["Thrust"], marker = '.', label = "0-D Thrust")
        filename = "ZeroAndOneDThrustProfileComparison.jpg"
        plt.grid()
        plt.legend()
        mng = plt.get_current_fig_manager()
        mng.full_screen_toggle()
        current_dir = os.getcwd()
        plt.savefig(current_dir + "/plots/" + filename, bbox_inches="tight")
        plt.close()
