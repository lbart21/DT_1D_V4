"""
Function:
Author: Luke Bartholomew
Edits:
"""
import os
import re
import numpy as np

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
                
            elif var == "molef":
                column_names = list(self.data_object.component_data.columns)
                molef_species_names = [species for species in column_names \
                                                if len(species) > 5 and species[:5] == "molef"]
                species_names = [species.split("_")[1] for species in molef_species_names]
                for index, species in enumerate(species_names):
                    split_name = re.findall('(\d+|[A-Za-z]+)', species)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    plt.scatter(self.data_object.component_data["pos_x"], \
                                self.data_object.component_data[molef_species_names[index]], \
                                marker = '.', label = r"$" + formatted_species_name + "$")
                plt.legend()

            elif len(var) > 5 and var[:5] == "molef": # Plot specific mass fractions
                species_name = var[6:] # Trim off starting molef_
                split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
                for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                SI_UNITS[var] = r'$-$'
                SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mole$' + ' ' + r'$Fraction$'
                plt.scatter(self.data_object.component_data["pos_x"], \
                            self.data_object.component_data[var], marker = '.')
            
            elif var == "conc":
                column_names = list(self.data_object.component_data.columns)
                conc_species_names = [species for species in column_names \
                                                if len(species) > 4 and species[:4] == "conc"]
                species_names = [species.split("_")[1] for species in conc_species_names]
                for index, species in enumerate(species_names):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    plt.scatter(self.data_object.component_data["pos_x"], \
                                self.data_object.component_data[conc_species_names[index]], \
                                marker = '.', label = r"$" + formatted_species_name + "$")

            else: # All other properties
                plt.scatter(self.data_object.component_data["pos_x"], self.data_object.component_data[var], marker = '.')
            plt.xlabel("Position (m)")
            plt.ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
            plt.title("Distribution of " + SYMBOLS[var] + " at t = " \
                                                    + formatted_title_time + r'$\mu$' + "s")
            
            filename = "Sim " + str(sim_number) + ' ' + var + " distribution at t = " + formatted_file_name_time + ".jpg"
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

            elif len(var) > 5 and var[:5] == "molef":
                species_name = var[6:] # Trim off starting molef_
                split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
                for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                SI_UNITS[var] = r'$-$'
                SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mole$' + ' ' + r'$Fraction$'

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
                
                elif var == "molef":
                    column_names = list(data_from_files[str(time)].columns)
                    molef_species_names = [species for species in column_names \
                                                if len(species) > 5 and species[:5] == "molef"]
                    species_names = [species.split("_")[1] for species in molef_species_names] 
                    for index, species in enumerate(species_names):
                        split_name = re.findall('(\d+|[A-Za-z]+)', species)
                        for ind, val in enumerate(split_name):
                            if val.isnumeric():
                                split_name[ind] = "_{" + split_name[ind] + "}"
                        formatted_species_name = ''.join(split_name)
                        plt.scatter(data_from_files[str(time)]["pos_x"], \
                                    data_from_files[str(time)][molef_species_names[index]], \
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
            
            elif len(var) > 5 and var[:5] == "molef":
                species_name = var[6:] # Trim off starting molef_
                split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
                for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                SI_UNITS[var] = r'$-$'
                SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mole$' + ' ' + r'$Fraction$'

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
            
            elif var == "molef":
                molef_names = [name for name in cell_data.cell_data.columns if "molef" in name]

                fig = plt.figure(figsize=(15, 5))
                plt.title("Transient Development of " + SYMBOLS[var] + " at Cell " + str(cell_data.cell_id))
                plt.ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                            rotation = "horizontal", ha = "right")
                plt.xlabel("Time (ms)")
                for species in molef_names:
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

class SingleSpatialPlotsWithMultipleYAxes():
    def __init__(self, data_file, plot_vars, visible_axes, plot_number) -> None:
        # plot_vars = [[var_1, var_2, ...], [var_n], etc]
        # visible_axes = [int, int, int, etc] -> len(plot_vars) == len(visible_axes)
        self.data_object = GenerateDataObject(data_file_name = data_file)

        if ["D"] in plot_vars:
            A_c = self.data_object.component_data["A_c"]
            self.data_object.component_data["D"] = (4.0 * A_c / np.pi) ** 0.5

        t_final = self.data_object.t_final
        sim_number = self.data_object.sim_number

        formatted_title_time = '{:.3f}'.format(t_final / 1e-6)
        formatted_file_name_time = '{:.9f}'.format(t_final)

        fig, ax = plt.subplots(figsize=(15, 5))
        axes = [ax] + [ax.twinx() for i in range(len(visible_axes) - 1)]

        # Make some space on the right side for the extra y-axis.
        num_active_axes = 0
        for i in visible_axes:
            num_active_axes += i
            
        fig.subplots_adjust(right=0.9 - 0.05 * (num_active_axes - 1))
        # Move the last y-axis spine over to the right by 20% of the width of the axes
        # To make the border of the right-most axis visible, we need to turn the frame
        # on. This hides the other plots, however, so we need to turn its fill off.
        new_axes_location = 1.0
        right_spine_taken = False
        for plot_num in range(1, len(plot_vars)): #Don't want to touch first variable
            
            if visible_axes[plot_num] == 1:
                if not right_spine_taken: # We're adding the first right-sided axis, don't need to do anything
                    right_spine_taken = True
                else: # Need to offset axes 
                    new_axes_location += 0.1
                    axes[plot_num].spines['right'].set_position(('axes', new_axes_location))
                    axes[plot_num].set_frame_on(True)
                    axes[plot_num].patch.set_visible(False)

            else:
                axes[plot_num].spines['right'].set_visible(False)
                axes[plot_num].get_yaxis().set_visible(False)

        plot_list = []
        marker_list = [".", "x", "+", "*", "|", "_", "^", "s", "o"]
        marker_location = 0
        for ind1, var_list in enumerate(plot_vars):
            if var_list == ["D"]:
                axes[ind1].set_ylim((0.0, 3.0 * np.max(self.data_object.component_data["D"])))
                p_ind, = axes[ind1].plot(self.data_object.component_data["pos_x"], \
                            self.data_object.component_data["D"], \
                            label = SYMBOLS["D"], color = "black")
                if visible_axes[ind1] == 1:
                    plot_list.append(p_ind)
                    axes[ind1].set_ylabel(SYMBOLS["D"] + " (" + SI_UNITS["D"] +")", \
                            ha = "right")
                    
            elif var_list == ["massf"]: # Plot all mass fractions on the one axis
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
                    p_ind = axes[ind1].scatter(self.data_object.component_data["pos_x"], \
                                self.data_object.component_data[massf_species_names[index]], \
                                marker = marker_list[marker_location], label = r"$f_{" + formatted_species_name + "}$")
                    plot_list.append(p_ind)
                    marker_location += 1
                    axes[ind1].set_ylabel(SYMBOLS[var_list[0]] + " (" + SI_UNITS[var_list[0]] +")", \
                                ha = "right")

            elif len(var_list[0]) > 5 and var_list[0][:5] == "massf": # Plot groups of mass fractions on one axis
                massf_label_list = []
                for var in var_list:
                    species_name = var[6:] # Trim off starting massf_
                    split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
                    for ind, val in enumerate(split_name):
                            if val.isnumeric():
                                split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    SI_UNITS[var] = r'$-$'
                    SYMBOLS[var] = r"$f_{" + formatted_species_name + r'}$'
                    p_ind = axes[ind1].scatter(self.data_object.component_data["pos_x"], \
                                self.data_object.component_data[var], marker = marker_list[marker_location], \
                                label = r"$f_{" + formatted_species_name + "}$")
                    plot_list.append(p_ind)
                    massf_label_list.append(SYMBOLS[var])
                    marker_location += 1
                axes[ind1].set_ylabel(', '.join(massf_label_list) + " (" + SI_UNITS[var_list[0]] +")", \
                                ha = "right")

            elif var_list == ["molef"]: # Plot all mole fractions on the one axis
                column_names = list(self.data_object.component_data.columns)
                molef_species_names = [species for species in column_names \
                                                if len(species) > 5 and species[:5] == "molef"]
                species_names = [species.split("_")[1] for species in molef_species_names] 
                for index, species in enumerate(species_names):
                    split_name = re.findall('(\d+|[A-Za-z]+)', species)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    p_ind = axes[ind1].scatter(self.data_object.component_data["pos_x"], \
                                self.data_object.component_data[molef_species_names[index]], \
                                marker = marker_list[marker_location], label = r"$f_{" + formatted_species_name + "}$")
                    plot_list.append(p_ind)
                    marker_location += 1
                    axes[ind1].set_ylabel(SYMBOLS[var_list[0]] + " (" + SI_UNITS[var_list[0]] +")", \
                                ha = "right")

            elif len(var_list[0]) > 5 and var_list[0][:5] == "molef": # Plot groups of mole fractions on one axis
                molef_label_list = []
                for var in var_list:
                    species_name = var[6:] # Trim off starting massf_
                    split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
                    for ind, val in enumerate(split_name):
                            if val.isnumeric():
                                split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    SI_UNITS[var] = r'$-$'
                    SYMBOLS[var] = r"$\chi_{" + formatted_species_name + r'}$'
                    p_ind = axes[ind1].scatter(self.data_object.component_data["pos_x"], \
                                self.data_object.component_data[var], marker = marker_list[marker_location], \
                                label = r"$\chi_{" + formatted_species_name + "}$")
                    plot_list.append(p_ind)
                    molef_label_list.append(SYMBOLS[var])
                    marker_location += 1
                axes[ind1].set_ylabel(', '.join(molef_label_list) + " (" + SI_UNITS[var_list[0]] +")", \
                                ha = "right")

            else:
                var_names = []
                for var in var_list:
                    var_names.append(var)
                    p_ind = axes[ind1].scatter(self.data_object.component_data["pos_x"], \
                            self.data_object.component_data[var], \
                            label = SYMBOLS[var], marker = marker_list[marker_location])
                    marker_location += 1
                    if visible_axes[ind1] == 1:
                        plot_list.append(p_ind)
                if visible_axes[ind1] == 1:
                    axes[ind1].set_ylabel(', '.join([SYMBOLS[var] for var in var_names])     + " (" + SI_UNITS[var_names[0]] +")", \
                                ha = "right")
            
            
        ax.legend(loc = 'upper center', bbox_to_anchor = (0.5, -0.12), handles = plot_list, ncol = len(plot_list))
        
        axes[0].set_xlabel("Position (m)")
        axes[0].grid(visible = True, axis = 'y', which = "both")
        axes[0].minorticks_on()
        plt.title("Distribution of Multiple Variables at t = " \
                                                    + formatted_title_time + r'$\mu$' + "s")
        filename = "Sim " + str(sim_number) + " plot " + str(plot_number) + \
            " multiple y axes spatial distributions at t = " + formatted_file_name_time + ".jpg"
        
        mng = plt.get_current_fig_manager()
        mng.full_screen_toggle()
        current_dir = os.getcwd()
        plt.savefig(current_dir + "/plots/" + filename, bbox_inches="tight")
        plt.close()

class SingleTemporalPlotsWithMultipleYAxes():
    def __init__(self, data_file, plot_vars, visible_axes, plot_number) -> None:
        cell_data = FormCellDataFromFile(data_file_name = data_file)
        sim_number = cell_data.sim_number
        

        fig, ax = plt.subplots(figsize=(15, 5))
        axes = [ax] + [ax.twinx() for i in range(len(visible_axes) - 1)]

        # Make some space on the right side for the extra y-axis.
        num_active_axes = 0
        for i in visible_axes:
            num_active_axes += i
            
        fig.subplots_adjust(right=0.9 - 0.05 * (num_active_axes - 1))

        # Move the last y-axis spine over to the right by 20% of the width of the axes
        # To make the border of the right-most axis visible, we need to turn the frame
        # on. This hides the other plots, however, so we need to turn its fill off.
        new_axes_location = 1.0
        right_spine_taken = False
        for plot_num in range(1, len(plot_vars)): #Don't want to touch first variable
            

            if visible_axes[plot_num] == 1:
                if not right_spine_taken: # We're adding the first right-sided axis, don't need to do anything
                    right_spine_taken = True
                else: # Need to offset axes 
                    new_axes_location += 0.1
                    axes[plot_num].spines['right'].set_position(('axes', new_axes_location))
                    axes[plot_num].set_frame_on(True)
                    axes[plot_num].patch.set_visible(False)

            else:
                axes[plot_num].spines['right'].set_visible(False)
                axes[plot_num].get_yaxis().set_visible(False)

        plot_list = []
        marker_list = [".", "x", "+", "*", "|", "_", "^", "s", "o"]
        marker_location = 0
        for ind1, var_list in enumerate(plot_vars):
            if var_list == ["massf"]: # Plot all mass fractions on the one axis
                column_names = list(cell_data.cell_data.columns)
                massf_species_names = [species for species in cell_data.cell_data \
                                                if len(species) > 5 and species[:5] == "massf"]
                species_names = [species.split("_")[1] for species in massf_species_names]
                for index, species in enumerate(species_names):
                    split_name = re.findall('(\d+|[A-Za-z]+)', species)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    p_ind = axes[ind1].scatter(cell_data.cell_data["time"] * 1e3, \
                                cell_data.cell_data[massf_species_names[index]], \
                                marker = marker_list[marker_location], label = r"$f_{" + formatted_species_name + "}$")
                    plot_list.append(p_ind)
                    marker_location += 1
                    axes[ind1].set_ylabel(SYMBOLS[var_list[0]] + " (" + SI_UNITS[var_list[0]] +")", \
                                ha = "right")
            
            elif len(var_list[0]) > 5 and var_list[0][:5] == "massf": # Plot groups of mass fractions on one axis
                massf_label_list = []
                for var in var_list:
                    species_name = var[6:] # Trim off starting massf_
                    split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
                    for ind, val in enumerate(split_name):
                            if val.isnumeric():
                                split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    SI_UNITS[var] = r'$-$'
                    SYMBOLS[var] = r"$f_{" + formatted_species_name + r'}$'
                    p_ind = axes[ind1].scatter(cell_data.cell_data["time"] * 1e3, \
                                cell_data.cell_data[var], marker = marker_list[marker_location], \
                                label = r"$f_{" + formatted_species_name + "}$")
                    plot_list.append(p_ind)
                    massf_label_list.append(SYMBOLS[var])
                    marker_location += 1
                axes[ind1].set_ylabel(', '.join(massf_label_list) + " (" + SI_UNITS[var_list[0]] +")", \
                                ha = "right")
            
            elif var_list == ["molef"]: # Plot all mole fractions on the one axis
                column_names = list(cell_data.cell_data.columns)
                molef_species_names = [species for species in column_names \
                                                if len(species) > 5 and species[:5] == "molef"]
                species_names = [species.split("_")[1] for species in molef_species_names] 
                for index, species in enumerate(species_names):
                    split_name = re.findall('(\d+|[A-Za-z]+)', species)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    p_ind = axes[ind1].scatter(cell_data.cell_data["time"] * 1e3, \
                                cell_data.cell_data[molef_species_names[index]], \
                                marker = marker_list[marker_location], label = r"$f_{" + formatted_species_name + "}$")
                    plot_list.append(p_ind)
                    marker_location += 1
                    axes[ind1].set_ylabel(SYMBOLS[var_list[0]] + " (" + SI_UNITS[var_list[0]] +")", \
                                ha = "right")

            elif len(var_list[0]) > 5 and var_list[0][:5] == "molef": # Plot groups of mole fractions on one axis
                molef_label_list = []
                for var in var_list:
                    species_name = var[6:] # Trim off starting massf_
                    split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
                    for ind, val in enumerate(split_name):
                            if val.isnumeric():
                                split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    SI_UNITS[var] = r'$-$'
                    SYMBOLS[var] = r"$\chi_{" + formatted_species_name + r'}$'
                    p_ind = axes[ind1].scatter(cell_data.cell_data["time"] * 1e3, \
                                cell_data.cell_data[var], marker = marker_list[marker_location], \
                                label = r"$\chi_{" + formatted_species_name + "}$")
                    plot_list.append(p_ind)
                    molef_label_list.append(SYMBOLS[var])
                    marker_location += 1
                axes[ind1].set_ylabel(', '.join(molef_label_list) + " (" + SI_UNITS[var_list[0]] +")", \
                                ha = "right")

            elif var_list == ["conc"]:
                column_names = list(cell_data.cell_data.columns)
                conc_species_names = [species for species in column_names \
                                                if len(species) > 4 and species[:4] == "conc"]
                
                species_names = [species.split("_")[1] for species in conc_species_names] 
                for index, species in enumerate(species_names):
                    split_name = re.findall('(\d+|[A-Za-z]+)', species)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    p_ind = axes[ind1].scatter(cell_data.cell_data["time"] * 1e3, \
                                cell_data.cell_data[conc_species_names[index]], \
                                marker = marker_list[marker_location], label = r"$f_{" + formatted_species_name + "}$")
                    plot_list.append(p_ind)
                    marker_location += 1
                    axes[ind1].set_ylabel(SYMBOLS[var_list[0]] + " (" + SI_UNITS[var_list[0]] +")", \
                                ha = "right")
            
            elif len(var_list[0]) > 4 and var_list[0][:4] == "conc": # Plot groups of concentrations on one axis
                conc_label_list = []
                for var in var_list:
                    species_name = var[5:] # Trim off starting conc_
                    split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
                    for ind, val in enumerate(split_name):
                            if val.isnumeric():
                                split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    SI_UNITS[var] = r'$moleL^{-1}$'
                    SYMBOLS[var] = r"$C_{" + formatted_species_name + r'}$'
                    p_ind = axes[ind1].scatter(cell_data.cell_data["time"] * 1e3, \
                                cell_data.cell_data[var], marker = marker_list[marker_location], \
                                label = r"$C_{" + formatted_species_name + "}$")
                    plot_list.append(p_ind)
                    conc_label_list.append(SYMBOLS[var])
                    marker_location += 1
                axes[ind1].set_ylabel(', '.join(conc_label_list) + " (" + SI_UNITS[var_list[0]] +")", \
                                ha = "right")

            else:
                var_names = []
                for var in var_list:
                    var_names.append(var)
                    p_ind = axes[ind1].scatter(cell_data.cell_data["time"] * 1e3, \
                            cell_data.cell_data[var], \
                            label = SYMBOLS[var], marker = marker_list[marker_location])
                    marker_location += 1
                    if visible_axes[ind1] == 1:
                        plot_list.append(p_ind)
                if visible_axes[ind1] == 1:
                    axes[ind1].set_ylabel(', '.join([SYMBOLS[var] for var in var_names])     + " (" + SI_UNITS[var_names[0]] +")", \
                                ha = "right")
                    
        ax.legend(loc = 'upper center', bbox_to_anchor = (0.5, -0.12), handles = plot_list, ncol = len(plot_list))
        
        axes[0].set_xlabel("time (ms)")
        axes[0].grid(visible = True, axis = 'y', which = "both")
        axes[0].minorticks_on()
        plt.title("Transient Development of Various Quantities at Cell " + str(cell_data.cell_id))
        filename = "Sim " + str(sim_number) + " plot " + str(plot_number) + \
            " multiple y axes transient development of various properties.jpg"
        
        mng = plt.get_current_fig_manager()
        mng.full_screen_toggle()
        current_dir = os.getcwd()
        plt.savefig(current_dir + "/plots/" + filename, bbox_inches="tight")
        plt.close()

class GenerateCustomSpatialPlotsFromMultipleSims():
    def __init__(self, spatial_data_files, custom_labels, plot_vars, convergence_parameter) -> None:
        data_from_files = [None] * len(spatial_data_files)

        for ind, file in enumerate(spatial_data_files):
            data = GenerateDataObject(data_file_name = file)
            data_from_files[ind] = data.component_data
            
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

            elif len(var) > 5 and var[:5] == "molef":
                species_name = var[6:] # Trim off starting molef_
                split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
                for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                SI_UNITS[var] = r'$-$'
                SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mole$' + ' ' + r'$Fraction$'

            for ind, data in enumerate(data_from_files):
                plt.scatter(data["pos_x"], data[var], \
                                label = custom_labels[ind], \
                                marker = ".")
            
            plt.title("Spatial convergence of " + SYMBOLS[var] + " against " + convergence_parameter)
            plt.ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            plt.legend()
            plt.grid()
            filename = "Spatial convergence of " + var + " against " + convergence_parameter + ".jpg"
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            current_dir = os.getcwd()
            plt.savefig(current_dir + "/plots/" + filename, bbox_inches="tight")
            plt.close()
        

class GenerateCustomCellTemporalPlotsFromMultipleSims():
    def __init__(self, cell_data_files, custom_labels, plot_vars, convergence_parameter) -> None:
        data_from_files = [None] * len(cell_data_files)

        for ind, file in enumerate(cell_data_files):
            data = FormCellDataFromFile(data_file_name = file)
            data_from_files[ind] = data
        
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

            elif len(var) > 5 and var[:5] == "molef":
                species_name = var[6:] # Trim off starting molef_
                split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
                for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                SI_UNITS[var] = r'$-$'
                SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mole$' + ' ' + r'$Fraction$'

            for ind, data in enumerate(data_from_files):
                plt.scatter(data.cell_data["time"], data.cell_data[var], \
                                label = custom_labels[ind], \
                                marker = ".")
            
            plt.title("Temporal convergence of " + SYMBOLS[var] + " against " + convergence_parameter)
            plt.ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            plt.legend()
            plt.grid()
            filename = "Temporal convergence of " + var + " against " + convergence_parameter + " at Cell" + str(data.cell_id) + ".jpg"
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            current_dir = os.getcwd()
            plt.savefig(current_dir + "/plots/" + filename, bbox_inches="tight")
            plt.close()


class GenerateCustomInterfaceTemporalPlotsFromMultipleSims():
    def __init__(self, interface_data_files, custom_labels, plot_vars, convergence_parameter) -> None:
        data_from_files = [None] * len(interface_data_files)

        for ind, file in enumerate(interface_data_files):
            data = FormInterfaceDataFromFile(data_file_name = file)
            data_from_files[ind] = data
        
        for var in plot_vars:
            fig = plt.figure(figsize=(15, 5))
            
            if var == "mass_flux":
                for data in data_from_files:
                    data.interface_data["mass_flux"] *= data.interface_data["A"]
            if var == "energy_flux":
                for data in data_from_files:
                    data.interface_data["energy_flux"] *= data.interface_data["A"]

            for ind, data in enumerate(data_from_files):
                plt.scatter(data.interface_data["time"], data.interface_data[var], \
                                label = custom_labels[ind], \
                                marker = ".")
            
            plt.title("Temporal convergence of " + SYMBOLS[var] + " against " + convergence_parameter)
            plt.ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            plt.legend()
            plt.grid()
            filename = "Temporal convergence of " + var + " against " + convergence_parameter + " at Interface" + str(data.interface_id) + ".jpg"
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            current_dir = os.getcwd()
            plt.savefig(current_dir + "/plots/" + filename, bbox_inches="tight")
            plt.close()
