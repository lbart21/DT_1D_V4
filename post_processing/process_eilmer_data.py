"""
Function:
Author: Luke Bartholomew
Edits:
"""
import os
import pandas

class ProcessEilmerData():
    def __init__(self, data_files) -> None:
        """
        dataFiles = list of raw data files
        """
        self.data = {
            "nCells"    : 0
        }
        
        self.t_final = None
        self.n_variables = None
        self.variable_names = []
        self.ReadEilmerFiles(data_files = data_files)

    def ReadEilmerFiles(self, data_files):
        first_file = True
        for file in data_files:
            cwd = os.getcwd()
            f = open(cwd + "/data/" + file)
            data_start = False
            current_file_data = {}
            if first_file:
                variable_flag = -1
                for row_ind, row in enumerate(f):
                    split_row = row.split(": ")
                    if split_row[0] == "sim_time":
                        self.t_final = float(split_row[1])
                    elif split_row[0] == "nicell":
                        self.data["nCells"] += int(split_row[1])
                        n_cells_in_row = int(split_row[1])
                        for name in self.variable_names:
                            self.data[name] = [] 
                            current_file_data[name] = [None] * n_cells_in_row
                    elif split_row[0] == "variables":
                        self.n_variables = int(split_row[1])
                        variable_flag = row_ind + 1
                    elif split_row[0] == "nkcell":
                        data_start = True
                        data_start_ind = row_ind + 1 
                        continue
                    if row_ind == variable_flag:
                        self.variable_names = row[1:-1].replace('"', "").split(" ")
                        self.variable_names = ["pos_x" if name == "pos.x" else name for name in self.variable_names]
                        self.variable_names = ["vel_x" if name == "vel.x" else name for name in self.variable_names]
                    if data_start:
                        if row_ind in range(data_start_ind, data_start_ind + n_cells_in_row):
                            for ind, value in enumerate(row[1:-1].split(" ")):
                                current_file_data[self.variable_names[ind]][row_ind - data_start_ind] = float(value)
                        else:
                            data_start = False
                first_file = False
                for name in self.variable_names:
                    self.data[name] += current_file_data[name]
            else:
                for row_ind, row in enumerate(f):
                    split_row = row.split(": ")
                    if split_row[0] == "nicell":
                        self.data["nCells"] += int(split_row[1])
                        n_cells_in_row = int(split_row[1])
                        for name in self.variable_names:
                            current_file_data[name] = [None] * n_cells_in_row
                    elif split_row[0] == "nkcell":
                        data_start = True
                        data_start_ind = row_ind + 1
                        continue
                    if data_start:
                        if row_ind in range(data_start_ind, data_start_ind + n_cells_in_row):
                            for ind, value in enumerate(row[1:-1].split(" ")):
                                current_file_data[self.variable_names[ind]][row_ind - data_start_ind] = float(value)
                        else:
                            data_start = False
                for name in self.variable_names:
                    self.data[name] += current_file_data[name]
            f.close()
            self.component_data = pandas.DataFrame(columns = self.variable_names, index = range(self.data["nCells"]))
            for cell in range(self.data["nCells"]):
                for var in self.variable_names:
                    self.component_data[var][cell] = self.data[var][cell]
            self.component_data.sort_values(by = ["pos_x"])