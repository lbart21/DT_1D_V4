"""
Function:
Author: Luke Bartholomew
Edits:
"""
import pandas
import os

class FormInterfaceDataFromFile():
    def __init__(self, data_file_name) -> None:
        cwd = os.getcwd()
        file_length = sum(1 for line in open(cwd + "/data/" + data_file_name))
        file = open(cwd + "/data/" + data_file_name)
        variable_flag = -1
        data_start = False
        for row_ind, row in enumerate(file):
            split_row = row.split(": ")
            if split_row[0] == "Variables":
                variable_flag = row_ind + 1
            elif split_row[0] == "ID":
                self.interface_id = int(split_row[1])
            elif split_row[0] == "Sim":
                self.sim_number = int(split_row[1])
            if row_ind == variable_flag:
                variable_names = [name for name in row[:-1].split(" ")]
                data_start_ind = row_ind + 1
                n_times = file_length - data_start_ind
                self.interface_data = pandas.DataFrame(columns = variable_names, index = range(n_times))
                data_start = True
                continue
            if data_start:
                for ind, name in enumerate(variable_names):
                    self.interface_data[name][row_ind - data_start_ind] = float(row.split(' ')[ind])
        file.close()
        self.interface_data = self.interface_data.sort_values(by = ["time"])