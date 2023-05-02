"""
Function:
Author: Luke Bartholomew
Edits:
"""
import os
def write_cell_data_to_file(cell, time, flow_property_variables, sim_number, simulation_description):
    cwd = os.getcwd()
    cell_id = cell.cell_id
    file_name = "Sim" + str(sim_number) + "DataAtCellID" + str(cell_id) + ".txt"
    if not os.path.exists(cwd + "/data/" + file_name):
        with open(cwd + "/data/" + file_name, "w") as file:
            cell_flow_data = {}

            if "vel_x" in flow_property_variables:
                cell_flow_data["vel_x"] = cell.flow_state.vel_x

            if "Ma" in flow_property_variables:
                cell_flow_data["Ma"] = cell.flow_state.vel_x / cell.flow_state.fluid_state.a

            if "p" in flow_property_variables:
                cell_flow_data["p"] = cell.flow_state.fluid_state.p

            if "T" in flow_property_variables:
                cell_flow_data["T"] = cell.flow_state.fluid_state.T

            if "rho" in flow_property_variables:
                cell_flow_data["rho"] = cell.flow_state.fluid_state.rho

            if "u" in flow_property_variables:
                cell_flow_data["u"] = cell.flow_state.fluid_state.u

            if "a" in flow_property_variables:
                cell_flow_data["a"] = cell.flow_state.fluid_state.a

            if "h" in flow_property_variables:
                cell_flow_data["h"] = cell.flow_state.fluid_state.enthalpy

            if "A_c" in flow_property_variables:
                cell_flow_data["A_c"] = cell.geo["A_c"]

            if "dV" in flow_property_variables:
                cell_flow_data["dV"] = cell.geo["dV"]

            if "massf" in flow_property_variables:
                species_names = cell.flow_state.fluid_state.gmodel.species_names
                mass_fraction_names = ["massf_" + name for name in species_names]
                for ind, species in enumerate(mass_fraction_names):
                    cell_flow_data[species] = cell.flow_state.fluid_state.massf[ind]

            if "conc" in flow_property_variables:
                gm_cell = cell.flow_state.fluid_state.gmodel
                species_names = gm_cell.species_names
                conc_names = ["conc_" + name for name in species_names]
                rho = cell.flow_state.fluid_state.rho
                for ind, species in enumerate(conc_names):
                    cell_flow_data[species] = cell.flow_state.fluid_state.massf[ind] * rho / gm_cell.mol_masses[ind]

            if "molef" in flow_property_variables:
                species_names = cell.flow_state.fluid_state.gmodel.species_names
                mole_fraction_names = ["molef_" + name for name in species_names]
                for ind, species in enumerate(mole_fraction_names):
                    cell_flow_data[species] = cell.flow_state.fluid_state.molef[ind]

            if "p_t" in flow_property_variables:
                h = cell.flow_state.fluid_state.enthalpy
                vel_x = cell.flow_state.vel_x
                h_t = h + 0.5 * vel_x ** 2.0

                gm_copy_object = cell.flow_state.fluid_state.gmodel.__class__
                gm_filename = cell.flow_state.fluid_state.gmodel.file_name
                gm_copy = gm_copy_object(gm_filename)

                gs_copy_object = cell.flow_state.fluid_state.__class__
                gs_copy = gs_copy_object(gm_copy)

                massf = cell.flow_state.fluid_state.massf
                gs_copy.massf = massf

                s = cell.flow_state.fluid_state.gmodel.entropy(gstate = cell.flow_state.fluid_state)
                gs_copy.update_thermo_from_hs(h = h_t, s = s) 

                p_t = gs_copy.p

                cell_flow_data["p_t"] = p_t

            if "T_t" in flow_property_variables:
                h = cell.flow_state.fluid_state.enthalpy
                vel_x = cell.flow_state.vel_x
                h_t = h + 0.5 * vel_x ** 2.0

                gm_copy_object = cell.flow_state.fluid_state.gmodel.__class__
                gm_filename = cell.flow_state.fluid_state.gmodel.file_name
                gm_copy = gm_copy_object(gm_filename)

                gs_copy_object = cell.flow_state.fluid_state.__class__
                gs_copy = gs_copy_object(gm_copy)

                massf = cell.flow_state.fluid_state.massf
                gs_copy.massf = massf

                s = cell.flow_state.fluid_state.gmodel.entropy(gstate = cell.flow_state.fluid_state)
                gs_copy.update_thermo_from_hs(h = h_t, s = s) 

                T_t = gs_copy.T

                cell_flow_data["T_t"] = T_t

            if "h_t" in flow_property_variables:
                h = cell.flow_state.fluid_state.enthalpy
                vel_x = cell.flow_state.vel_x
                cell_flow_data["h_t"] = h + 0.5 * vel_x ** 2.0

            variable_names = list(cell_flow_data.keys()) #Does not include time
            file.write("Sim: " + str(sim_number) + "\n")
            file.write("ID: " + str(cell_id) + "\n")
            file.write("Description: " + simulation_description + "\n")
            file.write("Variables: " + str(len(variable_names) + 1) + "\n")
            file.write("time " + " ".join(variable_names) + "\n")
            file.write(str(format(time, ".9f")) + " ")
            for name in variable_names:
                file.write(str(format(cell_flow_data[name], ".9f")) + " ")
            file.write("\n")
    else:
        with open(cwd + "/data/" + file_name, "a") as file:
            cell_flow_data = {}

            if "vel_x" in flow_property_variables:
                cell_flow_data["vel_x"] = cell.flow_state.vel_x

            if "Ma" in flow_property_variables:
                cell_flow_data["Ma"] = cell.flow_state.vel_x / cell.flow_state.fluid_state.a

            if "p" in flow_property_variables:
                cell_flow_data["p"] = cell.flow_state.fluid_state.p

            if "T" in flow_property_variables:
                cell_flow_data["T"] = cell.flow_state.fluid_state.T

            if "rho" in flow_property_variables:
                cell_flow_data["rho"] = cell.flow_state.fluid_state.rho

            if "u" in flow_property_variables:
                cell_flow_data["u"] = cell.flow_state.fluid_state.u

            if "a" in flow_property_variables:
                cell_flow_data["a"] = cell.flow_state.fluid_state.a

            if "h" in flow_property_variables:
                cell_flow_data["h"] = cell.flow_state.fluid_state.enthalpy

            if "A_c" in flow_property_variables:
                cell_flow_data["A_c"] = cell.geo["A_c"]

            if "dV" in flow_property_variables:
                cell_flow_data["dV"] = cell.geo["dV"]

            if "massf" in flow_property_variables:
                species_names = cell.flow_state.fluid_state.gmodel.species_names
                mass_fraction_names = ["massf_" + name for name in species_names]
                for ind, species in enumerate(mass_fraction_names):
                    cell_flow_data[species] = cell.flow_state.fluid_state.massf[ind]

            if "conc" in flow_property_variables:
                gm_cell = cell.flow_state.fluid_state.gmodel
                species_names = gm_cell.species_names
                conc_names = ["conc_" + name for name in species_names]
                rho = cell.flow_state.fluid_state.rho
                for ind, species in enumerate(conc_names):
                    cell_flow_data[species] = cell.flow_state.fluid_state.massf[ind] * rho / gm_cell.mol_masses[ind]

            if "molef" in flow_property_variables:
                species_names = cell.flow_state.fluid_state.gmodel.species_names
                mole_fraction_names = ["molef_" + name for name in species_names]
                for ind, species in enumerate(mole_fraction_names):
                    cell_flow_data[species] = cell.flow_state.fluid_state.molef[ind]

            if "p_t" in flow_property_variables:
                h = cell.flow_state.fluid_state.enthalpy
                vel_x = cell.flow_state.vel_x
                h_t = h + 0.5 * vel_x ** 2.0

                gm_copy_object = cell.flow_state.fluid_state.gmodel.__class__
                gm_filename = cell.flow_state.fluid_state.gmodel.file_name
                gm_copy = gm_copy_object(gm_filename)

                gs_copy_object = cell.flow_state.fluid_state.__class__
                gs_copy = gs_copy_object(gm_copy)

                massf = cell.flow_state.fluid_state.massf
                gs_copy.massf = massf

                s = cell.flow_state.fluid_state.gmodel.entropy(gstate = cell.flow_state.fluid_state)
                gs_copy.update_thermo_from_hs(h = h_t, s = s) 

                p_t = gs_copy.p

                cell_flow_data["p_t"] = p_t

            if "T_t" in flow_property_variables:
                h = cell.flow_state.fluid_state.enthalpy
                vel_x = cell.flow_state.vel_x
                h_t = h + 0.5 * vel_x ** 2.0

                gm_copy_object = cell.flow_state.fluid_state.gmodel.__class__
                gm_filename = cell.flow_state.fluid_state.gmodel.file_name
                gm_copy = gm_copy_object(gm_filename)

                gs_copy_object = cell.flow_state.fluid_state.__class__
                gs_copy = gs_copy_object(gm_copy)

                massf = cell.flow_state.fluid_state.massf
                gs_copy.massf = massf

                s = cell.flow_state.fluid_state.gmodel.entropy(gstate = cell.flow_state.fluid_state)
                gs_copy.update_thermo_from_hs(h = h_t, s = s) 

                T_t = gs_copy.T

                cell_flow_data["T_t"] = T_t
            
            if "h_t" in flow_property_variables:
                h = cell.flow_state.fluid_state.enthalpy
                vel_x = cell.flow_state.vel_x
                cell_flow_data["h_t"] = h + 0.5 * vel_x ** 2.0

            variable_names = list(cell_flow_data.keys()) #Does not include time
            file.write(str(format(time, ".9f")) + " ")
            for name in variable_names:
                file.write(str(format(cell_flow_data[name], ".9f")) + " ")
            file.write("\n")
    file.close()