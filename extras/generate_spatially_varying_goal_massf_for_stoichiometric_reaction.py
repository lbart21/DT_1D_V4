import numpy as np

def generate_spatially_varying_goal_massf_for_stoichiometric_reaction(profile, \
                    reaction, x, sp_name, n_cells):
    """
    profile = ["linear", x_pos_list, L, f_left, f_right]
    profile = ["constant", f_goal]
    profile = ["tanh", x_pos_list, L, f_left, f_right, beta, x_mid]
    profile = ["custom", goal_massf_list]
    """
    goal_massf_data = [None] * n_cells
    if profile[0] == "linear":
        [_, x_pos_list, L, f_left, f_right] = profile
        goal_massf_array = f_left + (f_right - f_left) * np.array(x_pos_list) / L
        goal_massf_list = goal_massf_array.tolist()
    
    elif profile[0] == "constant":
        [_, f_goal] = profile
        goal_massf_list = [f_goal] * n_cells
    
    elif profile[0] == "tanh":
        [_, x_pos_list, L, f_left, f_right, beta, x_mid] = profile
        goal_massf_array = f_left + 0.5 * (f_right - f_left) * (1.0 + np.tanh(beta * (np.array(x_pos_list) - x_mid) / L))
        goal_massf_list = goal_massf_array.tolist()

    elif profile[0] == "custom":
        [_, goal_massf_list] = profile

    for i in range(n_cells):
        goal_massf_data[i] = [sp_name, goal_massf_list[i], x, reaction]
    
    return goal_massf_data
