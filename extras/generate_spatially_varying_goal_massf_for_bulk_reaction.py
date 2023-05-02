import numpy as np
def generate_spatially_varying_goal_massf_for_bulk_reaction(profile, reaction, x, n_cells):
    """
    profile = ["constant", goal_massf_dict]
    profile = ["linear", goal_massf_dict_left, goal_massf_dict_right, x_pos_list, L]
    profile = ["tanh", goal_massf_dict_left, goal_massf_dict_right, x_pos_list, L, beta, x_mid]
    profile = ["custom", goal_massf_list]
    """

    goal_massf_data = [None] * n_cells

    if profile[0] == "constant":
        [_, goal_massf_dict] = profile
        goal_massf_list = [goal_massf_dict] * n_cells #list of dicts

    elif profile[0] == "linear":
        [_, goal_massf_dict_left, goal_massf_dict_right, x_pos_list, L] = profile
        goal_massf_lists = [None] * len(goal_massf_dict_left.keys())
        for ind, val in enumerate(goal_massf_dict_left.keys()):
            goal_massf_lists[ind] = goal_massf_dict_left[val] + (goal_massf_dict_right[val] - goal_massf_dict_left[val]) * np.array(x_pos_list) / L
        goal_massf_list = [{spcs:goal_massf_lists[ind][cell] for ind, spcs in enumerate(goal_massf_dict_left.keys())} for cell in range(n_cells)]

    elif profile[0] == "tanh":
        [_, goal_massf_dict_left, goal_massf_dict_right, x_pos_list, L, beta, x_mid] = profile
        goal_massf_lists = [None] * len(goal_massf_dict_left.keys())
        for ind, val in enumerate(goal_massf_dict_left.keys()):
            goal_massf_lists[ind] = goal_massf_dict_left[val] + 0.5 * (goal_massf_dict_right[val] - goal_massf_dict_left[val]) * (1.0 + np.tanh(beta * (np.array(x_pos_list) - x_mid) / L))

    elif profile[0] == "custom":
        [_, goal_massf_list] = profile

    for i in range(n_cells):
        goal_massf_data[i] = [x, goal_massf_list[i], x, reaction]
    
    return goal_massf_data