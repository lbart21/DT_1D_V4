import numpy as np

def bulk_reaction_fast_chemsitry(massf_current, x, massf_goal):
    for ind, val in enumerate(massf_goal):
        if massf_current[ind] + x * (massf_goal[ind] - massf_current[ind]) < 0.0:
            x = min(x, massf_current[ind] / (massf_current[ind] - massf_goal[ind]))
    print("x", x)
    massf_new = (np.array(massf_current) + x * (np.array(massf_goal) - np.array(massf_current))).tolist()

    return massf_new