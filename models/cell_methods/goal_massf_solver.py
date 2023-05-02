import numpy as np

def form_goal_massf_solver_matrices(massf_current, bulk_reaction_parameters, molar_masses, species_names):
    [_, specified_massf_dict, bulk_reaction] = bulk_reaction_parameters
    n_sp = len(molar_masses)
    n_constraints = len(bulk_reaction.elements)

    LHS_A = np.zeros((n_sp, n_sp))
    RHS_A = np.zeros((n_sp, n_sp))
    source = np.zeros((n_sp, 1))

    for row_ind, element in enumerate(bulk_reaction.elements):
        for column_ind, spcs in enumerate(species_names):
            C_idx = bulk_reaction.species_names.index(spcs)
            C_i_j = bulk_reaction.species[C_idx].base_representation[element]
            M_idx = species_names.index(spcs)
            M_i = molar_masses[M_idx]
            RHS_A[row_ind, column_ind] = C_i_j / M_i
            LHS_A[row_ind, column_ind] = C_i_j / M_i
    

    for row_ind_2, spcs_given in enumerate(specified_massf_dict.keys()):
        column_ind_2 = species_names.index(spcs_given)
        LHS_A[row_ind_2 + n_constraints, column_ind_2] = 1.0
        source[row_ind_2 + n_constraints, 0] = specified_massf_dict[spcs_given]

    return LHS_A, RHS_A, source

def goal_massf_solver(massf_current, bulk_reaction_parameters, species_names, \
                        LHS_A, RHS_A, source):

    RHS_total = np.matmul(RHS_A, np.transpose(np.array([massf_current]))) + source

    goal_massf = np.linalg.solve(LHS_A, RHS_total)
    goal_massf = np.transpose(goal_massf)[0].tolist()

    return goal_massf