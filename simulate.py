def growth_curve_diff_eqs(y, t, bioreactor, strain_list, species_key_list):
    """ Function to simulate bioreactor and strains
    compatible with scipy.odeint

    @param y
    @param t
    @param bioreactor object being simulated.
    @param strain_list is a list of strain objects being simulated.
    @species_key_list is a list of species defining the order they should be output.
    """

    # Create dictionary of species and values
    current_species_dict = dict(zip(species_key_list, y))

    # Update the object species values with the current y
    bioreactor.update_current_species_values(current_species_dict)
    for strain in strain_list:
        strain.update_current_species_values(current_species_dict)

    # Initialise solution dict
    sol_dict = {}

    # Calculate bioreactor differentials and update dictionary
    sol_dict.update(bioreactor.calculate_differentials(strain_list))

    # Calculate differentials for each strain and update dictionary
    for strain in strain_list:
        sol_dict.update(strain.calculate_differentials())


    # Load solutions into dictionary in correct order
    # defined by species_key_list
    y_sol = []
    for key in species_key_list:
        y_sol.append(sol_dict[key])

    return y_sol


def plasmid_loss_diff_eqs(y, t, bioreactor, strain_list, species_key_list, plasmid_bearing_strain, wt_strain, PCN):
    """ Function to simulate bioreactor and strains with plasmid drop
    compatible with scipy.odeint

    @param y
    @param t
    @param bioreactor object being simulated.
    @param strain_list is a list of strain objects being simulated.
    @param species_key_list is a list of species defining the order they should be output.
    @param plasmid_bearint_strain is an object of the plasmid bearing strain
    @param wt_strain is an object of the wt strain
    @param PCN is the plasmid copy number, used in the drop rate term.

    @return solution of the time step derivatives
    """
    
    # Create dictionary of species and values
    current_species_dict = dict(zip(species_key_list, y))

    # Update the object species values with the current y
    bioreactor.update_current_species_values(current_species_dict)
    for strain in strain_list:
        strain.update_current_species_values(current_species_dict)

    # Initialise solution dict
    sol_dict = {}
    sol_dict.update(bioreactor.calculate_differentials(strain_list))

    # Calculate differentials for each strain and update dictionary
    for strain in strain_list:
        sol_dict.update(strain.calculate_differentials())

    # Calculate plasmid drop rates
    plasmid_bearing_key = plasmid_bearing_strain.strain_prefix + '_N'
    wt_key = wt_strain.strain_prefix + '_N'

    # Calculate plasmid bearing strain growth rate
    P_mass = plasmid_bearing_strain.calculate_mass()
    P_gamma_a = plasmid_bearing_strain.calculate_gamma_a()
    lambd_plasmid_bearing = plasmid_bearing_strain.calculate_growth_rate(P_mass, P_gamma_a)

    # Calculate drop rate
    N_plasmid_bearing = sol_dict[plasmid_bearing_key]
    drop_rate = 2 ** (1 - PCN) * lambd_plasmid_bearing * N_plasmid_bearing

    # Add drop rate term to differential equations
    sol_dict[plasmid_bearing_key] = sol_dict[plasmid_bearing_key] - drop_rate
    sol_dict[wt_key] = sol_dict[wt_key] + drop_rate


    # Load solutions into dictionary in correct order
    # defined by species_key_list
    y_sol = []
    for key in species_key_list:
        y_sol.append(sol_dict[key])

    return y_sol

