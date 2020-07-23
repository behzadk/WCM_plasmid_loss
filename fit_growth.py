from scipy.integrate import ode
from scipy.integrate import odeint
from strain import Strain
from bioreactor import Bioreactor
import oyaml as  yaml

from simulate import growth_curve_diff_eqs
import numpy as np
import utils
import scipy.optimize as optimize
import time


def growth_optimize_func(free_params, bioreactor_config, strain_config, target_6hr_ratio=500, target_final_6hr_ratio=1):
    """ Runs optimises growth parameters
    Uses target ratios (hard coded) 

    @param free_params vector of params being fit, these should be unpacked
    @param target_6hr_ratio target cell density ratio for 6hr density / initial density
    @param target_final_6hr_ratio target cell density ratio for 24hr density / 6hr density
    """

    # Unpack parameters
    B_s = free_params[0] * float(10**10)
    n_s = free_params[1]

    # Initialise bioreactor object
    bioreactor = Bioreactor(bioreactor_config)
    
    # Initialise strain object
    Q_strain = Strain(bioreactor, strain_config)
    
    # Make list of strains
    strain_list = [Q_strain]
    species_keys, y0 = utils.generate_integrate_inputs(bioreactor, strain_list)

    # Set bioreactor initial S
    B_s_idx = species_keys.index('B_s')
    y0[B_s_idx] = float(B_s)

    # Set strain nutrient efficiency
    Q_strain.params['n_s'] = float(n_s)

    # Make time vector (minutes)
    t = np.arange(0, 1440, 0.01)

    # Simulate growth
    sol = odeint(growth_curve_diff_eqs, y0, t, 
        args=(bioreactor, strain_list, species_keys), mxstep=10**4)    

    species_str = 'Q_N'
    Q_N_idx = species_keys.index(species_str)

    # Get index of six hour time point
    six_hr_idx = np.argwhere(t==390)[0][0]

    # Get cell counts for timepoints of interest
    cell_count_0_hrs = sol[:, Q_N_idx][0]
    cell_count_6_hrs = sol[:, Q_N_idx][six_hr_idx]
    cell_count_24_hrs = sol[:, Q_N_idx][-1]

    # Get fold changes
    initial_6hr_ratio =  cell_count_6_hrs / cell_count_0_hrs
    final_six_hr_ratio = cell_count_6_hrs / cell_count_24_hrs

    objective_1_distance = (initial_6hr_ratio - target_6hr_ratio) ** 2
    objective_2_distance = (final_six_hr_ratio - target_final_6hr_ratio) ** 2

    sum_distance = objective_1_distance + objective_2_distance

    # print(cell_count_24_hrs)
    print("free_params: ", free_params)
    print("6hr / initial density distance: ", objective_1_distance)
    print("6hr / 24hr density distance: ", objective_2_distance)
    print("sum distance: ", sum_distance)
    print("")

    return sum_distance

def fit_growth_parameters():
    """ Fits parameters for growth mimicking batch plasmid loss criteria
    
    Bioreactor nutrient concentration (B_s) and strain nutrient efficiency (n_s)
    are fit to meet batch growth criteria. Scipy Basinhopping is used to optimise for 
    fold change in cell density at 6hr and 24hr time points
    """

    # Set initial guesses
    init_B_s = 9.75 # B_s is scaled up in the optimisation function B_s * 10 ** 10
    init_n_s = 300
    initial_guess = [init_B_s, init_n_s]

    # Set maximum iterations of minimisation routine
    max_iter = 50


    bioreactor_config_yaml = "./bioreactor_config.yaml"
    with open(bioreactor_config_yaml, 'r') as yaml_file:
        bioreactor_config = yaml.load(yaml_file, Loader=yaml.FullLoader)

    strain_config_yaml = "./Q_strain_config.yaml"
    with open(strain_config_yaml, 'r') as yaml_file:
        strain_config = yaml.load(yaml_file, Loader=yaml.FullLoader)

    res = optimize.minimize(
        growth_optimize_func, 
        initial_guess,
        method='Nelder-Mead',
        args=(bioreactor_config, strain_config),
        options={'maxiter': max_iter}
        )

    print("Finished: ")
    print(res)
    print("")
    print("")

    print("Fitted B_s value: ", res['x'][0], "E+10")
    print("Fitted n_s value: ", res['x'][1])



if __name__ == "__main__":
    fit_growth_parameters()

