from scipy.integrate import ode
from scipy.integrate import odeint
from strain import Strain
from bioreactor import Bioreactor
import yaml
import argparse
 
from simulate import growth_curve_diff_eqs
from simulate import plasmid_loss_diff_eqs

import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Runs a plasmid drop experiment for \
                                                a given plasmid copy number (--PCN) \
                                                and number of passages (--n_passages).\n\n \
                                                Parameters of the plasmid strain are defined in \
                                                P_strain_config.yaml. Parameters of the \
                                                plasmid free strain (wt strain) are defined \
                                                by Q_strain_config.yaml")

parser.add_argument('--n_passages', type=int, help="Number of passages to run")
parser.add_argument('--PCN', type=int, help="Plasmid copy number")
parser.add_argument('--time', type=int, help="Time per passage (mins)", default=600)

args = parser.parse_args()


def generate_integrate_inputs(bioreactor, strain_list):
    species_list = []
    y0 = []

    bioreactor.sample_initial_species()
    B_species_keys, B_y0 = bioreactor.get_initial_species()
    bioreactor.sample_parameters()

    species_list += B_species_keys
    y0 += B_y0

    for strain in strain_list:
        strain.sample_initial_species()
        strain.sample_parameters()
        strain.categorise_species()
        strain_species_keys, strain_y0 = strain.get_initial_species()
        species_list += strain_species_keys
        y0 += strain_y0

    return species_list, y0


def run_growth_curve(y0, t, bioreactor, strain_list, species_keys, plot_species, output_dir):
    """ Runs simple growth curve experiment.

    @param y0 initial species values
    @param t time vector
    @param bioreactor object
    @param strain_list list of strain objects to simulate
    @param species_keys list of species keys
    @param plot_species list of strings referring to which species to plot on the same line graph
    """
    sol = odeint(growth_curve_diff_eqs, y0, t, args=(bioreactor, strain_list, species_keys), mxstep=10**4)

    width_inches = 95*4 / 25.4
    height_inches = 51*4 / 25.4
    fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    
    # Plot selected species
    output_path = output_dir + "growth_curve.pdf"
    for species_str in plot_species:
        P_N_idx = species_keys.index(species_str)
        sns.lineplot(x = t, y = sol[:, P_N_idx], label=species_str, ax=ax)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    fig.tight_layout()
    plt.savefig(output_path, dpi=500)
    plt.close()



def run_plasmid_loss(y0, t, bioreactor, strain_list, species_keys, plasmid_bearing_strain, wt_strain, 
    heterologous_species_name, PCN, n_passages, output_dir):
    """ Runs plasmid drop experiment
    
    Similar to a normal growth experiment but incldes a plasmid drop term 
    defining the rate of transfer between plasmid bearing and wt strains, based on the PCN

    @param y0 initial species values
    @param t time vector
    @param bioreactor object
    @param strain_list list of strain objects to simulate
    @param species_keys list of species keys
    @param plasmid_bearing_strain is a strain object referring to the plasmid bearing strain
    @param wt_strain is a strain object referring to the wt strain
    @param heterologous_species_name is a string referring to heterologously expressed protein (do not include the prefix)
    @param PCN is the plasmid copy number, used to linearly increase the transcription rate of heterologous species
    @param n_passages is the number of passages to conduct
    """

    # Update transcription rate, 'w_x' by multiplication with PCN
    plasmid_bearing_strain.params['w_' + heterologous_species_name] =  float(plasmid_bearing_strain.params['w_' + heterologous_species_name] * PCN)
    
    # Make strain population keys
    plasmid_bearing_key = plasmid_bearing_strain.strain_prefix + '_N'
    wt_key = wt_strain.strain_prefix + '_N'

    plasmid_bearing_ratios = []

    for idx in range(n_passages):
        sol = odeint(plasmid_loss_diff_eqs, y0, t, args=(bioreactor, strain_list, species_keys, plasmid_bearing_strain, wt_strain, PCN), mxstep=10**5)

        P_N_idx = species_keys.index(plasmid_bearing_key)
        Q_N_idx = species_keys.index(wt_key)

        P_N_end = sol[:, P_N_idx][-1]
        Q_N_end = sol[:, Q_N_idx][-1]


        sum_cells = P_N_end + Q_N_end

        probability_plasmid_free = Q_N_end / sum_cells

        plasmid_bearing_ratios.append(P_N_end / sum_cells)


        # Reset initial values by sampling with probability being plasmid free
        y0[Q_N_idx] = 2e5 * probability_plasmid_free
        y0[P_N_idx] = 2e5 - y0[Q_N_idx]



    # Plot passage experiment
    width_inches = 95*4 / 25.4
    height_inches = 51*4 / 25.4
    fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    
    # Plot selected species
    output_name = "PCN_" + str(PCN) + "_loss.pdf"
    output_path = output_dir + output_name

    sns.lineplot(x = range(n_passages), y = plasmid_bearing_ratios, label="Plasmid bearing", ax=ax)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    ax.set_ylabel('Ratio plasmid bearing')
    ax.set_xlabel('N passages')
    ax.set_ylim([0.0, 1.0])
    fig.tight_layout()

    plt.savefig(output_path, dpi=500)


def main():
    # Load bioreactor config file and initialise bioreactor object
    bioreactor_config_yaml = "./bioreactor_config.yaml"
    with open(bioreactor_config_yaml, 'r') as yaml_file:
        bioreactor_config = yaml.load(yaml_file, Loader=yaml.FullLoader)
        bioreactor = Bioreactor(bioreactor_config)

    # Load strain config file and initialise strain object
    strain_config_yaml = "./P_strain_config.yaml"
    with open(strain_config_yaml, 'r') as yaml_file:
        strain_config = yaml.load(yaml_file, Loader=yaml.FullLoader)
        P_strain = Strain(bioreactor, strain_config)

    # Load strain config file and initialise strain object
    strain_config_yaml = "./Q_strain_config.yaml"
    with open(strain_config_yaml, 'r') as yaml_file:
        strain_config = yaml.load(yaml_file, Loader=yaml.FullLoader)
        Q_strain = Strain(bioreactor, strain_config)

    output_dir = "./output/"

    n_passages = args.n_passages
    PCN = args.PCN
    time = args.time

    # Make time vector (minutes)
    t = np.arange(0, time, 0.01)


    # Make list of strains
    strain_list = [P_strain]

    # Generate integration inputs and prepare objects
    species_keys, y0 = generate_integrate_inputs(bioreactor, strain_list)

    plot_species = ['P_N']

    # Run simple growth curve experiment
    run_growth_curve(y0, t, bioreactor, strain_list, species_keys, plot_species, output_dir)
    
    # Make list of strains to be simulated
    strain_list = [P_strain, Q_strain]

    # Generate integration inputs and prepare objects
    species_keys, y0 = generate_integrate_inputs(bioreactor, strain_list)


    # Run plasmid loss experiment 
    run_plasmid_loss(y0, t, bioreactor, strain_list, species_keys, 
        plasmid_bearing_strain=P_strain, wt_strain=Q_strain, 
        heterologous_species_name='h', PCN=PCN, n_passages=n_passages, output_dir=output_dir)


if __name__ == "__main__":
    main()