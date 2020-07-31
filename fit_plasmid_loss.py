from scipy.integrate import ode
from scipy.integrate import odeint
from WCM.strain import Strain
from WCM.bioreactor import Bioreactor
import oyaml as  yaml

from WCM.simulate import growth_curve_diff_eqs
from WCM.simulate import plasmid_loss_diff_eqs

import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import WCM.utils as utils
import scipy.optimize as optimize

import pandas as pd
import glob
import os

def run_plasmid_loss_passage_experiment(y0, species_keys, bioreactor, plasmid_bearing_strain, wt_strain, t, n_passages, PCN):
    print("Running plasmid loss experiment... ")
    strain_list = [plasmid_bearing_strain, wt_strain]
    
    # Make strain population keys
    plasmid_bearing_key = plasmid_bearing_strain.strain_prefix + '_N'
    wt_key = wt_strain.strain_prefix + '_N'

    plasmid_bearing_ratios = []

    for idx in range(n_passages):
        P_N_idx = species_keys.index(plasmid_bearing_key)
        Q_N_idx = species_keys.index(wt_key)

        sol = odeint(plasmid_loss_diff_eqs, y0, t, args=(bioreactor, strain_list, species_keys, plasmid_bearing_strain, wt_strain, PCN), mxstep=10**4)

        P_N_end = sol[:, P_N_idx][-1]
        Q_N_end = sol[:, Q_N_idx][-1]

        sum_cells = P_N_end + Q_N_end

        probability_plasmid_free = Q_N_end / sum_cells

        plasmid_bearing_ratios.append(P_N_end / sum_cells)

        # Reset initial values by sampling with probability being plasmid free
        y0[Q_N_idx] = 1000.0 * probability_plasmid_free
        y0[P_N_idx] = 1000.0 - y0[Q_N_idx]

    return plasmid_bearing_ratios


def data_preprocessing(df):
    event_num_samples_threshold = 1000

    # Drop sub threshold events
    df = df[df['num_samples'] > 5000]

    # Group biological replicates and get mean of technical replicates for
    # each passage
    df = df.groupby(['bio_rpt', 'passage'], as_index=False)['max_clust_prop'].mean()

    df = df.groupby(['passage'], as_index=False)['max_clust_prop'].mean()

    df.rename(columns={'max_clust_prop': 'mean_max_clust_prop'}, inplace=True)

    return df


def plot_simulation_bio_rpt_data_overlay(data_dict, sim_passage_data, heterologous_species_name, w_val, init_Q, PCN, output_dir="./output/"):
    data_path = data_dict['path']

    data_df = pd.read_csv(data_path)
    grouped = data_df.groupby(['bio_rpt', 'passage'], as_index=False).agg({'max_clust_prop': ['mean', 'std']}).reset_index()
    grouped['mean_max_clust_prop'] = grouped['max_clust_prop']['mean']
    grouped['std_max_clust_prop'] = grouped['max_clust_prop']['std']

    width_inches = 95*4 / 25.4
    height_inches = 51*4 / 25.4
    fig, ax = plt.subplots(figsize=(width_inches, height_inches))

    output_name = str(data_dict['plasmid_name']) + "_PCN_" + str(PCN) + "_w_" + str(w_val) + ".pdf"
    output_path = output_dir + "/overlay_" + output_name

    plot_title = "w_" + heterologous_species_name + ": " + str(w_val) + ", PCN: " + str(PCN) + ", init_Q: " + str(init_Q)
    sns.lineplot(data=grouped[grouped['bio_rpt'] == 1], x='passage', y='mean_max_clust_prop', label='bio_rpt_1')
    sns.lineplot(data=grouped[grouped['bio_rpt'] == 2], x='passage', y='mean_max_clust_prop', label='bio_rpt_2')
    sns.lineplot(data=grouped[grouped['bio_rpt'] == 2], x='passage', y=sim_passage_data, color='black', label='simulation')
    ax.errorbar(
        x=grouped[grouped['bio_rpt'] == 1]['passage'], 
        y=grouped[grouped['bio_rpt'] == 1]['mean_max_clust_prop'], 
        yerr=grouped[grouped['bio_rpt'] == 1]['std_max_clust_prop'], 
        fmt='-o', alpha=0.3
        )

    ax.errorbar(
        x=grouped[grouped['bio_rpt'] == 2]['passage'], 
        y=grouped[grouped['bio_rpt'] == 2]['mean_max_clust_prop'], 
        yerr=grouped[grouped['bio_rpt'] == 2]['std_max_clust_prop'], 
        fmt='-o', alpha=0.3
        )
    
    ax.set_ylabel('Ratio plasmid bearing')
    ax.set_xlabel('N passages')
    ax.set_title(plot_title)
    fig.tight_layout()

    plt.savefig(output_path, dpi=500)
    

def plasmid_loss_optim_func(free_params, bioreactor_config, plasmid_bearing_config, wt_config, heterologous_species_name, data_dict):
    # Unpack parameters
    w = free_params[0]

    t = np.arange(0, 1440, 0.5)

    distance = 0

    # Iterate all data sets calculating the distances from
    # simulations and the data
    for key in data_dict:
        PCN = data_dict[key]['PCN']
        data_df = data_dict[key]['processed_df']
        n_passages = np.shape(data_df['passage'])[0]

        # Initialise bioreactor object
        bioreactor = Bioreactor(bioreactor_config)
        
        # Initialise plasmid bearing strain
        plasmid_bearing_strain = Strain(bioreactor, plasmid_bearing_config)

        # Initialise wt strain
        wt_strain = Strain(bioreactor, wt_config)


        # Sample parameters and initial species from prior
        strain_list = [plasmid_bearing_strain, wt_strain]
        species_keys, y0 = utils.generate_integrate_inputs(bioreactor, strain_list)

        # Set heterologous gene transcription rate (w)
        # Update transcription rate, 'w_x' by multiplication with PCN
        plasmid_bearing_strain.params['w_' + heterologous_species_name] =  float(free_params[0] * PCN)#
        y0, species_keys, bioreactor, plasmid_bearing_strain, wt_strain, t, n_passages, PCN
        sim_plasmid_bearing_props = run_plasmid_loss_passage_experiment(y0=y0, species_keys=species_keys, bioreactor=bioreactor, 
            plasmid_bearing_strain=plasmid_bearing_strain, wt_strain=wt_strain, t=t, n_passages=n_passages, PCN=PCN)

        for real_data, sim_data in (zip(data_df['mean_max_clust_prop'].values, sim_plasmid_bearing_props)):
            distance += abs(real_data - sim_data)

        sim_data.append(sim_plasmid_bearing_props)

    return distance



def make_data_dict(data_dir="./data/"):
    SC_101_df = pd.read_csv(data_dir + "/SC101_crude_clusters.csv")
    SC_101_df = data_preprocessing(SC_101_df)

    pUC_df = pd.read_csv(data_dir + "/pUC_crude_clusters.csv")
    pUC_df = data_preprocessing(pUC_df)

    p15A_df = pd.read_csv(data_dir + "/p15A_crude_clusters.csv")
    p15A_df = data_preprocessing(p15A_df)

    data_sets = {
    'SC101': {'plasmid_name': 'SC101', 'prior_PCN': [5.0, 5.0], 'processed_df': SC_101_df, 'path': data_dir + "SC101_crude_clusters.csv"}, 
    'pUC': {'plasmid_name': 'pUC', 'prior_PCN': [13.0, 13.0], 'processed_df': pUC_df, 'path': data_dir + "pUC_crude_clusters.csv"},
    'p15A': {'plasmid_name': 'p15A', 'prior_PCN': [7.0, 7.0], 'processed_df': p15A_df, 'path': data_dir + "p15A_crude_clusters.csv"}
    }

    return data_sets

def plot_simulation_fitting_data_overlay(data_dict, bioreactor_config, plasmid_bearing_config, wt_config, heterologous_species_name, w_val, init_Q_N, PCN, output_dir="./output/"):
    data_df = data_dict['processed_df']
    n_passages = np.shape(data_df['passage'])[0]
    t = np.arange(0, 1440, 1.0)

    # Initialise bioreactor object
    bioreactor = Bioreactor(bioreactor_config)
    
    # Initialise plasmid bearing strain
    plasmid_bearing_strain = Strain(bioreactor, plasmid_bearing_config)

    # Initialise wt strain
    wt_strain = Strain(bioreactor, wt_config)


    strain_list = [plasmid_bearing_strain, wt_strain]

    # Sample parameters and initial species from prior
    # Generate integration inputs and prepare objects
    species_keys, y0 = utils.generate_integrate_inputs(bioreactor, strain_list)

    # Set plasmid bearing strain w value
    plasmid_bearing_strain.params['w_' + heterologous_species_name] =  float(w_val * PCN)
    
    # Set Q strain initial N
    Q_N_idx = species_keys.index('Q_N')
    y0[Q_N_idx] = float(init_Q_N)

    # Set P strain initial N 
    P_N_idx = species_keys.index('P_N')
    y0[P_N_idx] = 1000 - y0[Q_N_idx]

    initial_plasmid_bearing_prop = y0[P_N_idx] / (y0[P_N_idx] + y0[Q_N_idx])


    # Set heterologous gene transcription rate (w)
    # Update transcription rate, 'w_x' by multiplication with PCN

    sim_plasmid_bearing_props = run_plasmid_loss_passage_experiment(
        y0=y0[:], species_keys=species_keys, bioreactor=bioreactor, 
        plasmid_bearing_strain=plasmid_bearing_strain, wt_strain=wt_strain, 
        t=t, n_passages=n_passages, PCN=PCN
        )

    print(sim_plasmid_bearing_props)
    # Plot passage experiment
    width_inches = 95*4 / 25.4
    height_inches = 51*4 / 25.4
    fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    
    # Plot selected species
    output_name = str(data_dict['plasmid_name']) + "_PCN_" + str(PCN) + "_w_" + str(w_val) + ".pdf"
    output_path = output_dir + "/fitting_" + output_name

    plot_title = "w_" + heterologous_species_name + ": " + str(w_val) + ", PCN: " + str(PCN) + ", init_Q: " + str(init_Q_N)

    sns.lineplot(x=data_df['passage'].values, y=sim_plasmid_bearing_props, ax=ax)
    sns.scatterplot(data=data_df, x='passage', y='mean_max_clust_prop', ax=ax)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    ax.set_ylabel('Ratio plasmid bearing')
    ax.set_xlabel('N passages')
    ax.set_title(plot_title)
    fig.tight_layout()

    plt.savefig(output_path, dpi=500)

    return sim_plasmid_bearing_props


def grid_search_plasmid_drop(data_dict, bioreactor_config, plasmid_bearing_config, 
    wt_config, heterologous_species_name, w_values, init_Q_values, PCN_values, output_dir):
    output_data = {}
    output_data['plasmid_name'] = []
    output_data['w_val'] = []
    output_data['init_Q_N'] = []
    output_data['init_prop'] = []
    output_data['PCN'] = []
    output_data['distance'] = []
    output_data['sim_data'] = []

    for PCN in PCN_values:
        for w_val in w_values:
            for init_Q_N in init_Q_values:
                print("PCN: " + str(PCN))
                print("w_" + str(heterologous_species_name) + ": ", str(w_val))
                print("init_Q: " + str(init_Q_N))

                data_df = data_dict['processed_df']
                n_passages = np.shape(data_df['passage'])[0]

                t = np.arange(0, 1440, 1.0)

                # Initialise bioreactor object
                bioreactor = Bioreactor(bioreactor_config)
                
                # Initialise plasmid bearing strain
                plasmid_bearing_strain = Strain(bioreactor, plasmid_bearing_config)

                # Initialise wt strain
                wt_strain = Strain(bioreactor, wt_config)

                strain_list = [plasmid_bearing_strain, wt_strain]

                # Sample parameters and initial species from prior
                # Generate integration inputs and prepare objects
                species_keys, y0 = utils.generate_integrate_inputs(bioreactor, strain_list)

                # Set plasmid bearing strain w value
                plasmid_bearing_strain.params['w_' + heterologous_species_name] =  float(w_val * PCN)
                
                # Set Q strain initial N
                Q_N_idx = species_keys.index('Q_N')
                y0[Q_N_idx] = float(init_Q_N)

                # Set Q strain initial N
                Q_N_idx = species_keys.index('Q_N')
                y0[Q_N_idx] = float(init_Q_N)

                # Set P strain initial N 
                P_N_idx = species_keys.index('P_N')
                y0[P_N_idx] = 1000 - y0[Q_N_idx]

                initial_plasmid_bearing_prop = y0[P_N_idx] / (y0[P_N_idx] + y0[Q_N_idx])


                # Set heterologous gene transcription rate (w)
                # Update transcription rate, 'w_x' by multiplication with PCN
                sim_plasmid_bearing_props = run_plasmid_loss_passage_experiment(
                    y0=y0[:], species_keys=species_keys, bioreactor=bioreactor, 
                    plasmid_bearing_strain=plasmid_bearing_strain, wt_strain=wt_strain, 
                    t=t, n_passages=n_passages, PCN=PCN
                    )

                distance = 0
                for real_data, sim_data in (zip(data_df['mean_max_clust_prop'].values, sim_plasmid_bearing_props)):
                    distance += abs(real_data - sim_data)

                print("distance: ", str(distance))
                print("")

                # Add data to output dict
                output_data['plasmid_name'].append(data_dict['plasmid_name'])
                output_data['w_val'].append(w_val)
                output_data['init_Q_N'].append(y0[Q_N_idx])
                output_data['init_prop'].append(initial_plasmid_bearing_prop)
                output_data['PCN'].append(PCN)
                output_data['distance'].append(distance)
                output_data['sim_data'].append(sim_plasmid_bearing_props)

    # Make into df
    out_df = pd.DataFrame(output_data)
    out_path = output_dir + "/GS_plasmid_loss.csv"
    with open(out_path, 'a') as f:
        out_df.to_csv(f, mode='a', header=f.tell()==0)

    min_row = out_df.iloc[out_df['distance'].idxmin()]

    # plot lowest distance
    plot_simulation_fitting_data_overlay(data_dict, bioreactor_config, plasmid_bearing_config, 
        wt_config, heterologous_species_name, min_row['w_val'], min_row['init_Q_N'], min_row['PCN'])

    # plot lowest distance
    plot_simulation_bio_rpt_data_overlay(data_dict, min_row['sim_data'], heterologous_species_name, min_row['w_val'], min_row['init_Q_N'], min_row['PCN'])


def fit_plasmid_loss_params():
    plasmid_loss_dir = "./fit_plasmid_loss"
    plasmid_loss_data_dir = plasmid_loss_dir + "/data"
    output_dir = "./fit_plasmid_loss/output"

    try:
        os.mkdir(output_dir)
    except FileExistsError:
        pass

    # Load data and process for fitting
    data_dict = make_data_dict(plasmid_loss_dir + "/data")

    # Load config files
    bioreactor_config_yaml = plasmid_loss_data_dir + "/bioreactor_config.yaml"
    with open(bioreactor_config_yaml, 'r') as yaml_file:
        bioreactor_config = yaml.load(yaml_file, Loader=yaml.FullLoader)

    P_strain_config_yaml = plasmid_loss_data_dir + "/P_strain_config.yaml"
    with open(P_strain_config_yaml, 'r') as yaml_file:
        P_strain_config = yaml.load(yaml_file, Loader=yaml.FullLoader)

    Q_strain_config_yaml = plasmid_loss_data_dir + "/Q_strain_config.yaml"
    with open(Q_strain_config_yaml, 'r') as yaml_file:
        Q_strain_config = yaml.load(yaml_file, Loader=yaml.FullLoader)

    plasmid_bearing_config = P_strain_config
    wt_config = Q_strain_config

    # Setting values for w and init plasmid free to test
    w_values = [30, 35, 40, 45]
    init_Q_values = [1]

    # Perform grid search for each dataset
    for key in data_dict:
        print("Running grid search fit for ", key, " data...")
        print("")
        PCN_prior = data_dict[key]['prior_PCN']
        PCN_values = np.linspace(PCN_prior[0], PCN_prior[1], 10)
        PCN_values = [PCN_prior[0]]

        grid_search_plasmid_drop(
            data_dict=data_dict[key], bioreactor_config=bioreactor_config, 
            plasmid_bearing_config=plasmid_bearing_config, wt_config=wt_config, 
            heterologous_species_name='h', PCN_values=PCN_values, 
            w_values=w_values, init_Q_values=init_Q_values, 
            output_dir=output_dir
            )

    print("Finished")

if __name__ == "__main__":
    fit_plasmid_loss_params()

