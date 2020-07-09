# WCM plasmid loss
Whole cell model with plasmid loss experiments

## Requirements
Python 3.6+, see requirements.txt for python libraries used

## Config files
Configuration files are used to define the bioreactor environment and each strain (See example config files). 

Each parameter definition take two values for uniform sampling, `[lower_bound, upper_bound]` . Make `lower_bound == upper_bound` 
if the parameter is fixed

The prefix defined by in each config file is used to define the membership of paramters and species. For instance, if `prefix: P`, the species name in the 


`bioreactor_config.yaml` defines the bioreactor environment. For batch cultures efflux and influx are zerro, `k_in = d_s = 0`. Chemostat simulations can be made by setting efflux and influx to constant values, `k_in = d_s > 0`.


`Q_strain_config.yaml` defines a wild-type strain with no heterologous expression. 
All parameters and species here present are essential to the model, their values can be changed


`P_strain_config.yaml` defines an engineered strain, expressing a heterologous protein, `h`. We declare the necessary species
to the priors, the prefixes define the species type. `g_h` (protein), `m_h` (mRNA), `c_h` (mRNA-ribosome complex).
We also declare the necessary parameters for the heterologous expression. `w_h` (transcription rate), `n_h` (transcript length) and `theta_h` (non-ribosome transcription threshold)


