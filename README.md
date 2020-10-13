# Whole cell model for plasmid loss simulation and prediction

## Requirements
Python 3.6+, see requirements.txt for python libraries used

## Configuration files
Configuration files are used to define the bioreactor environment and each strain (See example config files). 

Each parameter definition take two values for uniform sampling, `[lower_bound, upper_bound]` . <br /> Make `lower_bound == upper_bound` 
if the parameter is fixed

The prefix set in each config file is used to define the membership of paramters and species.

### Bioreactor
`bioreactor_config.yaml` defines the bioreactor environment. For batch cultures efflux and influx are zerro,<br />
`k_in = d_s = 0`

Chemostat simulations can be made by setting efflux and influx to constant values, <br /> `k_in = d_s > 0`

### Strains

`Q_strain_config.yaml` defines a wild-type strain with no heterologous expression. 
All parameters and species here present are essential to the model, their values can be changed


`P_strain_config.yaml` defines an engineered strain, expressing a heterologous protein, `h`. We declare the necessary species
to the priors, the prefixes define the species type. <br /> `g_h` (protein), <br /> `m_h` (mRNA), <br /> `c_h` (mRNA-ribosome complex).


We also declare the necessary parameters for the heterologous expression. <br /> `w_h` (transcription rate), <br /> `n_h` (transcript length) and <br /> `theta_h` (non-ribosome transcription threshold)


