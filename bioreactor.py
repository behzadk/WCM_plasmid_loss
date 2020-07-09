import numpy as np

class Bioreactor:
    def __init__(self, bioreactor_config):
        self.prior_species_dict = bioreactor_config['prior_initial_species']
        self.prefix = bioreactor_config['bioreactor_prefix']
        self.prior_species_dict = bioreactor_config['prior_initial_species']
        self.prior_params_dict = bioreactor_config['prior_paramters']

    def sample_initial_species(self):
        self.initial_species = {}
        prior_dict = self.prior_species_dict

        for k in prior_dict.keys():
            sampled_value = np.random.uniform(np.float(prior_dict[k][0]), np.float(prior_dict[k][1]))

            key_w_prefix = self.prefix + '_' + k
            self.initial_species[key_w_prefix] = sampled_value

    def sample_parameters(self):
        """ Generates a dictionary containing parameters for the bioreactor
        """
        self.params = {}
        prior_dict = self.prior_params_dict

        for k in prior_dict.keys():
            if prior_dict[k][0] == prior_dict[k][1]:
                self.params[k] = float(prior_dict[k][0])

            else:
                sampled_value = np.random.uniform(float(prior_dict[k][0]), float(prior_dict[k][1]))
                self.params[k] = sampled_value


    """Returns a list of species keys and a list of species values
    """
    def get_initial_species(self):
        species_keys = []
        species_values = []

        for key, value in self.initial_species.items():
            species_keys.append(key)
            species_values.append(value)


        return species_keys, species_values


    def update_current_species_values(self, current_species_dict):
        self.species_values = {}

        for k in self.initial_species.keys():
            no_prefix_key = '_'.join(k.split('_')[1:])
            self.species_values[no_prefix_key] = current_species_dict[k]


    def calculate_differentials(self, strain_list):
        k_in = self.params['k_in']
        
        sum_nutrient_import = 0.0

        for strain in strain_list:
            nu_imp = strain.calculate_nu_imp()
            N = strain.species_values['N']

            sum_nutrient_import += nu_imp * N


        dot_s = k_in - sum_nutrient_import - self.params['d_s'] * self.species_values['s']

        output_dict = {}
        dict_key_prefix = self.prefix + '_'

        output_dict[dict_key_prefix + 's'] = dot_s


        return output_dict