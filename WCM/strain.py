import numpy as np
import oyaml as  yaml
from collections import OrderedDict
from . import functions as funcs


class Strain:
    def __init__(self, bioreactor, strain_config):
        """ Initialise strain

        @param bioreactor object in which the strain is growing
        @param strain_config dictionary containing priors and strain settings
        """
        self.prior_species_dict = strain_config['prior_initial_species']
        self.prior_params_dict = strain_config['prior_paramters']

        self.strain_prefix = strain_config['strain_prefix']
        self.bioreactor = bioreactor
        
        self.protein_species = []
        self.mrna_species = []
        self.c_mrna_species = []


        self.sample_initial_species()
        self.sample_parameters()
        self.categorise_species()


    def sample_initial_species(self):
        """ Generates a dictionary containing the initial species for this strain
        """
        self.initial_species = {}
        prior_dict = self.prior_species_dict

        for k in prior_dict.keys():
            key_with_prefix = self.strain_prefix + '_' + k

            if prior_dict[k][0] == prior_dict[k][1]:
                self.initial_species[key_with_prefix] = float(prior_dict[k][0])
            else:
                sampled_value = np.random.uniform(float(prior_dict[k][0]), float(prior_dict[k][1]))
                self.initial_species[key_with_prefix] = sampled_value

    def sample_parameters(self):
        """ Generates a dictionary containing parameters for this strain
        """
        self.params = {}
        prior_dict = self.prior_params_dict

        for k in prior_dict.keys():
            if prior_dict[k][0] == prior_dict[k][1]:
                self.params[k] = float(prior_dict[k][0])

            else:
                sampled_value = np.random.uniform(float(prior_dict[k][0]), float(prior_dict[k][1]))
                self.params[k] = float(sampled_value)

    def categorise_species(self):
        """ Adds species keys to the protein, mrna and c_mrna lists
        """
        for k in self.prior_species_dict.keys():
            species_tag = k.split('_')[0]

            if species_tag == 'g':
                self.protein_species.append(k)

            elif species_tag == 'm':
                self.mrna_species.append(k)

            elif species_tag == 'c':
                self.c_mrna_species.append(k)

            else:
                continue

   
    def get_initial_species(self):
        """ Returns a list of species keys and a list of species values
        """
        species_keys = []
        species_values = []

        for key, value in self.initial_species.items():
            species_keys.append(key)
            species_values.append(value)

        return species_keys, species_values


    def update_current_species_values(self, current_species_dict):
        """Creates a dictionary containing the current species values, 
        extracting species that are relavent to this particular strain.
        """

        self.species_values = {}

        for k in self.initial_species.keys():
            no_prefix_key = '_'.join(k.split('_')[1:])
            self.species_values[no_prefix_key] = current_species_dict[k]

    def calculate_mass(self):
        """ Mass is the sum of proteins (including bound ribosomes) Eq(12)
        """
        sum_protein_mass = 0

        # Sum of proteins
        for p in self.protein_species:
            protein_name = '_'.join(p.split('_')[1:])
            protein_length = self.params['n_' + protein_name]

            protein_mass = protein_length * self.species_values[p]
            sum_protein_mass += protein_mass

        # Sum of mrna bound ribosome
        for c in self.c_mrna_species:
            sum_protein_mass += self.species_values[c] * self.params['n_r']
        
        return sum_protein_mass
        

    def calculate_nu_imp(self):
        """ Rate of nutrient import Supp EQ (7)
        """
        e_t = self.species_values['g_e_t']
        s = self.bioreactor.species_values['s']
        v_t = self.params['v_t']
        K_t = self.params['K_t']

        nu_imp = funcs.nu_imp(e_t, v_t, s, K_t)

        return nu_imp

    def calculate_nu_cat(self):
        """ Rate of metabolic enzyme activity Supp EQ (7)
        """
        s_i = self.species_values['s_i']
        e_m = self.species_values['g_e_m']
        K_m = self.params['K_m']
        v_m = self.params['v_m']

        nu_cat = funcs.nu_cat(e_m, v_m, s_i, K_m)

        return nu_cat

    def calculate_gamma_a(self):
        """ Mass is the sum of proteins (including bound ribosomes) main text Eq(3)
        """

        gamma_max = self.params['gamma_max']
        a = self.species_values['a']
        K_gamma = self.params['K_gamma']

        return (gamma_max * a) / (K_gamma + a)


    def calculate_growth_rate(self, mass, gamma_a):
        """ Growth rate is a function of all translating ribosomes (c_x)
        """
        # Sum mrna bound ribosomes
        sum_c_mrna = 0
        for c in self.c_mrna_species:
            sum_c_mrna += self.species_values[c]


        lambd = (gamma_a / mass) * sum_c_mrna

        return lambd

    def dot_N(self, lambd):
        """ Differential for change in population
        """
        dot_N = (lambd * self.species_values['N']) - (self.params['d_N'] * self.species_values['N'])
        return dot_N

    def dot_s_i(self, nu_imp, nu_cat, lambd):
        """ Differential for intracellular nutrients Supp Eq(1)
        """
        s_i = self.species_values['s_i']
        a = self.species_values['a']

        dot_si = nu_imp - nu_cat - lambd * s_i

        return dot_si

    def dot_a(self, nu_cat, gamma_a, lambd):
        """ Differential for energy species Supp Eq (2)
        """
        n_s = self.params['n_s']
        a = self.species_values['a']

        sum_translation_costs = 0
        for c in self.c_mrna_species:
            c_x = self.species_values[c]
            n_x_param = 'n_' + '_'.join(c.split('_')[1:])
            n_x = self.params[n_x_param]

            sum_translation_costs += funcs.nu_x(c_x, n_x, gamma_a) * n_x

        dot_a = n_s * nu_cat - sum_translation_costs - lambd * a

        return dot_a

    def dot_g_r(self, gamma_a, lambd):
        """ Differential for ribosome protein species Supp Eq (3)
        """
        c_r = self.species_values['c_r']
        n_r = self.params['n_r']
        k_b = self.params['k_b']
        k_u = self.params['k_u']
        g_r = self.species_values['g_r']

        nu_r = funcs.nu_x(c_r, n_r, gamma_a)

        sum_c_mrna = 0

        for c in self.c_mrna_species:
            c_x = self.species_values[c]
            n_x_param = 'n_' + '_'.join(c.split('_')[1:])
            n_x = self.params[n_x_param]
            nu_x = funcs.nu_x(c_x, n_x, gamma_a)

            # mRNA species for x
            m_x = self.species_values['m_' + '_'.join(c.split('_')[1:])]

            sum_c_mrna += nu_x - (k_b * g_r * m_x) + (k_u * c_x)

        dot_r = nu_r - lambd * g_r + sum_c_mrna

        return dot_r

    def dot_g_x(self, g_x, c_x, n_x, gamma_a, lambd):
        """ Differential for translation of non ribosome protein species Supp Eq (4)
        """
        nu_x = funcs.nu_x(c_x, n_x, gamma_a)

        dot_g_x = nu_x - lambd * g_x

        return dot_g_x


    def dot_m_x(self, m_x, c_x, n_x, w_x, theta_x, gamma_a, lambd):
        """ Differential for transcription of an m_RNA species Supp Eq (5)
        """ 
        g_r = self.species_values['g_r']
        a = self.species_values['a']
        d_m = self.params['d_m']
        k_b = self.params['k_b']
        k_u = self.params['k_u']

        omega_x = funcs.omega_x(w_x, theta_x, a)
        nu_x = funcs.nu_x(c_x, n_x, gamma_a)

        dot_m_x = omega_x - (lambd + d_m) * m_x + nu_x - k_b * g_r * m_x + k_u * c_x

        return dot_m_x


    def dot_c_x(self, m_x, c_x, n_x, gamma_a, lambd):
        """ Differential for a mRNA transcript undergoing translation Supp Eq (6)

        @param m_x mRNA species value
        @param c_x mRNA-ribosome complex species value
        @param n_x protein length
        @param gamma_a 
        @param lambd
        
        @return dot_c_x is the derivative for a particular mRNA-ribosome complex species
        """
        g_r = self.species_values['g_r']
        k_b = self.params['k_b']
        k_u = self.params['k_u']
        nu_x = funcs.nu_x(c_x, n_x, gamma_a)

        dot_c_x = -(lambd * c_x) + k_b * g_r * m_x - k_u * c_x - nu_x

        return dot_c_x


    def calculate_all_g_x_diff(self, gamma_a, lambd):
        """ Iterates through all protein species of the strain and 
        calculates the differentials

        @param gamma_a
        @param lambd

        @return dot_g_x_dict containing differentials for each protein species
        """ 

        dot_g_x_dict = {}

        for g_x_str in self.protein_species:
            
            # Skip ribosome protein            
            if g_x_str == 'g_r':
                continue

            c_x_str = 'c_' + '_'.join(g_x_str.split('_')[1:])
            n_x_str = 'n_' + '_'.join(g_x_str.split('_')[1:])

            g_x = self.species_values[g_x_str]
            c_x = self.species_values[c_x_str]
            n_x = self.params[n_x_str]

            output_name = self.strain_prefix + '_' + g_x_str

            dot_g_x_dict[output_name] = self.dot_g_x(g_x, c_x, n_x, gamma_a, lambd)

        return dot_g_x_dict

    def calculate_all_m_x_diff(self, gamma_a, lambd):
        """ Iterates through all mRNA species of the strain and 
        calculates the differentials

        @param gamma_a
        @param lambd

        @return dot_m_x_dict containing differentials for each mRNA species
        """ 
        dot_m_x_dict = {}

        for m_x_str in self.mrna_species:
            c_x_str = 'c_' + '_'.join(m_x_str.split('_')[1:])
            n_x_str = 'n_' + '_'.join(m_x_str.split('_')[1:])
            w_x_str = 'w_' + '_'.join(m_x_str.split('_')[1:])
            theta_x_str = 'theta_' + '_'.join(m_x_str.split('_')[1:])

            m_x = self.species_values[m_x_str]
            c_x = self.species_values[c_x_str]

            n_x = self.params[n_x_str]
            w_x = self.params[w_x_str]
            theta_x = self.params[theta_x_str]

            dot_m_x = self.dot_m_x(m_x, c_x, n_x, w_x, theta_x, gamma_a, lambd)

            output_name = self.strain_prefix + '_' + m_x_str

            dot_m_x_dict[output_name] = dot_m_x

        return dot_m_x_dict

    def calculate_all_c_x_diff(self, gamma_a, lambd):
        """ Iterates through all mRNA-ribosome complex species of the strain and 
        calculates the differentials

        @param gamma_a
        @param lambd

        @return dot_c_x_dict containing differentials for each mRNA species
        """ 

        dot_c_x_dict = {}

        k_b = self.params['k_b']
        k_u = self.params['k_u']

        for c_x_str in self.c_mrna_species:
            m_x_str = 'm_' + '_'.join(c_x_str.split('_')[1:])
            n_x_str = 'n_' + '_'.join(c_x_str.split('_')[1:])

            c_x = self.species_values[c_x_str]
            m_x = self.species_values[m_x_str]
            n_x = self.params[n_x_str]

            dot_c_x = self.dot_c_x(m_x, c_x, n_x, gamma_a, lambd)
            
            output_name = self.strain_prefix + '_' + c_x_str

            dot_c_x_dict[output_name] = dot_c_x

        return dot_c_x_dict


    def calculate_differentials(self):
        """ Calculates differentials for a time step using current species values
        """

        output_dict = {}

        # Multuple use functions
        nu_imp = self.calculate_nu_imp()
        nu_cat = self.calculate_nu_cat()
        mass = self.calculate_mass()
        gamma_a = self.calculate_gamma_a()
        lambd = self.calculate_growth_rate(mass, gamma_a)

        dict_key_prefix = self.strain_prefix + '_'

        # Differential equations
        output_dict[dict_key_prefix + 'N'] = self.dot_N(lambd)
        output_dict[dict_key_prefix + 's_i'] = self.dot_s_i(nu_imp, nu_cat, lambd)
        output_dict[dict_key_prefix + 'a'] = self.dot_a(nu_cat, gamma_a, lambd)
        output_dict[dict_key_prefix + 'g_r'] = self.dot_g_r(gamma_a, lambd)

        # Diff eqs for protein, mrna dn c_mrna lists
        dot_g_x_dict = self.calculate_all_g_x_diff(gamma_a, lambd)
        dot_m_x_dict = self.calculate_all_m_x_diff(gamma_a, lambd)
        dot_c_x_dict = self.calculate_all_c_x_diff(gamma_a, lambd)

        output_dict.update(dot_g_x_dict)
        output_dict.update(dot_m_x_dict)
        output_dict.update(dot_c_x_dict)

        return output_dict
