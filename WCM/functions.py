""" Function on rate of nutrient import Supp Eq (7)
"""
def nu_imp(e_t, v_t, s, K_t):
	return e_t * (v_t * s) / (K_t + s)

""" Function for rate of catabolism Supp Eq (7)
"""
def nu_cat(e_m, v_m, s_i, K_m):
	return e_m * (v_m * s_i) / (K_m + s_i)

""" Function for translation rate of protien x Supp Eq (8)
"""
def nu_x(c_x, n_x, gamma_a):
	return c_x * (gamma_a / n_x)

""" Function for rate of transcription of a gene Supp Eq(10)
"""
def omega_x(w_x, theta_x, a):
	return w_x * (a) / (theta_x +a)

""" Function for rate of transcription of housekeeping genes Supp Eq (11)
"""
def omega_q(q, w_q, theta_x, a):
	I_q = 1 / (1 + (q / K_q)**h_q )
	return I_q * w_x * a / (theta_x +a)

