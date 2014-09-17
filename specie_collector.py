from __future__ import division  # normal division
from second_order_ode import *
from boundary_conditions import *
from bvp_ode import *
from bvp_pde import *
from coupled_pde import *
from specie_collector import *
import numexpr as ne


class SpecieCollector(object):

    """docstring for SpecieCollector"""

    def __init__(self):
        self.all = {}
        self.dict_of_conc_and_params = {}

    def add_specie(self, specie, D, w, dt, T, bc_x0_type, bc_x0_value, bc_xn_type, bc_xn_value, init_concentrations, x_min, x_max, num_x_nodes):
        self.all[specie] = self.create_single_container(D, w, dt, T, bc_x0_type, bc_x0_value, bc_xn_type, bc_xn_value, init_concentrations, x_min, x_max, num_x_nodes)
        self.add_to_dict_of_conc(specie)

    def create_single_container(self, D, w, dt, T, bc_x0_type, bc_x0_value, bc_xn_type, bc_xn_value, init_concentrations, x_min, x_max, num_x_nodes):
        """Create dict for specie"""
        container = {'D': D, 'w': w, 'dt': dt, 'T': T, 'x_min': x_min, 'x_max': x_max,
                     'num_x_nodes': num_x_nodes, 'bc_x0_type': bc_x0_type, 'bc_x0_value': bc_x0_value,
                     'bc_xn_type': bc_xn_type, 'bc_xn_value': bc_xn_value,
                     'init_concentrations': init_concentrations}
        container['pde'] = self.create_pde(container)
        return container

    def create_pde(self, specie):
        """Create pde from specie container(dict)
        specie dict {'D': D, 'w':w, 'dt':dt, 'T':T, 'bc_x0_type': bc_x0_type, 'bc_x0_value':bc_x0_value, 'bc_xn_type':bc_xn_type, 'bc_xn_value': bc_xn_value, 'init_concentrations': init_concentrations}"""

        def model_prob_1(x):
            return 0

        ode = SecondOrderOde1D(specie['D'], specie['w'], 0, model_prob_1, 0, specie['x_min'], specie['x_max'])
        bc = BoundaryConditions1D()
        if specie['bc_x0_type'] == 'Dirichlet':
            bc.set_x0_dirichlet_bc(specie['bc_x0_value'])
        if specie['bc_x0_type'] == 'Neumann':
            bc.set_x0_neumann_bc(specie['bc_x0_value'])
        if specie['bc_xn_type'] == 'Dirichlet':
            bc.set_xn_dirichlet_bc(specie['bc_xn_value'])
        if specie['bc_xn_type'] == 'Neumann':
            bc.set_xn_neumann_bc(specie['bc_xn_value'])
        return BvpPde1D(ode, bc, specie['dt'], 0, specie['T'], specie['num_x_nodes'], specie['init_concentrations'])

    def differentiate_all_1_TS(self):
        self.differentiate_transport_terms_all()
        self.differentiate_reaction_terms_all()

    def differentiate_transport_terms_all(self):
        for specie in self.all:
            self.all[specie]['pde'].differentiate_pde_1TS()
            self.update_dict_of_conc()

    def get_C_vector(self, specie):
        if self.all[specie]['pde'].Ut.shape == (self.all[specie]['pde'].num_x_nodes,):
            return self.all[specie]['pde'].Ut
        else:
            return self.all[specie]['pde'].Ut[-1]

    def add_to_dict_of_conc(self, specie):
        self.dict_of_conc_and_params[specie] = self.get_C_vector(specie)

    def update_dict_of_conc(self):
        for specie in self.dict_of_conc_and_params:
            self.dict_of_conc_and_params[specie] = self.get_C_vector(specie)

    def differentiate_current_reaction_term(self, specie):
        R_dt = self.all[specie]['dt'] * ne.evaluate(self.all[specie]['rate'], local_dict=self.dict_of_conc_and_params)
        self.all[specie]['pde'].Ut += R_dt

    def differentiate_reaction_terms_all(self):
        for specie in self.all:
            self.differentiate_current_reaction_term(specie)

    def create_variables_from_regexp(self):
        # re.sub(r'([A-z][A-z0-9]*)', r'%(\1)s', rate)
        pass
