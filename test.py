from __future__ import division  # normal division
import nose.tools as test
from nose.plugins.skip import SkipTest
from mock import *
from second_order_ode import *
from boundary_conditions import *
from bvp_ode import *
from bvp_pde import *
from coupled_pde import *
from specie_collector import *
from grid import *
import math


def model_prob_1(x):
    return 1.0


def model_prob_2(x):
    return 34.0 * math.sin(x)


def solution_prob_1(x):
    return 0.5 * x * (1 - x)


def solution_prob_2(x):
    return (4 * math.exp(x) + math.exp(-4 * x)) / (4 * math.exp(math.pi) + math.exp(-4 * math.pi)) - 5 * math.sin(x) - 3 * math.cos(x)


def init_cond(x):
    return 0.5 * x


class TestBoundaryConditions:

    def initilization_1D_test(cls):
        bc0 = BoundaryConditions1D()
        test.eq_(bc0.x0_is_dirichlet, False)
        test.eq_(bc0.xn_is_dirichlet, False)
        test.eq_(bc0.x0_is_neumann, False)
        test.eq_(bc0.xn_is_neumann, False)

    def dirichlet_1D_test(cls):
        bc1 = BoundaryConditions1D()
        bc1.set_x0_dirichlet_bc(0.15)
        bc1.set_xn_dirichlet_bc(0)
        test.eq_(bc1.x0_value, 0.15)
        test.eq_(bc1.xn_value, 0)
        test.eq_(bc1.x0_is_dirichlet, True)
        test.eq_(bc1.xn_is_dirichlet, True)
        test.eq_(bc1.x0_is_neumann, False)
        test.eq_(bc1.xn_is_neumann, False)

    def neumann_1D_test(cls):
        bc2 = BoundaryConditions1D()
        bc2.set_x0_neumann_bc(0.15)
        bc2.set_xn_neumann_bc(0)
        test.eq_(bc2.x0_value, 0.15)
        test.eq_(bc2.xn_value, 0)
        test.eq_(bc2.x0_is_neumann, True)
        test.eq_(bc2.xn_is_neumann, True)
        test.eq_(bc2.x0_is_dirichlet, False)
        test.eq_(bc2.xn_is_dirichlet, False)

    def initilization_2D_test(cls):
        bc0 = BoundaryConditions2D()
        test.eq_(bc0.y0_is_dirichlet, False)
        test.eq_(bc0.yn_is_dirichlet, False)
        test.eq_(bc0.y0_is_neumann, False)
        test.eq_(bc0.yn_is_neumann, False)

    def dirichlet_2D_test(cls):
        bc1 = BoundaryConditions2D()
        bc1.set_x0_dirichlet_bc(0.15)
        bc1.set_xn_dirichlet_bc(0)
        bc1.set_y0_dirichlet_bc(0.15)
        bc1.set_yn_dirichlet_bc(0)
        test.eq_(bc1.x0_value, 0.15)
        test.eq_(bc1.xn_value, 0)
        test.eq_(bc1.x0_is_dirichlet, True)
        test.eq_(bc1.xn_is_dirichlet, True)
        test.eq_(bc1.x0_is_neumann, False)
        test.eq_(bc1.xn_is_neumann, False)
        test.eq_(bc1.y0_is_dirichlet, True)
        test.eq_(bc1.yn_is_dirichlet, True)
        test.eq_(bc1.y0_is_neumann, False)
        test.eq_(bc1.yn_is_neumann, False)

    def neumann_2D_test(cls):
        bc1 = BoundaryConditions2D()
        bc1.set_x0_neumann_bc(0.15)
        bc1.set_xn_neumann_bc(0)
        bc1.set_y0_neumann_bc(0.15)
        bc1.set_yn_neumann_bc(0)
        test.eq_(bc1.x0_value, 0.15)
        test.eq_(bc1.xn_value, 0)
        test.eq_(bc1.x0_is_dirichlet, False)
        test.eq_(bc1.xn_is_dirichlet, False)
        test.eq_(bc1.x0_is_neumann, True)
        test.eq_(bc1.xn_is_neumann, True)
        test.eq_(bc1.y0_is_dirichlet, False)
        test.eq_(bc1.yn_is_dirichlet, False)
        test.eq_(bc1.y0_is_neumann, True)
        test.eq_(bc1.yn_is_neumann, True)


class TestODE:

    def initialization_test(cls):
        ode1 = SecondOrderOde1D(10, -10, 0, model_prob_1, 0, 0, 15)
        test.eq_(ode1.uxx, 10)
        test.eq_(ode1.ux, -10)
        test.eq_(ode1.u, 0)
        test.eq_(ode1.x_min, 0)
        test.eq_(ode1.x_max, 15)


class TestGrid:

    def initialization_of_1D_test(cls):
        ode1 = SecondOrderOde1D(10, -10, 0, model_prob_1, 0, 0, 15)
        bc1 = BoundaryConditions1D()
        bc1.set_x0_dirichlet_bc(0.15)
        bc1.set_xn_neumann_bc(0)
        bvp1 = BvpPde1D(ode1, bc1, 0.01, 0, 10, 120, init_cond)
        test.eq_((bvp1.grid_x == np.linspace(0, 15, 120)).all(), True)


class TestAccuracyOfSolution:

    def solution_variance_of_problem_1_test(cls):
        n = 128
        ode1 = SecondOrderOde1D(-1, 0.0, 0.0, model_prob_1, 0.0, 0, 1)
        bc1 = BoundaryConditions1D()
        bc1.set_x0_dirichlet_bc(0)
        bc1.set_xn_dirichlet_bc(0)
        bvp1 = BvpOde1D(ode1, bc1, n)
        bvp1.solve()
        sum = 0
        for i in xrange(0, 128):
            x = bvp1.grid_x[i]
            sum += math.pow((bvp1.U[i] - solution_prob_1(x)), 2)
        sum /= n - 1
        test.eq_((math.sqrt(sum) < math.pow(10, -10)), True)

    def solution_variance_of_problem_2_test(cls):
        n = 128
        ode1 = SecondOrderOde1D(1, 3.0, -4.0, model_prob_2, 0.0, 0, math.pi)
        bc1 = BoundaryConditions1D()
        bc1.set_x0_neumann_bc(-5.0)
        bc1.set_xn_dirichlet_bc(4.0)
        bvp1 = BvpOde1D(ode1, bc1, n)
        bvp1.solve()
        sum = 0
        for i in xrange(0, n):
            x = bvp1.grid_x[i]
            sum += math.pow((bvp1.U[i] - solution_prob_2(x)), 2)
        sum /= (n - 1)
        test.eq_((math.sqrt(sum) < math.pow(10, -2)), True)

    def with_custom_grid_test(cls):
        n = 128
        grid = np.hstack((np.linspace(0, 1, 64, endpoint=False), np.linspace(1, 2, 48, endpoint=False), np.linspace(2, math.pi, 16, endpoint=True)))
        ode1 = SecondOrderOde1D(1, 3.0, -4.0, model_prob_2, 0.0, 0, math.pi)
        bc1 = BoundaryConditions1D()
        bc1.set_x0_neumann_bc(-5.0)
        bc1.set_xn_dirichlet_bc(4.0)
        bvp1 = BvpOde1D(ode1, bc1, n, True, grid)
        bvp1.solve()
        sum = 0
        for i in xrange(0, n):
            x = bvp1.grid_x[i]
            sum += math.pow((bvp1.U[i] - solution_prob_2(x)), 2)
        sum /= (n - 1)
        test.eq_((math.sqrt(sum) < math.pow(10, -2)), True)


class TestSpecieCollector:

    def add_method_creates_container_with_specie_test(cls):
        species = SpecieCollector()
        species.create_pde = MagicMock()
        species.add_specie('ox', 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
        test.eq_(species.all['ox']['D'] == 1, True)
        test.eq_(species.all['ox']['w'] == 2, True)
        test.eq_(species.all['ox']['dt'] == 3, True)
        test.eq_(species.all['ox']['T'] == 4, True)
        test.eq_(species.all['ox']['bc_x0_type'] == 5, True)
        test.eq_(species.all['ox']['bc_x0_value'] == 6, True)
        test.eq_(species.all['ox']['bc_xn_type'] == 7, True)
        test.eq_(species.all['ox']['bc_xn_value'] == 8, True)
        test.eq_(species.all['ox']['init_concentrations'] == 9, True)
        test.eq_(species.all['ox']['x_min'] == 10, True)
        test.eq_(species.all['ox']['x_max'] == 11, True)
        test.eq_(species.all['ox']['num_x_nodes'] == 12, True)
        species.create_pde.assert_called_once_with(species.all['ox'])

    def creating_PDEs_automatically_test(cls):
        x_min = 0
        x_max = 20
        num_x_nodes = 12
        x = np.linspace(x_min, x_max, num=num_x_nodes)
        a = 'ox'
        b = 'om'
        d_ox = 500.0
        w_ox = -5
        dt = 0.01
        T = 1
        bc_x0_type = 'Dirichlet'
        bc_x0_value = 0.15
        bc_xn_type = 'Neumann'
        bc_xn_value = 0
        init_concentrations = init_cond
        container = SpecieCollector()
        container.add_specie('ox', d_ox, w_ox, dt, T, bc_x0_type, bc_x0_value, bc_xn_type, bc_xn_value, init_concentrations, x_min, x_max, num_x_nodes)
        assert isinstance(container.all['ox']['pde'], BvpPde1D)

        # Testing of updating of dictionary of species
        for x in xrange(1, 10):
            container.differentiate_transport_terms_all()
        test.eq_((container.get_C_vector('ox') == container.dict_of_conc_and_params['ox']).all(), True)

    def differentitate_1_time_step_for_each_specie_in_the_container_test(cls):
        pde_stub = MagicMock()
        species = SpecieCollector()
        species.dict_of_conc_and_params = MagicMock()
        species.all[1] = {'pde': pde_stub}
        species.all[2] = {'pde': pde_stub}
        species.all[3] = {'pde': pde_stub}
        species.differentiate_transport_terms_all()
        pde_stub.differentiate_pde_1TS.assert_called_with()

    def get_C_vector_test(cls):
        species = SpecieCollector()
        pde_stub = MagicMock()
        species.all['ox'] = {'pde': pde_stub}
        species.all['ox']['pde'].Ut = np.array([0.5, 0.5])
        test.eq_((species.get_C_vector('ox') == np.array([0.5, 0.5])).all(), True)

    def update_dictionary_of_concentrations_test(cls):
        species = SpecieCollector()
        pde_stub = MagicMock()
        species.all['ox'] = {'pde': pde_stub}
        species.all['ox']['pde'].Ut = np.array([0.5, 0.5])
        species.dict_of_conc_and_params['ox'] = 0
        species.update_dict_of_conc()
        test.eq_((species.dict_of_conc_and_params['ox'] == np.array([0.5, 0.5])).all(), True)

    def add_current_specie_to_dictionary_of_conctreations_test(cls):
        species = SpecieCollector()
        species.create_pde = MagicMock()
        species.add_to_dict_of_conc = MagicMock()
        species.add_specie('ox', 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
        species.add_to_dict_of_conc.assert_called_once_with('ox')

    def add_to_dictionary_of_concentration_test(cls):
        species = SpecieCollector()
        pde_stub = MagicMock()
        species.all['ox'] = {'pde': pde_stub}
        species.all['ox']['pde'].Ut = np.array([0.01, 0.01])
        species.add_to_dict_of_conc('ox')
        test.eq_((species.all['ox']['pde'].Ut == species.dict_of_conc_and_params['ox']).all(), True)

    def differentiate_reaction_terms_test(cls):
        species = SpecieCollector()
        species.all['ox'] = {}
        pde_stub = MagicMock()
        species.all['ox'] = {'pde': pde_stub}
        species.all['ox']['dt'] = 0.1
        species.all['ox']['pde'].Ut = np.array([0.1, 0.2, 0.1])
        species.all['ox']['rate'] = 'k*ox'
        species.dict_of_conc_and_params['ox'] = np.array([0.1, 0.2, 0.1])
        species.dict_of_conc_and_params['k'] = 10
        species.differentiate_reaction_terms_all()
        test.eq_((species.all['ox']['pde'].Ut == np.array([0.1, 0.2, 0.1]) * 10 * 0.1 + np.array([0.1, 0.2, 0.1])).all(), True)

    def differentiate_all_1_TS_test(cls):
        species = SpecieCollector()
        species.differentiate_reaction_terms_all = MagicMock()
        species.differentiate_transport_terms_all = MagicMock()
        species.differentiate_all_1_TS()
        species.differentiate_transport_terms_all.assert_called_once_with()
        species.differentiate_transport_terms_all.assert_called_once_with()
