from __future__ import division  # normal division
import nose.tools as test
from nose.plugins.skip import SkipTest
from mock import *
from species import *
from environment import *
from math_model import *
import pdb


class TestSpescieAndEnvironment:

    def function_test(cls):
        oxygen = Species('ox', 385, 'Fixed', 0.15)
        organic_m = Species('om', 5, 'Flux', -500)
        environment = Environment(0.85, -5)
        environment.add_species(oxygen)
        environment.add_species(organic_m)
        model = MathModel(environment)

    def add_species_test(cls):
        environment = Environment(0.85, -5)
        oxygen = MagicMock()
        environment.add_species(oxygen)
        test.eq_(environment.species[oxygen.name], oxygen)

    def viewver_test(cls):
        # raise SkipTest
        organic_m = Species('om', 5, 'Flux', -500)
        oxygen = Species('ox', 400, 'Fixed', 0.15)
        environment = Environment(0.85, -5)
        oxygen.reaction_term = '-om*ox'
        environment.add_species(organic_m)
        environment.add_species(oxygen)
        model = MathModel(environment)


class TestMathModel:

    def initialization_with_default_values_test(cls):
        env = MagicMock()
        model = MathModel(env)
        # pdb.set_trace()
        assert model.dt == 0.1

    def solution_of_ODE_using_Butcher5_test(cls):
        def test_rates(x): return x
        assert butcher5(4, test_rates, 1000) - 4 * math.exp(1) <= 1e-10

    def solution_of_ODE_using_RK4_test(cls):
        def test_rates(x): return x
        assert rk4(4, test_rates, 1000) - 4 * math.exp(1) <= 1e-10
