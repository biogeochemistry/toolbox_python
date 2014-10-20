from __future__ import division  # normal division
import nose.tools as test
from nose.plugins.skip import SkipTest
from mock import *
from specie import *
from environment import *
from math_model import *


class TestSpescieAndEnvironment:

    def function_test(cls):
        oxygen = Specie('ox', 385, 'Fixed', 0.15)
        organic_m = Specie('om', 5, 'Flux', -500)
        environment = Environment(0.85, -5)
        environment.add_specie(oxygen)
        environment.add_specie(organic_m)
        model = MathModel(environment)

    def add_specie_test(cls):
        environment = Environment(0.85, -5)
        oxygen = MagicMock()
        environment.add_specie(oxygen)
        test.eq_(environment.species[oxygen.name], oxygen)

    def viewver_test(cls):
        raise SkipTest
        organic_m = Specie('om', 5, 'Flux', -500)
        oxygen = Specie('ox', 400, 'Fixed', 0.15)
        environment = Environment(0.85, -5)
        oxygen.reaction_term = '-om*ox'
        environment.add_specie(organic_m)
        environment.add_specie(oxygen)
        model = MathModel(environment)


class TestMathModel:

    def initialization_with_default_values_test(cls):
        env = MagicMock()
        model = MathModel(env)
        assert model.dt == 0.001
