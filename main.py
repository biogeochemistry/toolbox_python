from __future__ import division  # normal division
from mock import *
from specie import *
from environment import *
from math_model import *


oxygen = Specie('ox', 385, 'Fixed', 0.15)
organic_m = Specie('om', 5, 'Flux', -500)
environment = Environment(0.85, -5)
oxygen.reaction_term = '-ox*om'

environment.add_specie(oxygen)
environment.add_specie(organic_m)
model = MathModel(environment)
model.run()

