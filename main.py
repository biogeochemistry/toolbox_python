from __future__ import division  # normal division
from mock import *
from specie import *
from environment import *
from math_model import *
import time
# from fipy.tools import parallel

t0 = time.time()

oxygen = Specie('ox', 385, 'Fixed', 0.15)
organic_m = Specie('om', 5, 'Flux', 500)
organic_m2 = Specie('om2', 5, 'Flux', 500)
organic_m3 = Specie('om3', 5, 'Flux', 500)
organic_m4 = Specie('om4', 5, 'Flux', 500)
organic_m5 = Specie('om5', 5, 'Flux', 500)
organic_m6 = Specie('om6', 5, 'Flux', 500)
organic_m7 = Specie('om7', 5, 'Flux', 500)
organic_m8 = Specie('om8', 5, 'Flux', 500)
environment = Environment(0.85, -5)
oxygen.reaction_term = '-ox*om'

environment.add_specie(oxygen)
environment.add_specie(organic_m)
# environment.add_specie(organic_m2)
# environment.add_specie(organic_m3)
# environment.add_specie(organic_m4)
# environment.add_specie(organic_m5)
# environment.add_specie(organic_m6)
# environment.add_specie(organic_m7)
# environment.add_specie(organic_m8)
model = MathModel(environment)
# model.run()
print DefaultSolver
print 'time:',parallel.procID,time.time() - t0

