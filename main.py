from __future__ import division  # normal division
from mock import *
from specie import *
from environment import *
from math_model import *
import time
from fipy.tools import parallel
import pdb
import numpy as np

t0 = time.time()

def rates(C0):
    dcdt = np.zeros(np.shape(np.array(C0)))
    dcdt[0] = -ox*om
    dcdt[1] = ox*om
    return dcdt

oxygen = Species('ox', 385, 'Fixed', 0.15)
organic_m = Species('om', 5, 'Flux', 5)
# organic_m2 = Species('om2', 5, 'Flux', 500)
# organic_m3 = Species('om3', 5, 'Flux', 500)
# organic_m4 = Species('om4', 5, 'Flux', 500)
# organic_m5 = Species('om5', 5, 'Flux', 500)
# organic_m6 = Species('om6', 5, 'Flux', 500)
# organic_m7 = Species('om7', 5, 'Flux', 500)
# organic_m8 = Species('om8', 5, 'Flux', 500)
environment = Environment(0, -5)
oxygen.reaction_term = '-ox*om'
oxygen.phreeqc_name = "O(-2)"

environment.add_species(oxygen)
environment.add_species(organic_m)
# environment.add_species(organic_m2)
# environment.add_species(organic_m3)
# environment.add_species(organic_m4)
# environment.add_species(organic_m5)
# environment.add_species(organic_m6)
# environment.add_species(organic_m7)
# environment.add_species(organic_m8)
model = MathModel(environment,nx=5)
model.run(Time=2)
# Tracer()()

# print locals()
print 'time:',parallel.procID,time.time() - t0

print rates(1,2)