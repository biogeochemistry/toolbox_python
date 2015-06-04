# from __future__ import print_function
from __future__ import division  # normal division
import numpy as np
from math_model import *
from species import *
from environment import *

MODE = 'dll'  # 'dll' or 'com'

if MODE == 'com':
    import phreeqpy.iphreeqc.phreeqc_com as phreeqc_mod
elif MODE == 'dll':
    import phreeqpy.iphreeqc.phreeqc_dll as phreeqc_mod
else:
    raise Exception('Mode "%s" is not defined use "com" or "dll".' % MODE)


def make_initial_conditions():
    """
    Specify initial conditions data blocks.

    Uniform initial conditions are assumed.
    """
    initial_conditions = """
    SOLUTION 0
        units            mmol/kgw
        temp             25.0
        Alkalinity 0.3 ## = [HCO3] + 2[CO3] + [OH] + 2[PO4] - [H+]
        C(+4)   0.1 ## = HCO3 + CO2 +CO3 + H2CO3 
        NO3-    0.1
        N(-3)   0.1 ## = NH3 + NH4
        N(+5)   0.1 ## = NO3
        S(-2)   0.1 ## = HS + H2S + FeS2
        S(+6)  0.1 ## = SO4
        Fe(+2)  0.1 ## = Fe2 + FeS2 + FeS
        Fe(+3)  0.1 ## = FeOH3 + FeOOH
        Ca  0.1 ## = Ca2 + Ca3PO4
        P 0.1 ## = PO4 + Ca3PO4 + PO4adsa + PO4adsb'
    END
        """
    return initial_conditions


def make_selected_output(components):
    """
    Build SELECTED_OUTPUT data block
    """
    headings = "-headings    cb    H    O    "
    for i in range(len(components)):
        headings += components[i] + "\t"
    selected_output = """
    SELECTED_OUTPUT
        -reset false
    USER_PUNCH
    """
    selected_output += headings + "\n"
    #
    # charge balance, H, and O
    #
    code = '10 w = TOT("water")\n'
    code += '20 PUNCH CHARGE_BALANCE, TOTMOLE("H"), TOTMOLE("O")\n'
    #
    # All other elements
    #
    lino = 30
    for component in components:
        code += '%d PUNCH w*TOT(\"%s\")\n' % (lino, component)
        lino += 10
    selected_output += code
    return selected_output


def get_selected_output(phreeqc):
    """Return calculation result as dict.

    Header entries are the keys and the columns
    are the values as lists of numbers.
    """
    output = phreeqc.get_selected_output_array()
    header = output[0]
    conc = {}
    for head in header:
        conc[head] = []
    for row in output[1:]:
        for col, head in enumerate(header):
            conc[head].append(row[col])
    return conc


def phreqc_run():
    phreeqc = phreeqc_mod.IPhreeqc()
    phreeqc.load_database(r"llnl.dat")
    initial_conditions = make_initial_conditions()
    phreeqc.run_string(initial_conditions)
    components = phreeqc.get_component_list()
    selected_output = make_selected_output(components)
    phreeqc.run_string(selected_output)
    phc_string = "RUN_CELLS; -cells 0\n"
    phc_string += 'USER_PUNCH\n'
    phc_string += '-headings C+4 N N+3 N+5 S-2 S+6 Fe+2 Fe+3 Fe Ca pH\n'
    phc_string += '1 PUNCH TOT("C(+4)")\n'
    phc_string += '2 PUNCH TOT("N")\n'
    phc_string += '3 PUNCH TOT("N(+3)")\n'
    phc_string += '4 PUNCH TOT("N(+5)")\n'
    phc_string += '5 PUNCH TOT("S(-2)")\n'
    phc_string += '6 PUNCH TOT("S(+6)")\n'
    phc_string += '7 PUNCH TOT("Fe(+2)")\n'
    phc_string += '8 PUNCH TOT("Fe(+3)")\n'
    phc_string += '9 PUNCH TOT("Fe")\n'
    phc_string += '10 PUNCH TOT("Ca")\n'
    phc_string += '11 PUNCH -LA("H+")\n -end'
    phreeqc.run_string(phc_string)
    conc = get_selected_output(phreeqc)
    # print(conc)
    initial = {}
    inflow = {}
    for name in conc:
        initial[name] = conc[name][0]

    for name in conc:
        value = [initial[name]] * len(conc[name])
        conc[name] = value

    all_names = conc.keys()
    names = [name for name in all_names if name not in ('cb', 'H', 'O')]

def create_function_from(string):
    pass


if __name__ == '__main__':
    # print rates(1)

    # print butcher5(4, test_rates, 100)
    # oxygen = Species('ox', 385, 'Fixed', 0.15, bc_xn_value=0)
    # organic_m = Species('om', 5, 'Flux', 5)

    # environment = Environment(0, -5)
    # oxygen.reaction_term = '-ox*om'
    # organic_m.reaction_term = 'ox*om'
    # environment.add_species(oxygen)
    # environment.add_species(organic_m)

    string = "-k*om*ox"
    k = 1
    om = 5
    ox = 1

    d = {}
    string = "def f(x): a=x[0] + x[1]; b;return a"
    exec string in d


    print d['f']([om,ox])

