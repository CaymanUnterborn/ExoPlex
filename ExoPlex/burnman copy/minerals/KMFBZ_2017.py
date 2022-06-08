# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
from Kurnosov et al.  Nature 2017
Kurnosov, A., et al. "Evidence for a Fe3+-rich pyrolitic lower mantle from (Al, Fe)-bearing bridgmanite elasticity data." 
Nature 543.7646 (2017): 543-546. doi:10.1038/nature21390
"""

from __future__ import absolute_import

from ..mineral import Mineral
from ..solidsolution import SolidSolution
from ..solutionmodel import *
from ..processchemistry import read_masses, dictionarize_formula, formula_mass

atomic_masses = read_masses()


"""
SOLID SOLUTIONS
"""


class bridgmanite(SolidSolution):

    def __init__(self, molar_fractions=None):
        self.name = 'bridgmanite/perovskite'
        self.endmembers = [[mg_perovskite(), '[Mg][Si]O3'],
            [fe_perovskite(), '[Fe][Si]O3'],
            [al_al_perovskite(), '[Al][Al]SiO3'],
            [fe_al_o3(), '[Fe][Al]SiO3'],]
        self.type = 'ideal'
        SolidSolution.__init__(self, molar_fractions)

class ferropericlase(SolidSolution):

    def __init__(self, molar_fractions=None):
        self.name = 'magnesiowustite/ferropericlase'
        self.type = 'symmetric'
        self.endmembers = [[periclase(), '[Mg]O'], [wuestite(), '[Fe]O']]
        self.energy_interaction = [[13.e3]]

        SolidSolution.__init__(self, molar_fractions)

"""
ENDMEMBERS
"""



class mg_perovskite (Mineral):
    
    def __init__(self):
        formula = 'MgSiO3'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'Mg_Perovskite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1368000.0,
            'V_0': 2.4445e-05,
            'K_0': 2.51e+11,
            'Kprime_0': 4.1,
            'Debye_0': 905.,
            'grueneisen_0': 1.57,
            'q_0': 1.1,
            'G_0': 1.73e+11,
            'Gprime_0': 1.7,
            'eta_s_0': 2.3,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

        Mineral.__init__(self)

class fe_perovskite (Mineral):

    def __init__(self):
        formula = 'FeSiO3'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'Fe_Perovskite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1043000.0,
            'V_0': 2.534e-05,
            'K_0': 2.72e+11,
            'Kprime_0': 4.1,
            'Debye_0': 871.,
            'grueneisen_0': 1.57,
            'q_0': 1.1,
            'G_0': 1.33+11,
            'Gprime_0': 1.4,
            'eta_s_0': 2.3,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

        Mineral.__init__(self)


class fe_al_o3 (Mineral):
    
    def __init__(self):
        formula = 'FeAlO3'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'FeAlO3',
            'formula': formula,
            'equation_of_state': 'slb3',
            'V_0': 2.69e-05,
            'K_0': 2.20e+11,
            'Kprime_0': 1.3,
            'Debye_0': 886.,
            'grueneisen_0': 1.57,
            'q_0': 1.1,
            'G_0': 96.0e+09,
            'Gprime_0': 3.4,
            'eta_s_0': 2.5,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}
    
        Mineral.__init__(self)


class al_al_perovskite(Mineral):
    
    def __init__(self):
        formula = 'AlAlO3'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'AlAl_perovskite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1533878.0,
            'V_0': 2.494e-05,
            'K_0': 2.58e+11,
            'Kprime_0': 4.1,
            'Debye_0': 886.,
            'grueneisen_0': 1.57,
            'q_0': 1.1,
            'G_0': 1.71e+11,
            'Gprime_0': 1.5,
            'eta_s_0': 2.5,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}
  
        Mineral.__init__(self)


class periclase (Mineral):

    def __init__(self):
        formula = 'MgO'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'Periclase',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -569000.,
            'V_0': 1.124e-05,
            'K_0': 1.602e+11,
            'Kprime_0': 3.99,
            'Debye_0': 767.0,
            'grueneisen_0': 1.36,
            'q_0': 1.7,
            'G_0': 1.31e+11,
            'Gprime_0': 2.1,
            'eta_s_0': 2.8,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

        Mineral.__init__(self)


class wuestite (Mineral):

    def __init__(self):
        formula = 'FeO'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'Wuestite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -242000.0,
            'V_0': 1.226e-05,
            'K_0': 1.49e+11,
            'Kprime_0': 3.6,
            'Debye_0': 454.,
            'grueneisen_0': 1.53,
            'q_0': 1.7,
            'G_0': 60.0e9,
            'Gprime_0': 1.8,
            'eta_s_0': 0.6,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

        Mineral.__init__(self)



class ca_perovskite (Mineral):
    
    def __init__(self):
        formula = 'CaSiO3'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'Ca_Perovskite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1463358.0,
            'V_0': 2.754e-05,
            'K_0': 2.36e+11,
            'Kprime_0': 3.9,
            'Debye_0': 802.,
            'grueneisen_0': 1.89,
            'q_0': 0.9,
            'G_0': 1.57e+11,
            'Gprime_0': 2.2,
            'eta_s_0': 1.3,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}
        

        Mineral.__init__(self)



class fe (Mineral):
    
    def __init__(self):
        formula = 'Fe'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'Fe',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1463358.0,
            'V_0': 0.673e-05,
            'K_0': 1.64e+11,
            'Kprime_0': 4.,
            'Debye_0': 422.,
            'grueneisen_0': 1.71,
            'q_0': 1.4,
            'G_0': .815e+11,
            'Gprime_0': 1.9,
            'eta_s_0': 7,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}
        
        
        Mineral.__init__(self)
'''
Mineral aliases
'''
perovskite=bridgmanite
