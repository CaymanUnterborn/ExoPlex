
# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.

"""
This example uses parallel processing to quickly calculate the best fit radius and pressure, temperature and density profiles
for a planet of a given mass, its uncertainty and a chosen composition.

The code begins by initializing the composition of the planet and retrieving the grids. In the main text code (at bottom)
one can set the number of samplings and the mass, radius, and uncertainties.

"""

import os
import sys
import multiprocessing as mp
import matplotlib.pyplot as plt

# hack to allow scripts to be placed in subdirectories next to exoplex:
import numpy as np

if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))
Pressure_range_mantle_UM = '1 1400000'
Temperature_range_mantle_UM = '1600 3500'

Pressure_range_mantle_LM = '1250000 40000000'
Temperature_range_mantle_LM = '1700 7000'

water_potential_temp = 300.

comp_keys = ['wt_frac_water','FeMg','SiMg','CaMg','AlMg','wt_frac_FeO_wanted','wt_frac_Si_core',
                          'wt_frac_O_core','wt_frac_S_core', 'combine_phases','use_grids','conserve_oxy']
struct_keys = ['Pressure_range_mantle_UM','Temperature_range_mantle_UM','resolution_UM',
                         'Pressure_range_mantle_LM', 'Temperature_range_mantle_LM', 'resolution_LM',
                         'Mantle_potential_temp','water_potential_temp']
combine_phases = True
use_grids = True

import ExoPlex as exo

from ExoPlex import run_perplex

# To have ExoPlex to give you compositional info and status of calculation set Verbose to TRUE.
# Note: setting this to True will slightly slow down the program
verbose = False

# create filename to store values

Output_filename = 'Filename'
# Next user must input the ratios by mole (Earth is Ca/Mg = .07, Si.Mg = 0.90, Al/Mg = 0.09, Fe/Mg = 0.9)
CaMg = 0.07
SiMg = 0.9
AlMg = 0.09
FeMg = 0.9

# How much water do you want in your planet? By mass fraction.
wt_frac_water = 0.0

# Don't forget that if you have water you need to add water layers
number_h2o_layers = 0

# The potential Temperature of Water, if present
water_potential_temp = 300.

# What fraction of the mantle would you like to be made of FeO? This Fe will be pulled from the core.
wt_frac_FeO_wanted = 0.  # by mass
conserve_oxy = False

# Now we can mix various elements into the core or mantle
wt_frac_Si_core = 0.  # by mass <1, note if you conserve oxygen this is calculated for you
wt_frac_O_core = 0.  # by mass
wt_frac_S_core = 0.  # by mass

# What potential temperature (in K) do you want to start your mantle adiabat?
Mantle_potential_temp = 1600.

# Input the resolution of your upper mantle and lower mantle composition, density grids
# These are input as number of T, P points. 50 50 = 2500 grid points, which takes about
# 5 minutes to calculate. Lower mantle resolution does not need to be higher since it's
# mostly ppv.
resolution_UM = '25 75'
resolution_LM = '75 75'

# lastly we need to decide how many layers to put in the planet. This is the resolution of
# the mass-radius sampling.
num_mantle_layers = 300
num_core_layers = 300

Output_radii = []
Output_mass = []

######### Initalize and run ExoPlex


compositional_params = dict(zip(comp_keys, [wt_frac_water, FeMg, SiMg, CaMg, AlMg, wt_frac_FeO_wanted, wt_frac_Si_core, \
                                            wt_frac_O_core, wt_frac_S_core, combine_phases, use_grids, conserve_oxy]))

if use_grids == True:
    filename = exo.functions.find_filename(compositional_params, verbose)
else:
    filename = ''

structure_params = dict(zip(struct_keys, [Pressure_range_mantle_UM, Temperature_range_mantle_UM, resolution_UM,
                                          Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                                          Mantle_potential_temp, water_potential_temp]))

layers = [num_mantle_layers, num_core_layers, number_h2o_layers]

Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac = exo.functions.get_percents(compositional_params, verbose)
Mantle_filename = exo.run_perplex.run_perplex(*[Mantle_wt_per,compositional_params,
                                                [structure_params.get('Pressure_range_mantle_UM'),structure_params.get('Temperature_range_mantle_UM'),
                                                structure_params.get('resolution_UM')],filename,verbose,combine_phases])
grids_low, names = exo.functions.make_mantle_grid(Mantle_filename,True,use_grids)
names.append('Fe')
if layers[-1] > 0:
    water_grid, water_phases = exo.functions.make_water_grid()
    for i in water_phases:
        names.append(i)
else:
    water_grid = []

Mantle_filename = exo.run_perplex.run_perplex(*[Mantle_wt_per,compositional_params,
                                            [structure_params.get('Pressure_range_mantle_LM'),structure_params.get('Temperature_range_mantle_LM'),
                                             structure_params.get('resolution_LM')],filename,verbose,False])
grids_high = exo.functions.make_mantle_grid(Mantle_filename,False,use_grids)[0]

core_grid = exo.functions.make_core_grid()

grids = [grids_low,grids_high,core_grid,water_grid]

def calc_planet(x):


    Planet = exo.functions.find_Planet_mass(x, core_mass_frac,structure_params, compositional_params, grids, Core_wt_per, layers,verbose)
    Planet['phase_names'] = names
    Planet['phases'],Planet['phase_names'] = exo.functions.get_phases(Planet, grids, layers,combine_phases)

    exo.functions.check(Planet)
    rad = Planet['radius'][-1]/6371e3
    mass = Planet['mass'][-1]/5.97e24
    CMF =  Planet['mass'][num_core_layers - 1] / Planet['mass'][-1]
    CRF =  Planet['radius'][num_core_layers - 1] / Planet['radius'][-1]
    CMB_P = Planet['pressure'][num_core_layers] / 1e4
    CMB_T = Planet['temperature'][num_core_layers]

    if number_h2o_layers > 0:
        WMB_P = Planet['pressure'][num_core_layers + num_mantle_layers] / 1e4
        WMB_T = Planet['temperature'][num_core_layers + num_mantle_layers]

    P_cen = Planet['pressure'][0] / 1e7
    T_cen = Planet['temperature'][0]

    if number_h2o_layers > 0:
        WMB_P = Planet['pressure'][num_core_layers + num_mantle_layers] / 1e4
        WMB_T = Planet['temperature'][num_core_layers + num_mantle_layers]
        keys = ['radius','mass','CMF','CRF','CMB_P','CMB_T','P_cen', 'T_cen','WMB_P','WMB_T']
        vals = [rad, mass, CMF, CRF, CMB_P, CMB_T, P_cen, T_cen, WMB_P, WMB_T]
        return(dict(zip(keys, vals)))

    else:
        keys = ['radius','mass','CMF','CRF','CMB_P','CMB_T','P_cen', 'T_cen']
        vals = [rad, mass, CMF, CRF, CMB_P, CMB_T, P_cen, T_cen]
        return(dict(zip(keys, vals)))


if __name__ == "__main__":
    num_pts = 5
    M = 2.49
    M_err = 0.425

    Mass_planet = np.random.normal(M, M_err, num_pts)

    pool = mp.Pool(processes=mp.cpu_count())

    Planets = pool.map_async(calc_planet,Mass_planet).get()

    pool.close()


    plt.scatter(Mass_planet, [Planets[i].get('radius') for i in range(num_pts)])
    plt.ylabel('Radius (Earth Radii)', size=20)
    plt.xlabel('Mass (Earth Masses)', size=20)
    plt.show()


