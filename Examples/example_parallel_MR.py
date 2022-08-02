
# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.


"""
This example uses parallel processing to quickly calculate the best fit Fe/Mg and CMF for a planet with a given
Mass, Radius and their respective uncertainties.

The code begins by initializing the composition of the planet and retrieving the grids. In the main text code (at bottom)
one can set the number of samplings and the mass, radius, and uncertainties.
"""

import os
import sys
from scipy.stats import norm
import matplotlib.pyplot as plt
import multiprocessing as mp
import statistics
import scipy.stats as sp
from scipy.optimize import root_scalar

# hack to allow scripts to be placed in subdirectories next to exoplex:
import numpy as np

if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))

from ExoPlex import functions
from ExoPlex import run_perplex as perp
from ExoPlex import make_grids

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


# To have ExoPlex to give you compositional info and status of calculation set Verbose to TRUE.
# Note: setting this to True will slightly slow down the program
verbose = False

# Next user must input the ratios by mole (Earth is Ca/Mg = .07, Si.Mg = 0.90, Al/Mg = 0.09, Fe/Mg = 0.9)
CaMg = 0.07
SiMg = 0.9
AlMg = 0.09
FeMg = 1.

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
num_mantle_layers = 400
num_core_layers = 500

Output_radii = []
Output_mass = []

######### Initalize and run ExoPlex


compositional_params = dict(zip(comp_keys, [wt_frac_water, FeMg, SiMg, CaMg, AlMg, wt_frac_FeO_wanted, wt_frac_Si_core, \
                                            wt_frac_O_core, wt_frac_S_core, combine_phases, use_grids, conserve_oxy]))

if use_grids == True:
    filename = functions.find_filename(compositional_params, verbose)
else:
    filename = ''

structure_params = dict(zip(struct_keys, [Pressure_range_mantle_UM, Temperature_range_mantle_UM, resolution_UM,
                                          Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                                          Mantle_potential_temp, water_potential_temp]))

layers = [num_mantle_layers, num_core_layers, number_h2o_layers]

Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac = functions.get_percents(compositional_params, verbose)
Mantle_filename = perp.run_perplex(*[Mantle_wt_per,compositional_params,structure_params,filename,verbose,combine_phases])
grids_low, names = make_grids.make_mantle_grid(Mantle_filename,Mantle_wt_per, True,use_grids)
names.append('Fe')
if layers[-1] > 0:
    water_grid, water_phases = make_grids.make_water_grid()
    for i in water_phases:
        names.append(i)
else:
    water_grid = []

grids_high = make_grids.make_mantle_grid(Mantle_filename,Mantle_wt_per, False,use_grids)[0]

core_grid = make_grids.make_core_grid()

grids = [grids_low,grids_high,core_grid,water_grid]


def run_planet(x, *args):
    Mass_planet, Rad_planet = args
    compositional_params['FeMg'] = x
    Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac = functions.get_percents(compositional_params, verbose)

    Planet = functions.find_Planet_mass(Mass_planet, core_mass_frac,structure_params, compositional_params, grids, Core_wt_per, layers,verbose)
    out = 1 - (Planet['radius'][-1] / 6371e3) / Rad_planet
    return (out)

def calc_planet(mass, radius):
    den = mass * 5.97e21 / ((4 * np.pi / 3) * pow((radius * 6371e3),3))
    min = 0.001
    if den < 10:
        max = 10
    else:
        max = 80

    try:
        FeMg = root_scalar(run_planet,bracket=[min,max] ,args=(mass, radius),x0 = 0.9,xtol=0.0001).root
    except:
        test = run_planet(min, *(mass, radius))
        if test < 0:
            #planet with smallest core produces radius too big
            #return very high FeMg

            return(95)

        if test > 0:
            #planet with largest core produces radius too big
            #return very low FeMg
            return(1e-6)
    else:
        return(FeMg)

if __name__ == "__main__":
    num_pts = 5

    filename = 'Earth'
    R=  1.
    R_err =0.05
    M = 1
    M_err = .1

    Mass_planet = np.random.normal(M, M_err, num_pts)

    Radius_planet = np.random.normal(R, R_err, num_pts)

    cov = [[pow(M_err,2),0],[0, pow(R_err,2)]]
    mean = [M, R]

    Mass_planet, Radius_planet = np.random.multivariate_normal(mean,cov, num_pts).T

    vals = zip(Mass_planet, Radius_planet)
    pool = mp.Pool(processes=mp.cpu_count())

    FeMg = pool.starmap_async(calc_planet,vals).get()

    pool.close()

    CMF = []

    for i in range(len(FeMg)):
        if FeMg[i] >0:
            compositional_params['FeMg'] = FeMg[i]
            CMF.append(functions.get_percents(compositional_params,verbose)[3])
        else:
            CMF.append(FeMg[i])

    mu= statistics.mean(FeMg)
    std = statistics.stdev(FeMg)
    print()
    print("Actual Fe/Mg:")
    print(round(mu, 2), "+/-", round(std, 2))
    print()

    #mu, std = sp.lognorm.fit(FeMg)
    shape, loc, scale = sp.lognorm.fit(FeMg)

    mu = np.exp(np.log(scale) + pow(shape,2)/2)
    std = mu*np.sqrt(pow(shape,2)-1)
    print()
    print("Log-Normal best fit Fe/Mg:")
    print(round(mu,2), "+/-",round(std,2))
    print()

    mu_CMF, std_CMF = sp.norm.fit(CMF)
    print()
    print("Normal best fit CMF:")
    print(round(mu_CMF,2), "+/-",round(std_CMF,2))
    print()

    output = list(zip(Mass_planet, Radius_planet,FeMg,CMF))

    header ='#Mass,Radius,FeMg,CMF'
    head = []
    names = header.split(',')

    for i in names:
        head.append(i+'_'+filename)
    header = ','.join(head)

    np.savetxt(filename+'.csv', output, header = header, comments='#', delimiter=',')

    fig, (ax1, ax2,ax3) = plt.subplots(1, 3, figsize=(10, 5))
    
    
    ax1.scatter(Mass_planet, Radius_planet)
    ax1.set_ylabel('Radius (Earth Radii)', size=20)
    ax1.set_xlabel('Mass (Earth Masses)', size=20)
    ax1.minorticks_on()
    ax1.set_xlim(M-5*M_err, M+5*M_err)
    ax1.set_ylim(R-5*R_err, R+5*R_err)

    ax2.hist(FeMg, bins=50, density=True, alpha=0.6, color='g')    #plt.show()
    xmin, xmax = ax2.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    p = sp.lognorm.pdf(x,shape, loc = loc, scale = scale)
    ax2.plot(x, p, 'k', linewidth=2)
    ax2.set_xlabel('Fe/Mg', size=20)
    ax2.set_xlim(0,max(FeMg))
    #plt.vlines(mu, 0, 2)

    ax3.hist(CMF, bins=50, density=True, alpha=0.6, color='g')    #plt.show()
    xmin, xmax = ax3.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu_CMF, std_CMF)
    ax3.plot(x, p, 'k', linewidth=2)
    ax3.set_xlabel('CMF', size=20)
    ax3.set_xlim(0,1)

    plt.savefig(filename+'.jpg', format = 'jpg')


