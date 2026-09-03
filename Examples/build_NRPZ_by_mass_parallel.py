
# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.

import os
import sys
import scipy.stats as sp
import requests
import multiprocessing as mp
import time

# hack to allow scripts to be placed in subdirectories next to exoplex:
import numpy as np
import pandas as pd


if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))
    
from ExoPlex import functions
from ExoPlex import run_perplex as perp
from ExoPlex import make_grids


#-----------------Output directory for rocky planet zone(s)-------------------#
output_path = './'
#-----------------------------------------------------------------------------#

#------------------------Mass ranges for the RZ-------------------------------#
mstart = 0.1 #minimum desired mass in Earth masses
mstop = 10 #maximum mass in Earth masses
num_mass_steps = 21
mass_array = 10**np.linspace(np.log10(mstart), np.log10(mstop), num_mass_steps)
#-----------------------------------------------------------------------------#

#---------------------Host star input data or NRPZ----------------------------#
# Input the star's name to generate a star-specific rocky planet zone. If you
# want to build a neaby galactic rocky planet zone, set the star_name = 'nearby_galactic'.  
# If you input a star's name, but the star does not have measured abundances 
# of Fe, Mg, AND Si, in Hypatia the code will default to building a nearby 
# galactic rocky planet zone. 

# Additionally, you need to specify the desired solar normalization from: 
# 'asplund05', 'lodders09', 'anders89', 'grevesse98', 'asplund09', 
# 'grevesse07', 'absolute', 'original'
# The Hypatia default is 'lodders09'. See Hypatia API page for more details.

star_name = 'HD 176981'
solar_norm_name = 'lodders09'
num_abundance_samples = 500 

# If you want to overwrite the default behavoir b/c, e.g., you have host abundances
# that are not yet listed on Hypatia, you can input them manually by uncommenting
# the lines below.

#FeH = -0.02; sig_FeH = 0.16
#MgH = -0.06; sig_MgH = 0.07
#SiH = -0.05; sig_SiH = 0.12

#-----------------------------------------------------------------------------#

#-------------------------------ExoPlex inputs--------------------------------#
# ExoPlex inputs if you feel the need to change them, but we suggest leaving
# them as is.

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

#-----------------------------------------------------------------------------#


def run_planet(*args):
    Mass_planet, compositional_params['FeMg'], compositional_params['SiMg'] = args
    Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac = functions.get_percents(compositional_params, verbose)
    try:
        Planet = functions.find_Planet_mass(Mass_planet, core_mass_frac,structure_params, compositional_params, grids, Core_wt_per, layers,verbose)
        return Planet['radius'][-1]/6371e3
    except:
        return np.nan

build_gal_rz = False

if star_name != 'nearby_galactic':
    # The following code first checks if the user has input their own abundances for
    # Fe, Mg, and Si. If no user input Fe, Mg, and Si values are present, 
    # then the code checks that star is listed in Hypatia AND has [Fe/H], [Mg/H], 
    # and [Si/H] abunds. If all four criteria are not satisfied, then it generates 
    # a galactic rocky planet zone, instead.
    
    if 'FeH' not in locals() or 'SiH' not in locals() or 'MgH' not in locals():

        'Checking Hypatia for desired star...'
        params = {"name": star_name, "element": ["fe"], "solarnorm": [solar_norm_name]}
        star_entry_fe = requests.get("https://hypatiacatalog.com/hypatia/api/v2/composition", params=params)
        star_entry_fe = star_entry_fe.json()[0]
        if star_entry_fe['name'] == 'not-found':
            print(f'{star_name} not found in Hypatia.')
            print('Building nearby galactic RZ instead...')
            build_gal_rz = True
        else:            
            FeH = star_entry_fe['median_value']
            
            params = {"name": star_name, "element": ["mg"], "solarnorm": [solar_norm_name]}
            star_entry_mg = requests.get("https://hypatiacatalog.com/hypatia/api/v2/composition", params=params)
            star_entry_mg = star_entry_mg.json()[0]
            MgH = star_entry_mg['median_value']
            
            params = {"name": star_name, "element": ["si"], "solarnorm": [solar_norm_name]}
            star_entry_si = requests.get("https://hypatiacatalog.com/hypatia/api/v2/composition", params=params)
            star_entry_si = star_entry_si.json()[0]
            SiH = star_entry_si['median_value']
        
            if FeH == None or MgH == None or SiH == None:
                print(f'{star_name} found in Hypatia, but it is missing one or more of the major rock-building abundances.')
                print('Building nearby galactic RZ instead...')
                build_gal_rz = True
            else:
                sig_FeH = star_entry_fe['plusminus']
                sig_MgH = star_entry_mg['plusminus']
                sig_SiH = star_entry_si['plusminus']

                print(f'Abundance values for {star_name} in Hypatia relative to {solar_norm_name}:')
                print(f"[Fe/H] = {FeH:.2f} +/- {sig_FeH:.2f}") 
                print(f"[Mg/H] = {MgH:.2f} +/- {sig_MgH:.2f}") 
                print(f"[Si/H] = {SiH:.2f} +/- {sig_SiH:.2f}") 
                build_gal_rz = False
            
    else:
        print(f'Building RZ for user input abundances with {solar_norm_name} normalization:')
        print(f"[Fe/H] = {FeH:.2f} +/- {sig_FeH:.2f}") 
        print(f"[Mg/H] = {MgH:.2f} +/- {sig_MgH:.2f}") 
        print(f"[Si/H] = {SiH:.2f} +/- {sig_SiH:.2f}") 
        
        build_gal_rz = False


if star_name == 'nearby_galactic' or build_gal_rz == True:
    star_name = 'nearby_galactic'
    hypatia_names = np.loadtxt('../hypatia_names.txt', dtype = str)
    inds = np.random.randint(0, len(hypatia_names), num_abundance_samples)
    hypatia_names = hypatia_names[inds]
    hypatia_names = list(hypatia_names)
    SiH = np.array([])
    MgH = np.array([])
    FeH = np.array([])
    for star in hypatia_names:
        params = {"name": star, "element": ["si"], "solarnorm": [solar_norm_name]}
        star_entry = requests.get("https://hypatiacatalog.com/hypatia/api/v2/composition", params=params)
        star_entry = star_entry.json()[0]
        SiH = np.append(SiH, star_entry['median_value'])

        params = {"name": star, "element": ["mg"], "solarnorm": [solar_norm_name]}
        star_entry = requests.get("https://hypatiacatalog.com/hypatia/api/v2/composition", params=params)
        star_entry = star_entry.json()[0]
        MgH = np.append(MgH, star_entry['median_value'])
        
        params = {"name": star, "element": ["fe"], "solarnorm": [solar_norm_name]}
        star_entry = requests.get("https://hypatiacatalog.com/hypatia/api/v2/composition", params=params)
        star_entry = star_entry.json()[0]
        FeH = np.append(FeH, star_entry['median_value'])
        
else:
    FeH = sp.norm.rvs(FeH, sig_FeH, num_abundance_samples)
    MgH = sp.norm.rvs(MgH, sig_MgH, num_abundance_samples)
    SiH = sp.norm.rvs(SiH, sig_SiH, num_abundance_samples)
    

# This block of code gets the A(X) values of each element X for the specified
# solar normalization name. See Hypatia API page for more details.
get_normalizations = requests.get("https://hypatiacatalog.com/hypatia/api/v2/solarnorm")
normalizations = get_normalizations.json()
solar_norm_list = np.array([normalizations[i]['id'] for i in range(0, len(normalizations))])
solar_norm = normalizations[np.where(solar_norm_list == solar_norm_name)[0][0]]
#print('Solar norm check. This should be the name of the specified normalization: ', solar_norm['id'])

sol_FeH = solar_norm['values']['Fe']
sol_SiH = solar_norm['values']['Si']
sol_MgH = solar_norm['values']['Mg']


# Finally calculate molar Fe/Mg and Si/Mg of the major rock-building elements
FeMg_array = 10**((FeH+sol_FeH)-(MgH+sol_MgH))
SiMg_array = 10**((SiH+sol_SiH)-(MgH+sol_MgH))

# EP grids go up to 0.1 < Si/Mg < 2.0 based off of the 2-sigma bounds of the Hypatia
# The following two lines of code shouldn't trigger much unless you are using poorly
# constrained Si/Mg values. While changing Si/Mg shouldn't ...
SiMg_array[np.where(SiMg_array >= 2.0)] = 1.99*np.ones(len(np.where(SiMg_array >= 2.0)[0]))
SiMg_array[np.where(SiMg_array <= 0.1)] = 0.101*np.ones(len(np.where(SiMg_array <= 0.1)[0]))



#-----------------------------------------------------------------------------#


#df = pd.DataFrame({'FeMg':FeMg_array, 'SiMg':SiMg_array})


if __name__ == "__main__":
    
        X = np.array([mass*np.ones(len(FeMg_array)) for mass in mass_array]).flatten()
        Y = np.array([FeMg_array for mass in mass_array]).flatten()
        Z = np.array([SiMg_array for mass in mass_array]).flatten()

        vals = zip(X, Y, Z)
        start = time.time()

        pool = mp.Pool(processes=mp.cpu_count())
        radius = pool.starmap(run_planet,vals)
        pool.close()
        
        end = time.time()
        length = round((end-start)/60.0,1)
        print(f"Done with RZ for {star_name}. Took {length} minutes.")
        print(f"Saving in {output_path}.")
        
        radius = np.array(radius)
        
        radius_lower_3sigma = np.zeros(len(mass_array))
        radius_median = np.zeros(len(mass_array))        
        radius_upper_3sigma = np.zeros(len(mass_array))
        
        for i in range(0, len(mass_array)):
            inds = np.where(X == mass_array[i])
            radius_lower_3sigma[i], radius_median[i], radius_upper_3sigma[i] = np.quantile(radius[inds], [0.5-(0.997/2), 0.5, 0.5+(0.997/2)])
            
        pd.DataFrame({'mass': mass_array, 
                      'radius_lower_3sigma': radius_lower_3sigma,
                      'radius_median': radius_median,
                      'radius_upper_3sigma': radius_upper_3sigma}).to_csv(output_path + star_name + '_RZ.csv', index = False)
            
        
            



