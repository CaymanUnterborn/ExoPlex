
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
mFe = 55.845
mMg = 24.306
mSi = 28.0867
mO = 15.9994
mS = 32.0650
mCa = 40.078
mAl = 26.981

if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))

from ExoPlex import functions
from ExoPlex import run_perplex as perp
from ExoPlex import make_grids
from ExoPlex import planet
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
AlMg = 0.09

# How much water do you want in your planet? By mass fraction.
wt_frac_water = 0.0

# Don't forget that if you have water you need to add water layers
number_h2o_layers = 0

# The potential Temperature of Water, if present
water_potential_temp = 300

# What fraction of the mantle would you like to be made of FeO? This Fe will be pulled from the core.
wt_frac_FeO_wanted = 0.  # by mass
conserve_oxy = False

# Now we can mix various elements into the core or mantle
wt_frac_Si_core = 0.  # by mass <1, note if you conserve oxygen this is calculated for you
wt_frac_O_core = 0.1  # by mass
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

core_grid = make_grids.make_core_grid()


structure_params = dict(zip(struct_keys, [Pressure_range_mantle_UM, Temperature_range_mantle_UM, resolution_UM,
                                          Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                                          Mantle_potential_temp, water_potential_temp]))

layers = [num_mantle_layers, num_core_layers, number_h2o_layers]
water_grid = make_grids.make_water_grid()[0]


Mass_planet = [round(.25 + .025 * i, 3) for i in range(31)]
Mass_planet = Mass_planet + [round(1.1 + .1 * i, 2) for i in range(40)]
Mass_planet = Mass_planet + [round(5.25 + .25 * i, 3) for i in range(34)]

Mass_planet = [round(Mass_planet[-1]+.25 + .25*i, 3) for i in range(4)]

def make_planet(femg, simg):
    compositional_params = dict(
        zip(comp_keys, [wt_frac_water, femg, simg, CaMg, AlMg, wt_frac_FeO_wanted, wt_frac_Si_core, \
                        wt_frac_O_core, wt_frac_S_core, combine_phases, use_grids, conserve_oxy]))

    filename = functions.find_filename(compositional_params, verbose)
    Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac = functions.get_percents(compositional_params, verbose)
    Mantle_filename = perp.run_perplex(
        *[Mantle_wt_per, compositional_params, structure_params, filename, verbose, combine_phases])
    grids_low, names = make_grids.make_mantle_grid(Mantle_filename, Mantle_wt_per, True, use_grids)

    grids_high = make_grids.make_mantle_grid(Mantle_filename, Mantle_wt_per, False, use_grids)[0]
    grids = [grids_low, grids_high, core_grid, water_grid]

    params = [grids, Core_wt_per,compositional_params,core_mass_frac]
    Rad = []
    for i in range(len(Mass_planet)):
        Planet = planet.initialize_by_mass(*[Mass_planet[i], structure_params, compositional_params, layers, core_mass_frac])

        try:
            Pl = calc_planet(params, Planet)
            functions.check(Pl)
        except:
            print("ere")
            Rad.append(np.nan)
        else:
            Rad.append(Pl['radius'][-1] / 6371e3)
    output = [[Mass_planet[i], femg, simg, Rad[i], 6.67e-11*Mass_planet[i]*5.97e24/pow(Rad[i]*6371e3,2), \
               Mass_planet[i]*5.97e24/((4*np.pi/3)*pow(Rad[i]*6371e3,3)), core_mass_frac] for i in range(len(Mass_planet))]
    return(output)


def calc_planet(params,Planet):
    grids = params[0]
    Core_wt_per = params[1]
    compositional_params = params[2]
    core_mass_frac = params[3]
    Planet = planet.compress_mass(*[Planet, grids, Core_wt_per, structure_params, compositional_params, core_mass_frac, layers,verbose])

    return(Planet)


if __name__ == "__main__":
    num_pts = 500
    hyp = open('hyp.csv', 'r')
    temp_file = hyp.readlines()

    data = temp_file[1:]
    num_rows = len(temp_file[1:])
    num_columns = len(temp_file[12].split(','))

    FeMg_get = []
    SiMg_get = []
    names = []
    for i in range(num_rows):
        columns = data[i].strip('\n').split(',')[1:]
        if columns[0] != '' and columns[1]!= '' and columns[2]!= '':
            FeMg_test = (pow(10,7.45-7.54)*pow(10,float(columns[0])-float(columns[1])))
            SiMg_test = ( pow(10,7.52-7.54)*pow(10,float(columns[2])-float(columns[1])))

            if FeMg_test <= 2 and SiMg_test <= 2:
                names.append(data[i].strip('\n').split(',')[0])
                FeMg_get.append(float(FeMg_test))
                SiMg_get.append(float(SiMg_test))

    dat = []
    for i in range(len(names)):
        dat.append([FeMg_get[i],SiMg_get[i]])

    #print(FeMg_get[0])
    #print(len(names))
    data = np.vstack((FeMg_get,SiMg_get)).T

    #head = 'femg,simg'
    #np.savetxt("hyp_sol.csv",data,header=head,fmt ='%f',delimiter=',',newline='\n',comments='#')

    np.stack((FeMg_get, SiMg_get))
    import scipy.stats as st
    data = np.vstack((FeMg_get,SiMg_get))
    cov = np.cov(data,bias= True)

    mean = [st.norm.fit(FeMg_get)[0],st.norm.fit(SiMg_get)[0] ]

    fe_new,si_new = np.random.default_rng().multivariate_normal(mean, cov, size=num_pts,check_valid='warn').T

    header = "Mass, FeMg, SiMg, Radius, Grav, Den, CMF"


    compositional_params = dict(
        zip(comp_keys,
            [wt_frac_water,1, 1, CaMg, AlMg, wt_frac_FeO_wanted, wt_frac_Si_core, \
             wt_frac_O_core, wt_frac_S_core, combine_phases, use_grids, conserve_oxy]))
    dat = []

    count = 0
    #file  = open('comps_water.csv','r')
    #temp_file = file.readlines()
    #data = temp_file[1:]

    #dats = np.asarray([(data[i].strip('\n').split(',')) for i in range(len(data))])

    #Fe_test = [float(i[0]) for i in range(len(data))]

    #fe_new  = [float(i) for i in dats[:,0]]
    #si_new = [float(i) for i in dats[:,1]]


    for i in range(len(si_new)):
        good = False
        while good == False:
            compositional_params['FeMg'] = fe_new[i]
            compositional_params['SiMg'] = si_new[i]

            #print(compositional_params)
            try:
                Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac = functions.get_percents(compositional_params, verbose)
            except:
                fe_new[i], si_new[i] = np.random.default_rng().multivariate_normal(mean, cov, 1)[0]
                count+=1
            else:
                SiMg_mant = ((Mantle_wt_per.get('SiO2') / Mantle_wt_per.get('MgO') * (40.30 / 60.08)))
                wt_frac_Si = Core_wt_per.get('Si')
                wt_frac_O = Core_wt_per.get('O')
                wt_frac_S = Core_wt_per.get('S')
                wt_frac_Fe = Core_wt_per.get('Fe')

                mol_total = ((wt_frac_Fe / mFe) + (wt_frac_O / mO) + (wt_frac_S / mS) + (wt_frac_Si / mSi)) / 100
                mol_frac_Fe = (wt_frac_Fe / mFe / 100) / mol_total
                mol_frac_Si = (wt_frac_Si / mSi / 100) / mol_total
                mol_frac_S = (wt_frac_S / mS / 100) / mol_total
                mol_frac_O = (wt_frac_O / mO / 100) / mol_total

                molar_weight_core = (mol_frac_Fe * mFe) + (mol_frac_Si * mSi) + (mol_frac_O * mO) + (mol_frac_S * mS)
                core_diff = molar_weight_core/mFe
                if si_new[i] > 2 or fe_new[i] <=0 or si_new[i] < 0.1 or SiMg_mant < 0.1:
                    fe_new[i], si_new[i] = np.random.default_rng().multivariate_normal(mean, cov, 1)[0]
                    print("bad")
                    count += 1

                else:
                    dat.append([fe_new[i], si_new[i], core_mass_frac,SiMg_mant,core_diff,wt_frac_Si])
                    good = True

    print(count,"bad")
    """
    from matplotlib import pyplot as plt
    plt.scatter(FeMg_get, SiMg_get,alpha=0.25)
    plt.scatter(fe_new,si_new,c='r')
    plt.xlim(0,2)
    plt.ylim(0,2)
    plt.show()
    """

    #head = 'femg,simg,CMF,SiMg_mant, core_rho_diff,Si_core'
    #np.savetxt("comps_lightel_2.csv",dat,header=head,fmt='%.10s',delimiter=',')
    #rewrite to make it chooose one composition, run every mass
    print(min(si_new))

    file = open('comps_lightel.csv')

    temp_file = file.readlines()

    num_rows = len(temp_file[1:])
    num_columns = len(temp_file[12].split(','))

    fe_new = []
    si_new = []
    for i in range(num_rows + 1)[1:501]:
        row = temp_file[i].strip('\n').split(',')
        fe_new.append(float(row[0]))
        si_new.append(float(row[1]))
    print(min(si_new))

    vals = zip(fe_new,si_new)
    pool = mp.Pool(processes=mp.cpu_count())

    Rad = np.vstack(pool.starmap_async(make_planet,vals).get())

    pool.close()

    np.savetxt("scan_cov_lightel_2.csv", Rad, delimiter=',', header=header, fmt='%.5f')
    #np.savetxt("d.csv", Rad, delimiter=',', header=header, fmt='%.5f')

