
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
from ExoPlex import exoplex_blackbox as epbb

#------------------- Inputs you'll likely want to change ---------------------#
#------------------------ Planet observables ---------------------------------#
R=  0.95
R_err = 0.01
M = 1.0
M_err = 0.05

planet = 'generic_super_mercury'
#-----------------------------------------------------------------------------#


#----------------- Host star compositional parameters ------------------------#
#----- Can input in terms of molar ratios or indivual abundances in dex ------#
#---- If the latter, make sure to double check your solar normalization! -----#
FeMg = [0.9, 0.1]
SiMg = 1.0
CaMg = 0.07
AlMg = 0.09

#FeH = [0.0, 0.05]
#MgH = [0.0, 0.05]
#SiH = 0.0
#CaH = 0.0
#AlH = 0.0

#A_Fe_sol = 7.5; A_Mg_sol = 7.6; A_Si_sol = 7.51; A_Ca_sol = 6.34; A_Al_sol = 6.45
#The above Solar values are from Asplund+09

#FeH = sp.norm.rvs(loc = FeH[0], scale = FeH[1], size = 1000)
#MgH = sp.norm.rvs(loc = MgH[0], scale = MgH[1], size = 1000)
#FeMg = 10**((FeH+A_Fe_sol)-(MgH+A_Mg_sol))
#FeMg = [sp.norm.fit(FeMg)[0], sp.norm.fit(FeMg)[1]]
#SiMg = 10**((SiH+A_Si_sol)-(MgH[0]+A_Mg_sol))
#AlMg = 10**((AlH+A_Al_sol)-(MgH[0]+A_Mg_sol))
#CaMg = 10**((CaH+A_Ca_sol)-(MgH[0]+A_Mg_sol))
#-----------------------------------------------------------------------------#

#------------------------ Other important inputs -----------------------------#
#filename = 'Earth'
#How many M-R and stellar abundance pairs do you want?
num_pts = 1000
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#


#------------------- Inputs you may want to change ---------------------------#
#-------------------- other comp/chemistry params ----------------------------#
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
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#



#------------------- Inputs probably won't want to change --------------------#
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
#-----------------------------------------------------------------------------#


#-------------------------- Under the hood -----------------------------------#
#Get the stellar Fe/Mg and CMF
FeMg_star = sp.norm.rvs(FeMg[0], FeMg[1], num_pts)

CMF_star = np.zeros(len(FeMg_star))
for i in range(0, len(FeMg_star)):
    stellar_compositional_params = dict(zip(comp_keys, [wt_frac_water, FeMg_star[i], SiMg, CaMg, AlMg, wt_frac_FeO_wanted, wt_frac_Si_core, \
                                            wt_frac_O_core, wt_frac_S_core, combine_phases, use_grids, conserve_oxy]))
    Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac = functions.get_percents(stellar_compositional_params, verbose)
    CMF_star[i] = core_mass_frac
    


######### Initalize and run ExoPlex
Output_radii = []
Output_mass = []

compositional_params = dict(zip(comp_keys, [wt_frac_water, FeMg[0], SiMg, CaMg, AlMg, wt_frac_FeO_wanted, wt_frac_Si_core, \
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
    g = 6.67e-11*Mass_planet*5.97e24/(pow(Planet['radius'][-1],2))
    g_act = 6.67e-11*Mass_planet*5.97e24/(pow(Rad_planet*6371e3,2))

    out = 1 - (g) / g_act
    return (out)

def calc_planet(mass, radius):
    try:
        den = mass * 5.97e21 / ((4 * np.pi / 3) * pow((radius * 6371e3),3))
        min = 1e-10
        if den < 10:
            max = 6
        else:
            max = 30
    
        try:
            FeMg = root_scalar(run_planet,bracket=[min,max] ,args=(mass, radius),x0 = 0.9,xtol=0.0001).root
        except:
            test = run_planet(min, *(mass, radius))
            if test < 0:
                # planet with smallest core produces radius too big
                # return very high FeMg
    
                return (1e-11)
    
            if test > 0:
                # planet with largest core produces radius too big
                # return very low FeMg
                return (95)
        else:
            return (FeMg)
    except:
        return (1e-13)

if __name__ == "__main__":
    Mass_planet = np.random.normal(M, M_err, num_pts)
    ind_keep = np.where(Mass_planet > 0)
    Mass_planet = Mass_planet[ind_keep]

    Radius_planet = np.random.normal(R, R_err, num_pts)
    Radius_planet = Radius_planet[ind_keep]

    cov = [[pow(M_err,2),0],[0, pow(R_err,2)]]
    mean = [M, R]

    Mass_planet, Radius_planet = np.random.multivariate_normal(mean,cov, num_pts).T

    vals = zip(Mass_planet, Radius_planet)
    pool = mp.Pool(processes=mp.cpu_count())

    FeMg = pool.starmap_async(calc_planet,vals).get()

    pool.close()

    CMF = []
    
    #CMF = np.zeros_like(FeMg)
    
    
    for i in range(len(FeMg)):
        if FeMg[i] >0:
            compositional_params['FeMg'] = FeMg[i]
            CMF.append(functions.get_percents(compositional_params,verbose)[3])
        else:
            CMF.append(FeMg[i])

    mu= statistics.mean(FeMg)
    std = statistics.stdev(FeMg)

    #mu, std = sp.lognorm.fit(FeMg)
    shape, loc, scale = sp.lognorm.fit(FeMg)
    shape_star, loc_star, scale_star = sp.lognorm.fit(FeMg_star)
    
    femg_rho_x1sig_low, femg_rho_xmed, femg_rho_x1sig_up = np.quantile(FeMg, [(1.0 - 0.68)/2.0, 0.5, 1.0 - (1.0 - 0.68)/2.0])
    femg_rho_y1sig_low, femg_rho_ymed, femg_rho_y1sig_up = sp.lognorm.pdf([femg_rho_x1sig_low, femg_rho_xmed, femg_rho_x1sig_up], shape, loc = loc, scale = scale)

    femg_star_x1sig_low, femg_star_xmed, femg_star_x1sig_up = np.quantile(FeMg_star, [(1.0 - 0.68)/2.0, 0.5, 1.0 - (1.0 - 0.68)/2.0])
    femg_star_y1sig_low, femg_star_ymed, femg_star_y1sig_up = sp.lognorm.pdf([femg_star_x1sig_low, femg_star_xmed, femg_star_x1sig_up], shape_star, loc = loc_star, scale = scale_star)

    slround = str(round(femg_rho_xmed - femg_rho_x1sig_low,2))
    suround = str(round(femg_rho_x1sig_up-femg_rho_xmed,2))
    medround = str(round(femg_rho_xmed, 2))
    print('Planet Fe/Mg = %s (+%s, -%s)' % (medround, suround, slround))

    slround = str(round(femg_star_xmed - femg_star_x1sig_low,2))
    suround = str(round(femg_star_x1sig_up-femg_star_xmed,2))
    medround = str(round(femg_star_xmed, 2))
    print('Star Fe/Mg = %s (+%s, -%s)' % (medround, suround, slround))
    print()
    
    
    mu_CMF_norm, std_CMF_norm = sp.norm.fit(CMF)
    shape_cmf_lognorm, loc_cmf_lognorm, scale_cmf_lognorm = sp.lognorm.fit(CMF)
    shape_cmf_star_lognorm, loc_cmf_star_lognorm, scale_cmf_star_lognorm = sp.lognorm.fit(CMF_star)

    cmfrho_x1sig_low, cmfrho_xmed, cmfrho_x1sig_up = np.quantile(CMF, [(1.0 - 0.68)/2.0, 0.5, 1.0 - (1.0 - 0.68)/2.0])
    cmfrho_y1sig_low, cmfrho_ymed, cmfrho_y1sig_up = sp.lognorm.pdf([cmfrho_x1sig_low, cmfrho_xmed, cmfrho_x1sig_up], shape_cmf_lognorm, loc = loc_cmf_lognorm, scale = scale_cmf_lognorm)

    cmfstar_x1sig_low, cmfstar_xmed, cmfstar_x1sig_up = np.quantile(CMF_star, [(1.0 - 0.68)/2.0, 0.5, 1.0 - (1.0 - 0.68)/2.0])
    cmfstar_y1sig_low, cmfstar_ymed, cmfstar_y1sig_up = sp.lognorm.pdf([cmfstar_x1sig_low, cmfstar_xmed, cmfstar_x1sig_up], shape_cmf_star_lognorm, loc = loc_cmf_star_lognorm, scale = scale_cmf_star_lognorm)

    slround = str(round(cmfrho_xmed - cmfrho_x1sig_low,2))
    suround = str(round(cmfrho_x1sig_up-cmfrho_xmed,2))
    medround = str(round(cmfrho_xmed, 2))
    print('Planet CMF = %s (+%s, -%s)' % (medround, suround, slround))

    slround = str(round(cmfstar_xmed - cmfstar_x1sig_low,2))
    suround = str(round(cmfstar_x1sig_up-cmfstar_xmed,2))
    medround = str(round(cmfstar_xmed, 2))
    print('Star CMF = %s (+%s, -%s)' % (medround, suround, slround))


    def numerator_function(x): 
        return sp.lognorm.pdf(x, shape_cmf_star_lognorm, loc = loc_cmf_star_lognorm, scale = scale_cmf_star_lognorm)*sp.lognorm.pdf(x, shape_cmf_lognorm, loc = loc_cmf_lognorm, scale = scale_cmf_lognorm)

    def denominator_function(x):
        xs = x + (cmfrho_xmed - cmfstar_xmed)
        return sp.lognorm.pdf(x, shape_cmf_star_lognorm, loc = loc_cmf_star_lognorm, scale = scale_cmf_star_lognorm)*sp.lognorm.pdf(xs, shape_cmf_lognorm, loc = loc_cmf_lognorm, scale = scale_cmf_lognorm)

    from scipy.integrate import quad
    numerator = quad(numerator_function, 0, 1)[0]
    denominator = quad(denominator_function, 0, 1)[0]
    ph0 = 100*numerator/denominator
    print('P(H0) = ', round(ph0,1))
        
        
    
    

    output = list(zip(Mass_planet, Radius_planet,FeMg,CMF))

    header ='#Mass,Radius,FeMg,CMF'
    head = []
    names = header.split(',')

    for i in names:
        head.append(i+'_'+filename)
    header = ','.join(head)

    #np.savetxt(planet + '_' + filename+'.csv', output, header = header, comments='#', delimiter=',')

    fig, (ax1, ax2,ax3) = plt.subplots(1, 3, figsize=(15, 5))
    
    import pandas as pd
    
    df = pd.DataFrame({'Mass': Mass_planet, 'Radius':Radius_planet, 'CMF': CMF, 'FeMg': FeMg}).to_csv(planet + '_ ' + filename + '.csv')
    df = pd.DataFrame({'FeMg': FeMg_star, 'CMF': CMF_star}).to_csv(planet + '_star_' + filename + '.csv')
    
    ax1.scatter(Mass_planet, Radius_planet, c = 'c', s = 5, alpha = 0.5)
    ax1.set_ylabel(r'Radius (R$_\oplus$)', size=20)
    ax1.set_xlabel('Mass (M$_\oplus$)', size=20)
    x = np.linspace(M-5*M_err, M+5*M_err, 10)
    femg_star_1sigl, femg_star_med, femg_star_1sigu = np.quantile(FeMg_star, [(1.0-0.68)/2.0, 0.5, 1.0 - (1.0-0.68)/2.0])
    ymed = np.zeros(len(x))
    ylow = np.zeros(len(x))
    yup = np.zeros(len(x))
    for i in range(0, len(x)):
        med = epbb.run_exoplex_mass(x[i], use_Earth_molar_ratios = False, FeMg = femg_star_med, SiMg = SiMg)
        ymed[i] = med['radius'][-1]/6371000.0
        low = epbb.run_exoplex_mass(x[i], use_Earth_molar_ratios = False, FeMg = femg_star_1sigl, SiMg = SiMg)
        ylow[i] = low['radius'][-1]/6371000.0
        up = epbb.run_exoplex_mass(x[i], use_Earth_molar_ratios = False, FeMg = femg_star_1sigu, SiMg = SiMg)
        yup[i] = up['radius'][-1]/6371000.0
    ax1.plot(x, ymed, 'k-', lw = 2, label = 'Median')
    ax1.plot(x, yup, 'k--', lw = 2, alpha = 0.5, label = r'1$\sigma$')
    ax1.plot(x, ylow, 'k--', lw = 2, alpha = 0.5, label = r'1$\sigma$')
    ax1.set_xlim(M-5*M_err, M+5*M_err)
    ax1.patch.set_linewidth(2)
    ax1.minorticks_on()
    ax1.tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, labelsize = 14)
    ax1.tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 5)
    ax1.patch.set_edgecolor('black')
    ax1.patch.set_linewidth(2)
    

    try:
        from labellines import labelLine, labelLines
        xv = (ax1.get_xlim()[1] + ax1.get_xlim()[0])/2.0
        labelLines(ax1.get_lines(), xvals = (xv, 0.8*xv, 1.2*xv))
    except:
        pass
    
    
    low = min(np.append(FeMg, FeMg_star))
    up = max(np.append(FeMg, FeMg_star))
    bins = np.linspace(low, up, int((up-low)/0.1) + 1)
    ax2.hist(FeMg, bins=bins, density=True, alpha=0.4, color='c', label = 'Planet')    #plt.show()
    ax2.hist(FeMg_star, bins=bins, density=True, alpha=0.4, color='k', label = 'Star')
    xmin, xmax = ax2.get_xlim()
    x = np.linspace(xmin, xmax, 1000)
    p = sp.lognorm.pdf(x,shape, loc = loc, scale = scale)
    ax2.plot([femg_rho_x1sig_low, femg_rho_x1sig_low], [0, femg_rho_y1sig_low], 'c--', lw = 2)
    ax2.plot([femg_rho_x1sig_up, femg_rho_x1sig_up], [0, femg_rho_y1sig_up], 'c--', lw = 2)
    ax2.plot([femg_rho_xmed, femg_rho_xmed], [0, femg_rho_ymed], 'c-', lw = 2)
    ax2.plot(x, p, 'c', linewidth=2)
    pstar = sp.lognorm.pdf(x, shape_star, loc = loc_star, scale = scale_star)
    ax2.plot([femg_star_x1sig_low, femg_star_x1sig_low], [0, femg_star_y1sig_low], 'k--', lw = 2)
    ax2.plot([femg_star_x1sig_up, femg_star_x1sig_up], [0, femg_star_y1sig_up], 'k--', lw = 2)
    ax2.plot([femg_star_xmed, femg_star_xmed], [0, femg_star_ymed], 'k-', lw = 2)
    ax2.plot(x, pstar, 'k', linewidth = 2)
    ax2.set_xlabel('Fe/Mg', size=20)
    if np.quantile(max(np.append(FeMg, FeMg_star)), [0.95])[0] > 5:
        ax2.set_xlim(0, 5)
    else:
        ax2.set_xlim(0, max(np.append(FeMg, FeMg_star)))
    ax2.legend(fontsize = 14, frameon = False)
    ax2.tick_params(which = 'both', direction = 'in', top = True, right = True)
    ax2.set_ylabel('Probability Density', fontsize = 20)
    ax2.patch.set_edgecolor('black')
    ax2.patch.set_linewidth(2)
    ax2.minorticks_on()
    ax2.tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, labelsize = 14)
    ax2.tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 5)

    low = 0; up = 1
    bins = np.linspace(low, up, int((up-low)/0.05) + 1)
    ax3.hist(CMF, bins=bins, density=True, alpha=0.4, color='c', label = 'Planet')    #plt.show()
    ax3.hist(CMF_star, bins=bins, density=True, alpha=0.4, color='k', label = 'Star')    #plt.show()
    xmin, xmax = ax3.get_xlim()
    x = np.linspace(xmin, xmax, 1000)
    p = sp.lognorm.pdf(x, shape_cmf_lognorm, loc = loc_cmf_lognorm, scale = scale_cmf_lognorm)
    ax3.plot([cmfrho_x1sig_low, cmfrho_x1sig_low], [0, cmfrho_y1sig_low], 'c--', lw = 2)
    ax3.plot([cmfrho_x1sig_up, cmfrho_x1sig_up], [0, cmfrho_y1sig_up], 'c--', lw = 2)
    ax3.plot([cmfrho_xmed, cmfrho_xmed], [0, cmfrho_ymed], 'c-', lw = 2)
    ax3.plot(x, p, 'c', linewidth=2)
    pstar = sp.lognorm.pdf(x, shape_cmf_star_lognorm, loc = loc_cmf_star_lognorm, scale = scale_cmf_star_lognorm)
    ax3.plot([cmfstar_x1sig_low, cmfstar_x1sig_low], [0, cmfstar_y1sig_low], 'k--', lw = 2)
    ax3.plot([cmfstar_x1sig_up, cmfstar_x1sig_up], [0, cmfstar_y1sig_up], 'k--', lw = 2)
    ax3.plot([cmfstar_xmed, cmfstar_xmed], [0, cmfstar_ymed], 'k-', lw = 2)
    ax3.set_ylim(0, 1.05*max(np.append(pstar, p)))
    ax3.plot(x, pstar, 'k', linewidth = 2)
    ax3.set_xlabel('CMF', size=20)
    ax3.legend(fontsize = 14, frameon = False)
    ax3.tick_params(which = 'both', direction = 'in', top = True, right = True)
    ax3.set_xlim(0,1)
    ax3.set_ylabel('Probability Density', fontsize = 20)
    ax3.patch.set_edgecolor('black')
    ax3.patch.set_linewidth(2)
    ax3.minorticks_on()
    ax3.tick_params(which = 'major', direction = 'in', top = True, right = True, length = 10, width = 2, labelsize = 14)
    ax3.tick_params(which = 'minor', direction = 'in', top = True, right = True, length = 5)
    plt.tight_layout()

    plt.savefig(planet+'.jpg', format = 'jpg')



