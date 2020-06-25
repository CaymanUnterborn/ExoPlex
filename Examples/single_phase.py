
# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.

import os
import sys
import numpy as np
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))
Pressure_range_mantle_UM = '1000 2000000'
Temperature_range_mantle_UM = '1600 3500'

Pressure_range_mantle_LM = '750000 10000000'
Temperature_range_mantle_LM = '2500 6000'

core_rad_frac_guess = 1.
water_rad_frac_guess = 0.1
water_potential_temp = 300.

combine_phases = True
use_grids = False
import ExoPlex as exo

if __name__ == "__main__":

    Mass_planet = 1 # in Earth masses
    Radius_planet =1
    #Need to give the run a name. This will be used as the name of the output files
    Star = 'Iron'

    #Next user must input the ratios by mole (Earth is Ca/Mg = .07, Si.Mg = 0.90, Al/Mg = 0.09, Fe/Mg = 0.9)
    CaMg = 0.07
    SiMg = 0.9
    AlMg = 0.09
    FeMg = .9


    #How much water do you want in your planet? By mass fraction.
    wt_frac_water = 0

    #Don't forget that if you have water you need to add water layers
    number_h2o_layers = 0

    #Now we can mix various elements into the core or mantle
    wt_frac_Si_core = 0. #by mass <1
    wt_frac_O_core = 0. #by mass
    wt_frac_S_core = 0. #by mass
    mol_frac_Fe_mantle =0.0#by mole

    #What potential temperature (in K) do you want to start your mantle adiabat?
    Mantle_potential_temp = 1600.

    #Input the resolution of your upper mantle and lower mantle composition, density grids
    #These are input as number of T, P points. 50 50 = 2500 grid points, which takes about
    #5 minutes to calculate. Lower mantle resolution does not need to be higher since it's
    #mostly ppv.
    resolution_UM = '35 35'
    resolution_LM = '50 50'

    #lastly we need to decide how many layers to put in the planet. This is the resolution of
    #the mass-radius sampling.
    num_mantle_layers = 0
    num_core_layers = 600

    #create filename to store values
    Output_filename = 'iron_some'

    number_of_runs = 1
    Output_radii = []
    Output_mass = []


    ######### Initalize and run ExoPlex

    mass = [.05,.1,.2,.3,.4,.5,.6,.7,.8,.9]
    #1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7,7.5, 8,8.5,9, 9.5,10
    Output_filename = Output_filename + 'MR'
    for i in mass:
        compositional_params = [wt_frac_water,FeMg,SiMg,CaMg,AlMg,mol_frac_Fe_mantle,wt_frac_Si_core, \
                              wt_frac_O_core,wt_frac_S_core,combine_phases,use_grids]

        if use_grids == True:
            filename = exo.functions.find_filename(compositional_params)
        else:
            filename=''

        structure_params =  [Pressure_range_mantle_UM,Temperature_range_mantle_UM,resolution_UM,
                             Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                             core_rad_frac_guess,Mantle_potential_temp,water_rad_frac_guess,water_potential_temp]


        layers = [num_mantle_layers,num_core_layers,number_h2o_layers]

        #This is where we actually run the planet. First PerPlex grids of mineralogy, density,
        #Cp and alpha are calculated and stored in the Solutions folder. If the file already exists
        #(in name, not necessarily in composition), then PerPlex is not run again.



        Planet = exo.run_planet_singlephase(i,compositional_params,structure_params,layers,filename)

        #Planet is a dictionary containing many parameters of interest:
        #Planet.get('radius') = list of the radial points from calculation (m)
        #Planet.get('mass') = list of the cumulative mass at each radius point from calculation (kg)
        #Planet.get('density') = list of densities from calculation (kg/m^3)
        #Planet.get('temperature') = list of temperature points from calculation (K)
        #Planet.get('gravity') = list of gravity points from calculation (SI)
        #Planet.get('pressure') = list of pressure points from calculation (bar)
        #Planet.get('alpha') = list of values of thermal expansivity points from calculation (1/K)
        #Planet.get('cp') = list of values of specific heat points from calculation (SI)
        #Planet.get('phases') = list of phases and their molar fractions
        print
        print "Mass = ", '%.3f'%(Planet['mass'][-1]/5.97e24), "Earth masses"
        print "Radius = ", '%.3f'%(Planet['radius'][-1]/6371e3), "Earth radii"

        Output_mass.append([Planet['mass'][-1]/5.97e24,Planet['radius'][-1]/6371e3])
    #If you'd like the full output, uncomment out these lines!


        np.savetxt(Output_filename + '.dat', Output_mass, '%.5f', "\t", newline='\n',header='M\tR', footer='', comments='# ')



    #Now let us plot
    import matplotlib.pyplot as plt

    figure = plt.figure(figsize=(12, 10))

    ax1 = plt.subplot2grid((6, 3), (0, 0), colspan=3, rowspan=3)
    ax2 = plt.subplot2grid((6, 3), (3, 0), colspan=3, rowspan=1)
    ax3 = plt.subplot2grid((6, 3), (4, 0), colspan=3, rowspan=1)
    ax4 = plt.subplot2grid((6, 3), (5, 0), colspan=3, rowspan=1)

    ax1.plot(Planet['radius'] / 1.e3, Planet['density'] / 1.e3, 'k', linewidth=2.)
    ax1.set_ylim(0., (max(Planet['density']) / 1.e3) + 1.)
    ax1.set_xlim(0., max(Planet['radius']) / 1.e3)
    ax1.set_ylabel("Density ( $\cdot 10^3$ kg/m$^3$)")
    ax1.minorticks_on()
    text = '%.3f' % (Planet['radius'][-1] / 6371e3) + ' Earth radii on last iteration'
    ax1.text(0.05, 0.95, text)

    # Make a subplot showing the calculated pressure profile
    ax2.plot(Planet['radius'] / 1.e3, Planet['pressure'] / 1.e4, 'b', linewidth=2.)
    ax2.set_ylim(0., (max(Planet['pressure']) / 1e4) + 10.)
    ax2.set_xlim(0., max(Planet['radius']) / 1.e3)
    ax2.set_ylabel("Pressure (GPa)")
    ax2.minorticks_on()

    # Make a subplot showing the calculated gravity profile
    ax3.plot(Planet['radius'] / 1.e3, Planet['gravity'], 'r', linewidth=2.)
    ax3.set_ylabel("Gravity (m/s$^2)$")
    ax3.set_xlim(0., max(Planet['radius']) / 1.e3)
    ax3.set_ylim(0., max(Planet['gravity']) + 0.5)
    ax3.minorticks_on()

    # Make a subplot showing the calculated temperature profile
    ax4.plot(Planet['radius'] / 1.e3, Planet['temperature'], 'g', linewidth=2.)
    ax4.set_ylabel("Temperature ($K$)")
    ax4.set_xlabel("Radius (km)")
    ax4.set_xlim(0., max(Planet['radius']) / 1.e3)
    ax4.set_ylim(0., max(Planet['temperature']) + 100)
    ax4.minorticks_on()

    plt.show()
