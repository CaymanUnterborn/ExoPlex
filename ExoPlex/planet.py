import sys
import os
import numpy as np


if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))
Earth_radius = 6371e3
Earth_mass = 5.97e24
ToPa = 100000.
ToBar = 1./ToPa
from ExoPlex import minphys as minphys

def initialize_by_mass(*args):

    """
   This module creates the dictionary of lists for each planetary parameter (e.g., density) for a planet of the mass
   input by user.

    Parameters
    ----------
    mass_planet: float
        input radius of planet in Earth radii

    structural_params: list
        Structural parameters of the planet; See example for description

    compositional_params: list
        Structural parameters of the planet; See example for description

    layers: list
        Number of layers for core, mantle and water
    Returns
    -------
    Planet: dictionary
        Dictionary of initial guess of pressure, temperature, expansivity, specific heat and phases for modeled planet
    """
    mass_planet = args[0]
    structural_params = args[1]
    compositional_params = args[2]

    wt_frac_water = compositional_params.get('wt_frac_water')
    num_layers = sum(args[3])

    num_mantle_layers, num_core_layers, number_h2o_layers = args[3]

    if wt_frac_water > 0 and number_h2o_layers ==0:
        print("You have water but no water layers, please add")
        sys.exit()
    Pressure_layers = np.zeros(num_layers)
    Temperature_layers = np.zeros(num_layers)

    core_mass_frac = args[4]
    water_mass = (wt_frac_water*mass_planet)*Earth_mass
    core_mass = ((mass_planet*(1.-wt_frac_water))* core_mass_frac *Earth_mass)

    mantle_mass = (mass_planet*Earth_mass)-water_mass-core_mass
    Mantle_potential_temp = structural_params.get('Mantle_potential_temp')
    water_potential_temp = structural_params.get('water_potential_temp')
    if mass_planet > 5:
        Radius_planet_guess = 1.5
    elif mass_planet >8:
        Radius_planet_guess = 1.75
    else:
        Radius_planet_guess = 1.

    CMB_T_guess = (4180.*(Radius_planet_guess-0)-2764.*pow(Radius_planet_guess-0,2.)+1219.*pow(Radius_planet_guess-0,3.)) \
                  + (Mantle_potential_temp-1600)*(0.82+pow(Radius_planet_guess-0,1.81))
    if wt_frac_water > 0:
        WMB_pres = 1e5 #in bar
    else:
        WMB_pres = 1 #in bar
    if CMB_T_guess < Mantle_potential_temp:
        CMB_T_guess = Mantle_potential_temp + 1000
    if CMB_T_guess > 7000:
        CMB_T_guess = 6800
    CMB_P_guess = 10000*(262.*(Radius_planet_guess)-550.*pow(Radius_planet_guess,2.) + 432.*pow(Radius_planet_guess,3.))
    dP_dr = (CMB_P_guess)/(num_mantle_layers)
    dT_dr = (CMB_T_guess-Mantle_potential_temp)/(num_mantle_layers)

    for i in range(num_layers):

            if i<number_h2o_layers:
                Pressure_layers[i] = 1 + .5*WMB_pres*(i/number_h2o_layers)
                Temperature_layers[i] = water_potential_temp+200*(i/number_h2o_layers)

            elif i <= number_h2o_layers+num_mantle_layers-1:
                Pressure_layers[i] = WMB_pres + 1*dP_dr*(i-number_h2o_layers)
                Temperature_layers[i] = Mantle_potential_temp+1*dT_dr*(i-number_h2o_layers)

            else:
                Pressure_layers[i] = Pressure_layers[number_h2o_layers+num_mantle_layers-1] + ((i-number_h2o_layers-num_mantle_layers)/num_core_layers)*1e7
                Temperature_layers[i] = Temperature_layers[i-1]+2*(i-number_h2o_layers-num_mantle_layers)/num_core_layers

    val = core_mass/num_core_layers
    mass_layers_core = [val for k in range(num_core_layers)]
    val =mantle_mass / num_mantle_layers

    mass_layers_mantle = [val for k in range(num_mantle_layers)]


    if number_h2o_layers > 0:
        val = water_mass / number_h2o_layers
        mass_layers_water = [val for k in range(number_h2o_layers)]
        mass_layers = np.concatenate((mass_layers_core, mass_layers_mantle,mass_layers_water))
        mass_layers = [(sum(mass_layers[:i + 1])) for i in range(len(mass_layers))]


    else:
        mass_layers = np.concatenate((mass_layers_core,mass_layers_mantle))
        mass_layers = [(sum(mass_layers[:i + 1])) for i in range(len(mass_layers))]

    Pressure_layers = Pressure_layers[::-1]
    Temperature_layers = Temperature_layers[::-1]
    keys = ['mass', 'temperature', 'pressure']
    return dict(zip(keys, [mass_layers, Temperature_layers, Pressure_layers]))


def compress_mass(*args):
    """
   This module iterates the density, mass within a sphere, adiabatic temperature and gravity integrals for a planet of Mass M
   until convergence is reached. Convergence is defined as the change from the previous run to the current is the
   difference in the density of all layers is <1e-6.

    Parameters
    ----------
    Planet: dictionary
        Dictionary of initial guess of pressure, temperature, expansivity, specific heat and phases for modeled planet

     grids: list of lists
        UM and LM grids containing pressure, temperature, density, expansivity, specific heat and phases

    Core_wt_per: float
        Composition of the Core

    structural_params: list
        Structural parameters of the planet; See example for description

    layers: list
        Number of layers for core, mantle and water

    verbose: bool
    Returns
    -------
    Planet: dictionary
        Dictionary of initial guess of pressure, temperature, expansivity, specific heat and phases for modeled planet
    """
    Planet = args[0]
    grids = args[1]
    Core_wt_per = args[2]
    structural_params= args[3]
    compositional_params = args[4]
    core_mass_frac = args[5]
    layers= args[6]
    verbose = args[7]
    n_iterations = 1
    wt_frac_water = compositional_params.get('wt_frac_water')
    max_iterations = 100

    old_r = [10  for i in range(len(Planet['mass']))]
    converge = False

    while n_iterations <= max_iterations and converge == False:
        if verbose == True:
            print ("iteration #",n_iterations)
        if n_iterations>1:
            converge,old_r = minphys.check_convergence(Planet['density'],old_r)

        Planet['density'] = minphys.get_rho(Planet,grids,Core_wt_per,layers)
        Planet['radius'] = minphys.get_radius(Planet,wt_frac_water,core_mass_frac,layers)
        Planet['gravity'] = minphys.get_gravity(Planet,layers)
        Planet['pressure'] = minphys.get_pressure(Planet,layers)
        Planet['temperature'] = minphys.get_temperature(Planet, grids, structural_params, layers)
        n_iterations+=1
    return Planet


