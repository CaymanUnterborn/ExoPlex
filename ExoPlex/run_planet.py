# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.

import os
import sys
import ExoPlex.minphys as minphys
import ExoPlex.functions as functions
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))


import ExoPlex.run_perplex as run_perplex


def run_planet_radius(radius_planet, compositional_params, structure_params, layers,filename, verbose):
    """
   This module creates the Planet dictionary for a planet of defined radius R

    Parameters
    ----------
    radius_planet: float
        input radius of planet in Earth radii

    compositional_params: list
        Structural parameters of the planet; See example for description

    structural_params: list
        Structural parameters of the planet; See example for description

    layers: list
        Number of layers for core, mantle and water

    filename: string
       chosen filename for output file

    verbose: boolean
        enable verbose mode
    Returns
    -------
    Planet: dictionary
        Dictionary of final pressure, temperature, expansivity, specific heat and phases for modeled planet
        keys = 'radius','density','temperature','gravity','pressure', 'alpha','cp','Vphi''Vp','Vs','K'
    """

    Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac = functions.get_percents(compositional_params,verbose)
    get_phases = compositional_params[9]

    use_grids = compositional_params[-1]
    #Run fine mesh grid
    Mantle_filename = run_perplex.run_perplex(*[Mantle_wt_per,compositional_params,[structure_params[0],structure_params[1],structure_params[2]],filename,verbose,get_phases])
    grids_low, names = functions.make_mantle_grid(Mantle_filename,True,use_grids)
    names.append('Fe')
    if layers[-1] > 0:
        names.append('water')
        names.append('ice_Ih')
        names.append('ice_II')
        names.append('ice_III')
        names.append('ice_V')
        names.append('ice_VI')
        names.append('ice_VII')

    Mantle_filename = run_perplex.run_perplex(*[Mantle_wt_per,compositional_params,[structure_params[3],structure_params[4],structure_params[5]],filename,verbose,False])
    grids_high = functions.make_mantle_grid(Mantle_filename,False,use_grids)[0]

    grids = [grids_low,grids_high]
    Planet = functions.find_Planet_radius(radius_planet, core_mass_frac,structure_params, compositional_params, grids, Core_wt_per, layers,verbose)

    Planet['mass'] = minphys.get_mass(Planet,layers)
    Planet['phase_names'] = names
    Planet['phases'], Planet['phase_names'] = functions.get_phases(Planet, grids, layers, get_phases)

    return Planet

def run_planet_mass(mass_planet, compositional_params, structure_params, layers,filename,verbose):
    """
   This module creates the Planet dictionary for a planet of defined mass M

    Parameters
    ----------
    mass_planet: float
        input mass of planet in Earth masses

    compositional_params: list
        Structural parameters of the planet; See example for description

    structural_params: list
        Structural parameters of the planet; See example for description

    layers: list
        Number of layers for core, mantle and water

    filename: string
       chosen filename for output file

    verbose: boolean
        enable verbose mode
    Returns
    -------
    Planet: dictionary
        Dictionary of final pressure, temperature, expansivity, specific heat and phases for modeled planet
        keys = 'radius','density','temperature','gravity','pressure', 'alpha','cp','Vphi''Vp','Vs','K'
    """
    Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac = functions.get_percents(compositional_params,verbose)

    #Run fine mesh grid
    use_grids = compositional_params[10]
    Mantle_filename = run_perplex.run_perplex(*[Mantle_wt_per,compositional_params,[structure_params[0],structure_params[1],structure_params[2]],filename,verbose,True])
    grids_low, names = functions.make_mantle_grid(Mantle_filename,True,use_grids)
    names.append('Fe')
    if layers[-1] > 0:
        names.append('water')
        names.append('ice_Ih')
        names.append('ice_II')
        names.append('ice_III')
        names.append('ice_V')
        names.append('ice_VI')
        names.append('ice_VII')


    Mantle_filename = run_perplex.run_perplex(*[Mantle_wt_per,compositional_params,[structure_params[3],structure_params[4],structure_params[5]],filename,verbose,False])
    grids_high = functions.make_mantle_grid(Mantle_filename,False,use_grids)[0]

    grids = [grids_low,grids_high]
    Planet = functions.find_Planet_mass(mass_planet, core_mass_frac,structure_params, compositional_params, grids, Core_wt_per, layers,verbose)
    Planet['phase_names'] = names

    get_phases = compositional_params[9]
    Planet['phases'],Planet['phase_names'] = functions.get_phases(Planet, grids, layers,get_phases)



    return Planet

