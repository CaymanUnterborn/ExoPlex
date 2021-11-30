
import numpy as np
import ExoPlex.burnman as burnman
from scipy import interpolate
import math
from scipy.integrate import odeint
import sys
import seafreeze as sf
import ExoPlex.functions as functions


ToPa = 100000.
ToBar = 1./ToPa
G = 6.67408e-11

def get_rho(Planet,grids,Core_wt_per,layers):
    """
   This module stitches together the lower and upper mantle grids and interpolates within them to determine the density
   of each material in each layer in the planet. It also calculates the density of the core for those pressures and temperatures
   within the core layers.

    Parameters
    ----------
    Planet: dictionary
        Dictionary of pressure, temperature, expansivity, specific heat and phases for modeled planet
    grids: list of lists
        UM and LM grids containing pressure, temperature, density, expansivity, specific heat and phases
    Core_wt_per: float
        Composition of the Core
    layers: list
        Number of layers for core, mantle and water

    Returns
    -------
    rho_layers: list
        list of densities for water, mantle and core layers [kg/m^3']

    """
    Pressure_layers = Planet.get('pressure')
    Temperature_layers = Planet.get('temperature')
    rho_layers = Planet.get('density')
    num_mantle_layers, num_core_layers, number_h2o_layers = layers

    if num_mantle_layers > 0:
        if num_core_layers > 0:
            for i in range(num_core_layers):
                if i < num_core_layers:
                    rho_layers[i] = get_core_rho(Pressure_layers[i],Temperature_layers[i],Core_wt_per)
            if number_h2o_layers > 0:
                P_points_water = Pressure_layers[(num_mantle_layers+num_core_layers-1):]
                T_points_water = Temperature_layers[(num_mantle_layers+num_core_layers-1):]
                water_rho = get_water_rho(P_points_water, T_points_water)
                rho_layers[num_core_layers+num_mantle_layers-1:] = water_rho

            P_points_UM = []
            T_points_UM = []
            P_points_LM = []
            T_points_LM = []


            for i in range(num_mantle_layers):
                if Pressure_layers[i+num_core_layers] >=1150000:
                    P_points_LM.append(Pressure_layers[i+num_core_layers])
                    T_points_LM.append(Temperature_layers[i+num_core_layers])
                else:
                    P_points_UM.append(Pressure_layers[i+num_core_layers])
                    T_points_UM.append(Temperature_layers[i+num_core_layers])


            P_points_LM = np.array(P_points_LM)
            T_points_LM = np.array(T_points_LM)
            P_points_UM = np.array(P_points_UM)
            T_points_UM = np.array(T_points_UM)

            UM_data = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                                   grids[0]['density'],(P_points_UM, T_points_UM), method='linear')


            LM_data = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                                   grids[1]['density'],(P_points_LM, T_points_LM), method='linear')

            to_switch_P = []
            to_switch_T = []
            to_switch_index = []

            for i in range(len(UM_data)):
                if np.isnan(UM_data[i]) == True:
                    #try lower mantle:
                    to_switch_P.append(P_points_UM[i])
                    to_switch_T.append(T_points_UM[i])
                    to_switch_index.append(i)

            if len(to_switch_P) > 0:
                test = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                            grids[1]['density'], (to_switch_P, to_switch_T), method='linear')

                for i in range(len(test)):

                    if np.isnan(test[i]) == True:
                        print ("UM Rho Outside of range! ")
                        sys.exit()
                    else:
                        UM_data[to_switch_index[i]] = test[i]

            to_switch_P = []
            to_switch_T = []
            to_switch_index = []

            for i in range(len(LM_data)):
                if np.isnan(LM_data[i]) == True:
                    # try upper mantle:
                    to_switch_P.append(P_points_LM[i])
                    to_switch_T.append(T_points_LM[i])
                    to_switch_index.append(i)

            if len(to_switch_P) > 0:
                test = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                            grids[0]['density'], (to_switch_P, to_switch_T), method='linear')
                for i in range(len(test)):
                    if np.isnan(test[i]) == True:
                        print ("LM Rho Outside of range!")
                        sys.exit()
                    else:
                        LM_data[to_switch_index[i]] = test[i]

            mantle_data = np.append(LM_data,UM_data)

            rho_layers[num_core_layers:num_core_layers+num_mantle_layers] = mantle_data

            return rho_layers
        else:
            P_points_UM = []
            T_points_UM = []
            P_points_LM = []
            T_points_LM = []

            for i in range(num_mantle_layers):
                if Pressure_layers[i + num_core_layers] >= 1150000:
                    P_points_LM.append(Pressure_layers[i + num_core_layers])
                    T_points_LM.append(Temperature_layers[i + num_core_layers])
                else:
                    P_points_UM.append(Pressure_layers[i + num_core_layers])
                    T_points_UM.append(Temperature_layers[i + num_core_layers])

            P_points_LM = np.array(P_points_LM)
            T_points_LM = np.array(T_points_LM)
            P_points_UM = np.array(P_points_UM)
            T_points_UM = np.array(T_points_UM)


            UM_data = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                           grids[0]['density'], (P_points_UM, T_points_UM), method='linear')

            LM_data = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                           grids[1]['density'], (P_points_LM, T_points_LM), method='linear')

            to_switch_P = []
            to_switch_T = []
            to_switch_index = []

            for i in range(len(UM_data)):
                if np.isnan(UM_data[i]) == True:
                    # try lower mantle:
                    to_switch_P.append(P_points_UM[i])
                    to_switch_T.append(T_points_UM[i])
                    to_switch_index.append(i)

            if len(to_switch_P) > 0:
                test = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                            grids[1]['density'], (to_switch_P, to_switch_T), method='linear')

                for i in range(len(test)):

                    if np.isnan(test[i]) == True:
                        print ("UM Rho Outside of range!")
                        sys.exit()
                    else:
                        UM_data[to_switch_index[i]] = test[i]

            to_switch_P = []
            to_switch_T = []
            to_switch_index = []

            for i in range(len(LM_data)):
                if np.isnan(LM_data[i]) == True:
                    # try upper mantle:
                    to_switch_P.append(P_points_LM[i])
                    to_switch_T.append(T_points_LM[i])
                    to_switch_index.append(i)

            if len(to_switch_P) > 0:
                test = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                            grids[0]['density'], (to_switch_P, to_switch_T), method='linear')
                for i in range(len(test)):
                    if np.isnan(test[i]) == True:
                        print ("LM Rho Outside of range!")
                        sys.exit()
                    else:
                        LM_data[to_switch_index[i]] = test[i]

            mantle_data = np.append(LM_data, UM_data)
            return mantle_data
    elif number_h2o_layers > 0:
        P_points_water = Pressure_layers
        T_points_water = Temperature_layers
        water_rho = get_water_rho(P_points_water, T_points_water)
        rho_layers[num_core_layers + num_mantle_layers:] = water_rho

        return rho_layers
    else:
        for i in range(num_core_layers):
            if i < num_core_layers:
                rho_layers[i] = get_core_rho(Pressure_layers[i],Temperature_layers[i],Core_wt_per)
        return rho_layers
def get_water_rho(Pressure,Temperature):
    """
   This module calculates the density of water (either as an ice or liquid) given a list of pressures and temperatures

    Parameters
    ----------
    pressure: list
        List of pressures to evaluate for density [bar]
    temperature: list
        List of temperatures to evaluate for density [bar]

    Returns
    -------
    density: list
        list of calculated density of water [kg/m^3]

    """
    pressure_sf = []
    temperature_sf = []
    phase_sf = []
    density_sf=[]
    density_other =[]
    PT = np.empty((len(Pressure),), dtype=object)
    for i in range(len(Pressure)):
        PT[i] = Pressure[i]*ToPa/1.e6, Temperature[i]

    phase_num = sf.whichphase(PT)


    for i in range(len(phase_num)):
        if (phase_num[i] < 7):

            pressure_sf.append(Pressure[i]*ToPa/1.e6)
            temperature_sf.append(Temperature[i])
            phase_sf.append(phase_num[i])

        else:
        #if (math.isnan(phase_num[i])== True or sf.whichphase(PT)[0] == 7.0):

            class Ice_VII(burnman.Mineral):

                def __init__(self):
                    self.params = {
                        'name': 'ice_VII',
                        'equation_of_state': 'bm2',
                        'V_0':12.49e-6,
                        'K_0': 20.15e9,
                        'Kprime_0':4.,
                        'molar_mass': 0.01801528,
                    }
                    burnman.Mineral.__init__(self)

            rock = Ice_VII()
            density_other.append(rock.evaluate(['density'], Pressure[i]*ToPa, Temperature[i])[0])


    PT = np.empty((len(pressure_sf),), np.object)
    for i in range(len(pressure_sf)):
        PT[i] = pressure_sf[i], temperature_sf[i]

    phase_names = [sf.phasenum2phase[pn] for pn in phase_sf]

    for i in range(len(phase_sf)):
        new = sf.seafreeze(np.array(PT[i]), phase_names[i])
        density_sf.append(new.rho[0][0])

    density = np.concatenate((density_other,density_sf),axis=0)

    return density

def get_water_Cp(Pressure, Temperature):
    """
   This module calculates the specific heat at constant pressure of water (either as an ice or liquid) given a
   single pressure and temperature

    Parameters
    ----------
    pressure: float
        pressure to evaluate for Cp [bar]
    temperature: list
        temperature to evaluate for Cp [bar]

    Returns
    -------
    Cp: float
        Calculated specific heat [J/(kg*K)]

    """

    pressure_sf = []
    temperature_sf = []
    phase_sf = []
    cp_sf=[]
    cp_other =[]

    PT = np.empty((len(Pressure),), np.object)
    for i in range(len(Pressure)):
        PT[i] = Pressure[i]*ToPa/1.e6, Temperature[i]

    phase_num = sf.whichphase(PT)

    for i in range(len(phase_num)):
        if (phase_num[i] < 7):
            pressure_sf.append(Pressure[i]*ToPa/1.e6)
            temperature_sf.append(Temperature[i])
            phase_sf.append(phase_num[i])
        else:

            cp_other.append(1000* 3.3 + 22.1 * np.exp(-0.058 * Pressure[i]*ToPa/1e9) ) # Asahara 2010

    PT = np.empty((len(pressure_sf),), np.object)
    for i in range(len(pressure_sf)):
        PT[i] = pressure_sf[i], temperature_sf[i]

    phase_names = [sf.phasenum2phase[pn] for pn in phase_sf]
    for i in range(len(phase_sf)):
        new = sf.seafreeze(np.array(PT[i]), phase_names[i])
        cp_sf.append(new.Cp[0][0])

    Cp = np.concatenate((cp_other, cp_sf), axis=0)
    return Cp

def get_water_alpha(Pressure,Temperature):
    """
   This module calculates the thermal expansivity of water (either as an ice or liquid)
   given a single pressure and temperature

    Parameters
    ----------
    pressure: float
        pressure to evaluate for Cp [bar]
    temperature: list
        temperature to evaluate for Cp [bar]

    Returns
    -------
    alpha: float
        Calculated thermal expansivity [1/K]

    """

    pressure_sf = []
    temperature_sf = []
    phase_sf = []
    alpha_sf=[]
    alpha_other =[]

    PT = np.empty((len(Pressure),), np.object)
    for i in range(len(Pressure)):
        PT[i] = Pressure[i]*ToPa/1.e6, Temperature[i]

    phase_num = sf.whichphase(PT)

    for i in range(len(phase_num)):
        if (phase_num[i] < 7):
            pressure_sf.append(Pressure[i]*ToPa/1.e6)
            temperature_sf.append(Temperature[i])
            phase_sf.append(phase_num[i])
        else:
            Ks = 23.7
            Ksp = 4.15
            a0 = -3.9e-4
            a1 = 1.5e-6
            at = a0 + a1 * Temperature[i]
            alpha_other.append( at * (1 + (Ksp / Ks) * Pressure[i]) ** (-0.9))  # fei 1993)

    PT = np.empty((len(pressure_sf),), np.object)
    for i in range(len(pressure_sf)):
        PT[i] = pressure_sf[i], temperature_sf[i]

    phase_names = [sf.phasenum2phase[pn] for pn in phase_sf]
    for i in range(len(phase_sf)):
        new = sf.seafreeze(np.array(PT[i]), phase_names[i])
        alpha_sf.append(new.alpha[0][0])
    alpha = np.concatenate((alpha_other, alpha_sf), axis=0)
    return alpha




def get_core_rho(Pressure,Temperature,Core_wt_per):
    """
   This module calculates the density of the core given a list of pressures and temperatures and the light element
   composition of the core. See documentation for description of how light elements are treated.

    Parameters
    ----------
    pressure: list
        List of pressures to evaluate for density [bar]
    temperature: list
        List of temperatures to evaluate for density [bar]

    Returns
    -------
    density: list
        list of calculated density of the core [kg/m^3]

    """

    wt_frac_Si = Core_wt_per.get('Si')
    wt_frac_O = Core_wt_per.get('O')
    wt_frac_S = Core_wt_per.get('S')
    wt_frac_Fe = Core_wt_per.get('Fe')

    mFe = 55.845 #molar weights
    mSi = 28.0867
    mO = 15.9994
    mS = 32.0650

    mol_total = (wt_frac_Fe/mFe)+(wt_frac_O/mO)+(wt_frac_S/mS)+(wt_frac_Si/mSi)
    mol_frac_Fe = (wt_frac_Fe/mFe) / mol_total

    mol_frac_Si = (wt_frac_Si/mSi) / mol_total
    mol_frac_S = (wt_frac_S/mS) / mol_total
    mol_frac_O = (wt_frac_O/mO) / mol_total

    molar_weight_core = (mol_frac_Fe*mFe) + (mol_frac_Si * mSi) + (mol_frac_O*mO) + (mol_frac_S*mS)


    class iron(burnman.Mineral):

        #Anderson&Ahrens 1994
        def __init__(self):
            self.params = {
                'name': 'iron',
                'equation_of_state': 'bm4',
                'V_0': 7.95626e-6,
                'K_0': 109.7e9,
                'Kprime_0': 4.66,
                'Kprime_prime_0': -0.043e-9,
                'molar_mass': molar_weight_core/1000.,
            }
            burnman.Mineral.__init__(self)

    rock = iron()

    density = rock.evaluate(['density'], 1.e9*(Pressure/10000.), Temperature)[0]
    return density


def get_gravity(Planet,layers):
    """
    This module calculates the gravity profile of the planet
    using Gauss' Law of gravity given the density and radius lists.
    Parameters
    ----------
    Planet: dictionary
        Dictionary of pressure, temperature, expansivity, specific heat and phases for modeled planet
    layers: list
        Number of layers for core, mantle and water

    Returns
    -------
    gravity: list
        list of gravities in each shell for water, mantle and core layers [kg/m^2]

    """
    radii = Planet.get('radius')
    density = Planet.get('density')
    num_mantle_layers, num_core_layers, number_h2o_layers = layers


    if num_mantle_layers > 0:
        if num_core_layers >0:
            radii_core = radii[:num_core_layers]
            density_core = density[:num_core_layers]

            radii_mantle = radii[num_core_layers:(num_core_layers + num_mantle_layers)]
            density_mantle = density[num_core_layers:(num_core_layers + num_mantle_layers)]

            radii_water = radii[(num_core_layers + num_mantle_layers):]
            density_water = density[(num_core_layers + num_mantle_layers):]

            rhofunc_core = interpolate.UnivariateSpline(radii_core, density_core)
            rhofunc_mantle = interpolate.UnivariateSpline(radii_mantle, density_mantle)

            # Create a spline fit of density as a function of radius

            # Numerically integrate Poisson's equation
            poisson_core = lambda p, x: 4.0 * np.pi * G * rhofunc_core(x) * x * x
            poisson_mantle = lambda p, x: 4.0 * np.pi * G * rhofunc_mantle(x) * x * x

            gravity_layers_core = np.ravel(odeint(poisson_core, 0., radii_core))
            gravity_layers_mantle = np.ravel(odeint(poisson_mantle,gravity_layers_core[-1],radii_mantle))

            if number_h2o_layers>0:
                rhofunc_water = interpolate.UnivariateSpline(radii_water, density_water)
                poisson_water = lambda p, x: 4.0 * np.pi * G * rhofunc_water(x) * x * x
                gravity_layers_water = np.ravel(odeint(poisson_water,gravity_layers_mantle[-1],radii_water))
                gravity_layers = np.concatenate((gravity_layers_core,gravity_layers_mantle,gravity_layers_water),axis=0)

            else:
                gravity_layers = np.concatenate((gravity_layers_core, gravity_layers_mantle), axis=0)

            gravity_layers[1:] = gravity_layers[1:]/radii[1:]/radii[1:]
            gravity_layers[0] = 0

            return gravity_layers
        else:
            radii_mantle = radii[num_core_layers-1:(num_core_layers + num_mantle_layers)]
            density_mantle = density[num_core_layers-1:(num_core_layers + num_mantle_layers)]
            rhofunc_mantle = interpolate.UnivariateSpline(radii_mantle, density_mantle)
            poisson_mantle = lambda p, x: 4.0 * np.pi * G * rhofunc_mantle(x) * x * x
            gravity_layers = np.ravel(odeint(poisson_mantle,0.,radii_mantle))
            gravity_layers[1:] = gravity_layers[1:]/radii[1:]/radii[1:]
            gravity_layers[0] = 0
            return gravity_layers


    elif number_h2o_layers > 0:
        radii_water = radii[(num_core_layers + num_mantle_layers):]
        density_water = density[(num_core_layers + num_mantle_layers):]
        rhofunc_water = interpolate.UnivariateSpline(radii_water, density_water)
        poisson_water = lambda p, x: 4.0 * np.pi * G * rhofunc_water(x) * x * x
        gravity_layers_water = np.ravel(odeint(poisson_water, 0., radii_water))
        gravity_layers = np.array(gravity_layers_water)
        gravity_layers[1:] = gravity_layers[1:] / radii[1:] / radii[1:]
        gravity_layers[0] = 0
        return gravity_layers
    else:
        radii_core = radii
        density_core = density


        rhofunc_core = interpolate.UnivariateSpline(radii_core, density_core)
        poisson_core = lambda p, x: 4.0 * np.pi * G * rhofunc_core(x) * x * x
        gravity_layers = np.ravel(odeint(poisson_core, 0., radii_core))
        gravity_layers[1:] = gravity_layers[1:] / radii[1:] / radii[1:]
        gravity_layers[0] = 0

        return gravity_layers

def get_pressure(Planet,layers):
    """
        This module calculates the pressure profile of the planet
        using the equation for hydrostatic equilibrium given the density and radius lists.
        Parameters
        ----------
        Planet: dictionary
            Dictionary of pressure, temperature, expansivity, specific heat and phases for modeled planet
        layers: list
            Number of layers for core, mantle and water

        Returns
        -------
        pressure: list
            list of pressures in each shell for water, mantle and core layers [bar]

    """
    radii = Planet.get('radius')
    density = Planet.get('density')
    gravity = Planet.get('gravity')
    num_mantle_layers, num_core_layers, number_h2o_layers = layers


    # convert radii to depths
    depths = radii[-1] - radii
    depths_core = depths[:num_core_layers]
    gravity_core = gravity[:num_core_layers]
    density_core = density[:num_core_layers]

    depths_mant = depths[num_core_layers:(num_core_layers+num_mantle_layers)]
    gravity_mant = gravity[num_core_layers:(num_core_layers+num_mantle_layers)]
    density_mant = density[num_core_layers:(num_core_layers+num_mantle_layers)]


    if num_mantle_layers > 0:
        if num_core_layers > 0:
            # Make a spline fit of density as a function of depth
            rhofunc_mant = interpolate.UnivariateSpline(depths_mant[::-1], density_mant[::-1])
            # Make a spline fit of gravity as a function of depth
            gfunc_mant = interpolate.UnivariateSpline(depths_mant[::-1], gravity_mant[::-1])

            rhofunc_core = interpolate.UnivariateSpline(depths_core[::-1], density_core[::-1])
            gfunc_core = interpolate.UnivariateSpline(depths_core[::-1], gravity_core[::-1])

            if number_h2o_layers >0:
                depths_water = depths[(num_core_layers + num_mantle_layers):]
                gravity_water = gravity[(num_core_layers + num_mantle_layers):]
                density_water = density[(num_core_layers + num_mantle_layers):]

                rhofunc_water = interpolate.UnivariateSpline(depths_water[::-1], density_water[::-1])
                gfunc_water = interpolate.UnivariateSpline(depths_water[::-1], gravity_water[::-1])
                #integrate from 1 bar
                pressure_water = np.ravel(odeint((lambda p, x: gfunc_water(x) * rhofunc_water(x)),(1./10000.)*1.e9, depths_water[::-1]))
                WMB_pres = pressure_water[-1]


            else:
                WMB_pres = 3.e8

            # integrate the hydrostatic equation
            pressure_mant = np.ravel(odeint((lambda p, x: gfunc_mant(x) * rhofunc_mant(x)),WMB_pres, depths_mant[::-1]))
            CMB_pres = pressure_mant[-1]

            pressure_core = np.ravel(odeint((lambda p, x: gfunc_core(x) * rhofunc_core(x)),CMB_pres, depths_core[::-1]))

            if number_h2o_layers > 0:
                pressure = np.concatenate((pressure_water,pressure_mant,pressure_core),axis=0)
                pressure = [((i/1.e9)*10000.) for i in pressure]

                return np.array(pressure[::-1])
            else:
                pressure = np.concatenate((pressure_mant,pressure_core),axis=0)
                pressure = [((i/1.e9)*10000.) for i in pressure]
                return np.array(pressure[::-1])
        else:
            # Make a spline fit of density as a function of depth
            depths_mant = depths
            density_mant = density
            gravity_mant = gravity
            rhofunc_mant = interpolate.UnivariateSpline(depths_mant[::-1], density_mant[::-1])
            # Make a spline fit of gravity as a function of depth
            gfunc_mant = interpolate.UnivariateSpline(depths_mant[::-1], gravity_mant[::-1])


            # integrate the hydrostatic equation
            pressure_mant = np.ravel(
                odeint((lambda p, x: gfunc_mant(x) * rhofunc_mant(x)), 0.1e9, depths_mant[::-1]))

            pressure = np.array(pressure_mant)
            pressure = [((i / 1.e9) * 10000.) for i in pressure]
            return np.array(pressure[::-1])

    elif number_h2o_layers > 0:
        depths_water = depths[(num_core_layers + num_mantle_layers):]
        gravity_water = gravity[(num_core_layers + num_mantle_layers):]
        density_water = density[(num_core_layers + num_mantle_layers):]

        rhofunc_water = interpolate.UnivariateSpline(depths_water[::-1], density_water[::-1])
        gfunc_water = interpolate.UnivariateSpline(depths_water[::-1], gravity_water[::-1])
        # integrate from 1 bar
        pressure = np.ravel(
            odeint((lambda p, x: gfunc_water(x) * rhofunc_water(x)), (1. / 10000.) * 1.e9, depths_water[::-1]))
        pressure = [((i/1.e9)*10000.) for i in pressure]
        return np.array(pressure[::-1])
    else:
        rhofunc_core = interpolate.UnivariateSpline(depths_core[::-1], density_core[::-1])
        gfunc_core = interpolate.UnivariateSpline(depths_core[::-1], gravity_core[::-1])
        CMB_pres = 1e8

        pressure = np.ravel(odeint((lambda p, x: gfunc_core(x) * rhofunc_core(x)), CMB_pres, depths_core[::-1]))
        pressure = [((i/1.e9)*10000.) for i in pressure]

        return np.array(pressure[::-1])


def get_mass(Planet,layers):
    """
        This module calculates the mass of each individual shell profile of the planet
        using the equation for mass within a sphere given the density and radius lists.
        Parameters
        ----------
        Planet: dictionary
            Dictionary of pressure, temperature, expansivity, specific heat and phases for modeled planet
        layers: list
            Number of layers for core, mantle and water

        Returns
        -------
        mass: list
            list of masses of each shell for water, mantle and core layers [kg/m^2]

    """
    radii = Planet.get('radius')
    density = Planet.get('density')
    num_mantle_layers, num_core_layers, number_h2o_layers = layers

    if num_mantle_layers >0:
        if num_core_layers >0:
            if len(radii)< num_core_layers:

                rhofunc = interpolate.UnivariateSpline(radii, density)
                mass_in_sphere = lambda p, x: 4.0 * np.pi * rhofunc(x) * x * x

                mass = np.ravel(odeint(mass_in_sphere,0.,radii))
                return mass

            elif len(radii) > (num_core_layers+num_mantle_layers):

                radii_core = radii[:num_core_layers]
                density_core = density[:num_core_layers]
                rhofunc_core = interpolate.UnivariateSpline(radii_core, density_core)
                mass_in_core = lambda p, x: 4.0 * np.pi * rhofunc_core(x) * x * x

                mass_core = np.ravel(odeint(mass_in_core, 0., radii_core))

                radii_mantle = radii[num_core_layers:num_core_layers + num_mantle_layers]
                density_mantle = density[num_core_layers:num_core_layers + num_mantle_layers]
                rhofunc_mantle = interpolate.UnivariateSpline(radii_mantle, density_mantle)

                mass_in_mantle = lambda p, x: 4.0 * np.pi * rhofunc_mantle(x) * x * x

                mass_mantle = np.ravel(odeint(mass_in_mantle,mass_core[-1],radii_mantle))

                radii_water = radii[num_core_layers + num_mantle_layers:]
                density_water = density[num_core_layers + num_mantle_layers:]
                rhofunc_water = interpolate.UnivariateSpline(radii_water, density_water)

                mass_in_water = lambda p, x: 4.0 * np.pi * rhofunc_water(x) * x * x

                mass_water= np.ravel(odeint(mass_in_water, mass_mantle[-1], radii_water))

                return np.concatenate((mass_core,mass_mantle,mass_water),axis=0)
            else:
                radii_core = radii[:num_core_layers]
                density_core = density[:num_core_layers]
                rhofunc_core = interpolate.UnivariateSpline(radii_core, density_core)
                mass_in_core = lambda p, x: 4.0 * np.pi * rhofunc_core(x) * x * x

                mass_core = np.ravel(odeint(mass_in_core, 0., radii_core))
                radii_mantle = radii[num_core_layers:(num_core_layers + num_mantle_layers)]
                density_mantle = density[num_core_layers:(num_core_layers + num_mantle_layers)]
                rhofunc_mantle = interpolate.UnivariateSpline(radii_mantle, density_mantle)

                mass_in_mantle = lambda p, x: 4.0 * np.pi * rhofunc_mantle(x) * x * x

                mass_mantle = np.ravel(odeint(mass_in_mantle,mass_core[-1],radii_mantle))

                return np.concatenate((mass_core,mass_mantle),axis=0)
        else:
            radii_mantle = radii[num_core_layers:num_core_layers + num_mantle_layers]
            density_mantle = density[num_core_layers:num_core_layers + num_mantle_layers]
            rhofunc_mantle = interpolate.UnivariateSpline(radii_mantle, density_mantle)

            mass_in_mantle = lambda p, x: 4.0 * np.pi * rhofunc_mantle(x) * x * x

            mass_mantle = np.ravel(odeint(mass_in_mantle, 0., radii_mantle))

            return np.array(mass_mantle)

def get_radius(Planet,layers):
    """
        This module calculates the radius of each individual shell profile of the planet
        using the equation for density given the density and mass lists.
        Parameters
        ----------
        Planet: dictionary
            Dictionary of pressure, temperature, expansivity, specific heat and phases for modeled planet
        layers: list
            Number of layers for core, mantle and water

        Returns
        -------
        radius: list
            list of radius of each shell for water, mantle and core layers [kg/m^2]

    """
    mass = Planet['mass']
    density = Planet['density']
    radius = np.zeros(len(mass))

    val = 4. * np.pi / 3.
    #Not getting radius for fiinal
    radius[1] = pow(mass[1]/density[1]/val,1./3)
    for i in range(len(radius))[1:]:
        mass_lay = mass[i]-mass[i-1]
        radius[i] = np.cbrt(mass_lay/density[i]/val +pow(radius[i-1],3.))

    #print(radius[-1]/1000.)
    return radius

def get_temperature(Planet,grids,structural_parameters,layers):
    """
    This module calculates the adiabatic temperature profile of the planet
    using the equation for an adiabat given the density, gravity, pressure and radius lists as well as the
    grids of expansivity and specific heat at constant pressure.

    Parameters
    ----------
    Planet: dictionary
        Dictionary of pressure, temperature, expansivity, specific heat and phases for modeled planet

    grids: list of lists
        UM and LM grids containing pressure, temperature, density, expansivity, specific heat and phases

    structural_params: list
        Structural parameters of the planet; See example for description

    layers: list
        Number of layers for core, mantle and water

    Returns
    -------
    temperature: list
        list of temperatures in each shell for water, mantle and core layers [kg/m^2]

    """
    radii = Planet.get('radius')
    gravity = Planet.get('gravity')
    temperature = Planet.get('temperature')
    pressure = Planet.get('pressure')

    num_mantle_layers, num_core_layers, number_h2o_layers = layers



    if num_mantle_layers > 0:
        if num_core_layers > 0:
            Mantle_potential_temp = structural_parameters[7]
            Water_potential_temp = structural_parameters[9]

            radii = radii[num_core_layers:]
            gravity = gravity[num_core_layers:]

            pressure = pressure[num_core_layers:]
            temperature = temperature[num_core_layers:]

            depths = radii[-1] - radii

            P_points_UM = []
            T_points_UM = []
            P_points_LM = []
            T_points_LM = []

            P_points_water = pressure[num_mantle_layers:]
            T_points_water = temperature[num_mantle_layers:]

            spec_heat_water = []
            alpha_water = []

            if len(P_points_water) > 0:
                spec_heat_water = get_water_Cp(P_points_water,T_points_water)
                alpha_water= get_water_alpha(P_points_water,T_points_water)

            for i in range(num_mantle_layers):
                if pressure[i] >=1250000:
                    P_points_LM.append(pressure[i])
                    T_points_LM.append(temperature[i])
                else:
                    P_points_UM.append(pressure[i])
                    T_points_UM.append(temperature[i])

            depths_mantle = depths[:num_mantle_layers]
            gravity_mantle = gravity[:num_mantle_layers]



            depths_water = depths[num_mantle_layers:]
            gravity_water = gravity[num_mantle_layers:]

            UM_cp_data = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                                   grids[0]['cp'],(P_points_UM, T_points_UM), method='linear')

            LM_cp_data = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                                   grids[1]['cp'],(P_points_LM, T_points_LM), method='linear')

            to_switch_P = []
            to_switch_T = []
            to_switch_index = []
            for i in range(len(UM_cp_data)):
                if np.isnan(UM_cp_data[i]) == True:
                    # try lower mantle:
                    to_switch_P.append(P_points_UM[i])
                    to_switch_T.append(T_points_UM[i])
                    to_switch_index.append(i)

            if len(to_switch_P) > 0:
                test = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                            grids[1]['cp'], (to_switch_P, to_switch_T), method='linear')

                for i in range(len(test)):

                    if np.isnan(test[i]) == True:
                        print (to_switch_P[i] / 1e5, to_switch_T[i])
                        print ("UM Cp Outside of range!")
                        sys.exit()
                    else:
                        UM_cp_data[to_switch_index[i]] = test[i]

            to_switch_P = []
            to_switch_T = []
            to_switch_index = []

            for i in range(len(LM_cp_data)):
                if np.isnan(LM_cp_data[i]) == True:
                    # try lower mantle:
                    to_switch_P.append(P_points_LM[i])
                    to_switch_T.append(T_points_LM[i])
                    to_switch_index.append(i)

            if len(to_switch_P) > 0:
                test = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                            grids[0]['cp'], (to_switch_P, to_switch_T), method='linear')

                for i in range(len(test)):
                    if np.isnan(test[i]) == True:
                        print (to_switch_P[i], to_switch_T[i])
                        print ("LM Cp Outside of range!")
                        sys.exit()
                    else:
                        LM_cp_data[to_switch_index[i]] = test[i]

            spec_heat_mantle = np.append(LM_cp_data,UM_cp_data)

            UM_alpha_data = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                                   grids[0]['alpha'],(P_points_UM, T_points_UM), method='linear')

            LM_alpha_data = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                                   grids[1]['alpha'],(P_points_LM, T_points_LM), method='linear')

            to_switch_P = []
            to_switch_T = []
            to_switch_index = []
            for i in range(len(UM_alpha_data)):
                if np.isnan(UM_alpha_data[i]) == True:
                    # try lower mantle:
                    to_switch_P.append(P_points_UM[i])
                    to_switch_T.append(T_points_UM[i])
                    to_switch_index.append(i)

            if len(to_switch_P) > 0:
                test = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                            grids[1]['alpha'], (to_switch_P, to_switch_T), method='linear')

                for i in range(len(test)):

                    if np.isnan(test[i]) == True:
                        print ("UM Alpha Outside of range!")
                        sys.exit()
                    else:
                        UM_alpha_data[to_switch_index[i]] = test[i]

            to_switch_P = []
            to_switch_T = []
            to_switch_index = []
            for i in range(len(LM_alpha_data)):
                if np.isnan(LM_alpha_data[i]) == True:
                    # try lower mantle:
                    to_switch_P.append(P_points_LM[i])
                    to_switch_T.append(T_points_LM[i])
                    to_switch_index.append(i)

            if len(to_switch_P) > 0:
                test = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                            grids[0]['alpha'], (to_switch_P, to_switch_T), method='linear')

                for i in range(len(test)):

                    if np.isnan(test[i]) == True:
                        print (to_switch_P[i] / 1e5, to_switch_T[i])
                        print ("LM Alpha Outside of range!")
                        sys.exit()
                    else:
                        LM_alpha_data[to_switch_index[i]] = test[i]

            alpha_mantle = np.append(LM_alpha_data,UM_alpha_data)

            grav_func = interpolate.UnivariateSpline(depths_mantle[::-1],gravity_mantle[::-1])
            spec_heat_func = interpolate.UnivariateSpline(depths_mantle[::-1],spec_heat_mantle[::-1])
            alpha_func = interpolate.UnivariateSpline(depths_mantle[::-1],alpha_mantle[::-1])


            adiabat_mantle = lambda p,x:  alpha_func(x)*grav_func(x) / spec_heat_func(x)

            gradient_mantle = np.ravel(odeint(adiabat_mantle, 0.,depths_mantle[::-1]))

            if number_h2o_layers > 0:
                mantle_temperatures = [math.exp(k)*(Mantle_potential_temp+13.659*P_points_water[0]/1e4) for k in gradient_mantle][::-1]

            else:
                mantle_temperatures = [math.exp(k) * (Mantle_potential_temp) for k in
                                       gradient_mantle][::-1]
            core_temperatures = Planet['temperature'][:num_core_layers]

            if number_h2o_layers > 0:
                grav_func = interpolate.UnivariateSpline(depths_water[::-1], gravity_water[::-1])
                spec_heat_func = interpolate.UnivariateSpline(depths_water[::-1], spec_heat_water[::-1])
                alpha_func = interpolate.UnivariateSpline(depths_water[::-1], alpha_water[::-1])

                adiabat_water = lambda p, x: alpha_func(x) * grav_func(x) / spec_heat_func(x)

                gradient_water = np.ravel(odeint(adiabat_water, 0., depths_water[::-1]))
                if gradient_water[-1] > 1.:
                    gradient_water = [0 for i in gradient_water]
                water_temperatures = [math.exp(i)*Water_potential_temp for i in gradient_water][::-1]

                return np.concatenate((core_temperatures,mantle_temperatures,water_temperatures),axis=0)

            else:
                return np.concatenate((core_temperatures, mantle_temperatures),axis=0)
        else:
            Mantle_potential_temp = structural_parameters[7]
            depths = radii[-1] - radii

            P_points_UM = []
            T_points_UM = []
            P_points_LM = []
            T_points_LM = []

            for i in range(num_mantle_layers):
                if pressure[i] >= 60000:
                    P_points_LM.append(pressure[i])
                    T_points_LM.append(temperature[i])
                else:
                    P_points_UM.append(pressure[i])
                    T_points_UM.append(temperature[i])

            depths_mantle = depths[:num_mantle_layers]
            gravity_mantle = gravity[:num_mantle_layers]

            depths_water = depths[num_mantle_layers:]
            gravity_water = gravity[num_mantle_layers:]

            UM_cp_data = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                              grids[0]['cp'], (P_points_UM, T_points_UM), method='linear')

            LM_cp_data = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                              grids[1]['cp'], (P_points_LM, T_points_LM), method='linear')

            to_switch_P = []
            to_switch_T = []
            to_switch_index = []
            for i in range(len(UM_cp_data)):
                if np.isnan(UM_cp_data[i]) == True:
                    # try lower mantle:
                    to_switch_P.append(P_points_UM[i])
                    to_switch_T.append(T_points_UM[i])
                    to_switch_index.append(i)

            if len(to_switch_P) > 0:
                test = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                            grids[1]['cp'], (to_switch_P, to_switch_T), method='linear')

                for i in range(len(test)):

                    if np.isnan(test[i]) == True:
                        print (to_switch_P[i] / 1e5, to_switch_T[i])
                        print ("UM Cp Outside of range!")
                        sys.exit()
                    else:
                        UM_cp_data[to_switch_index[i]] = test[i]

            to_switch_P = []
            to_switch_T = []
            to_switch_index = []

            for i in range(len(LM_cp_data)):
                if np.isnan(LM_cp_data[i]) == True:
                    # try lower mantle:
                    to_switch_P.append(P_points_LM[i])
                    to_switch_T.append(T_points_LM[i])
                    to_switch_index.append(i)

            if len(to_switch_P) > 0:
                test = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                            grids[0]['cp'], (to_switch_P, to_switch_T), method='linear')

                for i in range(len(test)):
                    if np.isnan(test[i]) == True:
                        print (to_switch_P[i], to_switch_T[i])
                        print ("LM Cp Outside of range!")
                        sys.exit()
                    else:
                        LM_cp_data[to_switch_index[i]] = test[i]

            spec_heat_mantle = np.append(LM_cp_data, UM_cp_data)

            UM_alpha_data = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                                 grids[0]['alpha'], (P_points_UM, T_points_UM), method='linear')

            LM_alpha_data = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                                 grids[1]['alpha'], (P_points_LM, T_points_LM), method='linear')

            to_switch_P = []
            to_switch_T = []
            to_switch_index = []
            for i in range(len(UM_alpha_data)):
                if np.isnan(UM_alpha_data[i]) == True:
                    # try lower mantle:
                    to_switch_P.append(P_points_UM[i])
                    to_switch_T.append(T_points_UM[i])
                    to_switch_index.append(i)

            if len(to_switch_P) > 0:
                test = interpolate.griddata((grids[1]['pressure'], grids[1]['temperature']),
                                            grids[1]['alpha'], (to_switch_P, to_switch_T), method='linear')

                for i in range(len(test)):

                    if np.isnan(test[i]) == True:
                        print (to_switch_P[i] / 1e5, to_switch_T[i])
                        print ("UM Alpha Outside of range!")
                        sys.exit()
                    else:
                        UM_alpha_data[to_switch_index[i]] = test[i]

            to_switch_P = []
            to_switch_T = []
            to_switch_index = []
            for i in range(len(LM_alpha_data)):
                if np.isnan(LM_alpha_data[i]) == True:
                    # try lower mantle:
                    to_switch_P.append(P_points_LM[i])
                    to_switch_T.append(T_points_LM[i])
                    to_switch_index.append(i)

            if len(to_switch_P) > 0:
                test = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                            grids[0]['alpha'], (to_switch_P, to_switch_T), method='linear')

                for i in range(len(test)):

                    if np.isnan(test[i]) == True:
                        print (to_switch_P[i] / 1e5, to_switch_T[i])
                        print ("LM Alpha Outside of range!")
                        sys.exit()
                    else:
                        LM_alpha_data[to_switch_index[i]] = test[i]

            alpha_mantle = np.append(LM_alpha_data, UM_alpha_data)

            grav_func = interpolate.UnivariateSpline(depths_mantle[::-1], gravity_mantle[::-1])
            spec_heat_func = interpolate.UnivariateSpline(depths_mantle[::-1], spec_heat_mantle[::-1])
            alpha_func = interpolate.UnivariateSpline(depths_mantle[::-1], alpha_mantle[::-1])


            adiabat_mantle = lambda p, x: alpha_func(x) * grav_func(x) / spec_heat_func(x)

            gradient_mantle = np.ravel(odeint(adiabat_mantle, 0., depths_mantle[::-1]))

            for i in range(len(gradient_mantle)):
                if gradient_mantle[i] < 1e-6 and i != 0:
                    gradient_mantle[i] = gradient_mantle[i-1] + .01

            mantle_temperatures = [math.exp(k) * (Mantle_potential_temp) for k in gradient_mantle][::-1]

            return np.array(mantle_temperatures)
    elif number_h2o_layers > 0:
        """
        Water_potential_temp = structural_parameters[9]
        depths_water = radii[-1] - radii
        temperature = temperature
        gravity_water = gravity

        P_points_water = pressure[num_mantle_layers:]
        T_points_water = temperature[num_mantle_layers:]

        spec_heat_water = []
        alpha_water = []
        for i in range(len(P_points_water)):
            spec_heat_water.append(get_water_Cp(P_points_water[i],T_points_water[i]))
            alpha_water.append(get_water_alpha(P_points_water[i],T_points_water[i]))

        grav_func = interpolate.UnivariateSpline(depths_water[::-1], gravity_water[::-1])
        spec_heat_func = interpolate.UnivariateSpline(depths_water[::-1], spec_heat_water[::-1])
        alpha_func = interpolate.UnivariateSpline(depths_water[::-1], alpha_water[::-1])

        adiabat_water = lambda p, x: alpha_func(x) * grav_func(x) / spec_heat_func(x)
        gradient_water = np.ravel(odeint(adiabat_water, 0., depths_water[::-1]))

        if gradient_water[-1] > 1.:
            gradient_water = [0 for i in gradient_water]

    
        water_temperatures = [math.exp(i) * Water_potential_temp for i in gradient_water][::-1]
        """
        return Planet['temperature']

    else:
        temperatures = Planet['temperature']
        return temperatures
def check_convergence(new_rho,old_rho):
    """
         This module checks convergence by comparing the density layers of the previous and current runs

         Parameters
         ----------
         new_rho: list
             densities calculated in current iteration
         old_rho: list
             densities calculated in previous iteration
         Returns
         -------
         converged: boolean
            True if converged
         new_rho: list
            densities of current iteration
         """

    delta = ([(1.-(old_rho[i]/new_rho[i])) for i in range(len(new_rho))[1:]])
    new_rho = [i for i in new_rho]

    for i in range(len(delta)):
        if abs(delta[i]) >= pow(10,-3):
            return False,new_rho

    return True, new_rho

