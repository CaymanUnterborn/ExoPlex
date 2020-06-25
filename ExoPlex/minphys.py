
import numpy as np
import ExoPlex.burnman as burnman
from scipy import interpolate
import math
from scipy.integrate import odeint
import sys
import ExoPlex.functions


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
    Planet: dictionary
        Dictionary of pressure, temperature, expansivity, specific heat and phases for modeled planet

    """
    Pressure_layers = Planet.get('pressure')
    Temperature_layers = Planet.get('temperature')
    rho_layers = Planet.get('density')
    num_mantle_layers, num_core_layers, number_h2o_layers = layers

    num_layers = num_mantle_layers+num_core_layers+number_h2o_layers

    if num_mantle_layers > 0:
        if num_core_layers > 0:
            for i in range(num_core_layers):
                if i <= num_core_layers:
                    rho_layers[i] = get_core_rho(Pressure_layers[i],Temperature_layers[i],Core_wt_per)

            if number_h2o_layers >0:
                P_points_water = Pressure_layers[(num_mantle_layers+num_core_layers):]
                T_points_water = Temperature_layers[(num_mantle_layers+num_core_layers):]
                water_rho = get_water_rho(P_points_water, T_points_water)
                rho_layers[num_core_layers+num_mantle_layers:] = water_rho

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
                        print (to_switch_P[i] / 1e4, to_switch_T[i])

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
                        print (to_switch_P[i]/10000, to_switch_T[i])
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
                        print (to_switch_P[i] / 1e4, to_switch_T[i])

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
                        print (to_switch_P[i] / 10000, to_switch_T[i])
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
            if i <= num_core_layers:
                rho_layers[i] = get_core_rho(Pressure_layers[i],Temperature_layers[i],Core_wt_per)
        return rho_layers
def get_water_rho(Pressure,Temperature):
    phase = []
    P_water = []
    T_water = []
    P_ice = []
    T_ice = []
    density_ice = []

    for i in range(len(Pressure)):
        if (functions.find_water_phase(Pressure[i],Temperature[i]) == 'Water'):
            P_water.append(Pressure[i]*(1.e9/10000.)-((73.5e-5)*2.06e9*(Temperature[i]-373.)))
            T_water.append(Temperature[i])

        else:
            phase.append(functions.find_water_phase(Pressure[i],Temperature[i]))
            P_ice.append(Pressure[i])
            T_ice.append(Temperature[i])


    for i in range(len(phase)):

        if phase[i] == 'Ice_VII':
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

            density_ice.append(rock.evaluate(['density'], 1.e9 * (P_ice[i] / 10000.), T_ice[i])[0])

        if phase[i] == 'Ice_VI':
            class Ice_VI(burnman.Mineral):
                def __init__(self):
                    self.params = {
                        'name': 'ice_VI',
                        'equation_of_state': 'bm2',
                        'V_0': 14.17e-6,
                        'K_0': 14.01e9,
                        'Kprime_0': 4.,
                        'molar_mass': 0.01801528,
                    }
                    burnman.Mineral.__init__(self)

            rock = Ice_VI()
            density_ice.append(rock.evaluate(['density'], 1.e9 * (P_ice[i] / 10000.), T_ice[i])[0])

        if phase[i] == 'Ice_Ih':
            class Ice_Ih(burnman.Mineral):

                def __init__(self):
                    self.params = {
                        'name': 'Ice_Ih',
                        'equation_of_state': 'bm3',
                        'V_0': 1./(916.72/0.01801528),
                        'K_0': 9.2e9,
                        'Kprime_0': 5.5,
                        'molar_mass': 0.01801528,
                    }
                    burnman.Mineral.__init__(self)

            rock = Ice_Ih()
            print ("uh oh")
            density_ice.append(rock.evaluate(['density'], 1.e9 * (P_ice[i] / 10000.), T_ice[i])[0])

    rock = burnman.minerals.other.water()
    density_water = rock.evaluate(['density'],P_water,T_water)[0]

    density = np.concatenate((density_ice,density_water),axis=0)
    test = np.concatenate((density_ice,P_ice),axis=0)

    return density

def get_water_Cp(Pressure, Temperature):
    phase = functions.find_water_phase(Pressure,Temperature)
    Pressure = Pressure/10000.
    if phase == 'Water':
        return 4.184e3

    if phase == 'Ice_VII' or phase=='Ice_VI':
        cp = 3.3 + 22.1 * np.exp(-0.058 * Pressure)  # Asahara 2010
        return 1000.*cp

    if phase == 'Ice_Ih':
        return 4.184e3


def get_water_alpha(Pressure,Temperature):
    phase = functions.find_water_phase(Pressure, Temperature)
    Pressure = Pressure / 10000.
    if phase == 'Water':
        return 214.e-6

    if phase == 'Ice_VII' or phase=='Ice_VI':
        Ks = 23.7
        Ksp = 4.15
        a0 = -3.9e-4
        a1 = 1.5e-6
        at = a0 + a1 * Temperature
        alpha = at * (1 + (Ksp / Ks) * Pressure) ** (-0.9)  # fei 1993
        return alpha

    if phase == 'Ice_Ih':
        return 214.e-6

    return 0

def get_core_rho(Pressure,Temperature,Core_wt_per):
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

def get_core_speeds(Pressure,Temperature,Core_wt_per):
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

    Pressure = [i*((1.e9)/10000.) for i in Pressure]
    Temperature = [i for i in Temperature]

    rock = iron()
    speeds = rock.evaluate(['v_phi', 'v_p', 'v_s'], Pressure, Temperature)

    return speeds


def get_gravity(Planet,layers):
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
            radii_mantle = radii[num_core_layers:(num_core_layers + num_mantle_layers)]
            density_mantle = density[num_core_layers:(num_core_layers + num_mantle_layers)]
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
    radii = Planet.get('radius')
    density = Planet.get('density')
    num_mantle_layers, num_core_layers, number_h2o_layers = layers

    if num_mantle_layers >0:
        if num_core_layers >0:
            if len(radii)<=num_core_layers:

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
    mass = Planet['mass']
    density = Planet['density']
    radius = np.zeros(len(mass))

    for i in range(len(radius))[1:]:
        Vol = (mass[i]-mass[i-1])/density[i]
        val = 4.*np.pi/3.
        radius[i] = pow((Vol/val)+pow(radius[i-1],3.),1./3.)

    return radius

def get_temperature(Planet,grids,structural_parameters,layers):
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

            for i in range(len(P_points_water)):
                spec_heat_water.append(get_water_Cp(P_points_water[i],T_points_water[i]))
                alpha_water.append(get_water_alpha(P_points_water[i],T_points_water[i]))

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

    delta = ([(1.-(old_rho[i]/new_rho[i])) for i in range(len(new_rho))[1:]])
    new_rho = [i for i in new_rho]

    for i in range(len(delta)):

        if delta[i] >= 1.e-6:
            return False,new_rho

    return True, new_rho

