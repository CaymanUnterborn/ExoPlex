
import numpy as np
from scipy import interpolate
import math
from scipy.integrate import odeint
import sys

ToPa = 100000
ToBar = 1/ToPa
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

    if number_h2o_layers > 0:
        P_points_water = Pressure_layers[(num_mantle_layers+num_core_layers):]
        T_points_water = Temperature_layers[(num_mantle_layers+num_core_layers):]
        water_rho = get_water_rho(P_points_water, T_points_water,grids)

    P_core = Pressure_layers[:num_core_layers]
    T_core = Temperature_layers[:num_core_layers]
    core_data = get_core_rho(grids[2], Core_wt_per, P_core, T_core)

    for i in range(num_core_layers):
        if i < num_core_layers:
            rho_layers[i] = core_data[i]

    if number_h2o_layers > 0:
        P_points_water = Pressure_layers[(num_mantle_layers + num_core_layers ):]
        T_points_water = Temperature_layers[(num_mantle_layers + num_core_layers):]
        water_rho =  get_water_rho(P_points_water, T_points_water,grids)

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
                print("Density: Pressure and/or Temperature Exceeds Mantle Grids ")

                print (to_switch_P[i]/10/1000, "GPa", to_switch_T[i] ,"K")
                print("Grid max:", max(grids[1]['pressure'])/10/1000, "GPa",  max(grids[1]['temperature']),"K")
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
                print("Density: Pressure and/or Temperature Exceeds Mantle Grids ")

                print (to_switch_P[i]/10/1000, "GPa", to_switch_T[i] ,"K")
                print("Grid max:", max(grids[1]['pressure'])/10/1000, "GPa",  max(grids[1]['temperature']),"K")
                sys.exit()
            else:
                LM_data[to_switch_index[i]] = test[i]

    mantle_data = np.append(LM_data, UM_data)

    if number_h2o_layers > 0:
        rho_layers = np.concatenate((core_data,mantle_data,water_rho))
    else:
        rho_layers= np.append(core_data,mantle_data)
    for i in range(sum(layers)-1)[::-1]:
        if rho_layers[i+1] > rho_layers[i]:
            drhodP = (rho_layers[i+1]-rho_layers[i+2])/(Pressure_layers[i+1]-Pressure_layers[i+2])
            rho_layers[i] =rho_layers[i+1] + (Pressure_layers[i]-Pressure_layers[i+1])* drhodP

    return rho_layers

def get_core_rho(grid,Core_wt_per,Pressure,Temperature):
    wt_frac_Si = Core_wt_per.get('Si')
    wt_frac_O = Core_wt_per.get('O')
    wt_frac_S = Core_wt_per.get('S')
    wt_frac_Fe = Core_wt_per.get('Fe')

    mFe = 55.845  # molar weights
    mSi = 28.0867
    mO = 15.9994
    mS = 32.0650

    mol_total = ((wt_frac_Fe / mFe) + (wt_frac_O / mO) + (wt_frac_S / mS) + (wt_frac_Si / mSi)) / 100
    mol_frac_Fe = (wt_frac_Fe / mFe / 100) / mol_total
    mol_frac_Si = (wt_frac_Si / mSi / 100) / mol_total
    mol_frac_S = (wt_frac_S / mS / 100) / mol_total
    mol_frac_O = (wt_frac_O / mO / 100) / mol_total

    molar_weight_core = (mol_frac_Fe * mFe) + (mol_frac_Si * mSi) + (mol_frac_O * mO) + (mol_frac_S * mS)

    core_rho = interpolate.griddata((grid['pressure'], grid['temperature']),
                                     grid['density'], (Pressure, Temperature), method='linear')


    core_rho = (molar_weight_core/mFe)*core_rho
    return core_rho

def get_water_rho(Pressure,Temperature,grids):
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
    from ExoPlex import burnman
    class Ice_VII(burnman.Mineral):

        def __init__(self):
            self.params = {
                'name': 'ice_VII',
                'equation_of_state': 'bm2',
                'V_0': 12.49e-6,
                'K_0': 20.15e9,
                'Kprime_0': 4.,
                'molar_mass': 0.01801528,
            }
            burnman.Mineral.__init__(self)
    rock = Ice_VII()
    Water_density = interpolate.griddata((grids[3]['pressure'], grids[3]['temperature']), grids[3]['density'],
    (Pressure, Temperature), method = 'linear')

    for i in range(len(Water_density)):
        if np.isnan(Water_density[i]) == True and Pressure[i] >= min(grids[3]['pressure']) and Temperature[i] >= min(grids[3]['temperature']):
            den_noT = rock.evaluate(['density'], Pressure[i]*ToPa, 300)[0]
            corr = np.exp((Temperature[i]-300)*11.58e-5)
            Water_density[i] = den_noT/corr

        elif np.isnan(Water_density[i]) == True:
            print("outside of water grid")
            print(Water_density[i])
            print(Pressure[i]/10/1000, "GPa",Temperature[i],"K")
            print(min((grids[3]['temperature'])),min(grids[3]['pressure']/10/1000))
            print(max((grids[3]['temperature'])),max(grids[3]['pressure']/10/1000))

            sys.exit()


    return Water_density

def get_water_Cp(Pressure, Temperature,grids):
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

    Water_Cp = interpolate.griddata((grids[3]['pressure'], grids[3]['temperature']), grids[3]['cp'],
    (Pressure, Temperature), method = 'linear')

    for i in range(len(Water_Cp)):
        if np.isnan(Water_Cp[i]) == True:
            Water_Cp[i] = 1000 * (3.3 + 22.1 * np.exp(-0.058 * Pressure[i]/10/1000)) #wants GPA

    return Water_Cp

def get_water_alpha(Pressure,Temperature,grids):
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
    Ks = 23.9
    Ksp = 4.15
    a0 = -3.9e-4
    a1 = 1.5e-6

    Water_alpha = interpolate.griddata((grids[3]['pressure'], grids[3]['temperature']), grids[3]['alpha'],
    (Pressure, Temperature), method = 'linear')

    for i in range(len(Water_alpha)):
        if np.isnan(Water_alpha[i]) == True:
            at = a0 + (a1 * Temperature[i])
            Water_alpha[i] = (at * pow(1 + ((Ksp / Ks) * Pressure[i]/10/1000), -0.9))  # Fei 94, pressure in Gpa

    return Water_alpha


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

    if number_h2o_layers>0:

        radii_core = radii[:num_core_layers]
        density_core = density[:num_core_layers]
        rhofunc_core = interpolate.InterpolatedUnivariateSpline(radii_core, density_core,k=3)
        poisson_core = lambda p, x: 4.0 * np.pi * G * rhofunc_core(x) * x * x
        gravity_layers_core = np.ravel(odeint(poisson_core, 0., radii_core))

        radii_rock = radii[num_core_layers:num_core_layers+num_mantle_layers]
        density_rock = density[num_core_layers:num_core_layers+num_mantle_layers]
        rhofunc_rock = interpolate.InterpolatedUnivariateSpline(radii_rock, density_rock,k=3)
        poisson_rock = lambda p, x: 4.0 * np.pi * G * rhofunc_rock(x) * x * x
        gravity_layers_rock = np.ravel(odeint(poisson_rock, gravity_layers_core[-1], radii_rock))

        radii_water = radii[num_core_layers+num_mantle_layers:]
        density_water = density[num_core_layers+num_mantle_layers:]
        rhofunc_water = interpolate.InterpolatedUnivariateSpline(radii_water, density_water,k=3)
        poisson_water = lambda p, x: 4.0 * np.pi * G * rhofunc_water(x) * x * x
        gravity_layers_water = np.ravel(odeint(poisson_water,gravity_layers_rock[-1],radii_water))

        gravity_layers = np.concatenate((gravity_layers_core,gravity_layers_rock,gravity_layers_water),axis=0)

        gravity_layers[1:] = gravity_layers[1:] / radii[1:] / radii[1:]
        gravity_layers[0] = 0

        return(gravity_layers)
    else:

        radii_core = radii[:num_core_layers]
        density_core = density[:num_core_layers]

        rhofunc_core = interpolate.InterpolatedUnivariateSpline(radii_core, density_core,k=4)
        poisson_core = lambda p, x: 4.0 * np.pi * G * rhofunc_core(x) * x * x
        gravity_layers_core = np.ravel(odeint(poisson_core, 0., radii_core))

        radii_rock = radii[num_core_layers:num_core_layers+num_mantle_layers]
        density_rock = density[num_core_layers:num_core_layers+num_mantle_layers]


        rhofunc_rock = interpolate.InterpolatedUnivariateSpline(radii_rock, density_rock,k=4)
        poisson_rock = lambda p, x: 4.0 * np.pi * G * rhofunc_rock(x) * x * x
        gravity_layers_rock = np.ravel(odeint(poisson_rock, gravity_layers_core[-1], radii_rock))

        gravity_layers = np.concatenate((gravity_layers_core,gravity_layers_rock),axis=0)


    gravity_layers[1:] = gravity_layers[1:]/radii[1:]/radii[1:]
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

    depths = radii[-1] - radii


    depths_mant = depths[num_core_layers:(num_core_layers+num_mantle_layers)]
    gravity_mant = gravity[num_core_layers:(num_core_layers+num_mantle_layers)]
    density_mant = density[num_core_layers:(num_core_layers+num_mantle_layers)]

    depths_core = depths[:num_core_layers]
    gravity_core = gravity[:num_core_layers]
    density_core = density[:num_core_layers]
    if number_h2o_layers > 0:


        depths_water = depths[num_core_layers + num_mantle_layers:]
        gravity_water = gravity[num_core_layers + num_mantle_layers:]
        density_water = density[num_core_layers + num_mantle_layers:]

        rhofunc_water = interpolate.InterpolatedUnivariateSpline(depths_water[::-1], density_water[::-1], k=4)
        gfunc_water = interpolate.InterpolatedUnivariateSpline(depths_water[::-1], gravity_water[::-1], k=4)
        p_func = lambda p, x: gfunc_water(x) * rhofunc_water(x)

        pressure_water = np.ravel(odeint(p_func, 1e5, depths_water[::-1]))

        WMB_pres = pressure_water[-1]

        rhofunc_mant = interpolate.InterpolatedUnivariateSpline(depths_mant[::-1], density_mant[::-1])
        gfunc_mant = interpolate.InterpolatedUnivariateSpline(depths_mant[::-1], gravity_mant[::-1])
        p_func = lambda p, x: gfunc_mant(x) * rhofunc_mant(x)

        pressure_mant = np.ravel(odeint(p_func, WMB_pres, depths_mant[::-1]))

        rhofunc_core = interpolate.InterpolatedUnivariateSpline(depths_core[::-1], density_core[::-1], k=5)
        gfunc_core = interpolate.InterpolatedUnivariateSpline(depths_core[::-1], gravity_core[::-1], k=5)

        p_func = lambda p, x: gfunc_core(x) * rhofunc_core(x)
        pressure_core = np.ravel(odeint(p_func, pressure_mant[-1] , depths_core[::-1]))

        pressure= np.concatenate((pressure_water,pressure_mant,pressure_core),axis=0) #in pascals
        pressure = np.asarray([(i*ToBar) for i in pressure]) #to bar

        return pressure[::-1]
    else:

        rhofunc_mant = interpolate.InterpolatedUnivariateSpline(depths_mant[::-1], density_mant[::-1],k=4)
        gfunc_mant = interpolate.InterpolatedUnivariateSpline(depths_mant[::-1], gravity_mant[::-1],k=4)
        p_func = lambda p, x: gfunc_mant(x) * rhofunc_mant(x)

        pressure_mant = np.ravel(odeint(p_func, 1e5, depths_mant[::-1]))

        rhofunc_core = interpolate.InterpolatedUnivariateSpline(depths_core[::-1], density_core[::-1],k=4)
        gfunc_core = interpolate.InterpolatedUnivariateSpline(depths_core[::-1], gravity_core[::-1],k=4)
        p_func = lambda p, x: gfunc_core(x) * rhofunc_core(x)

        pressure_core = np.ravel(odeint(p_func, pressure_mant[-1], depths_core[::-1]))

        pressure= np.concatenate((pressure_mant,pressure_core),axis=0)

        pressure = np.asarray([(i*ToBar) for i in pressure])
        return pressure[::-1]


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
    radius =np.zeros(len(mass))
    val = 4. * np.pi / 3.


    for i in range(len(radius))[1:]:
        mass_lay = mass[i]-mass[i-1]
        radius[i] = np.cbrt((mass_lay/((density[i]+density[i-1])/2)/val) +pow(radius[i-1],3.))

    return radius

def get_mantle_values_T(pressure, temperature, layers,grids):
    P_points_UM = []
    T_points_UM = []
    P_points_LM = []
    T_points_LM = []

    num_mantle_layers = layers[0]
    num_core_layers = layers[1]
    for i in range(num_mantle_layers):
        if pressure[i+num_core_layers] >=125e9*ToBar:
            P_points_LM.append(pressure[i+num_core_layers])
            T_points_LM.append(temperature[i+num_core_layers])
        else:
            P_points_UM.append(pressure[i+num_core_layers])
            T_points_UM.append(temperature[i+num_core_layers])

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


    if len(to_switch_P) > 0: #try in upper mantle
        test = interpolate.griddata((grids[0]['pressure'], grids[0]['temperature']),
                                    grids[0]['cp'], (to_switch_P, to_switch_T), method='linear')


        for i in range(len(test)):
            if np.isnan(test[i]) == True:
                if to_switch_P[i] < (2500e9)*ToBar:
                    print("Specific Heat: Pressure and/or Temperature Exceeds Mantle Grids ")
                    print("%.2f" % (to_switch_P[i] / 10 / 1000), "GPa", "%.2f" % (to_switch_T[i]), "K")
                    print("Grid max:", max(grids[1]['pressure']) / 10 / 1000, "GPa", max(grids[1]['temperature']), "K")
                    sys.exit()
                else:
                    LM_cp_data[to_switch_index[i]] = LM_cp_data[i+1]
            else:
                LM_cp_data[to_switch_index[i]] = test[i]

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
                if to_switch_P[i] < (2500e9 * ToBar):
                    print("Thermal Expans: Pressure and/or Temperature Exceeds Mantle Grids ")
                    print("%.2f" % (to_switch_P[i] / 10 / 1000), "GPa", "%.2f" % (to_switch_T[i]), "K")
                    print("Grid max:", max(grids[1]['pressure']) / 10 / 1000, "GPa", max(grids[1]['temperature']), "K")
                    sys.exit()
                else:
                    LM_alpha_data[to_switch_index[i]] = LM_alpha_data[i + 1]

            else:
                LM_alpha_data[to_switch_index[i]] = test[i]

    spec_heat_mantle = np.concatenate((LM_cp_data, UM_cp_data), axis=0)
    alpha_mant = np.concatenate((LM_alpha_data, UM_alpha_data), axis=0)

    return(alpha_mant,spec_heat_mantle)

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

    P_core = pressure[:num_core_layers]
    T_core = temperature[:num_core_layers]
    core_alpha = interpolate.griddata((grids[2]['pressure'], grids[2]['temperature']),
                                           grids[2]['alpha'],(P_core, T_core), method='linear')
    core_cp = interpolate.griddata((grids[2]['pressure'], grids[2]['temperature']),
                                           grids[2]['cp'],(P_core, T_core), method='linear')

    depths = radii[-1] - radii

    Mantle_potential_temp = structural_parameters.get('Mantle_potential_temp')
    Water_potential_temp = structural_parameters.get('Mantle_potential_temp')
    depths_mantle = depths[num_core_layers:num_core_layers+num_mantle_layers]
    gravity_mantle = gravity[num_core_layers:num_core_layers+num_mantle_layers]
    depths_mantle = depths[num_core_layers:num_core_layers + num_mantle_layers]
    gravity_mantle = gravity[num_core_layers:num_core_layers + num_mantle_layers]
    alpha_mant, spec_heat_mantle = get_mantle_values_T(pressure, temperature, layers, grids)

    depths_core = depths[:num_core_layers]
    gravity_core = gravity[:num_core_layers]



    if number_h2o_layers > 0:

        P_points_water = pressure[num_core_layers+num_mantle_layers:]
        T_points_water = temperature[num_core_layers+num_mantle_layers:]

        depths_water = depths[num_core_layers+num_mantle_layers:]
        gravity_water = gravity[num_core_layers+num_mantle_layers:]

        spec_heat_water = get_water_Cp(P_points_water,T_points_water,grids)
        alpha_water= get_water_alpha(P_points_water,T_points_water,grids)



        grav_func_water = interpolate.InterpolatedUnivariateSpline(depths_water[::-1], gravity_water[::-1],k=3)
        spec_heat_func_water = interpolate.InterpolatedUnivariateSpline(depths_water[::-1], spec_heat_water[::-1],k=3)
        alpha_func_water = interpolate.InterpolatedUnivariateSpline(depths_water[::-1], alpha_water[::-1],k=3)



        adiabat_water = lambda p, x: alpha_func_water(x) * grav_func_water(x) / spec_heat_func_water(x)

        gradient_water = np.ravel(odeint(adiabat_water,0., depths_water[::-1]))

        temperatures_water = [math.exp(i)*Water_potential_temp for i in gradient_water][::-1]

        for i in range(len(alpha_mant))[::-1]:
            if i < len(alpha_mant)-1:
                if alpha_mant[i] > 2*alpha_mant[i+1]:
                    dalpha_dr =  (alpha_mant[i+1]- alpha_mant[i+2])/(depths_mantle[i+1]-depths_mantle[i+2])
                    alpha_mant[i] =alpha_mant[i+1] + (depths_mantle[i]-depths_mantle[i+1])*dalpha_dr

        grav_func_mant = interpolate.InterpolatedUnivariateSpline(depths_mantle[::-1], gravity_mantle[::-1], k=5)
        spec_heat_func_mant = interpolate.InterpolatedUnivariateSpline(depths_mantle[::-1], spec_heat_mantle[::-1], k=5)
        alpha_func_mant = interpolate.InterpolatedUnivariateSpline(depths_mantle[::-1], alpha_mant[::-1], k=5)


        wmb_pressure = pressure[num_core_layers+num_mantle_layers]/10/1000
        starting_grad = 0.0041453 + 0.00221*(wmb_pressure) + -6.6523e-06 * pow(wmb_pressure,2.)
        adiabat_mant = lambda p, x: alpha_func_mant(x) * grav_func_mant(x) / spec_heat_func_mant(x)

        gradient_mant = np.ravel(odeint(adiabat_mant, starting_grad, depths_mantle[::-1]))
        temperatures_mant = [math.exp(i) * Mantle_potential_temp for i in gradient_mant][::-1]

        grav_func_core = interpolate.InterpolatedUnivariateSpline(depths_core[::-1], gravity_core[::-1], k=5)
        spec_heat_func_core = interpolate.InterpolatedUnivariateSpline(depths_core[::-1], core_cp[::-1], k=5)
        alpha_func_core = interpolate.InterpolatedUnivariateSpline(depths_core[::-1], core_alpha[::-1], k=5)

        adiabat_core = lambda p, x: alpha_func_core(x) * grav_func_core(x) / spec_heat_func_core(x)

        gradient_core = np.ravel(odeint(adiabat_core, 0., depths_core[::-1]))

        temperatures_core = [math.exp( i) * temperatures_mant[0] for i in gradient_core][::-1]

        temperatures = np.concatenate((temperatures_core,temperatures_mant,temperatures_water),axis=0)

        return temperatures

    else:

        for i in range(len(alpha_mant))[::-1]:
            if i < len(alpha_mant)-1:
                if alpha_mant[i] > 3*alpha_mant[i+1]:
                    dalpha_dr =  (alpha_mant[i+1]- alpha_mant[i+2])/(depths_mantle[i+1]-depths_mantle[i+2])
                    alpha_mant[i] =alpha_mant[i+1] + (depths_mantle[i]-depths_mantle[i+1])*dalpha_dr

        grav_func_mant = interpolate.InterpolatedUnivariateSpline(depths_mantle[::-1], gravity_mantle[::-1],k=4)
        spec_heat_func_mant = interpolate.InterpolatedUnivariateSpline(depths_mantle[::-1], spec_heat_mantle[::-1],k=4)
        alpha_func_mant = interpolate.InterpolatedUnivariateSpline(depths_mantle[::-1], alpha_mant[::-1],k=4)

        adiabat_mant = lambda p, x: alpha_func_mant(x) * grav_func_mant(x) / spec_heat_func_mant(x)

        gradient_mant = np.ravel(odeint(adiabat_mant, 0., depths_mantle[::-1]))

        temperatures_mant = [math.exp(k) * Mantle_potential_temp for k in gradient_mant][::-1]

        grav_func_core = interpolate.InterpolatedUnivariateSpline(depths_core[::-1], gravity_core[::-1],k=3)
        spec_heat_func_core = interpolate.InterpolatedUnivariateSpline(depths_core[::-1], core_cp[::-1],k=3)
        alpha_func_core = interpolate.InterpolatedUnivariateSpline(depths_core[::-1], core_alpha[::-1],k=3)

        adiabat_core = lambda p, x: alpha_func_core(x) * grav_func_core(x) / spec_heat_func_core(x)

        gradient_core = np.ravel(odeint(adiabat_core, 0., depths_core[::-1]))
        temperatures_core = [math.exp(k) * temperatures_mant[0] for k in gradient_core][::-1]

        temperatures = np.concatenate((temperatures_core,temperatures_mant),axis=0)


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

