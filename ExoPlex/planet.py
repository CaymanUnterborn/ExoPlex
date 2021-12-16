import sys
import numpy as np
Earth_radius = 6371e3
Earth_mass = 5.97e24
ToPa = 100000.
ToBar = 1./ToPa
from ExoPlex import minphys as minphys
from scipy import interpolate


def initialize_by_radius(*args):
    from ExoPlex import burnman as bm

    """
   This module creates the dictionary of lists for each planetary parameter (e.g., density) for a planet of the radius
   input by user.

    Parameters
    ----------
    radius_planet: float
        input radius of planet in Earth radii

    structural_params: list
        Structural parameters of the planet; See example for description

    compositional_params: list
        Structural parameters of the planet; See example for description

    Core_wt_per: dict
        composition of the core
    layers: list
        Number of layers for core, mantle and water
    Returns
    -------
    Planet: dictionary
        Dictionary of initial guess of pressure, temperature, expansivity, specific heat and phases for modeled planet
        keys = 'radius','density','temperature','gravity','pressure', 'alpha','cp','Vphi''Vp','Vs','K'
    """

    rock = bm.Composite([bm.minerals.SLB_2011.mg_perovskite(),
                              bm.minerals.SLB_2011.periclase()], \
                             [0.8, 0.2])
    radius_planet = args[0]
    structural_params = args[1]
    compositional_params = args[2]
    Core_wt_per = args[3]

    num_mantle_layers, num_core_layers, number_h2o_layers = args[4]
    core_grid = args[5][2]




    wt_frac_water, FeMg, SiMg, CaMg, AlMg, mol_frac_Fe_mantle, wt_frac_Si_core, \
     wt_frac_O_core, wt_frac_S_core,combine_phases,use_grids = compositional_params
    core_mass_frac = compositional_params[4]

    core_rad_frac = structural_params[6]
    Mantle_potential_temp = structural_params[7]
    water_rad_frac = structural_params[8]
    water_potential_temp = structural_params[9]


    # if there is a water layer, the imput temperature is lowered because that temperature is for the crustal layer
    # also 50 shells are used for the water layer hence the nh20 vaiable

    if wt_frac_water == 0. and number_h2o_layers > 0:
       print ("You have layers of water but no water!")
       number_h2o_layers = 0

    if number_h2o_layers >0:
        water_rad_frac = 0.1

    num_layers = num_core_layers+num_mantle_layers + number_h2o_layers # add 50 shells if there is an h2O layer
    # arrays to be used

    # these grids are initialized in this function and are then a passed around through the other routines and edited each iteration
    radius_layers = np.zeros(num_layers)
    density_layers = np.zeros(num_layers)
    volume_layers = np.zeros(num_layers)
    mass_layers = np.zeros(num_layers)
    cumulative_mass = np.zeros(num_layers)
    Temperature_layers = np.zeros(num_layers)
    # used for compressiofh2on funciton

    gravity_layers = np.zeros(num_layers)
    Pressure_layers = np.zeros(num_layers)
    alpha = np.zeros(num_layers)
    cp = np.zeros(num_layers)
    Vphi = np.zeros(num_layers )

    Vs = np.zeros(num_layers )
    Vp = np.zeros(num_layers)
    K =  np.zeros(num_layers)

    mass_planet_guess = 6000. * (4*np.pi/3.)*pow(radius_planet*6371e3,3)
    CMF_guess = 25000* (4*np.pi/3.)*pow(core_rad_frac* radius_planet* 6371e3,3)/mass_planet_guess

    planet_radius_guess = radius_planet*Earth_radius
    water_thickness_guess = water_rad_frac*planet_radius_guess
    core_thickness_guess = core_rad_frac * (planet_radius_guess-water_thickness_guess)
    mantle_thickness_guess = planet_radius_guess - water_thickness_guess - core_thickness_guess



    water_mass = (wt_frac_water * mass_planet_guess)
    core_mass = (CMF_guess * (mass_planet_guess - water_mass))
    mantle_mass = (mass_planet_guess ) - water_mass - core_mass

    CMB_T_guess = (4180.*(radius_planet-0)-2764.*pow(radius_planet-0,2.)+1219.*pow(radius_planet-0,3.)) + (Mantle_potential_temp-1600)*(0.82+pow(radius_planet-0,1.81))
    if CMB_T_guess > 7000:
        CMB_T_guess = 6900

    CMB_P_guess = 10000*(262.*(radius_planet)-550.*pow(radius_planet,2.) + 432.*pow(radius_planet,3.))

    dP_dr = (CMB_P_guess-5000)/(num_mantle_layers)
    dT_dr = (CMB_T_guess-Mantle_potential_temp)/(num_mantle_layers)

    for i in range(num_layers):

            if i<number_h2o_layers:
                Pressure_layers[i] =  1 + i*(wt_frac_water*10000)
                Temperature_layers[i] = water_potential_temp+i/10.

            elif i <= number_h2o_layers+num_mantle_layers-1:

                Pressure_layers[i] = 5000. + dP_dr*(i-number_h2o_layers)
                Temperature_layers[i] = Mantle_potential_temp+dT_dr*(i-number_h2o_layers)

            else:
                Pressure_layers[i] = Pressure_layers[number_h2o_layers+num_mantle_layers-1]  + 3.5*(dP_dr*num_mantle_layers/num_core_layers)*(i-number_h2o_layers-num_mantle_layers)
                Temperature_layers[i] = 1900+i


    Pressure_layers = Pressure_layers[::-1]
    Temperature_layers= Temperature_layers[::-1]


    P_core = Pressure_layers[:num_core_layers]
    T_core = Temperature_layers[:num_core_layers]
    rho_core = interpolate.griddata((core_grid['pressure'], core_grid['temperature']),
                                            core_grid['density'], (P_core, T_core), method='linear')

    for i in range(num_layers):
        if i < num_core_layers:
            if i == 0:
                radius_layers[i]= 1
                density_layers[i] = rho_core[i]
                mass_layers[i] = density_layers[i]*(4*np.pi/3.) *pow(radius_layers[i],3)

            else:
                radius_layers[i] =((float(i) / num_core_layers) * core_thickness_guess)
                density_layers[i] = rho_core[i]
                mass_layers[i] = density_layers[i] * (4 * np.pi / 3.) * (pow(radius_layers[i], 3)-pow(radius_layers[i-1],3) )


        elif i <= (num_core_layers + num_mantle_layers):
            radius_layers[i] = core_thickness_guess+ ((((i - num_core_layers) / num_mantle_layers) * mantle_thickness_guess))
            density_layers[i]=(rock.evaluate(['density'], Pressure_layers[i]*ToPa, 300))
            mass_layers[i] = density_layers[i] * (4 * np.pi / 3.) * (
                        pow(radius_layers[i], 3) - pow(radius_layers[i - 1], 3))

        else:
            radius_layers[i] = core_thickness_guess + mantle_thickness_guess + \
                               ((float(
                                   i - num_core_layers - num_mantle_layers) / number_h2o_layers) * water_thickness_guess)
            mass_layers[i] = (water_mass / num_mantle_layers)
            Volume_out = 4. * np.pi / 3 * pow(radius_layers[i] , 3)
            Volume_in = 4. * np.pi / 3 * pow(radius_layers[i - 1] , 3)
            density_layers[i] = 1000+i

    mass_update = np.zeros(num_layers)

    radius_layers[-1] = planet_radius_guess

    for i in range(len(mass_layers)):
        mass_update[i] = (sum(mass_layers[:i + 1]))

    mass_layers = mass_update
    mass_layers[-1] = mass_planet_guess


    keys = ['radius','mass','density','temperature','gravity','pressure',\
            'alpha','cp','Vphi''Vp','Vs','K']



    return dict(zip(keys,[radius_layers, mass_layers,density_layers,Temperature_layers,gravity_layers, Pressure_layers,
                          alpha, cp,Vphi,Vp,Vs,K]))

def initialize_by_mass(*args):
    from ExoPlex import burnman as bm

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
    rock = bm.Composite([bm.minerals.SLB_2011.mg_perovskite(),
                         bm.minerals.SLB_2011.periclase()], \
                        [0.8, 0.2])
    mass_planet = args[0]
    structural_params = args[1]
    compositional_params = args[2]

    num_mantle_layers, num_core_layers, number_h2o_layers = args[3]
    core_mass_frac = args[4]
    core_grid = args[5][2]

    wt_frac_water, FeMg, SiMg, CaMg, AlMg, mol_frac_Fe_mantle, wt_frac_Si_core, \
     wt_frac_O_core, wt_frac_S_core,combine_phases,use_grids = compositional_params

    core_rad_frac = structural_params[6]
    Mantle_potential_temp = structural_params[7]
    water_rad_frac = structural_params[8]
    water_potential_temp = structural_params[9]

    Radius_planet_guess = pow(mass_planet*5.97e24 / 5500. / (4*np.pi/3.),1/3.)/6371e3
    if Radius_planet_guess > 2:
        Radius_planet_guess = 2.

    if number_h2o_layers > 0:
        water_thickness_guess = water_rad_frac*Radius_planet_guess*6371e3
        core_thickness_guess = core_rad_frac * (Radius_planet_guess*6371e3-water_thickness_guess)
        mantle_thickness_guess = Radius_planet_guess*6371e3 - water_thickness_guess - core_thickness_guess
    else:
        core_thickness_guess = core_rad_frac * (Radius_planet_guess*6371e3)
        mantle_thickness_guess = Radius_planet_guess*6371e3 - core_thickness_guess

    if wt_frac_water == 0. and number_h2o_layers > 0:
       print ("You have layers of water but no water!")
       number_h2o_layers = 0

    num_layers = num_core_layers+num_mantle_layers + number_h2o_layers # add 50 shells if there is an h2O layer
    # arrays to be used

    # these grids are initialized in this function and are then a passed around through the other routines and edited each iteration
    radius_layers = np.zeros(num_layers)
    density_layers = np.zeros(num_layers)
    volume_layers = np.zeros(num_layers)
    mass_layers = np.zeros(num_layers)
    cumulative_mass = np.zeros(num_layers)
    Temperature_layers = np.zeros(num_layers)
    radius_layers = np.zeros(num_layers)

    # used for compressiofh2on funciton
    gravity_layers = np.zeros(num_layers)
    Pressure_layers = np.zeros(num_layers)
    alpha = np.zeros(num_layers)
    cp = np.zeros(num_layers)
    Vphi = np.zeros(num_layers )

    Vs = np.zeros(num_layers )
    Vp = np.zeros(num_layers)
    K =  np.zeros(num_layers)

    # 15 mineral phases + 2ice + liquid water #phasechange
    water_mass = (wt_frac_water*mass_planet)*Earth_mass
    core_mass = (core_mass_frac * (mass_planet*Earth_mass-water_mass))
    mantle_mass = (mass_planet*Earth_mass)-water_mass-core_mass

    #Update to use BurnMan
    radius_planet = Radius_planet_guess
    CMB_T_guess = (4180.*(radius_planet-0)-2764.*pow(radius_planet-0,2.)+1219.*pow(radius_planet-0,3.)) + (Mantle_potential_temp-1600)*(0.82+pow(radius_planet-0,1.81))
    if CMB_T_guess > 7000:
        CMB_T_guess = 6800
    #print("CMBT guess",CMB_T_guess)
    CMB_P_guess = 10000*(262.*(radius_planet-0)-550.*pow(radius_planet-0,2.) + 432.*pow(radius_planet-0,3.))
    #print(radius_planet,CMB_T_guess)

    dP_dr = (CMB_P_guess-5000)/(num_mantle_layers)
    dT_dr = (CMB_T_guess-Mantle_potential_temp)/(num_mantle_layers)

    for i in range(num_layers):

            if i<number_h2o_layers:
                Pressure_layers[i] = 1 + i*(wt_frac_water*10000)
                Temperature_layers[i] = water_potential_temp+i/10.

            elif i <= number_h2o_layers+num_mantle_layers-1:
                Pressure_layers[i] = 5000. + dP_dr*(i-number_h2o_layers)
                Temperature_layers[i] = Mantle_potential_temp+dT_dr*(i-number_h2o_layers)

            else:
                Pressure_layers[i] = Pressure_layers[number_h2o_layers+num_mantle_layers-1] + (i-number_h2o_layers-num_mantle_layers)*5000
                #Pressure_layers[number_h2o_layers+num_mantle_layers-1]  + 3.5*(dP_dr*num_mantle_layers/num_core_layers)*(i-number_h2o_layers-num_mantle_layers)
                Temperature_layers[i] = 1900+i


    Pressure_layers = Pressure_layers[::-1]
    Temperature_layers= Temperature_layers[::-1]

    P_core = Pressure_layers[:num_core_layers]
    T_core = Temperature_layers[:num_core_layers]
    rho_core = interpolate.griddata((core_grid['pressure'], core_grid['temperature']),
                                            core_grid['density'], (P_core, T_core), method='linear')
    for i in range(num_layers):

        if i < num_core_layers:
                radius_layers[i] =((float(i) / num_core_layers) * core_thickness_guess)
                density_layers[i] = rho_core[i]
                mass_layers[i] = core_mass/num_core_layers

                #print(i, radius_layers[i] / 6371e3, Volume_out, Volume_in, Volume_out - Volume_in, mass_layers[i],density_layers[i] / 1000)


        elif i <= (num_core_layers + num_mantle_layers):

            density_layers[i]=(rock.evaluate(['density'], Pressure_layers[i]*ToPa, 300.))

            mass_layers[i] = mantle_mass/num_mantle_layers

            radius_layers[i] =  core_thickness_guess+ ((((i - num_core_layers) / num_mantle_layers) * mantle_thickness_guess))



        else:
            radius_layers[i] = core_thickness_guess + mantle_thickness_guess + \
                               ((float(
                                   i - num_core_layers - num_mantle_layers) / number_h2o_layers) * water_thickness_guess)
            mass_layers[i] = (water_mass / number_h2o_layers)
            Volume_out = 4. * np.pi / 3 * pow(radius_layers[i] , 3)
            Volume_in = 4. * np.pi / 3 * pow(radius_layers[i - 1] , 3)
            density_layers[i] = 1000+i
            # print(i,radius_layers[i],Volume_out-Volume_in,mass_layers[i],density_layers[i])

    mass_update = np.zeros(num_layers)

    for i in range(len(mass_layers)):
        mass_update[i] = (sum(mass_layers[:i + 1]))

    mass_layers = mass_update

    #mass_layers[-1] = mass_planet * Earth_mass
    radius_layers[-1] = radius_planet*6371e3


    keys = ['mass','radius', 'density', 'temperature', 'gravity', 'pressure', \
            'alpha', 'cp', 'Vphi''Vp', 'Vs', 'K']


    return dict(zip(keys, [mass_layers, radius_layers,density_layers, Temperature_layers, gravity_layers, Pressure_layers,
                           alpha, cp, Vphi, Vp, Vs, K]))
def compress_radius(*args):
    """
   This module iterates the density, mass within a sphere, adiabatic temperature and gravity integrals for a planet of radius R
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
    Returns
    -------
    Planet: dictionary
        Dictionary of initial guess of pressure, temperature, expansivity, specific heat and phases for modeled planet
    """

    Planet = args[0]
    grids = args[1]
    Core_wt_per = args[2]

    structural_params= args[3]
    layers= args[4]
    verbose = args[5]
    n_iterations = 1
    max_iterations = 100

    old_rho = [10  for i in range(len(Planet['density']))]
    converge = False
    while n_iterations <= max_iterations and converge == False:
        if verbose  ==  True:
            print ("iteration #",n_iterations)

        for i in range(len(Planet['density'])):

            if np.isnan(Planet['density'][i]) == True:
                print ("Density has a nan")
                print (i, Planet['pressure'][i]/10/1000,Planet['temperature'][i])
                print

                sys.exit()

        Planet['density'] = minphys.get_rho(Planet,grids,Core_wt_per,layers)

        Planet['gravity'] = minphys.get_gravity(Planet,layers)

        Planet['pressure'] = minphys.get_pressure(Planet,layers)

        if n_iterations > 1:

            converge, old_rho = minphys.check_convergence(Planet['density'], old_rho)

            Planet['temperature'] = minphys.get_temperature(Planet, grids, structural_params, layers)

            for i in range(len(Planet['temperature'])):

                if Planet['temperature'][i] >= 6900 and i >layers[1]:
                    import matplotlib.pyplot as plt
                    print(i)
                    #plt.plot(Planet['radius'], Planet['temperature'])
                    #plt.ylabel('tempreature')
                    #()
                    #print("in first loop")
                    #print((Planet['tempeerature'][i]))
                    #Planet['temperature'][i-1] = Planet['temperature'][i]
                    for k in range(len(Planet['temperature']))[::-1]:

                        if k > layers[0]:
                            #print("k", k)
                            #if Planet['temperature'][k]>6000:
                                #print(k,Planet['temperature'][k])

                            Planet['temperature'][k] = 1600 + 9.*((sum(layers))-k)
                            #print(sum(layers),k,Planet['temperature'][k])
                    break


            for i in range(sum(layers)-1)[::-1]:
                if Planet['density'][i] < Planet['density'][i+1] and i >0:
                   #print("here",i,Planet['density'][i + 1],Planet['density'][i])
                    Planet['density'][i] = Planet['density'][i+1]

        n_iterations+=1



    return Planet

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
    layers= args[4]
    verbose = args[5]
    n_iterations = 1
    max_iterations = 100

    old_r = [10  for i in range(len(Planet['mass']))]
    converge = False

    while n_iterations <= max_iterations and converge == False:
        if verbose == True:
            print ("iteration #",n_iterations)
        if n_iterations>1:
            converge,old_r = minphys.check_convergence(Planet['density'],old_r)

        for i in range(len(Planet['density'])):
            if np.isnan(Planet['density'][i]) == True:
                print ("Density has a nan at P (GPa), T (K):")
                print (i,Planet['pressure'][i]*0.0001,Planet['temperature'][i])
                print
                sys.exit()

        #import matplotlib.pyplot as plt
        #plt.plot(Planet['radius']/1000,Planet['pressure']/10/1000)
        #plt.show()
        Planet['density'] = minphys.get_rho(Planet,grids,Core_wt_per,layers)

        Planet['radius'] = minphys.get_radius(Planet, layers)

        #for i in range(len(Planet['radius'])-1)[::-1]:
        #    if Planet['density'][i]<Planet['density'][i+1] or Planet['radius'][i+1] < Planet['radius'][i]:
        #        print("Bad",i,Planet['radius'][i]/6371e3, Planet['radius'][i+1]/6371e3)


        Planet['gravity'] = minphys.get_gravity(Planet,layers)

        Planet['temperature'] = minphys.get_temperature(Planet, grids, structural_params, layers)

        Planet['pressure'] = minphys.get_pressure(Planet,layers)

        n_iterations+=1

    return Planet


