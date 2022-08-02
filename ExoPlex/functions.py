import numpy as np
import sys

from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator

from ExoPlex import planet as planet
ToPa = 100000.

#constants, atomic masses
mFe = 55.845
mMg = 24.306
mSi = 28.0867
mO = 15.9994
mS = 32.0650
mCa = 40.078
mAl = 26.981
range_FeO = np.array([0., .02, .04, .06, .08, .1, .15, .20])


def get_FeO(wt_frac_FeO_wanted,FeMg,SiMg,AlMg,CaMg):

    mol_frac_Fe_mantle = ((wt_frac_FeO_wanted / ((mFe + mO)  * FeMg)) * (
                SiMg * mSi + AlMg * mAl + CaMg * mCa + mMg + mO * (
                    2 * SiMg + 1.5 * AlMg + CaMg + 1))) / (1 - wt_frac_FeO_wanted)

    return mol_frac_Fe_mantle

def get_Si_core_w_FeO(compositional_params, verbose):


    FeMg = compositional_params.get('FeMg')
    SiMg = compositional_params.get('SiMg')
    CaMg = compositional_params.get('CaMg')
    AlMg = compositional_params.get('AlMg')
    wt_fe_m = compositional_params.get('wt_frac_FeO_wanted')
    wt_frac_Si_core = compositional_params.get('wt_frac_Si_core')

    conserve_oxy = compositional_params.get('conserve_oxy')

    if conserve_oxy == True:
        if verbose ==True:
            print("Please note this will over-write your choice of wt% Si in core since O is conserved")
        m_feo = mO+mFe
        OTR = mMg + CaMg*mCa + AlMg*mAl + mO*(1+CaMg + 1.5*AlMg)

        mols_fe_m = wt_fe_m*(OTR + SiMg*(mSi+2*mO))/(m_feo - wt_fe_m*(mFe-.5*mSi))
        wt_frac_Si_core = .5*mols_fe_m*mSi / (FeMg*mFe + mols_fe_m*(.5*mSi - mFe))

        mol_frac_Fe_mantle = mols_fe_m/FeMg

        return(wt_frac_Si_core,mol_frac_Fe_mantle)

    else:
        if verbose ==True:
            if wt_fe_m > 0 or wt_frac_Si_core> 0:
                print("Please note that you have FeO or Si in the core but are not conserving oxygen")
        if wt_fe_m > 0:
            mol_frac_Fe_mantle = get_FeO(wt_fe_m,FeMg, SiMg, AlMg, CaMg)

            return(wt_frac_Si_core, mol_frac_Fe_mantle)
        else:
            return(wt_frac_Si_core,0)


def get_percents(compositional_params,verbose):

    """
    This module calculates the bulk composition of the mantle given the elemental ratios provided in the input files

    Parameters
    ----------
    compositional_params : list of floats
        List containing molar ratios of planets, fraction of Fe in core, core composition and flags for writing
        individual phase and whether to use grids

    Returns
    -------
    Core_wt_per : dictionary
        composition of core for individual elements:math:`[wt\%]`
    Mantle_wt_per : dictionary
        composition of mantle in oxides :math:`[wt\%]`
    Core_mol_per : dictionary
        composition of core in mole percent
    core_mass_frac: float
        Mass of core divided by mass of planet
    """
    FeMg = compositional_params.get('FeMg')
    SiMg = compositional_params.get('SiMg')
    CaMg = compositional_params.get('CaMg')
    AlMg = compositional_params.get('AlMg')
    wt_fe_m = compositional_params.get('wt_frac_FeO_wanted')
    wt_frac_Si_core = compositional_params.get('wt_frac_Si_core')
    wt_frac_O_core = compositional_params.get('wt_frac_O_core')
    wt_frac_S_core = compositional_params.get('wt_frac_S_core')
    conserve_oxy = compositional_params.get('conserve_oxy')
    use_grids = compositional_params.get('use_grids')

    MgSi = 1./SiMg
    FeSi = FeMg*MgSi
    CaSi = CaMg*MgSi
    AlSi = AlMg*MgSi


    # system of equation which when solved gives molar values of elements in mantle and core
    # x = [1 nFec,2 nSic, 3 nOc, 4 nSc | ,5 nFem, 6 nMgm,7 nSim,8 nOm, 9 nCam, 10 nAlm]
    #Sum of all masses = 100 g
    b = np.array([0., 0. , 0. , 0. ,  0. , 0. , 0 , 0.,0., 100.])

    wt_frac_Si_core, mol_frac_Fe_mantle = get_Si_core_w_FeO(compositional_params,verbose)

############################################################################################


    #Return only moles of elements, no molecules
    A = np.array([ [1., -1. * FeSi, 0., 0., 1., 0., -1. * FeSi, 0., 0., 0.] ,
                   [wt_frac_Si_core * mFe, (wt_frac_Si_core - 1.) * mSi, wt_frac_Si_core * mO, wt_frac_Si_core * mS ,\
                    0., 0., 0., 0., 0., 0.],
                   [wt_frac_O_core * mFe, (wt_frac_O_core) * mSi, (wt_frac_O_core - 1.) * mO, wt_frac_O_core * mS
                       , 0., 0., 0., 0., 0., 0.],
                   [wt_frac_S_core * mFe, (wt_frac_S_core) * mSi, (wt_frac_S_core) * mO, (wt_frac_S_core - 1.) * mS
                       , 0., 0., 0., 0., 0., 0.],
                   [mol_frac_Fe_mantle, 0., 0., 0., (mol_frac_Fe_mantle - 1.), 0., 0., 0., 0., 0.],
                   [0., -1. * MgSi, 0., 0., 0., 1., -1. * MgSi, 0., 0., 0.],
                   [0., -1. * CaSi, 0., 0., 0., 0., -1. * CaSi, 0., 1., 0.],
                   [0., -1. * AlSi, 0., 0., 0., 0., -1. * AlSi, 0., 0., 1.],
                   [0., 0., 0. , 0. , -1. ,-1., -2. , 1. , -1., -1.5] ,
          [mFe , mSi , mO, mS , mFe , mMg , mSi ,mO, mCa, mAl]])



    #returns number of moles of each element between core and mantle assuming 100 moles total.


    Num_moles = np.linalg.solve(A,b)

    ## find masses and wt% below for perplex ##

    #Splitting up into lists
    #in order Fe, Si, O, S in core
    Core_moles = Num_moles[:4]


    #in order Fe, Mg, Si, O, Ca, Al in mantle
    Mantle_moles = Num_moles[4:]

    for i in Core_moles:
        if i <-1*np.finfo(float).eps:
            print("You have a negative core composition, therefore your chosen composition is not valid stoichiometrically" )
            print(Core_moles)
            sys.exit()
    for i in Mantle_moles:
        if i <-1*np.finfo(float).eps:
            print("You have a negative mantle composition, therefore your chosen composition is not valid stoichiometrically" )
            sys.exit()
    tot_moles_core = sum(Core_moles)


    mass_of_Core = (mFe*(Core_moles[0])+(mSi*Core_moles[1])\
                    +((mO)*Core_moles[2])+(mS*Core_moles[3]))

    mass_of_Mantle = (mFe*Mantle_moles[0])+(mMg*Mantle_moles[1])\
                   +(mSi*Mantle_moles[2])+(mO*Mantle_moles[3])\
                   +(mCa*Mantle_moles[4])+(mAl*Mantle_moles[5])

    # Mass of 100 mole planet in g
    Mtot= mass_of_Core+mass_of_Mantle

    #Core mass fraction of planet
    core_mass_frac = mass_of_Core/Mtot

    if core_mass_frac <=0:
        print("Your chosen composition has a negative core mass fraction")
        print("This is likely due to too low of an Fe/Mg and too high of a mantle FeO content")
        sys.exit()


    #Weight percents of mantle oxides
    #Weight percents assuming FeO, MgO, SiO2, CaO, Al2O3

    Fe_moles_mant = Mantle_moles[0]
    Mg_moles_mant = Mantle_moles[1]
    Si_moles_mant = Mantle_moles[2]
    Ca_moles_mant = Mantle_moles[4]
    Al_moles_mant = Mantle_moles[5]


    if Si_moles_mant/Mg_moles_mant <0.5 or Si_moles_mant/Mg_moles_mant >2:
        if use_grids == True:
            print("Using grids and have a mantle Si/Mg out of range of ExoPlex grids")
            print("which are valid 0.5 <= Si/Mg <= 2")
            print(Si_moles_mant/Mg_moles_mant)
            print(core_mass_frac)
            sys.exit()

    FeO_mant_wt = Fe_moles_mant*(mFe+mO)/mass_of_Mantle
    MgO_mant_wt = Mg_moles_mant*(mMg+mO)/mass_of_Mantle
    SiO2_mant_wt = Si_moles_mant*(mSi+(2.*mO))/mass_of_Mantle
    CaO_mant_wt = Ca_moles_mant*(mCa+mO)/mass_of_Mantle
    Al2O3_mant_wt = (Al_moles_mant/2.)*(2.*mAl+3.*mO)/mass_of_Mantle

    #Throw exception, not if statement
    total = float(FeO_mant_wt+MgO_mant_wt+SiO2_mant_wt+CaO_mant_wt+Al2O3_mant_wt)

    assert (total < 1.+(5.*np.finfo(float).eps) and total > 1.-(5.*np.finfo(float).eps)) == True, str(total) + "total mantle wt% do not add up to 1. Total is " + str(total)

    #Fix for PerPlex to input and round to the 8th decimal point
    #converts to percentages

    FeO_mant_wt = abs(round(FeO_mant_wt * 100., 8))
    SiO2_mant_wt = abs(round(SiO2_mant_wt * 100., 8))
    MgO_mant_wt  = abs(round(MgO_mant_wt * 100., 8))
    CaO_mant_wt  = abs(round(CaO_mant_wt * 100., 8))
    Al2O3_mant_wt  = abs(round(Al2O3_mant_wt * 100., 8))
    Mantle_wt_per = {'FeO': FeO_mant_wt, 'SiO2': SiO2_mant_wt, 'MgO': MgO_mant_wt, \
                     'CaO': CaO_mant_wt,'Al2O3':Al2O3_mant_wt}

    # this is the wt% of elements in the core

    Fe_core_wt   = (Core_moles[0])*(mFe)/mass_of_Core
    Si_core_wt   = Core_moles[1]*(mSi)/mass_of_Core
    O_core_wt   = Core_moles[2]*(mO)/mass_of_Core
    S_core_wt   = Core_moles[3]*(mS)/mass_of_Core

    total_core = (S_core_wt+O_core_wt+Si_core_wt+Fe_core_wt)
    assert (total_core < 1.+(5.*np.finfo(float).eps) and total_core > 1.-(5.*np.finfo(float).eps)) == True, str(total_core) + \
                    "total core wt% do not add up to 1. Total is " + str(total_core)

    #make inequality not, absolute if. Use machine precision

    Fe_core_wt = abs(round(Fe_core_wt*100.,8))
    Si_core_wt  = abs(round(Si_core_wt*100.,8))
    O_core_wt = abs(round(O_core_wt*100.,8))
    S_core_wt  = abs(round(S_core_wt*100.,8))

    Core_wt_per = {'Fe':Fe_core_wt,'Si':Si_core_wt,'O':O_core_wt,'S':S_core_wt}

    if verbose == True:
        print()
        print("Core composition: ",Core_wt_per)
        print("Mantle composition: ", Mantle_wt_per)

        print("Mantle Fe#", Mantle_moles[0]/(Mantle_moles[0]+Mantle_moles[1]))
        print("Core Mass Percent = ", '%.3f'%(core_mass_frac*100.))
        print()


    return(Core_wt_per,Mantle_wt_per,core_mass_frac)


def get_phases(Planet,grids,layers,combine_phases):
    """
    This module creates the output file of Pressure, Temperature and phase fractions

    Parameters
    ----------
    Planet: dictionary
        Dictionary of pressure, temperature, expansivity, specific heat and phases for modeled planets
    grids: list of lists
        Upper and lower mantle phase diagrams
    layers: list
        number of core layers, mantle layers and water layers
    combine_phases: boolean
        True if you'd like individual endmembers combined into constituent rocks, False if you'd like all endmembers

    Returns
    -------
    Phases : list
        List of all phases present
    new_names : list
        composition of mantle in oxides :math:`[wt\%]`
    """
    num_mantle_layers, num_core_layers, number_h2o_layers = layers


    mantle_pressures = Planet['pressure'][num_core_layers:(num_core_layers+num_mantle_layers)]
    mantle_temperatures = Planet['temperature'][num_core_layers:(num_core_layers+num_mantle_layers)]

    P_points_UM = []
    T_points_UM = []
    P_points_LM = []
    T_points_LM = []
    for i in range(len(mantle_pressures)):
        if mantle_pressures[i]<=1250000.:
            P_points_UM.append(mantle_pressures[i])
            T_points_UM.append(mantle_temperatures[i])
        else:
            P_points_LM.append(mantle_pressures[i])
            T_points_LM.append(mantle_temperatures[i])

    interpolator_phase_UM = grids[0]['phases']
    interpolator_phase_LM = grids[1]['phases']

    mesh_UM = np.vstack((P_points_UM, T_points_UM)).T
    Mantle_phases_UM = interpolator_phase_UM(mesh_UM)

    mesh_LM = np.vstack((P_points_LM, T_points_LM)).T
    Mantle_phases_LM = interpolator_phase_LM(mesh_LM)

    mesh_LM = np.vstack((P_points_LM, T_points_LM)).T
    Mantle_phases_LM = interpolator_phase_LM(mesh_LM)

    if (num_mantle_layers-len(P_points_UM)) > 0:
        Mantle_phases = np.concatenate((Mantle_phases_LM,Mantle_phases_UM),axis=0)

    else:
        Mantle_phases = Mantle_phases_UM

    for i in range(len(Mantle_phases)):
        entry = np.zeros(len(Mantle_phases[i]))
        tot = sum(Mantle_phases[i])
        for j in range(len(Mantle_phases[i])):
            entry[j] = 100.*Mantle_phases[i][j]/tot
        Mantle_phases[i] = entry


    if number_h2o_layers>0:
        interpolator_phase_water = grids[3]['phases']

        water_P = Planet['pressure']
        water_T = Planet['temperature']

        mesh_water = np.vstack((water_P, water_T)).T
        Water_phases = interpolator_phase_water(mesh_water)


        for i in range(len(Water_phases)):
            for j in range(len(Water_phases[i])):
                if np.isnan(Water_phases[i][j])==True:
                    Water_phases[i][j]=0

        Core_phases = np.asarray([[100] for i in range(num_core_layers)])
        rest = np.asarray([[0] for i in range((num_mantle_layers+number_h2o_layers))])
        Core_phases = np.vstack((Core_phases,rest))

        man_core = np.asarray([np.zeros(len(Mantle_phases[0])) for i in range(num_core_layers)])
        man_wat = np.asarray([np.zeros(len(Mantle_phases[0])) for i in range(number_h2o_layers)])

        Mantle_phases = np.vstack((man_core, Mantle_phases, man_wat))

        Phases = np.zeros(shape=(sum(layers),len(Mantle_phases[0])+len(Core_phases[0])+len(Water_phases[0])))
    else:
        Core_phases = np.asarray([[100] for i in range(num_core_layers)])
        rest = np.asarray([[0] for i in range((num_mantle_layers))])
        Core_phases = np.vstack((Core_phases, rest))

        man_core = np.asarray([np.zeros(len(Mantle_phases[0])) for i in range(num_core_layers)])

        Mantle_phases = np.vstack((man_core, Mantle_phases))
        Phases = np.zeros(shape=(sum(layers),len(Mantle_phases[0])+len(Core_phases[0])))



    for i in range(sum(layers)):
        if number_h2o_layers > 0:
            Phases[i] = np.hstack((Mantle_phases[i],Core_phases[i],Water_phases[i]))
        else:
            Phases[i] = np.hstack((Mantle_phases[i],Core_phases[i]))

    phase_names = Planet['phase_names']

    if combine_phases == True:
        solutions = {'c2/c':'c2c_ss', 'fc2/c':'c2c_ss', 'per':'ferropericlase', 'wus': 'ferropericlase', 'aperov':'perovskite', 'fperov':'perovskite','perov':'perovskite',
        'ab':'plagioclase','an':'plagioclase','sp':'spinel','herc':'spinel', 'fo':'olivine', 'fa':'olivine', 'wad':'wadslyite','fwad':'wadslyite',
        'ring':'ringwoodite','fring':'ringwoodite','odi':'orthopyroxene','en':'orthopyroxene', 'fs':'orthopyroxene','ts':'orthopyroxene',
        'jd':'clinopyroxene','di':'clinopyroxene','hed':'clinopyroxene','cen':'clinopyroxene','cts':'clinopyroxene', 'cor':'akimotoite',
        'aki':'akimotoite','faki':'akimotoite','gr':'garnet','alm':'garnet','maj':'garnet','cmaj':'garnet','py':'garnet','jmaj':'garnet','appv':'postperovskite',
        'ppv':'postperovskite','fppv':'postperovskite','mfer':'cf','ffer':'cf','nfer':'cf', 'ca-pv':'ca-perovskite','cfs':'cfs',
        'ky':'kyanite','neph':'nephaline','coe':'coesite','seif':'seiferite','q':'quartz','s':'stishovite'}

        phase_amounts = {'c2c_ss':0.,'ferropericlase':0.,'perovskite':0.,'plagioclase':0.,'spinel':0.,'olivine':0.,'wadslyite':0.,
                         'ringwoodite':0.,'orthopyroxene':0.,'clinopyroxene':0.,'akimotoite':0.,'garnet':0.,'postperovskite':0.,
                         'kyanite':0.,'nephaline':0.,'coesite':0.,'seiferite':0.,'quartz':0.,'stishovite':0.,'cf':0.,'ca-perovskite':0.,'cfs':0.}


        for j in phase_names:
            if j not in solutions.keys():
                phase_amounts.update({j:0.})
                solutions.update({j:j})

        new_names = phase_amounts.keys()
        final_output = []
        for i in range(len(Phases)):
            phase_amounts = dict.fromkeys(phase_amounts,0)
            line_output = []
            phase_line = Phases[i]

            for j in range(len(phase_names)):
                phase_name = phase_names[j]
                solution_to_put = solutions.get(phase_name)
                phase_amounts.update({solution_to_put: phase_amounts.get(solution_to_put) + phase_line[j]})


            for k in new_names:
               line_output.append(phase_amounts.get(k))

            final_output.append(line_output)

        Phases = np.array(final_output)
    else:
        new_names = phase_names

    return Phases, new_names

def write(Planet,filename):
    """
    This module creates the output file for reading later

    Parameters
    ----------
    Planet: dictionary
        Dictionary of pressure, temperature, expansivity, specific heat and phases for modeled planets
    filename: string
       chosen filename for output file

    Returns
    -------
    None
    """

    output = []
    for i in range(len(Planet['pressure'])):
        line_item = [(Planet['radius'][-1]-Planet['radius'][i])/1000.,Planet['radius'][i]/1000., Planet['mass'][i]/5.97e24,
                     Planet['density'][i]/1000.,Planet['pressure'][i]/10000.,Planet['temperature'][i],Planet['gravity'][i]]

        for j in range(len(Planet['phases'][i])):
            line_item.append(Planet['phases'][i][j])

        output.append(line_item)
    line_name = []
    line_name.append('Depth_km')
    line_name.append('Radius_km')
    line_name.append('Mass_Me')
    line_name.append('Density_g_cc')
    line_name.append('Pressure_GPa')
    line_name.append('Temperature_K')
    line_name.append('Gravity_m_s-2')

    for i in Planet['phase_names']:
        line_name.append(str(i))

    string_element = ','.join(line_name)
    filename = 'Data/'+filename

    np.savetxt(filename+'.csv', output, '%.5f', ",", newline='\n',
                header=string_element, footer='', comments='# ')

    print()
    print("Detailed data file written to:", filename+'.txt')
    print()
    return 0

def calc_planet_radius(mass_guess,args):
    radius_planet = args[0]
    structure_params = args[1]
    compositional_params = args[2]
    layers = args[3]
    grids = args[4]
    Core_wt_per = args[5]
    core_mass_frac = args[6]
    verbose = args[7]

    if verbose == True:
        print("Trying mass: ", round(mass_guess,5), "Earth masses")
    Planet = planet.initialize_by_mass(
        *[mass_guess, structure_params, compositional_params, layers, core_mass_frac])

    Planet = planet.compress_mass(*[Planet, grids, Core_wt_per, structure_params, layers, verbose])
    if verbose == True:
        print()
        print("Difference between actual and guess radius = ", radius_planet - (Planet['radius'][-1]/6371e3))
        print()
    return (1.-(radius_planet/(Planet['radius'][-1]/6371e3)))

def find_Planet_radius(radius_planet, core_mass_frac, structure_params, compositional_params, grids, Core_wt_per, layers,verbose):
    """
    This module contains functions to determine the a planet's mass when given radius and composition. Because we conserve the mass ratios of the planet \
    we must iterate over core mass fraction to match the input composition.

    Parameters
    ----------
    radius_planet: float
        Radius of planet input by user :math:`[m]`
    core_mass_frac: float
        Core mass fraction from composition :math:`[wt\%]`
    structure_params: list
        Structural parameters of the planet; See example for description
    compositional_params: list
        Structural parameters of the planet; See example for description
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

    from scipy.optimize import brentq
    args = [radius_planet, structure_params, compositional_params, layers, grids, Core_wt_per, core_mass_frac, verbose]
    den_guess = 1000 * (2.43 + 3.39 * radius_planet)  # From Weiss and Marcy, 2013

    #Mass = brentq(calc_planet_radius, 0.5*mass_guess,  1.3*mass_guess, args=args, xtol=1e-5)
    Mass = brentq(calc_planet_radius, .2,  13., args=args, xtol=1e-5)

    Planet = planet.initialize_by_mass(
        *[Mass, structure_params, compositional_params, layers, core_mass_frac])

    Planet = planet.compress_mass(*[Planet, grids, Core_wt_per, structure_params, layers, verbose])

    return Planet

def find_Planet_mass(mass_planet, core_mass_frac, structure_params, compositional_params, grids, Core_wt_per, layers,verbose):
    """
    This module contains functions to determine the a planet's radius when given mass and composition.

    Parameters
    ----------
    radius_planet: float
        Radius of planet input by user :math:`[m]`
    core_mass_frac: float
        Core mass fraction from composition :math:`[wt\%]`
    structure_params: list
        Structural parameters of the planet; See example for description
    compositional_params: list
        Structural parameters of the planet; See example for description
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
    Planet = planet.initialize_by_mass(*[mass_planet, structure_params, compositional_params, layers,core_mass_frac])

    Planet = planet.compress_mass(*[Planet, grids, Core_wt_per, structure_params, layers,verbose])


    return Planet


def find_filename(compositional_params,verbose):
    """
    This module determines the closest compositional grid to pull from for a given composition

    Parameters
    ----------
    compositional_params: list
        Structural parameters of the planet; See example for description
    verbose: bool


    Returns
    -------
    filename: string
        Name of file containing the grids based on this composition
    """
    mFe      = 55.845
    mMg      = 24.306
    mSi      = 28.0867
    mO       = 15.9994
    mS       = 32.0650
    mCa      = 40.078
    mAl      = 26.981

    mantle_wt_percents = get_percents(compositional_params,verbose)[1]


    mol_Mg = mantle_wt_percents.get('MgO')/(mMg + mO)
    SiMg = round((mantle_wt_percents.get('SiO2')/(mSi+2*mO))/mol_Mg,4)
    CaMg = round((mantle_wt_percents.get('CaO')/(mCa+mO))/mol_Mg,4)
    AlMg =round((2.*mantle_wt_percents.get('Al2O3')/(2.*mAl+3.*mO))/mol_Mg,4)
    FeO = mantle_wt_percents.get('FeO')/100.


    range_SiMg = np.array([0.5 + .1*x for x in range(16)])
    range_CaMg = np.array([0.02 + 0.01*x for x in range(8)])

    range_AlMg = np.array([0.04 + 0.01*x for x in range(8)])
    range_FeO = np.array([0.,.02,.04,.06,.08,.1,.15,.20])

    assert SiMg <= max(range_SiMg) ,"Si/Mg is too high and outside of range of grids. max = "+ str(max(range_SiMg))+ " your value = "+ str(SiMg)
    assert  SiMg >= min(range_SiMg),"Si/Mg is outside of range of grids, try using PerPlex and set use_grids = False"
    assert CaMg <= max(range_CaMg), "Ca/Mg is outside of range of grids, try using PerPlex and set use_grids = False"
    assert CaMg >= min(range_CaMg), "Ca/Mg is outside of range of grids, try using PerPlex and set use_grids = False"
    assert AlMg <= max(range_AlMg), "Al/Mg is outside of range of grids, try using PerPlex and set use_grids = False"
    assert AlMg >= min(range_AlMg), "Al/Mg is outside of range of grids, try using PerPlex and set use_grids = False"
    assert FeO <= max(range_FeO), "FeO is outside of range of grids, try using PerPlex and set use_grids = False"
    assert FeO >= min(range_FeO), "FeO is outside of range of grids, try using PerPlex and set use_grids = False"

    SiMg_file = range_SiMg[(np.abs(range_SiMg - SiMg)).argmin()]
    CaMg_file = range_CaMg[(np.abs(range_CaMg - CaMg)).argmin()]
    AlMg_file = range_AlMg[(np.abs(range_AlMg - AlMg)).argmin()]
    FeO_file = str('%.02f'% FeO)

    if SiMg_file >= 1.1:
        SiMg_file = '%.1f'%(SiMg_file)

    FeMg_file  = '0.00'
    NaMg_file = 0.0

    filename = str(CaMg_file)+"CaMg_"+str(FeMg_file)+"FeMg_"+str(AlMg_file)+"AlMg_"+str(SiMg_file)+"SiMg_"+str(NaMg_file)+"NaMg_"+FeO_file+"Fe"

    if verbose ==True:
        print("Your closest grid filename is: ", filename+"_U/LM_results.txt")
        print()

    return filename

def check(Planet):
    def monotonic(x):
        dx = np.diff(x)
        return np.all(dx <= 0) or np.all(dx >= 0)
    T_test = monotonic(Planet['temperature'])
    P_test = monotonic(Planet['pressure'])
    rho_test = monotonic(Planet['density'])
    plot = False
    if T_test == False:
        print("Temperature is not monotonic, thus an error. Will plot.")
        plot = True
    if P_test == False:
        print("Pressure is not monotonic, thus an error. Will plot.")
        plot = True
    if rho_test == False:
        print("density is not monotonic, thus an error. Will plot.")
        plot = True

    if plot == True:
        import matplotlib.pyplot as plt
        figure = plt.figure(figsize=(10, 9))

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
        ax4.set_ylim(0., max(Planet['temperature']) + 300)
        ax4.minorticks_on()

        plt.show()
        sys.exit()


