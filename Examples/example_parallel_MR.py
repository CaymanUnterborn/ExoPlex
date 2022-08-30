
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
from ExoPlex import make_grids
from ExoPlex import run_perplex as perp
from ExoPlex import make_grids

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
SiMg = 0.9
AlMg = 0.09
FeMg = 1.

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


compositional_params = dict(zip(comp_keys, [wt_frac_water, FeMg, SiMg, CaMg, AlMg, wt_frac_FeO_wanted, wt_frac_Si_core, \
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
Mantle_filename = perp.run_perplex(*[Mantle_wt_per,compositional_params,
                                                [structure_params.get('Pressure_range_mantle_UM'),structure_params.get('Temperature_range_mantle_UM'),
                                                structure_params.get('resolution_UM')],filename,verbose,combine_phases])
grids_low, names = make_grids.make_mantle_grid(Mantle_filename,Mantle_wt_per,True, use_grids)
names.append('Fe')
if layers[-1] > 0:
    water_grid, water_phases = make_grids.make_water_grid()
    for i in water_phases:
        names.append(i)
else:
    water_grid = []

Mantle_filename = perp.run_perplex(*[Mantle_wt_per,compositional_params,
                                            [structure_params.get('Pressure_range_mantle_LM'),structure_params.get('Temperature_range_mantle_LM'),
                                             structure_params.get('resolution_LM')],filename,verbose,False])
grids_high = make_grids.make_mantle_grid(Mantle_filename,Mantle_wt_per,False,use_grids)[0]

core_grid = make_grids.make_core_grid()

grids = [grids_low,grids_high,core_grid,water_grid]



'''JGS: I added this b/c some more massive planets (~5 Me) planets were getting a heat capacity error during the root solving.
Wendy suggested lowering the mantle potential temperature. Haven't tested it a ton so IDK if this is really needed
or not. Anyway, I replaced the sys.exit() associated with this error in MinPhys with a RuntimeError (not sure if this is the best error to raise).
If the code hits this RuntimeError, it iteratively lowers the max Fe/Mg value given to the root solver
until it no longer hits the error, and then asks if it can still find a solution within this new range.
This is called within calc_planet.
'''
def iter_root_scalar(minimum, maximum, mass, radius, step = -10):
    #Created to solve issue that it would break when encoutered SMs.
    try:
        FeMg = root_scalar(run_planet,bracket=[minimum, maximum] ,args=(mass, radius),x0 = 0.9,xtol=0.0001).root


    except ValueError:
        test = run_planet(minimum, *(mass, radius))

        if test < 0:
            #Even largest CMF produces radius too large, return very large Fe/Mg
            return 10**6, maximum

        if test > 0:
            #Even smallest CMF produces radius too small, return very small/negative Fe/Mg
            return -10**(6), maximum

    except RuntimeError:
        maximum = maximum + step
        FeMg, maximum = iter_root_scalar(minimum, maximum, mass, radius, step = step)

    return FeMg, maximum



def run_planet(x, *args):
    Mass_planet, Rad_planet = args

    #print("MR",Mass_planet,Rad_planet)
    grav = 6.67e-11*Mass_planet*5.97e24/(pow(Rad_planet*6371e3,2))
    compositional_params['FeMg'] = x
    Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac = functions.get_percents(compositional_params, verbose)

    Planet = functions.find_Planet_mass(Mass_planet, core_mass_frac,structure_params, compositional_params, grids, Core_wt_per, layers,verbose)
    #print(grav,Planet['gravity'][-1])

    out = 1. - (Planet['gravity'][-1] / grav)
    return out

def calc_planet(mass, radius):
    #print(mass,radius)
    den = mass * 5.97e21 / ((4 * np.pi / 3) * pow((radius * 6371e3),3))
    pure_rock = 3.012 + 1.208*pow(radius,1.6496)
    if den < pure_rock:
        #print("lower than pure rock")
        return(np.nan)

    if den < 8:
        minimum = 1e-9
        maximum = 20
    else:
        minimum = 1
        maximum = 70

    """
    JGS: The changes below call iter_root_scalar. If it finds a solution within the minimum and
    maximum Fe/Mg bounds, then the code runs exactly the same as it would have before my changes.
    If it does not find a solution and it gets the RuntimeError in MinPhys then it starts by
    lowering the max Fe/Mg by 10 until it doesn't get the error, then it increases Fe/Mg by 1 until
    it either reaches the previous max Fe/Mg that broke it or it finds a solution. This will certainly
    increase the run time for some planets, but will hopefully increase the number of useable samples.

    This is a very recent addition, so it certainly needs more testing.



    FeMg, maximum = iter_root_scalar(minimum, maximum, mass, radius, step = -10)
    mhold = maximum
    maximum = maximum + 10
    flag = 0

    while (FeMg == 10**6) and (maximum > mhold) and (flag == 0):
        m1 = maximum
        FeMg, maximum = iter_root_scalar(minimum, maximum, mass, radius, step = -1)
        m2 = maximum
        if m2-m1 == 0:
            flag = 1

    return FeMg
    """

    try:
        FeMg = root_scalar(run_planet,bracket=[minimum,maximum] ,args=(mass, radius),x0 = 0.9,xtol=0.0001).root
    except:
        test = run_planet(minimum, *(mass, radius))
        if test < 0:
            # planet with smallest core produces radius too big
            # return very high FeMg

            return (-10**(6))

        if test > 0:
            # planet with largest core produces radius too big
            # return very low FeMg
            return (10**6)
    else:
        return (FeMg)


'''JGS: just a general function for fitting for SkewNorm parameters. Not much to see here.'''
def sn(x, a, loc, scale):
    return sp.skewnorm.cdf(x, a, loc, scale)



'''JGS: This function can handle assymmetric mass and radius errors. The input for err can be a single scalar or 2D vector
of the form err = [upper_err, lower_err] -- aka [pl_radeerr1, pl_radeerr2] or [pl_bmasseerr1, pl_bmasseerr2] if
you are running samples from an NEA file.

If the errors are assymmetric, then it fits to a skewnorm CDF to get the skewnorm parameters
where x = [mu - 2 sigma lower, mu - 1 sigma lower, mu, mu + 1 sigma upper, mu + 2 sigma upper] and
y = what the corresponding CDF values should be = [0.0228, 0.1587, 0.5, 0.8413, 0.9772]. Otherwise, it just
samples from a normal dist.

I left some plots (commented out) so you can see what I'm doing if my description doesn't make sense.
'''
def sample_norm_or_skewnorm(mu, err, num_pts):
    if (type(err) == np.ndarray or type(err) == list):
        if (len(err)>1) and err[0] != err[1]:
            mx = np.array([mu - 2.0*abs(err[1]), mu - abs(err[1]), mu, mu + err[0], mu + 2.0*err[0]])
            y = np.array([0.5 - (0.3413+0.1359), 0.5 - 0.3413, 0.5, 0.5 + 0.3413, 0.5 + (0.3413+0.1359)])
            mparams, mcov = so.curve_fit(sn, mx, y, p0 = [(abs(err[0]) - abs(err[1]))/abs(abs(err[0]) - abs(err[1])), mu, err[0]])
            dist = sp.skewnorm.rvs(a = mparams[0], loc = mparams[1], scale = mparams[2], size = num_pts)

            #fig, ax = plt.subplots(1,2, figsize = (10,5))
            #ax[0].plot(mx, y, 'ko')
            #x = np.linspace(mu - 3.0*err[0], mu + 3.0*err[1], 1000)
            #ax[0].plot(x, sp.skewnorm.cdf(x, a = mparams[0], loc = mparams[1], scale = mparams[2]))
            #ax[0].set_ylabel('CDF')
            #ax[0].set_xlabel('x')
            #plt.show()
        else:
            dist = sp.norm.rvs(loc = mu, scale = abs(err[0]), size = num_pts)
            #plt.hist(dist, bins = 50)
            #plt.show()
    else:
        dist = sp.norm.rvs(loc = mu, scale = err, size = num_pts)
        #plt.hist(dist, bins = 50)
        #plt.show()


    return dist



'''JGS: Calculates the probability that the planet and host Fe/Mg distributions are actually the same given
measurement uncertainty. I broke the probability calculation in to seperate functions for Fe/Mg and CMF
b/c I found it to be much easier to work in log space for Fe/Mg. But, I can probably combine
them at some point.

Anyway, in both cases, I opted to just do a data integral instead of working with KDEs or trying to
fit some distribution to the data. For well-characterized super-Earths,
the KDE or distribution fitting should work just fine, But, neither work well for super-Mercuries, LDSPs, or planets
with huge M-R uncertainties -- i.e., planets with more than a handful of M-R samples that are more dense than a pure Fe
planet, less dense than pure sil, or break ExoPlex for another, unkown reason.
'''
def calc_PH0_FeMg(FeMg, FeMg_star):


    '''JGS: The first chuck of this code is function is just finding the bins to use for the data integral (and plotting),
    using the Freedman-Diaconis rule -- which will hopefully be more robust than just bins = some number. I calculate
    the bin width only using the 'good' samples -- i.e., those that fall between pure Fe and pure sil -- but I still
    account for the bad samples in the Fe/Mg PDF.
    '''

    indbl = np.where(FeMg == -10**6)
    indgood = np.where((FeMg != -10**(6)) & (FeMg != 10**(6)))

    FeMg[indbl] = 10**(-10)
    lFeMg = np.log10(FeMg)
    lFeMg_star = np.log10(FeMg_star)

    #Freedman–Diaconis rule
    a, b = np.percentile(lFeMg[indgood], [25, 75])
    IQR = b-a
    bw = 2.0*IQR/(len(lFeMg)**(1.0/3.0))

    hrange = [min(lFeMg[indgood]), max(lFeMg[indgood])]
    nbins = math.ceil((math.ceil(hrange[1]) - math.floor(hrange[0]))/bw)
    bins = np.array([math.floor(hrange[0]) + i*bw for i in range(0, nbins)])


    '''JGS: I get the bin width using the Feedman-Diaconis rule using the 'good' samples only first.
    Then I get the counts for every bin and normalize these using the entire sample to get lFeMg_planet_dens which
    is just the pdf histogram of lFeMg_planet accounting for the bad samples.'''
    counts, bins = np.histogram(lFeMg, bins = bins, range = hrange)
    lFeMg_planet_dens = counts/(num_pts*np.diff(bins))

    counts, bins = np.histogram(lFeMg_star, bins = bins)
    lFeMg_star_dens = counts/(num_pts*np.diff(bins))

    x = np.array([(bins[i] + bins[i+1])/2.0 for i in range(0, len(bins)-1)])

    '''JGS: this if/else statement just ensures I don't get negative Fe/Mg values. If, planet Fe/Mg
    is greater than host Fe/Mg, shift the stellar Fe/Mg dist up to the planet Fe/Mg and do the probability calculation.
    If the reverse is true, shift the planet Fe/Mg dist up to the stellar Fe/Mg dist and do the prob calc there.'''
    #change sv to difference in teh averages of the distributions
    planet = calc_planet(M,R)
    if planet < 0:
        planet = 0
    if planet > sFeMg:
        sv = planet - sFeMg
        FeMg_star_shift = FeMg_star + sv
        lFeMg_star_shift = np.log10(FeMg_star_shift)
        counts, bins = np.histogram(lFeMg_star_shift, bins = bins)
        lFeMg_star_shift_dens = counts/(num_pts*np.diff(bins))

        ph0 = si.simpson(lFeMg_planet_dens*lFeMg_star_dens, x)/si.simpson(lFeMg_star_shift_dens*lFeMg_planet_dens, x)

    else:
        sv = sFeMg - planet
        FeMg_planet_shift = FeMg + sv
        lFeMg_planet_shift = np.log10(FeMg_planet_shift)
        counts, bins = np.histogram(lFeMg_planet_shift, bins = bins)
        lFeMg_planet_shift_dens = counts/(num_pts*np.diff(bins))

        #print(lFeMg_planet_shift_dens)
        ph0 = si.simpson(lFeMg_planet_dens*lFeMg_star_dens, x)/si.simpson(lFeMg_star_dens*lFeMg_planet_shift_dens, x)



    return ph0, lFeMg_planet_dens, lFeMg_star_dens, bins, hrange



'''JGS: This is much the same as calc_PH0_FeMg. See comments above. The big difference is I don't need an if/else statement
to shift the distributions since I'm not working in log space for CMF. I can shift either distribution to the other to compute
the denominator of the probability calculation since negative values are allowed (well not physically, but mathematically it
doesn't make a difference).
'''
def calc_PH0_CMF(CMF, CMF_star):

    #Freedman–Diaconis rule
    ind = np.where((CMF != -10**(6)) & (CMF != 10**(6)))
    a, b = np.percentile(CMF[ind], [25, 75])
    IQR = b-a
    bw = 2.0*IQR/(len(CMF)**(1.0/3.0))

    ind = np.where(abs(CMF) != 10**(6))
    hrange = [min(CMF[ind]), max(CMF[ind])]
    nbins = math.ceil((math.ceil(hrange[1]) - math.floor(hrange[0]))/bw)
    bins = np.array([math.floor(hrange[0]) + i*bw for i in range(0, nbins+10)])

    counts, bins = np.histogram(CMF, bins = bins)
    CMF_planet_dens = counts/(num_pts*np.diff(bins))

    counts, bins = np.histogram(CMF_star, bins = bins)
    CMF_star_dens = counts/(num_pts*np.diff(bins))


    compositional_params['FeMg'] = calc_planet(M, R)
    CMF_cent = functions.get_percents(compositional_params,verbose)[3]

    compositional_params['FeMg'] = sFeMg
    CMF_cent_star = functions.get_percents(compositional_params,verbose)[3]

    sv = CMF_cent - CMF_cent_star
    CMF_star_shift = np.array(CMF_star) + sv

    counts, bins = np.histogram(CMF_star_shift, bins = bins)
    CMF_star_shift_dens = counts/(num_pts*np.diff(bins))

    x = np.array([(bins[i] + bins[i+1])/2.0 for i in range(0, len(bins)-1)])


    ph0 = si.simpson(CMF_star_dens*CMF_planet_dens, x)/si.simpson(CMF_star_shift_dens*CMF_planet_dens, x)



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

'''JGS: This function just finds the 2sigma and 1sigma data ranges for Fe/Mg and CMF'''
def find_drange(data):
    cent = 50
    sigl, sigu = cent - (34.13), cent + 34.13
    tsigl, tsigu = cent - (34.13+13.59), cent + (34.13+13.59)

    dranges = np.percentile(data, [tsigl, sigl, cent, sigu, tsigu])


    return dranges

if __name__ == "__main__":
    import multiprocessing as mp
    num_pts = 100

    filename = 'TOI_ww'
    R=  1.67
    R_err = [0.07, 0.07]
    M = 4.8
    M_err = [1.3, 1.3]

    #stellar Fe/Mg
    sFeMg =  1.48*(24.305/55.845)
    sFeMg_err = [.36*(24.305/55.845),.29*(24.305/55.845)]

    Mass_planet = sample_norm_or_skewnorm(M, M_err, num_pts)
    Radius_planet = sample_norm_or_skewnorm(R, R_err, num_pts)


    #FeMg_star = sp.norm.rvs(loc = sFeMg, scale = sFeMg_err, size = num_pts)
    FeMg_star = sample_norm_or_skewnorm(sFeMg, sFeMg_err, num_pts)
    indup = np.where(FeMg_star < 0)[0]
    if len(indup) > 0:
        FeMg_star = np.delete(FeMg_star, indup)

    test = calc_planet(Mass_planet[0],Radius_planet[0])

    vals = zip(Mass_planet, Radius_planet)
    pool = mp.Pool(processes=mp.cpu_count())

    #pool = mp.Pool(processes=1)

    FeMg = pool.starmap_async(calc_planet,vals).get()

    #for i in range(len(FeMg)):
    #    print(Mass_planet[i], Radius_planet[i],FeMg[i])
    pool.close()

    CMF = []

    for i in range(len(FeMg)):
        if FeMg[i] >0:
            compositional_params['FeMg'] = FeMg[i]
            CMF.append(functions.get_percents(compositional_params,verbose)[3])
        else:
            CMF.append(FeMg[i])

    mu= statistics.mean(FeMg)
    std = statistics.stdev(FeMg)
    print()
    print("Actual Fe/Mg:")
    print(round(mu, 2), "+/-", round(std, 2))
    print()

    #mu, std = sp.lognorm.fit(FeMg)
    shape, loc, scale = sp.lognorm.fit(FeMg)

    mu = np.exp(np.log(scale) + pow(shape,2)/2)
    std = mu*np.sqrt(pow(shape,2)-1)
    print()
    print("Log-Normal best fit Fe/Mg:")
    print(round(mu,2), "+/-",round(std,2))
    print()

    mu_CMF, std_CMF = sp.norm.fit(CMF)
    print()
    print("Normal best fit CMF:")
    print(round(mu_CMF,2), "+/-",round(std_CMF,2))
    print()

    output = list(zip(Mass_planet, Radius_planet,FeMg,CMF))

    header ='#Mass,Radius,FeMg,CMF'
    head = []
    names = header.split(',')

    for i in names:
        head.append(i+'_'+filename)
    header = ','.join(head)

    np.savetxt(filename+'.csv', output, header = header, comments='#', delimiter=',')

    fig, (ax1, ax2,ax3) = plt.subplots(1, 3, figsize=(10, 5))


    ax1.scatter(Mass_planet, Radius_planet)
    ax1.set_ylabel('Radius (Earth Radii)', size=20)
    ax1.set_xlabel('Mass (Earth Masses)', size=20)


    #abs doesn't work for lists. try to build a work around cause right now if lower uncertainty is larger and input as a negative number
    #the axis limits may be a bit wonky. Pretty low priority item i think

    x = np.linspace(M-5*max(M_err), M+5*max(M_err), num = 10)
    R_high = -.48+1.54*pow(x,0.2)

    R_low = -.35+1.32*pow(x,0.21)

    ax1.fill_between(x,R_low,R_high, alpha=0.4, label= 'NRPZ')
    ax1.legend(loc='lower right',fontsize = 'large')
    ax1.set_xlim(M-5*max(M_err), M+5*max(M_err))
    ax1.set_ylim(R-5*max(R_err), R+5*max(R_err))
    ax1.minorticks_on()

    FeMg = np.array(FeMg)

    #calc ph0 using Fe/Mg
    ph0_FeMg, lFeMg_planet_dens, lFeMg_star_dens, bins, hrange  = calc_PH0_FeMg(FeMg, FeMg_star)
    FeMg_range = find_drange(FeMg)


    sigFeMgu = FeMg_range[3] - FeMg_range[2]
    sigFeMgl = FeMg_range[2] - FeMg_range[1]

    print('Fe/Mg: ' + str(round(FeMg_range[2],2)) + ' +' + str(round(sigFeMgu,2)) + ' -' + str(round(sigFeMgl,2)))
    print('-- 2-sig data range: [' + str(round(FeMg_range[0],2)) + ' , ' + str(round(FeMg_range[4],2)) + ']')
    print('-- 1-sig data range: [' + str(round(FeMg_range[1],2)) + ' , ' + str(round(FeMg_range[3],2)) + ']')
    print('P(H0) from Fe/Mg: ', round(100*ph0_FeMg, 2))



    x = np.linspace(hrange[0], hrange[1], 1000)

    ax2.bar(bins[:-1], height = lFeMg_planet_dens, width = np.diff(bins), align = 'edge', alpha = 0.6, color = 'c')
    ax2.bar(bins[:-1], height = lFeMg_star_dens, width = np.diff(bins), align = 'edge', alpha = 0.6, color = 'mediumorchid')

    ax2.set_xlim(hrange[0], hrange[1])

    ax2.set_xlabel(r'log$_{10}$ Fe/Mg', size=20)
    """
    CMF = np.array(CMF)

    #calc ph0 using Fe/Mg
    ph0_CMF, CMF_planet_dens, CMF_star_dens, bins, hrange  = calc_PH0_CMF(CMF, CMF_star)
    CMF_range = find_drange(CMF)


    sigCMFu = CMF_range[3] - CMF_range[2]
    sigCMFl = CMF_range[2] - CMF_range[1]

    print('CMF: ' + str(round(CMF_range[2],2)) + ' +' + str(round(sigCMFu,2)) + ' -' + str(round(sigCMFl,2)))
    print('-- 2-sig data range: [' + str(round(CMF_range[0],2)) + ' , ' + str(round(CMF_range[4],2)) + ']')
    print('-- 1-sig data range: [' + str(round(CMF_range[1],2)) + ' , ' + str(round(CMF_range[3],2)) + ']')
    print('P(H0) from CMF: ', round(100*ph0_CMF, 2))




    ax3.bar(bins[:-1], height = CMF_planet_dens, width = np.diff(bins), align = 'edge', alpha = 0.6, color = 'c')
    ax3.bar(bins[:-1], height = CMF_star_dens, width = np.diff(bins), align = 'edge', alpha = 0.6, color = 'mediumorchid')

    ax3.set_xlim(0, 1)
    

    ax3.set_xlabel(r'Core Mass Fraction', size=20)
    """




    plt.savefig(filename+'.jpg', format = 'jpg')
