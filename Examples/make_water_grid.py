import os
import sys
import numpy as np
# hack to allow scripts to be placed in subdirectories next to exoplex:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))

import seafreeze as sf
from ExoPlex import burnman
def ice_VII_vals(P,T):
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
    den_noT = rock.evaluate(['density'], P *1e6, 300)[0]
    corr = np.exp((T- 300) * 11.58e-5)
    density = den_noT/corr
    Cp = 1000*(3.3 + 22.1 * np.exp(-0.058 * P/1000.))  # Asahara 2010
    Ks = 23.9
    Ksp = 4.2
    a0 = -3.9e-4
    a1 = 1.5e-6
    at = a0 + a1 * T
    alpha = at * (1. + (Ksp / Ks) * P/1000) ** (-0.9)  # fei 1993)
    return (density, Cp, alpha)

if __name__ == "__main__":

    max_P = 50*1000. #50 GPa
    min_P = 0.1 #1 bar

    max_T = 500. #K
    min_T = 300.
    num_P_points = 200
    num_T_points = 100


    output = []
    P = np.arange(min_P,max_P,(max_P-min_P)/num_P_points)
    T = np.arange(min_T,max_T,(max_T-min_T)/num_T_points)
    PT = np.empty(((len(P)*len(T)),), dtype=object)

    counter = 0
    for i in range(len(P)):
        for j in range(len(T)):
            PT[counter] = P[i],T[j]
            counter+=1
    out = sf.whichphase(PT)

    print("out",out)
    counter = 0
    from collections import Counter
    outs = []

    outs = list(Counter(out).keys())
    names = []
    for i in outs:
        if i < 7:
            names.append(sf.phasenum2phase[i])
        else:
            names.append('ice_VII')

    name_for_index = list(Counter(names).keys())

    names_index = [i for i in range(len(name_for_index))]
    len_append = len(name_for_index)

    print(names_index)
    print("name",name_for_index)

    appends = []
    for i in range(len(name_for_index)):
        dummy = [0 for i in range(len(name_for_index))]
        dummy[names_index[i]] = 1
        appends.append(dummy)

    print("dummy",appends)
    for i in range(len(P)):
        for j in range(len(T)):
            if out[counter] < 7:
                phase = sf.phasenum2phase[out[counter]]
                for k in range(len(name_for_index)):
                    if phase == name_for_index[k]:
                        phase_append = appends[k]

                        break
                new = sf.seafreeze(np.array(PT[counter]), phase)

                to_go = [P[i]*10,T[j],new.rho[0][0],new.Cp[0][0],np.log10(new.alpha[0][0])]
                for l in phase_append:
                    to_go.append(l)
                output.append(to_go)
                counter+=1

            else:
                phase = 'Ice_VII'
                for k in range(len(name_for_index)):
                    if phase == name_for_index[k]:
                        phase_append = appends[k]
                        break
                density, Cp, alpha = ice_VII_vals(P[i],T[j])
                to_go = [P[i] * 10, T[j], density, Cp, np.log10(alpha)]
                for l in phase_append:
                    to_go.append(l)
                output.append(to_go)
                counter+=1

    for i in range(len(name_for_index)):
        if name_for_index[i] == "water1":
            name_for_index[i] = 'liq_water'
        elif name_for_index[i] == 'VI':
            name_for_index[i] = 'ice_VI'
        elif name_for_index[i] != 'ice_VII':
            print("new phase"), name_for_index[i]
            sys.exit()

    head = ['P','T','density','Cp','alpha']
    for i in name_for_index:
        head.append(i)
    header = ','.join(head)
    format = '%.3f,%.1f,%.5f,%.5f,%.5f,%i,%i,%i'
    #format = '%s'
    np.savetxt('water_grid.dat', output, fmt=format, delimiter=",", newline='\n',
                header=header, footer='', comments='# ')