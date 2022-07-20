import os
import sys
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))
import numpy as np
import ExoPlex.burnman as bm
import math
if __name__ == "__main__":

    filename_out = 'new_liquid_iron_grid'
    header = 'P_bar,T_K,rho_kgm3,alpha_1_K,Cp_J_K_mol'
    rock = bm.minerals.other.liquid_iron()
    output = []
    counter = 0
    """
    Pressure = [(10+i*5.)*1e9 for i in range(50)]
    Temperature = [1600+40*i for i in range(50)]


    output = []
    counter = 0
    counter_full = 2* len(Pressure)
    for i in Pressure:
        print("on P",counter," of ",counter_full)
        for j in Temperature:
            rho, alpha, cp = rock.evaluate(['density','thermal_expansivity','molar_heat_capacity_p'],i,j)
            output.append([i*1e-5,j,rho,math.log10(alpha),cp*(1000./55.845)])


        if (counter+1)%5 ==0:
            np.savetxt(filename_out + '.dat', output, delimiter=",", newline='\n',
                   header=header, footer='', comments='# ')

        counter+=1
    
    
    Pressure = [(1+i*1)*1e9 for i in range(50)]
    Temperature = [1600+40*i for i in range(125)]

    counter = 0
    counter_full = len(Pressure)
    for i in Pressure:
        print("on P",counter," of ",counter_full)
        for j in Temperature:

            rho, alpha, cp = rock.evaluate(['density','thermal_expansivity','molar_heat_capacity_p'],i,j)
            output.append([i*1.e-5,j,rho,math.log10(alpha),cp*(1000./55.845)])


        if (counter+1)%5 ==0:
            np.savetxt(filename_out + '.dat', output, delimiter=",", newline='\n',
                   header=header, footer='', comments='# ')

        counter+=1
    """
    Pressure = [(10+i*75)*1e9 for i in range(200)]
    Temperature = [1700+100*i for i in range(88)]

    counter = 0


    counter_full = len(Pressure)
    for i in Pressure:
        print("on P",counter," of ",counter_full)
        for j in Temperature:
            #print(i/1e9,j)
            rho, alpha, cp = rock.evaluate(['density','thermal_expansivity','molar_heat_capacity_p'],i,j)
            output.append([i*1.e-5,j,rho,math.log10(alpha),cp*(1000./55.845)])


        if (counter+1)%5 ==0:
            np.savetxt(filename_out + '.dat', output, delimiter=",", newline='\n', fmt = '%.5e',
                   header=header, footer='', comments='# ')

        counter+=1
