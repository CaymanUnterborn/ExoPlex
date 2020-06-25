import os
import sys
import pexpect as pe

#PerPlex_path = os.path.dirname(os.path.realpath(__file__))+"/PerPlex"

# hack to allow scripts to be placed in subdirectories next to ExoPlex:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))

def run_perplex(*args):
    PerPlex_path = os.path.dirname(os.path.abspath(__file__))+ '/PerPlex'

    Mantle_wt_per = args[0]

    Mantle_wt_per['Al2O3'] = 3.96
    Mantle_wt_per['SiO2'] = 46.3
    Mantle_wt_per['MgO'] = 38.4
    Mantle_wt_per['CaO'] = 3.2
    Mantle_wt_per['FeO'] = 7.7




    FeMg = args[1][1]
    SiMg = args[1][2]
    CaMg = args[1][3]
    AlMg = args[1][4]
    mol_frac_Fe_mantle = args[1][5]
    mol_frac_Fe_mantle = args[1][5]
    use_grids = args[1][-1]

    wt_frac_Si_core = args[1][6]

    Pressure_range_mantle = args[2][0]
    Temperature_range_mantle = args[2][1]
    resolution = args[2][2]

    filename = args[3]
    UMLM = args[4]

    plxMan = str(Mantle_wt_per.get('MgO')) + ' ' + str(Mantle_wt_per.get('SiO2')) + ' ' \
             + str(Mantle_wt_per.get('FeO')) + ' ' + str(Mantle_wt_per.get('CaO')) \
             + ' ' + str(Mantle_wt_per.get('Al2O3'))+ ' ' + str(0.) #last value included for Na

    solfileparamsString0 = '_' + str(round(SiMg, 3)) + '_' + str(round(FeMg, 3)) + '_' + str(
        round(CaMg, 3)) + '_' + str(round(AlMg, 3)) \
                           + '_' + str(round(mol_frac_Fe_mantle, 3)) + '_' + str(round(wt_frac_Si_core, 3))

    # changes periods to commas
    solfileparamsString = solfileparamsString0.replace('.', ',')

    solutionFileNameMan = 'SiMg_FeMg_CaMg_AlMg_XFeO_fSic' + solfileparamsString + '_MANTLE'
    solutionFileNameMan_short = list(solutionFileNameMan)
    solutionFileNameMan_short[0:30] = []

    solutionFileNameMan = "".join(solutionFileNameMan_short)

    if use_grids == False:
        filename = solutionFileNameMan
    #print "file",filename+'_UM_results.txt'


    if os.path.isfile('../Solutions_Small/'+filename+'_UM_results.txt') and UMLM == True and use_grids==True:
        print ('The Upper mantle .tab already exists, please wait briefly for solution\n')
        return '../Solutions_Small/' + filename

    if os.path.isfile('../Solutions_Small/'+filename+'_LM_results.txt') and UMLM == False and use_grids==True:
        print ('The Lower mantle .tab already exists, please wait briefly for solution\n')
        return '../Solutions_Small/' + filename

    else:


        if os.path.isfile('../Calc_Solutions/'+solutionFileNameMan+'_UM_results.txt') == True and UMLM == True:
            filename = solutionFileNameMan
            print ('The Upper mantle .tab already exists, please wait briefly for solution\n')
            return '../Calc_Solutions/' + filename

        if os.path.isfile('../Calc_Solutions/'+solutionFileNameMan+'_LM_results.txt') == True and UMLM == False:
            filename = solutionFileNameMan

            print ('The Lower mantle .tab already exists, please wait briefly for solution\n')
            return '../Calc_Solutions/' + filename
        else:


            if  UMLM == True:
                print ('Making upper mantle PerPlex phase file. \n This will be stored in: ../Calc_Solutions/'+ filename+'_UM_results.txt')
            else:

                print ('Making lower mantle PerPlex phase file. \n This will be stored in: ../Calc_Solutions/'+ filename+'_LM_results.txt')


        #we need to shorten the file name for PerPlex to accept it


    # define perplex inputs in terms of components, this is for legibility
    #print "If you see this error, something has gone wrong and PerPlex needs to run. Let's avoid that by exiting."
    #sys.exit()
    component1 = 'MGO'
    component2 = 'SIO2'
    component3 = 'FEO'
    component4 = 'CAO'
    component5 = 'AL2O3'
    component6 = 'NA2O'
    if os.path.isfile(solutionFileNameMan+'.arf')==True:
        os.remove(solutionFileNameMan+'.arf')
        os.remove(solutionFileNameMan+'.blk')
        os.remove(solutionFileNameMan+'.dat')
        os.remove(solutionFileNameMan+'.plt')
        os.remove(solutionFileNameMan+'.tof')

        os.remove(solutionFileNameMan + '_VERTEX_options.txt')
        os.remove(solutionFileNameMan + '_WERAMI_options.txt')
        os.remove(solutionFileNameMan + '_auto_refine.txt')

    p = pe.spawn(PerPlex_path+"/./build")


    p.sendline(solutionFileNameMan)
    #p.sendline('../ExoPlex_dist/PerPlex/stx11ver.dat')
    p.sendline(PerPlex_path + '/stx11ver.dat')

    p.sendline(PerPlex_path + '/perplex_options.dat')
    # Transform them (Y/N)?
    p.sendline('N')
    # Calculations with saturated components (Y/N)?
    p.sendline('N')
    # Use chemical potentials, activities or fugacities as independent variables (Y/N)?
    p.sendline('N')
    # Select thermodynamic components from the set:
    p.sendline(component1)  # MGO
    p.sendline(component2)  # SIO2
    p.sendline(component3)  # FEO
    p.sendline(component4)  # CAO
    p.sendline(component5)  # AL2O3
    p.sendline(component6)  # NA2O
    p.sendline('')
    # Specify computational mode:

    p.sendline('2')

    # Make one dependent on the other, e.g., as along a geothermal gradient (y/n)?
    p.sendline('N')

    # Select x-axis variable:
    p.sendline('2')
    # Enter minimum and maximum values, respectively, for: T(K)
    p.sendline(Temperature_range_mantle)
    # Enter minimum and maximum values, respectively, for: P(bar)
    # P(Pa) = P(bar)*100000
    p.sendline(Pressure_range_mantle)
    # Specify component amounts by weight (Y/N)?
    p.sendline('Y')
    # Enter weight amounts of the components:
    # MGO SIO2 FEO CAO AL2O3
    # for the bulk composition of interest:
    # NOTE*** This is a wt%
    p.sendline(plxMan)
    # Output a print file (Y/N)?
    p.sendline('N')
    # Exclude pure and/or endmember phases (Y/N)?
    p.sendline('N')

    # Include solution models (Y/N)?
    p.sendline('Y')
    p.sendline(PerPlex_path + '/stx11_solution_model.dat')
    p.sendline('C2/c') #C2C Phase of clinopyroxene
    p.sendline('Wus')
    p.sendline('Pv')
    p.sendline('Pl')
    p.sendline('Sp')
    p.sendline('O')
    p.sendline('Wad')
    p.sendline('Ring')
    p.sendline('Opx')
    p.sendline('Cpx')
    p.sendline('Aki')
    p.sendline('Gt_maj') #kitchen sink
    p.sendline('Ppv')
    p.sendline('CF')
    p.sendline('')
    # Enter calculation title:

    p.sendline(solutionFileNameMan + 'calc')

    p.logfile = open('build.log','wb')
    p.read()
    p.wait()

    print ("Done with Build, moving on to Vertex")

    # Spawn Vertex ----------------#
    # Enter the project name (the name assigned in BUILD) [default = my_project]:
    p = pe.spawn(PerPlex_path+"/./vertex",timeout=2400)

    p.sendline(solutionFileNameMan)

    #p.expect('$$$$',timeout=None)
    p.logfile = open('vertex.log', 'wb')
    p.read()
    p.wait()

    print ('Finished with Vertex, beginning Werami')

    p = pe.spawn(PerPlex_path+"/./werami",timeout=None)


    p.sendline(solutionFileNameMan)
    # select 2D grid
    p.sendline('2')
    # Below, select parameters density, alpha, cp.
    # Ns for no calculating individual phase properties
    p.sendline('2')
    p.sendline('N')
    p.sendline('12')
    p.sendline('N')
    p.sendline('13')
    p.sendline('N')
    p.sendline('14')
    p.sendline('N')
    p.sendline('4')
    p.sendline('N')
    p.sendline('19')
    p.sendline('N')
    ####### the next lines will pass requests to perplex to print phases and their proportions into the .tab file

    # 21 species, in all for Fe-Si-Mg-O regime
    p.sendline('7')
    p.sendline('C2/c')  # 0
    p.sendline('7')
    p.sendline('Wus')  # 1
    p.sendline('7')
    p.sendline('Pv')  # 2
    p.sendline('7')
    p.sendline('an')  # 3
    p.sendline('7')
    p.sendline('Sp')  #4
    p.sendline('7')
    p.sendline('O')  # 4
    p.sendline('7')
    p.sendline('Wad')  # 5
    p.sendline('7')
    p.sendline('Ring')  # 6  #if statement about no FeO or some shit
    p.sendline('7')
    p.sendline('Opx')  # 7
    p.sendline('7')
    p.sendline('Cpx')  # 8
    p.sendline('7')
    p.sendline('Aki')  # 9
    p.sendline('7')
    p.sendline('Gt_maj')  # 10
    p.sendline('7')
    p.sendline('Ppv')  # 11
    p.sendline('7')
    p.sendline('CF')   # 12
    p.sendline('7')
    p.sendline('st')  # 12
    p.sendline('7')
    p.sendline('q')  # 13
    p.sendline('7')
    p.sendline('ca-pv')  # 14
    p.sendline('7')
    p.sendline('cfs')  # 15
    p.sendline('7')
    p.sendline('coe')  # 16
    p.sendline('7')
    p.sendline('ky')  # 17
    p.sendline('7')
    p.sendline('seif')  # 18

    # exit parameter choosing

    p.sendline('0')
    # Change default variable range (y/n)?
    p.sendline('N')

    # Enter number of nodes in the T(K)     and P(bar)   directions:

    p.sendline(resolution)
    p.logfile = open('werami.log','wb')
    p.expect('    0 - EXIT', timeout=None)
    p.sendline('0')
    p.read()
    p.terminate()
    print ("Done with PerPlex")

    if UMLM == True:
        os.rename(solutionFileNameMan+'_1.tab','../Calc_Solutions/'+filename+'_UM_results.txt')
    else:
        os.rename(solutionFileNameMan + '_1.tab', '../Calc_Solutions/' + filename + '_LM_results.txt')

    os.remove(solutionFileNameMan+'.arf')
    os.remove(solutionFileNameMan+'.blk')
    os.remove(solutionFileNameMan+'.dat')
    os.remove(solutionFileNameMan+'.plt')
    os.remove(solutionFileNameMan+'.tof')

    os.remove(solutionFileNameMan+'_VERTEX_options.txt')
    os.remove(solutionFileNameMan+'_WERAMI_options.txt')
    os.remove(solutionFileNameMan+'_auto_refine.txt')

    filename = '../Calc_Solutions/'+filename

    return filename

