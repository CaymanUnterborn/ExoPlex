import os
import sys
import pexpect as pe


# hack to allow scripts to be placed in subdirectories next to ExoPlex:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))
PerPlex_path = os.path.dirname(os.path.abspath(__file__))+ '/PerPlex'

def run_perplex(*args):
    """
   This module runs PerPlex to produce mantle phase diagrams for a custom composition
    if the user opts to not use the premade grids.

    Parameters
    ----------
    Mantle_wt_per : dictionary
        composition of mantle in oxides :math:`[wt\%]`

    compositional_params: list
        Structural parameters of the planet; See example for description

    structural_params: list
        Structural parameters of the planet; See example for description

    filename: string
       chosen filename for output file

    UMLM: boolean
        True if creating grid for upper mantle, false if lower mantle

    Returns
    -------
    Phase diagram: file
        Mantle phase diagram for the chosen composition. Contains P, T, expansivity, density and specific heat.
        Stored in /Calc_Solutions/'+filename, where filename is the chosen user name
    """
    compositional_params = args[1]
    verbose = args[4]
    UMLM = args[5]
    use_grids = compositional_params.get('use_grids')

    FeMg = compositional_params.get('FeMg')
    SiMg = compositional_params.get('SiMg')
    CaMg = compositional_params.get('CaMg')
    AlMg = compositional_params.get('AlMg')
    wt_frac_FeO_wanted = compositional_params.get('wt_frac_FeO_wanted')

    solfileparamsString0 = '_' + str(round(SiMg, 3)) + '_' + str(round(FeMg, 3)) + '_' + str(
        round(CaMg, 3)) + '_' + str(round(AlMg, 3)) \
                           + '_' + str(round(wt_frac_FeO_wanted, 3))

    # changes periods to commas
    solfileparamsString = solfileparamsString0.replace('.', ',')

    solutionFileNameMan = 'SiMg_FeMg_CaMg_AlMg_XFeO_fSic' + solfileparamsString + '_MANTLE'
    solutionFileNameMan_short = list(solutionFileNameMan)
    solutionFileNameMan_short[0:30] = []

    solutionFileNameMan = "".join(solutionFileNameMan_short)

    if use_grids == False:
        if os.path.isfile('../Calc_Solutions/' + solutionFileNameMan + '_UM_results.txt') == True and UMLM == True:
            if verbose == True:
                print('The Upper mantle .tab already exists, please wait briefly for solution\n')
            return '../Calc_Solutions/' + solutionFileNameMan

        if os.path.isfile('../Calc_Solutions/' + solutionFileNameMan + '_LM_results.txt') == True and UMLM == False:
            if verbose == True:
                print('The Lower mantle .tab already exists, please wait briefly for solution\n')
            return '../Calc_Solutions/' + solutionFileNameMan

    else:
        filename = args[3]

        check_FeO = float(filename.split('_')[-1].split('Fe')[0])
        if check_FeO <= 0.2 and check_FeO >0:
            test = filename.split('_')
            test[-1] = '0.15Fe'
            filename = '_'.join(test)

        if os.path.isfile('../Solutions_Small/'+filename+'_UM_results.txt') and UMLM == True:
            if verbose == True:
                print ('The Upper mantle .tab already exists\n')
            return '../Solutions_Small/' + filename

        if os.path.isfile('../Solutions_Small/'+filename+'_LM_results.txt') and UMLM == False:
            if verbose == True:
                print ('The Lower mantle .tab already exists\n')
            return '../Solutions_Small/' + filename

        else:
            print("You have set use_grids to True but your composition does not exist within the grids.")
            print("Try again with use_grids = False")
            sys.exit()

    Mantle_wt_per = args[0]
    structure_params = args[2]
    if UMLM == True:
        print("Running PerPlex for Upper Mantle:")
        Pressure_range_mantle = structure_params.get('Pressure_range_mantle_UM')
        Temperature_range_mantle = structure_params.get('Temperature_range_mantle_UM')
        resolution = structure_params.get('resolution_UM')
    else:
        print("Running PerPlex for Lower Mantle:")
        Pressure_range_mantle = structure_params.get('Pressure_range_mantle_LM')
        Temperature_range_mantle = structure_params.get('Temperature_range_mantle_LM')
        resolution = structure_params.get('resolution_LM')




    plxMan = str(Mantle_wt_per.get('MgO')) + ' ' + str(Mantle_wt_per.get('SiO2')) + ' ' \
             + str(Mantle_wt_per.get('FeO')) + ' ' + str(Mantle_wt_per.get('CaO')) \
             + ' ' + str(Mantle_wt_per.get('Al2O3')) + ' ' + str(0.)  # last value included for Na

    component1 = 'MGO'
    component2 = 'SIO2'
    component3 = 'FEO'
    component4 = 'CAO'
    component5 = 'AL2O3'
    component6 = 'NA2O'

    if os.path.isfile(solutionFileNameMan+'.arf')==True:
        os.remove(solutionFileNameMan + '.arf')
        os.remove(solutionFileNameMan + '.blk')
        os.remove(solutionFileNameMan + '.dat')
        os.remove(solutionFileNameMan + '.plt')
        os.remove(solutionFileNameMan + '.tof')
        os.remove(solutionFileNameMan + '_seismic_data.txt')
        os.remove(solutionFileNameMan + '_auto_refine.txt')

    p = pe.spawn(PerPlex_path+"/./build")


    p.sendline(solutionFileNameMan)
    p.sendline(PerPlex_path + '/stx11ver.dat')

    p.sendline(PerPlex_path + '/perplex_options.dat')
    # Transform them (Y/N)?
    p.sendline('N')
    #Specify computational mode:
    p.sendline('2')

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
    if verbose == True:
        print ("Done with Build, moving on to Vertex")

    # Spawn Vertex ----------------#
    # Enter the project name (the name assigned in BUILD) [default = my_project]:
    p = pe.spawn(PerPlex_path+"/./vertex",timeout=2400)

    p.sendline(solutionFileNameMan)

    #p.expect('$$$$',timeout=None)
    p.logfile = open('vertex.log', 'wb')
    p.read()
    p.wait()

    if verbose == True:
        print ('Finished with Vertex, beginning Werami')

    p = pe.spawn(PerPlex_path+"/./werami",timeout=None)



    p.sendline(solutionFileNameMan)
    # select 2D grid
    p.sendline('2')
    # Below, select parameters density, alpha, cp.
    # Ns for no calculating individual phase properties
    p.sendline('2')
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
    p.sendline('Ring')  # 6
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
    p.sendline('Y')
    p.sendline(Temperature_range_mantle)
    p.sendline(Pressure_range_mantle)
    # Enter number of nodes in the T(K)     and P(bar)   directions:

    p.sendline(resolution)
    p.logfile = open('werami.log','wb')
    p.expect('    0 - EXIT', timeout=None)
    p.sendline('0')
    p.read()
    p.terminate()
    if verbose == True:
        print ("Done with PerPlex")

    if UMLM == True:
        os.rename(solutionFileNameMan+'_1.tab','../Calc_Solutions/'+solutionFileNameMan+'_UM_results.txt')
    else:
        os.rename(solutionFileNameMan + '_1.tab', '../Calc_Solutions/' + solutionFileNameMan + '_LM_results.txt')

    os.remove(solutionFileNameMan+'.arf')
    os.remove(solutionFileNameMan+'.blk')
    os.remove(solutionFileNameMan+'.dat')
    os.remove(solutionFileNameMan+'.plt')
    os.remove(solutionFileNameMan+'.tof')
    os.remove(solutionFileNameMan+'_seismic_data.txt')
    os.remove(solutionFileNameMan+'_auto_refine.txt')

    filename = '../Calc_Solutions/'+solutionFileNameMan

    return filename

