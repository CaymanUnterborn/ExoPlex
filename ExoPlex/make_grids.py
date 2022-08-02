import sys

import numpy as np
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator


#constants, atomic masses
mFe = 55.845
mMg = 24.306
mSi = 28.0867
mO = 15.9994
mS = 32.0650
mCa = 40.078
mAl = 26.981
range_FeO = np.array([0., .02, .04, .06, .08, .1, .15, .20])


def make_core_grid():

    file = open('../Solutions_Small/liquid_iron_grid.dat')

    temp_file = file.readlines()
    num_rows = len(temp_file[1:])
    num_columns = len(temp_file[12].split(','))

    for i in temp_file[1:]:
        if i[0] == '#':
            num_rows = num_rows-1
            start+=1

    start = 1

    data = temp_file[start:]
    grid = np.zeros((num_rows, num_columns))

    for i in range(num_rows):
        # for j in range(num_columns):
        columns = data[i].strip('\n').split(',')
        grid[i] = [float(j) for j in columns]


    pressure_grid = np.array([row[0] for row in grid])
    temperature_grid = np.array([row[1] for row in grid])
    density_grid = np.array([row[2] for row in grid])
    alpha_grid = [pow(10,row[3]) for row in grid]
    cp_grid = [row[4] for row in grid]

    PT = np.vstack((pressure_grid, temperature_grid)).T
    tri_PT = Delaunay(PT)  # Compute the triangulation
    interpolator_rho = LinearNDInterpolator(tri_PT, np.array(density_grid))
    interpolator_alpha = LinearNDInterpolator(tri_PT, np.array(alpha_grid))
    interpolator_CP = LinearNDInterpolator(tri_PT, np.array(cp_grid))

    keys = ['density','alpha','cp']
    return dict(zip(keys,[interpolator_rho,interpolator_alpha,interpolator_CP]))

def make_water_grid():

    file = open('../Solutions_Small/water_grid.dat')
    temp_file = file.readlines()
    num_rows = len(temp_file[1:])
    num_columns = len(temp_file[12].split(','))

    for i in temp_file[1:]:
        if i[0] == '#':
            num_rows = num_rows-1
            start+=1

    start = 1
    phase_grid = []


    data = temp_file[start:]
    header = temp_file[0].strip('\n').split(',')
    Phases = header[5:]


    grid = np.zeros((num_rows, num_columns))
    for i in range(num_rows):
        # for j in range(num_columns):
        columns = data[i].strip('\n').split(',')
        grid[i] = [float(j) for j in columns]

    pressure_grid = np.array([row[0] for row in grid])
    temperature_grid = np.array([row[1] for row in grid])
    density_grid = np.array([row[2] for row in grid])
    cp_grid = np.array([row[3] for row in grid])
    alpha_grid = np.array([row[4] for row in grid])
    phase_grid =  np.array([100*row[5:] for row in grid])

    PT = np.vstack((pressure_grid, temperature_grid)).T
    tri_PT = Delaunay(PT)  # Compute the triangulation
    interpolator_rho = LinearNDInterpolator(tri_PT, np.array(density_grid))
    interpolator_alpha = LinearNDInterpolator(tri_PT, np.array(alpha_grid))
    interpolator_CP = LinearNDInterpolator(tri_PT, np.array(cp_grid))
    interpolator_phases = LinearNDInterpolator(tri_PT, np.array(phase_grid))

    keys = ['density', 'alpha', 'cp', 'phases']

    return dict(zip(keys, [interpolator_rho, interpolator_alpha, interpolator_CP, interpolator_phases])),Phases

def make_mantle_feo_grid(Mantle_filename,Mantle_wt_per,UMLM):

    test = Mantle_filename.split('_')
    Ca = float(test[1].split('/')[1].split('Ca')[0])
    Al = float(test[3].split('Al')[0])
    Si = float(test[4].split('Si')[0])

    mu_bar = Ca * (mCa + mO) + Al * (mAl + 1.5 * mO) + Si * (mSi + 2 * mO) + (mMg + mO)

    FeO_act = Mantle_wt_per.get('FeO') / 100
    mol_Fe_act = (mu_bar / (mFe + mO)) * (-1 + 1 / (1 - FeO_act))

    FeO_file_1 = float('%.02f' % (range_FeO[(np.abs(range_FeO - FeO_act)).argmin()]))
    FeO_file_1_id = int(np.where(range_FeO == FeO_file_1)[0])
    if FeO_act > FeO_file_1:

        FeO_file_2_id = FeO_file_1_id + 1
        test[-1] = str(format(range_FeO[FeO_file_1_id], '.2f')) + 'Fe'
        filename_down = '_'.join(test)
        test[-1] = str(format(range_FeO[FeO_file_2_id], '.2f')) + 'Fe'
        filename_up = '_'.join(test)
        X_wt_up = (FeO_act - range_FeO[FeO_file_1_id]) / (range_FeO[FeO_file_2_id] - range_FeO[FeO_file_1_id])
        mol_Fe_down = (mu_bar / (mO + mFe)) * (-1 + (1 / (1 - range_FeO[FeO_file_1_id])))
        mol_Fe_up = (mu_bar / (mO + mFe)) * (-1 + (1 / (1 - range_FeO[FeO_file_2_id])))

    else:
        FeO_file_2_id = FeO_file_1_id - 1
        test[-1] = str(format(range_FeO[FeO_file_1_id], '.2f')) + 'Fe'
        filename_up = '_'.join(test)
        test[-1] = str(format(range_FeO[FeO_file_2_id], '.2f')) + 'Fe'
        filename_down = '_'.join(test)
        X_wt_up = (FeO_act - range_FeO[FeO_file_2_id]) / (range_FeO[FeO_file_1_id] - range_FeO[FeO_file_2_id])
        mol_Fe_down = (mu_bar / (mO + mFe)) * (-1 + (1 / (1 - range_FeO[FeO_file_2_id])))
        mol_Fe_up = (mu_bar / (mO + mFe)) * (-1 + (1 / (1 - range_FeO[FeO_file_1_id])))

    FeO_mol_per_up = (mol_Fe_up) / (1 + mol_Fe_up)
    FeO_mol_per_down = (mol_Fe_down) / (1 + mol_Fe_down)
    FeO_mol_per_act = (mol_Fe_act) / (1 + mol_Fe_act)

    X_mol_up = (FeO_mol_per_act - FeO_mol_per_down) / (FeO_mol_per_up - FeO_mol_per_down)

    assert X_mol_up >= 0 and X_wt_up > 0, "Problem"

    if UMLM == True:
        file_1 = open(filename_up + '_UM_results.txt', 'r')
        file_2 = open(filename_down + '_UM_results.txt', 'r')
        P_up = 1390000
        P_down = 1
        T_up = 3300
        T_down = 1500
    else:
        file_1 = open(filename_up + '_LM_results.txt', 'r')
        file_2 = open(filename_down + '_LM_results.txt', 'r')
        P_up = 27000000.0
        P_down = 1250000.0
        T_up = 6800
        T_down = 1750

    temp_file = file_1.readlines()
    num_rows = len(temp_file[1:])
    num_columns = len(temp_file[12].split(','))
    start = 1

    for i in temp_file[1:]:
        if i[0] == '#':
            num_rows = num_rows - 1
            start += 1

    header = temp_file[0].strip('\n').split(',')

    Phases = header[5:-1]
    for i in range(len(Phases)):
        Phases[i] = Phases[i].strip()

    data = temp_file[start:]
    grid = np.zeros((num_rows, num_columns))

    for i in range(num_rows):
        # for j in range(num_columns):
        columns = data[i].strip('\n').split(',')
        grid[i] = [float(j) for j in columns]

    ##
    data = temp_file[start:]
    grid = np.zeros((num_rows, num_columns))

    for i in range(num_rows):
        # for j in range(num_columns):
        columns = data[i].strip('\n').split(',')
        grid[i] = [float(j) for j in columns]

    num_phases = len(grid[0][5:]) - 1

    temperature_grid_up = np.array([row[1] for row in grid])
    pressure_grid_up = np.array([row[0] for row in grid])
    density_grid_up = np.array([X_wt_up * row[2] * 1000 for row in grid])
    alpha_grid_up = [X_mol_up * pow(10, row[3]) for row in grid]
    cp_grid_up = [X_wt_up * row[4] for row in grid]
    phase_grid_up = [X_mol_up * row[5:-1] for row in grid]

    PT = np.vstack((pressure_grid_up, temperature_grid_up)).T
    tri_PT = Delaunay(PT)  # Compute the triangulation
    interpolator_rho_up = LinearNDInterpolator(tri_PT, np.array(density_grid_up))
    interpolator_alpha_up = LinearNDInterpolator(tri_PT, np.array(alpha_grid_up))
    interpolator_CP_up = LinearNDInterpolator(tri_PT, np.array(cp_grid_up))
    interpolator_phases_up = LinearNDInterpolator(tri_PT, np.array(phase_grid_up))

    ##
    temp_file = file_2.readlines()
    num_rows = len(temp_file[1:])
    num_columns = len(temp_file[12].split(','))
    start = 1

    for i in temp_file[1:]:
        if i[0] == '#':
            num_rows = num_rows - 1
            start += 1


    data = temp_file[start:]
    grid = np.zeros((num_rows, num_columns))

    for i in range(num_rows):
        # for j in range(num_columns):
        columns = data[i].strip('\n').split(',')
        grid[i] = [float(j) for j in columns]

        ##
    data = temp_file[start:]
    grid = np.zeros((num_rows, num_columns))

    for i in range(num_rows):
        # for j in range(num_columns):
        columns = data[i].strip('\n').split(',')
        grid[i] = [float(j) for j in columns]

    temperature_grid_down = np.array([row[1] for row in grid])
    pressure_grid_down = np.array([row[0] for row in grid])
    density_grid_down = np.array([(1 - X_wt_up) * row[2] * 1000 for row in grid])
    alpha_grid_down = [(1 - X_mol_up) * pow(10, row[3]) for row in grid]
    cp_grid_down = [(1 - X_wt_up) * row[4] for row in grid]
    phase_grid_down = [(1 - X_mol_up) * row[5:-1] for row in grid]

    PT = np.vstack((pressure_grid_down, temperature_grid_down)).T
    tri_PT = Delaunay(PT)  # Compute the triangulation
    interpolator_rho_down = LinearNDInterpolator(tri_PT, np.array(density_grid_down))
    interpolator_alpha_down = LinearNDInterpolator(tri_PT, np.array(alpha_grid_down))
    interpolator_CP_down = LinearNDInterpolator(tri_PT, np.array(cp_grid_down))
    interpolator_phases_down = LinearNDInterpolator(tri_PT, np.array(phase_grid_down))

    n_pts = 55
    dP = P_up - P_down
    dT = T_up - T_down

    P_new = np.hstack([[P_down + (dP / (n_pts - 1)) * j for i in range(n_pts)] for j in range(n_pts)])

    T_new = np.hstack([T_down + (dT / (n_pts - 1)) * j for i in range(n_pts) for j in range(n_pts)])

    mesh = np.vstack((P_new, T_new)).T
    rho = interpolator_rho_up(mesh) + interpolator_rho_down(mesh)

    alpha = interpolator_alpha_up(mesh) + interpolator_alpha_down(mesh)
    Cp = interpolator_CP_up(mesh) + interpolator_CP_down(mesh)
    phases = interpolator_phases_up(mesh) + interpolator_phases_down(mesh)
    phases = np.asarray([[phases[j][i] / sum(phases[j]) for i in range(num_phases)] for j in range(len(phases))])

    PT = np.vstack((P_new, T_new)).T
    tri_PT = Delaunay(PT)  # Compute the triangulation
    interpolator_rho = LinearNDInterpolator(tri_PT, rho)
    interpolator_alpha = LinearNDInterpolator(tri_PT, alpha)
    interpolator_CP = LinearNDInterpolator(tri_PT, Cp)
    interpolator_phases = LinearNDInterpolator(tri_PT, phases)

    keys = ['density', 'alpha', 'cp', 'phases']

    return dict(zip(keys, [interpolator_rho, interpolator_alpha, interpolator_CP, interpolator_phases])), Phases

def make_mantle_grid(Mantle_filename,Mantle_wt_per,UMLM,use_grids):
    """
    This module converts the PerPlex or premade grids into a dictionary of individual lists (e.g., pressure) for use
    by ExoPlex integrators

    Parameters
    ----------
    Mantle_filename: string
        name of file either from PerPlex or premade grids

    UMLM: boolean
        True for upper mantle grids, False for lower mantle grids

    use_grids: boolean
        True is user is using premade grids, false if using perplex-derived grids

    Returns
    -------
    grid_dictionary: dictionary of lists
        dictionary of individual parameters taken from the phase diagram.
        Keys include: 'temperature','pressure','density','alpha','cp','phases'

    """
    #Use ExoPlex pre-made grid
    if use_grids==True:
        in_grid = False
        if np.where(range_FeO==Mantle_wt_per['FeO']/100)[0] > -1:
            in_grid = True

        if Mantle_wt_per.get('FeO') > 0 and in_grid == False:
            return(make_mantle_feo_grid(Mantle_filename,Mantle_wt_per,UMLM))

        if UMLM == True:
            file = open(Mantle_filename+'_UM_results.txt','r')
        else:
            file = open(Mantle_filename+'_LM_results.txt','r')
        temp_file = file.readlines()
        num_rows = len(temp_file[1:])
        num_columns = len(temp_file[12].split(','))
        start = 1

        for i in temp_file[1:]:
            if i[0] == '#':
                num_rows = num_rows-1
                start+=1


        header = temp_file[0].strip('\n').split(',')

        Phases = header[5:-1]
        for i in range(len(Phases)):
            Phases[i] = Phases[i].strip()


        #calculate number of rows getting rid of #'s

        data = temp_file[start:]
        grid = np.zeros((num_rows,num_columns))

        for i in range(num_rows):
            #for j in range(num_columns):
            columns = data[i].strip('\n').split(',')
            grid[i] = [float(j) for j in columns]

        num_phases = len(grid[0][5:])-1

        temperature_grid = np.array([row[1] for row in grid])

        pressure_grid = np.array([row[0] for row in grid])
        density_grid = np.array([row[2]*1000 for row in grid])
        alpha_grid = [pow(10,row[3]) for row in grid]
        cp_grid = [row[4] for row in grid]
        phase_grid = [row[5:-1] for row in grid]

        PT = np.vstack((pressure_grid, temperature_grid)).T
        tri_PT = Delaunay(PT)  # Compute the triangulation
        interpolator_rho = LinearNDInterpolator(tri_PT, np.array(density_grid))
        interpolator_alpha = LinearNDInterpolator(tri_PT, np.array(alpha_grid))
        interpolator_CP = LinearNDInterpolator(tri_PT, np.array(cp_grid))
        interpolator_phases = LinearNDInterpolator(tri_PT, np.array(phase_grid))

        keys = ['density','alpha','cp','phases']

        return dict(zip(keys,[interpolator_rho,interpolator_alpha,interpolator_CP,interpolator_phases])),Phases

    else:
        #Use PerPlex derived grid
        if UMLM == True:
            file = open(Mantle_filename + '_UM_results.txt', 'r')
        else:
            file = open(Mantle_filename + '_LM_results.txt', 'r')

        temp_file = file.readlines()
        num_rows = len(temp_file[13:])
        num_columns = len(temp_file[12].split())

        header = temp_file[12].strip('\n').split()
        Phases = header[5:]

        for i in range(len(Phases)):
            Phases[i] = Phases[i].strip(",mo%")

        data = temp_file[13:]
        grid = np.zeros((num_rows, num_columns))

        for i in range(num_rows):
            columns = data[i].strip('\n').split()
            grid[i] = [float(j) for j in columns]

        num_phases = len(grid[0][5:])
        phases_grid = np.zeros((num_rows, num_phases))
        for i in range(num_rows):
            phases_grid[i] = grid[i][5:]

        temperature_grid = [row[0] for row in grid]
        pressure_grid = [row[1] for row in grid]
        density_grid = [row[2] for row in grid]
        alpha_grid = [row[3] for row in grid]
        cp_grid = [row[4] for row in grid]
        phase_grid = [row[5:] for row in grid]

        PT = np.vstack((pressure_grid, temperature_grid)).T
        tri_PT = Delaunay(PT)  # Compute the triangulation
        interpolator_rho = LinearNDInterpolator(tri_PT, np.array(density_grid))
        interpolator_alpha = LinearNDInterpolator(tri_PT, np.array(alpha_grid))
        interpolator_CP = LinearNDInterpolator(tri_PT, np.array(cp_grid))
        interpolator_phases = LinearNDInterpolator(tri_PT, np.array(phase_grid))

        keys = ['density', 'alpha', 'cp', 'phases']

        return dict(zip(keys, [interpolator_rho, interpolator_alpha, interpolator_CP, interpolator_phases])), Phases


