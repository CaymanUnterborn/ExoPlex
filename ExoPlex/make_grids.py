import numpy as np
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
import pandas as pd

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

    filename = '../Solutions_Small/liquid_iron_grid.dat'
    df = pd.read_csv(filename)
    df = df.iloc[::10]
    df.rename(columns={'# P_bar': 'P_bar'},inplace=True)
    pressure_grid = np.array(df['P_bar'])
    temperature_grid = np.array(df['T_K'])
    density_grid = np.array(df['rho_kgm3'])
    alpha_grid = np.array(pow(10,df['alpha_1_K']))
    cp_grid = np.array(df['Cp_J_K_mol'])

    PT = np.vstack((pressure_grid, temperature_grid)).T
    tri_PT = Delaunay(PT)  # Compute the triangulation
    interpolator_rho = LinearNDInterpolator(tri_PT, density_grid)
    interpolator_alpha = LinearNDInterpolator(tri_PT, alpha_grid)
    interpolator_CP = LinearNDInterpolator(tri_PT, cp_grid)

    keys = ['density','alpha','cp']
    return dict(zip(keys,[interpolator_rho,interpolator_alpha,interpolator_CP]))

def make_water_grid():
    filename = '../Solutions_Small/water_grid.dat'
    df = pd.read_csv(filename)
    df.rename(columns={'# P': 'P'}, inplace=True)
    pressure_grid = np.array(df['P'])
    temperature_grid = np.array(df['T'])
    density_grid = np.array(df['density'])
    alpha_grid = np.array(pow(10, df['alpha']))
    cp_grid = np.array(df['Cp'])
    phase_names = (np.array(df.columns[df.columns.get_loc("alpha") + 1:])).tolist()
    phase_grid = df[phase_names]
    df['sum_phases'] = phase_grid.sum(axis=1)
    phase_grid = np.asarray(100 * df[phase_names].div(df["sum_phases"], axis=0))



    PT = np.vstack((pressure_grid, temperature_grid)).T
    tri_PT = Delaunay(PT)  # Compute the triangulation
    interpolator_rho = LinearNDInterpolator(tri_PT, np.array(density_grid))
    interpolator_alpha = LinearNDInterpolator(tri_PT, np.array(alpha_grid))
    interpolator_CP = LinearNDInterpolator(tri_PT, np.array(cp_grid))
    interpolator_phases = LinearNDInterpolator(tri_PT, np.array(phase_grid))

    keys = ['density', 'alpha', 'cp', 'phases']

    return dict(zip(keys, [interpolator_rho, interpolator_alpha, interpolator_CP, interpolator_phases])),phase_names

def make_mantle_feo_grid(Mantle_filename,Mantle_wt_per,UMLM):

    test = Mantle_filename.split('_')
    Ca = float(test[1].split('/')[1].split('Ca')[0])
    Al = float(test[3].split('Al')[0])
    Si = float(test[4].split('Si')[0])

    mu_bar = Ca * (mCa + mO) + Al * (mAl + 1.5 * mO) + Si * (mSi + 2 * mO) + (mMg + mO)

    FeO_act = Mantle_wt_per.get('FeO') / 100
    mol_Fe_act = (mu_bar / (mFe + mO)) * (-1 + 1 / (1 - FeO_act))

    FeO_file_1 = float('%.02f' % (range_FeO[(np.abs(range_FeO - FeO_act)).argmin()]))
    FeO_file_1_id = int((np.where(range_FeO == FeO_file_1))[0][0])
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
        file_open_1 = filename_up + '_UM_results.txt'
        file_open_2 = filename_down + '_UM_results.txt'
        P_up = 1390000
        P_down = 1
        T_up = 3300
        T_down = 1500
    else:
        file_open_1 = filename_up + '_LM_results.txt'
        file_open_2 = filename_down + '_LM_results.txt'
        P_up = 27000000.0
        P_down = 1250000.0
        T_up = 6800
        T_down = 1750

    df = pd.read_csv(file_open_1)
    df = df.iloc[::3]

    df.rename(columns={'#P[bar]': 'P[bar]'}, inplace=True)
    df.columns = df.columns.str.lstrip()
    pressure_grid_up = np.array(df['P[bar]'])
    temperature_grid_up = np.array(df['T[K]'])
    density_grid_up = np.array(X_wt_up *1000 * df['rho[g/cm3]'])
    alpha_grid_up = np.array(X_mol_up*pow(10, df['log10(alpha)[1/K]']))
    cp_grid_up = np.array(X_wt_up*df['cp[J/(kg*K)]'])
    phase_names_up = (np.array(df.columns[df.columns.get_loc("cp[J/(kg*K)]") + 1:-1])).tolist()
    phase_grid_up = X_wt_up*df[phase_names_up]
    num_phases = len(phase_names_up)

    PT = np.vstack((pressure_grid_up, temperature_grid_up)).T
    tri_PT = Delaunay(PT)  # Compute the triangulation
    interpolator_rho_up = LinearNDInterpolator(tri_PT, np.array(density_grid_up))
    interpolator_alpha_up = LinearNDInterpolator(tri_PT, np.array(alpha_grid_up))
    interpolator_CP_up = LinearNDInterpolator(tri_PT, np.array(cp_grid_up))
    interpolator_phases_up = LinearNDInterpolator(tri_PT, np.array(phase_grid_up))

    ##
    df = pd.read_csv(file_open_2)
    df = df.iloc[::3]
    df.rename(columns={'#P[bar]': 'P[bar]'}, inplace=True)
    df.columns = df.columns.str.lstrip()
    pressure_grid_down= np.array(df['P[bar]'])
    temperature_grid_down = np.array(df['T[K]'])
    density_grid_down = np.array((1.-X_wt_up) *1000 * df['rho[g/cm3]'])
    alpha_grid_down = np.array((1.-X_mol_up)*pow(10, df['log10(alpha)[1/K]']))
    cp_grid_down = np.array((1.-X_wt_up)*df['cp[J/(kg*K)]'])
    phase_names_down = (np.array(df.columns[df.columns.get_loc("cp[J/(kg*K)]") + 1:-1])).tolist()
    phase_grid_down = (1.-X_wt_up)*df[phase_names_down]

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

    return dict(zip(keys, [interpolator_rho, interpolator_alpha, interpolator_CP, interpolator_phases])), phase_names_up

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
        test = np.where(range_FeO == Mantle_wt_per['FeO'] / 100)[0]
        if len(test) > 0:
            in_grid = True

        if Mantle_wt_per.get('FeO') / 100 > 0 and Mantle_wt_per['FeO'] / 100 <= 0.2 and in_grid == False:
            return (make_mantle_feo_grid(Mantle_filename, Mantle_wt_per, UMLM))


        if UMLM == True:
            filename = Mantle_filename+'_UM_results.txt'

        else:
            filename = Mantle_filename+'_LM_results.txt'

        df = pd.read_csv(filename)
        df = df.iloc[::2]

        df.rename(columns={'#P[bar]': 'P[bar]'}, inplace=True)
        df.columns = df.columns.str.lstrip()
        pressure_grid = np.array(df['P[bar]'])
        temperature_grid = np.array(df['T[K]'])
        density_grid = np.array(1000*df['rho[g/cm3]'])
        alpha_grid = np.array(pow(10, df['log10(alpha)[1/K]']))
        cp_grid = np.array(df['cp[J/(kg*K)]'])
        phase_names = (np.array(df.columns[df.columns.get_loc("cp[J/(kg*K)]") + 1:])).tolist()
        phase_grid = df[phase_names]
        df['sum_phases'] = phase_grid.sum(axis=1)
        phase_grid = np.asarray(100*df[phase_names].div(df["sum_phases"], axis=0))

        PT = np.vstack((pressure_grid, temperature_grid)).T
        tri_PT = Delaunay(PT)  # Compute the triangulation
        interpolator_rho = LinearNDInterpolator(tri_PT, density_grid)
        interpolator_alpha = LinearNDInterpolator(tri_PT, alpha_grid)
        interpolator_CP = LinearNDInterpolator(tri_PT, cp_grid)
        interpolator_phases = LinearNDInterpolator(tri_PT, phase_grid)

        keys = ['density','alpha','cp','phases']

        return dict(zip(keys,[interpolator_rho,interpolator_alpha,interpolator_CP,interpolator_phases])),phase_names

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


