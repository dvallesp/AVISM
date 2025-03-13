"""
read_voids module
Provides a function to read the void catologue and a function to read
the 3D void structure outputs

Created by Óscar Monllor and Mónica Hernández
"""

# GENERAL PURPOSE AND SPECIFIC LIBRARIES USED IN THIS MODULE
import numpy as np
import os
from scipy.io import FortranFile as FF

# FUNCTIONS DEFINED IN THIS MODULE:

def read_void_catalogue(it, path='', output_format = 'dictionaries', ret_details=True, ret_catalogue=True):
    """
    Reads the voidXXXXX files, containing the void catalogue of the MASCLET void finder.
    
    Args:
        - it: iteration number (int)
        - path: path of the stellar_haloes (str)
        - output_format: 'dictionaries' or 'arrays'

    Returns:
        - if ret_details=True: dictionary with voidfinder details
        - if ret_catalogue=True: list of dictionaries, one per void level
            - if output_format='dictionaries': list of dictionaries, one per void
            - if output_format='arrays': dictionary of arrays, one array per void property
        - if both are true: dictionary with voidfinder details and list of dictionaries, one per void level
    """
    
    #Check output format
    if output_format not in ['dictionaries', 'arrays']:
        raise ValueError('output_format must be "dictionaries" or "arrays"')

    with open(os.path.join(path, 'voids{:05d}'.format(it)), 'r') as void_catalogue:
        # Level data
        void_details = {}
        header = void_catalogue.readline().split()
        nlev = int(header[0])
        levmin = int(header[1])
        levmax = int(header[2])
        nmax = int(header[3]) #l = 0 grid size !!!
        nmay = int(header[4])
        nmaz = int(header[5])
        L0 = float(header[6])
        void_details['nlev'] = nlev
        void_details['levmin'] = levmin
        void_details['levmax'] = levmax
        void_details['nmax'] = nmax
        void_details['nmay'] = nmay
        void_details['nmaz'] = nmaz
        void_details['L0'] = L0

        if ret_catalogue:
            levels = []
            for ilev in range(nlev):
                this_level = {}
                lev, ncubes, nvoids, nparents, FF = void_catalogue.readline().split()
                lev = int(lev)
                ncubes = int(ncubes)
                nvoids = int(nvoids)
                nparents = int(nparents)
                FF = float(FF)
                this_level['level'] = lev
                this_level['ncubes'] = ncubes
                this_level['nvoids'] = nvoids
                this_level['nparents'] = nparents
                this_level['FF'] = FF
                
                # Void data
                voids=[]
                for iv in range(nvoids):
                    void = {}
                    data_line = np.array(void_catalogue.readline().split()).astype('f4')
                    void['id'] = int(data_line[0])
                    void['xc'] = data_line[1]
                    void['yc'] = data_line[2]
                    void['zc'] = data_line[3]
                    void['gxc'] = data_line[4]
                    void['gyc'] = data_line[5]
                    void['gzc'] = data_line[6]
                    void['vol'] = data_line[7]
                    void['R'] = data_line[8]
                    void['mean_overdensity'] = data_line[9]
                    void['ellipticity'] = data_line[10]
                    void['inv_porosity'] = data_line[11]
                    void['id_father'] = int(data_line[12])
                    void['R_father'] = data_line[13]
                    void['M'] = data_line[14]
                    voids.append(void)

                if output_format=='dictionaries':
                    this_level['voids'] = voids
                elif output_format=='arrays':
                    this_level['voids'] = {k: np.array([v[k] for v in voids]) for k in voids[0].keys()}

                levels.append(this_level)


    if ret_details and ret_catalogue:
        return void_details, levels
    elif ret_catalogue:
        return levels
    elif ret_details:
        return void_details


def read_cubes(it, path='', output_format = 'dictionaries'):
    """
    Reads the cubesXXXXX files, containing all cubes making up the final void catalogue

    Args:
        it: iteration number (int)
        path: path of the grids file in the system (str)
        output_format: 'dictionaries' or 'arrays'   

    Returns:
        List of dictionaries, one per level, containing the cube data:
        - 'level': level number (int)
        - 'ncubes': number of cubes (int)
        - 'nmains': number of main voids (int)
        - 'cubes': list of dictionaries, one per cube, containing the following
            - 'id': cube ID (int)
            - 'uvoid': -1 (it's main), 0 (discarted), >0 main void ID to which it belongs
            - 'ini_ix', 'end_ix': initial and final x indices (int) on the grid
            - 'ini_iy', 'end_iy': initial and final y indices (int) on the grid
            - 'ini_iz', 'end_iz': initial and final z indices (int) on the grid
            - 'ini_rx', 'end_rx': initial and final x coordinates (float)
            - 'ini_ry', 'end_ry': initial and final y coordinates (float)
            - 'ini_rz', 'end_rz': initial and final z coordinates (float)
            - 
    """

    fname = 'cubes{:05d}'.format(it)

    with open(os.path.join(path, fname), 'r') as cube_catalogue:
        # Level data
        header = cube_catalogue.readline().split()
        nlev = int(header[0])
        levmin = int(header[1])
        levmax = int(header[2])

        levels = []
        for ilev in range(nlev):
            this_level = {}
            lev, ncubes, nmains = cube_catalogue.readline().split()
            lev = int(lev)
            ncubes = int(ncubes)
            nmains = int(nmains)
            this_level['level'] = lev
            this_level['ncubes'] = ncubes
            this_level['nmains'] = nmains
    
            # cube data
            cubevoids=[]
            for iv in range(ncubes):
                C = {}
                data_line = np.array(cube_catalogue.readline().split()).astype('f4')
                C['id'] = int(data_line[0])
                C['uvoid'] = int(data_line[1])
                C['ini_ix'] = int(data_line[2])
                C['end_ix'] = int(data_line[3])
                C['ini_iy'] = int(data_line[4])
                C['end_iy'] = int(data_line[5])
                C['ini_iz'] = int(data_line[6])
                C['end_iz'] = int(data_line[7])
                C['ini_rx'] = data_line[8]
                C['end_rx'] = data_line[9]
                C['ini_ry'] = data_line[10]
                C['end_ry'] = data_line[11]
                C['ini_rz'] = data_line[12]
                C['end_rz'] = data_line[13]

                cubevoids.append(C)

            if output_format=='dictionaries':
                this_level['cubes'] = cubevoids
            elif output_format=='arrays':
                this_level['cubes'] = {k: np.array([v[k] for v in cubevoids]) for k in cubevoids[0].keys()}

            levels.append(this_level)

    return levels


def read_void_map(it, path='', output_marca=True, output_deltaTot=True, output_div=False):
    """
    Reads the 3D arrays (map) file

    Args:
        it: iteration number (int)
        path: path of the grids file in the system (str)
        output_marca: boolean, whether to output marca (default=True)
        output_deltaTot: boolean, whether to output deltaTot (default=True)
        output_div: boolean, whether to output div (default=False)
    Returns:
        List of arrays (one for each level); marca[:,:,:], deltaTot[:,:,:], deltaGas[:,:,:]
    """
    # First, call read_void_catalogue to get the voidfinder details
    void_details = read_void_catalogue(it, path, ret_catalogue=False)
    nlev = void_details['nlev']
    levmin = void_details['levmin']
    levmax = void_details['levmax']
    nmax = void_details['nmax'] #l = 0 grid size !!!
    nmay = void_details['nmay']
    nmaz = void_details['nmaz']
    L0 = void_details['L0']

    if output_marca:
        marca = []
    if output_deltaTot:
        deltaTot = []
    if output_div:
        div = []

    with FF(os.path.join(path, f'map{it:05d}')) as f:
        for ilev in range(levmin, levmax+1):
            nx_lev = int(nmax * 2**(ilev-levmin))
            ny_lev = int(nmay * 2**(ilev-levmin))
            nz_lev = int(nmaz * 2**(ilev-levmin))
            if output_marca:
                marca.append(f.read_ints('i4').reshape((nx_lev, ny_lev, nz_lev)).T)
            else:
                f.read_ints('i4')

            if output_deltaTot:
                deltaTot.append(f.read_reals('f4').reshape((nx_lev, ny_lev, nz_lev)).T)
            else:
                f.read_reals('f4')


            if output_div:
                div.append(f.read_reals('f4').reshape((nx_lev, ny_lev, nz_lev)).T)
            else:
                f.read_reals('f4')
                

    output = []
    if output_marca:
        output.append(marca)
    if output_deltaTot:
        output.append(deltaTot)
    if output_div:
        output.append(div)

    return output
