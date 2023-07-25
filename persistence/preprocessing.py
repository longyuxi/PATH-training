from pathlib import Path
import pickle
import types

import pandas as pd
from Bio.PDB import PDBParser
from biopandas.mol2 import PandasMol2
from biopandas.pdb import PandasPdb
import numpy as np

def load_pdbbind_data_index(index_filename: str) -> pd.DataFrame:
    index = pd.read_csv(index_filename, delim_whitespace=True, skiprows=6, names=['PDB code', "resolution", "release year", "-logKd/Ki", "Kd/Ki", "slashes", "reference", "ligand name"])

    index.drop(columns='slashes', inplace=True)
    index['ligand name'] = index.apply(lambda row:  row['ligand name'][1:][:-1], axis=1)

    return index

def get_pdb_coordinates(file) -> np.ndarray:
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = parser.get_structure(file, file)
    # Generate a list of the protein's atoms' R^3 coordinates
    coords = []
    for atom in structure.get_atoms():
        coords.append(list(atom.get_vector()))
    coords = np.array(coords)

    return coords

def get_pdb_coordinates_heavy(file) -> np.ndarray:
    # Similar to get_pdb_coordinates, but only gets the C, N, O, S coordinates
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = parser.get_structure(file, file)

    coords = []
    for atom in structure.get_atoms():
        if atom.element in ['C', 'N', 'O', 'S']:
            coords.append(list(atom.get_vector()))

    coords = np.array(coords)
    coords = np.round(coords, 5)

    return coords

def get_pdb_coordinates_by_element(file, element) -> np.ndarray:
    # Similar to get_pdb_coordinates, but only gets the coordinate matching the input element
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = parser.get_structure(file, file)

    coords = []
    for atom in structure.get_atoms():
        if atom.element == element:
            coords.append(list(atom.get_vector()))

    if len(coords) == 0:
        return np.array([]).reshape(0, 3)

    coords = np.array(coords)
    coords = np.round(coords, 5)

    return coords

def get_mol2_coordinates(file):
    # Gets coordinates from a mol2 file
    file = str(file)
    if '.mol2' in file:
        pmol = PandasMol2().read_mol2(file)
    elif '.pdb' in file:
        pmol = types.SimpleNamespace()
        pmol.df = pd.concat([PandasPdb().read_pdb(file).df['ATOM'], PandasPdb().read_pdb(file).df['HETATM']])
        pmol.df['atom_type'] = pmol.df['atom_name'].str[0]
        pmol.df['x'] = pmol.df['x_coord']
        pmol.df['y'] = pmol.df['y_coord']
        pmol.df['z'] = pmol.df['z_coord']
    else:
        raise Exception('Invalid file type.')

    return pmol.df[['x', 'y', 'z']].to_numpy()

def get_mol2_coordinates_heavy(file) -> np.ndarray:
    # Similar to get_mol2_coordinates, but only selects the following heavy atoms as detailed Cang and Wei: {C; N; O; S; P; F; Cl; Br}
    file = str(file)
    if '.mol2' in file:
        pmol = PandasMol2().read_mol2(file)
    elif '.pdb' in file:
        pmol = types.SimpleNamespace()
        pmol.df = pd.concat([PandasPdb().read_pdb(file).df['ATOM'], PandasPdb().read_pdb(file).df['HETATM']])
        pmol.df['atom_type'] = pmol.df['atom_name'].str[0]
        pmol.df['x'] = pmol.df['x_coord']
        pmol.df['y'] = pmol.df['y_coord']
        pmol.df['z'] = pmol.df['z_coord']
    else:
        raise Exception('Invalid file type.')

    coords = []

    for idx in range(len(pmol.df)):
        sybyl_atom_type = pmol.df.loc[idx, 'atom_type']
        if sybyl_atom_type.split('.')[0] in ['C', 'N', 'O', 'S', 'F', 'P', 'Cl', 'Br', 'I']:
            coords.append(pmol.df.loc[idx, ['x', 'y', 'z']].to_numpy())

    if len(coords) == 0:
        raise Exception('No heavy element.')
    coords = np.array(coords)
    return coords

def get_mol2_coordinates_by_element(file, element) -> np.ndarray:
    # Similar to get_mol2_coordinates, but only selects the following heavy atoms as detailed Cang and Wei: {C; N; O; S; P; F; Cl; Br}
    file = str(file)
    if '.mol2' in file:
        pmol = PandasMol2().read_mol2(file)
    elif '.pdb' in file:
        pmol = types.SimpleNamespace()
        pmol.df = pd.concat([PandasPdb().read_pdb(file).df['ATOM'], PandasPdb().read_pdb(file).df['HETATM']])
        pmol.df['atom_type'] = pmol.df['atom_name'].str[0]
        pmol.df['x'] = pmol.df['x_coord']
        pmol.df['y'] = pmol.df['y_coord']
        pmol.df['z'] = pmol.df['z_coord']
    else:
        raise Exception('Invalid file type.')

    coords = []

    for idx in range(len(pmol.df)):
        sybyl_atom_type = pmol.df.loc[idx, 'atom_type']
        if sybyl_atom_type.split('.')[0] == element:
            coords.append(pmol.df.loc[idx, ['x', 'y', 'z']].to_numpy())

    if len(coords) == 0:
        return np.array([]).reshape(0, 3)

    coords = np.array(coords)
    return coords
