## This file is part of PATH, which is part of OSPREY 3.0
##
## OSPREY Protein Redesign Software Version 3.0
## Copyright (C) 2001-2023 Bruce Donald Lab, Duke University
##
## OSPREY is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License version 2
## as published by the Free Software Foundation.
##
## You should have received a copy of the GNU General Public License
## along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
##
## OSPREY relies on grants for its development, and since visibility
## in the scientific literature is essential for our success, we
## ask that users of OSPREY cite our papers. See the CITING_OSPREY
## document in this distribution for more information.
##
## Contact Info:
##    Bruce Donald
##    Duke University
##    Department of Computer Science
##    Levine Science Research Center (LSRC)
##    Durham
##    NC 27708-0129
##    USA
##    e-mail: www.cs.duke.edu/brd/
##
## <signature of Bruce Donald>, Mar 1, 2023
## Bruce Donald, Professor of Computer Science

# Initial prototype: Calculate just the atomic persistence diagram
# Later: Move on to use pairiwise opposition homology as detailed in Wei's paper

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent / '..' / 'gbr'  / 'gbr_inference'))

import numpy as np
from scipy.spatial.distance import cdist
from gtda.homology import VietorisRipsPersistence

from preprocessing import load_pdbbind_data_index, get_mol2_coordinates_by_element, get_pdb_coordinates_by_element

################

platform = 'CS'  # Cluster specification

################

def atom_persistence_homology(coords):
    # Rewrite this (to-do)
    # Track connected components, loops, and voids
    homology_dimensions = [0, 1, 2]

    # Collapse edges to speed up H2 persistence calculation!
    persistence = VietorisRipsPersistence(
        homology_dimensions=homology_dimensions,
        collapse_edges=True,
        max_edge_length=200
    )

    diagrams_basic = persistence.fit_transform(coords[None, :, :])

    return diagrams_basic

def opposition_homology(protein_coords, ligand_coords):
    def opposition_distance_metric(vec1, vec2):
        if np.abs(vec1[-1] - vec2[-1]) > 2:  # If the two atoms are not of the same type
            return np.linalg.norm(vec1[:3] - vec2[:3])
        else:
            return np.Inf

    # Append each coordinate with 1 for protein and 2 for ligand
    protein_coords = np.concatenate((protein_coords, np.ones((len(protein_coords), 1))), axis=1)
    ligand_coords = np.concatenate((ligand_coords, 4 * np.ones((len(ligand_coords), 1))), axis=1)

    combined_coords = np.concatenate((protein_coords, ligand_coords), axis=0)

    if combined_coords.shape[0] == 0:
        return None

    if protein_coords is None or ligand_coords is None:
        return None

    # Track connected components, loops, and voids
    homology_dimensions = [0, 1, 2]

    # Collapse edges to speed up H2 persistence calculation!
    persistence = VietorisRipsPersistence(
        metric=opposition_distance_metric,
        homology_dimensions=homology_dimensions,
        collapse_edges=True,
        max_edge_length=200
    )

    diagrams_basic = persistence.fit_transform(combined_coords[None, :, :])

    return diagrams_basic

def get_pairwise_opposition_persistence_diagrams(pdb_file, mol2_file):
    protein_heavy_elements = ['C', 'N', 'O', 'S']
    ligand_heavy_elements = ['C', 'N', 'O', 'S', 'F', 'P', 'Cl', 'Br', 'I']
    diagrams = []

    for pe in protein_heavy_elements:
        for le in ligand_heavy_elements:
            # calculate homology (pe, le)
            # opposition_homologies.append(...)
            protein_coords = get_pdb_coordinates_by_element(pdb_file, pe)
            ligand_coords = get_mol2_coordinates_by_element(mol2_file, le)
            diagram = opposition_homology(protein_coords, ligand_coords)
            diagrams.append(diagram)

    return diagrams

def get_2345_persistence_diagrams(pdb_file, mol2_file):

    homologies = [] # this is used to store all of the calculated persistence diagrams

    def concatenate_coordinates(list_of_coordinates):
        # input: list of ndarray of size (*, 3)
        output = None
        for i in range(len(list_of_coordinates) - 1):
            if i == 0:
                output = np.concatenate((list_of_coordinates[i], list_of_coordinates[i+1]))
            else:
                output = np.concatenate((output, list_of_coordinates[i+1]))

        return output

    protein_heavy_elements = ['C', 'N', 'O', 'S']
    ligand_heavy_elements = ['C', 'N', 'O', 'S', 'F', 'P', 'Cl', 'Br', 'I']

    # 2: all heavy atoms of protein
    protein_heavy_atom_coords = []
    for pe in protein_heavy_elements:
        protein_coords = get_pdb_coordinates_by_element(pdb_file, pe)
        protein_heavy_atom_coords.append(protein_coords)

    protein_heavy_atom_coords = concatenate_coordinates(protein_heavy_atom_coords)

    homologies.append(atom_persistence_homology(protein_heavy_atom_coords))

    # 3: all heavy atoms of protein and all heavy atoms of ligand
    ligand_heavy_atom_coords = []
    for le in ligand_heavy_elements:
        ligand_coords = get_mol2_coordinates_by_element(mol2_file, le)
        ligand_heavy_atom_coords.append(ligand_coords)

    ligand_heavy_atom_coords = concatenate_coordinates(ligand_heavy_atom_coords)
    all_heavy_atom_coords = np.concatenate((protein_heavy_atom_coords, ligand_heavy_atom_coords))

    homologies.append(atom_persistence_homology(all_heavy_atom_coords))

    # 4: all carbon atoms of protein
    protein_carbon_coords = get_pdb_coordinates_by_element(pdb_file, 'C')
    homologies.append(atom_persistence_homology(protein_carbon_coords))

    # 5: all carbon atoms of protein and all heavy atoms of ligand
    ligand_carbon_coords = get_mol2_coordinates_by_element(mol2_file, 'C')
    all_carbon_coords = np.concatenate((protein_carbon_coords, ligand_carbon_coords))
    homologies.append(atom_persistence_homology(all_carbon_coords))

    return homologies
