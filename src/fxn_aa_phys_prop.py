import pandas as pd
import numpy as np


# This file contains a function to return the physical properties of an individual amino acid
# And compare properties of the pairs of amino acids
NaN = np.nan
aa_name_dict = {'A': 'Alanine',
                'R': 'Arginine',
                'N': 'Asparagine',
                'D': 'Aspartic Acid',
                'C': 'Cysteine',
                'E': 'Glutamic acid',
                'Q': 'Glutamine',
                'G': 'Glycine',
                'H': 'Histidine',
                'O': 'Hydroxoyproline',
                'I': 'Isoleucine',
                'L': 'Leucine',
                'K': 'Lysine',
                'M': 'Methionine',
                'F': 'Phenylalanine',
                'P': 'Proline',
                'U': 'Pyroglutamic',
                'S': 'Serine',
                'T': 'Threonine',
                'W': 'Tryptophan',
                'Y': 'Tyrosine',
                'V': 'Valine'}

# pKa of side chain residue
aa_pKx_dict = {'A': NaN,
               'R': 12.38,
               'N': NaN,
               'D': 3.65,
               'C': 8.18,
               'E': 4.25,
               'Q': NaN,
               'G': NaN,
               'H': 6.0,
               'O': NaN,
               'I': NaN,
               'L': 'NA',
               'K': 10.53,
               'M': NaN,
               'F': NaN,
               'P': NaN,
               'U': NaN,
               'S': NaN,
               'T': NaN,
               'W': NaN,
               'Y': 10.07,
               'V': NaN}

# hydrophobicity scaled from 100 to -100
# taken from sigma-aldrich website,
# 0 is glycine, 100 is most hydrophobic, -100 is least hydrophilic
aa_ph7_hydrophob_dict = {'A': 41,
                         'R': -14,
                         'N': -28,
                         'D': -55,
                         'C': 49,
                         'E': -31,
                         'Q': -10,
                         'G': 0,
                         'H': 8,
                         'O': NaN,
                         'I': 99,
                         'L': 97,
                         'K': -23,
                         'M': 74,
                         'F': 100,
                         'P': -46,
                         'U': NaN,
                         'S': -5,
                         'T': 13,
                         'W': 97,
                         'Y': 63,
                         'V': 76}

aa_binary_props = {  'A': ['hydrophobic'],
                     'R': ['hydrophilic', 'ionic', 'plus', 'donor', 'acceptor'],
                     'N': ['hydrophilic', 'acceptor', 'donor'],
                     'D': ['hydrophilic', 'ionic', 'minus'],
                     'C': ['hydrophilic', 'donor'],
                     'E': ['hydrophilic', 'donor', 'minus', 'ionic'],
                     'Q': ['hydrophilic', 'donor', 'minus', 'ionic'],
                     'G': ['donor'],
                     'H': ['hydrophilic', 'donor', 'acceptor'],
                     'O': ['hydrophilic', 'ionic', 'minus'],
                     'I': ['hydrophobic'],
                     'L': ['hydrophobic'],
                     'K': ['hydrophobic', 'ionic', 'plus'],
                     'M': ['hydrophobic', 'acceptor'],
                     'F': ['hydrophobic'],
                     'P': ['hydrophobic'],
                     'U': ['hydrophilic', 'ionic', 'minus'],
                     'S': ['hydrophilic', 'donor', 'acceptor'],
                     'T': ['hydrophilic', 'donor', 'acceptor'],
                     'W': ['hydrophobic', 'donor'],
                     'Y': ['hydrophilic', 'donor', 'acceptor'],
                     'V': ['hydrophobic']}


aa_phys_prop_df = pd.DataFrame([aa_name_dict, aa_pKx_dict, aa_ph7_hydrophob_dict, aa_binary_props],
                               index=['name', 'pKx', 'hydrophobicity', 'binaryProp']).T

# Define some initial stuff:
AA_key = ['A', 'G', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S',
          'P', 'V', 'I', 'C', 'Y', 'H', 'R', 'N', 'D', 'T']
# We can also instead look at some properties like charge or kidera factors
# Kidera factors reference can be found in:
# Kidera, A., Konishi, Y., Oka, M. et al. J Protein Chem (1985) 4: 23. https://doi.org/10.1007/BF01025492
AA_num_key = np.arange(20)+1
# electrostatic Charges
charge_key = [0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, -0.048, 0, 0.091, 1, 0, -1, 0]
# Helix/bend preference
factor1 = [-1.56, 1.46, -1.04, -1.4, -0.21, 0.3, -0.34, -0.47, -1.45, 0.81,
           2.06, -0.74, -0.73, 0.12, 1.38, -0.41, 0.22, 1.14, 0.58, 0.26]
# Side-chain size
factor2 = [-1.67, -1.96, 0, 0.18, 0.98, 2.1, 0.82, 0.24, 0.19, -1.08,
           -0.33, -0.71, -0.16, -0.89, 1.48, 0.52, 1.27, -0.07, -0.22, -0.7]
# Extended structure preference
factor3 = [-0.97, -0.23, -0.24, -0.42, -0.36, -0.72, -0.23, 0.07, -1.61, 0.16,
           -1.15, 2.04, 1.79, 0.45, 0.8, -0.28, 1.37, -0.12, -1.58, 1.21]
# Hydrophobicity
factor4 = [-0.27, -0.16, -1.1, -0.73, -1.43, -1.57, 1.7, 1.1, 1.17, 0.42,
           -0.75, -0.4, -0.77, -1.05, -0.56, 0.28, 1.87, 0.81, 0.81, 0.63]
# Double-bend preference
factor5 = [-0.93, 0.1, -0.55, 2, 0.22, -1.16, 1.54, 1.1, -1.31, -0.21,
           0.88, 0.5, -0.54, -0.71, 0, 1.61, -1.7, 0.18, -0.92, -0.1]
# Partial specific volume
factor6 = [-0.78, -0.11, -2.05, 1.52, -0.81, 0.57, -1.62, 0.59, 0.4, -0.43,
           -0.45, -0.81, 0.03, 2.41, -0.68, 1.01, 0.46, 0.37, 0.15, 0.21]
# Flat extended preference
factor7 = [-0.2, 1.32, 0.96, 0.26, 0.67, -0.48, 1.15, 0.84, 0.04, -1.89,
           0.3, -1.07, -0.83, 1.52, -0.31, -1.85, 0.92, -0.09, -1.52, 0.24]
# Occurrence in alpha region
factor8 = [-0.08, 2.36, -0.76, 0.11, 1.1, -0.4, -0.08, -0.71, 0.38, -1.15,
           -2.3, 0.06, 0.51, -0.69, 1.03, 0.47, -0.39, 1.23, 0.47, -1.15]
# pK-C
factor9 = [0.21, -1.66, 0.45, -1.27, 1.71, -2.3, -0.48, -0.03, -0.35, -0.97,
           0.74, -0.46, 0.66, 1.13, -0.05, 1.13, 0.23, 1.1, 0.76, -0.56]
# Surrounding hydrophobicity
factor10 = [-0.48, 0.46, 0.93, 0.27, -0.44, -0.6, 0.6, -2.33, -0.12, -0.23,
            -0.28, 0.65, -1.78, 1.1, 0.53, 1.63, 0.93, -1.73, 0.7, 0.19]
props = np.vstack([charge_key, factor1, factor2, factor3, factor4,
                   factor5, factor6, factor7, factor8, factor9, factor10])
labels = ['charge', 'KF1', 'KF2', 'KF3', 'KF4', 'KF5', 'KF6', 'KF7', 'KF8', 'KF9', 'KF10']

kf_df = pd.DataFrame(props, columns=AA_key, index=labels).T

aa_phys_prop_df = aa_phys_prop_df.join(kf_df, how='outer')


def get_phys_prop_AA(aa_code, prop_code):
    """ Input the amino acid single letter code, and the prop code from the list below
        This will return a physical property for the amino acid of interest. \n
        Prop codes are: name, pKx, hydrophobicity,  charge,
        KF1, KF2, KF3, KF4, KF5, KF6, KF7, KF8, KF9,  KF10.\n
        The interpretation of the KF codes are as follows:\n
        1: Helix/bend preference. \n
        2: Side-chain size. \n
        3: Extended structure preference. \n
        4: Hydrophobicity. \n
        5: Double-bend preference. \n
        6: Partial specific volume. \n
        7: Flat extended preference. \n
        8: Occurrence in alpha region. \n
        9: pK-C\n
        10: Surrounding hydrophobicity. """
    df_row_index = aa_code
    column_index = prop_code
    return aa_phys_prop_df.loc[df_row_index, column_index]

# TODO: BLOSUM scores, block matrix for substitution and alignments similarity score
