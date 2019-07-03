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

aa_pKx_dict = {'A': NaN,
               'R': 12.38,
               'N': NaN,
               'D': 3.65,
               'C': 8.18,
               'E': 4.25,
               'Q': NaN,
               'G': 'NA',
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

aa_phys_prop_df = pd.DataFrame([aa_name_dict, aa_pKx_dict, aa_ph7_hydrophob_dict],
                               index=['name', 'pKx', 'hydrophobicity']).T

print(aa_phys_prop_df)
