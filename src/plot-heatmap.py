import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import preprocessing
import os, sys

#df = pd.read_table('top-10-MI-cognate-combined.txt', delim_whitespace=True, names=('Dpr','DIP'),
#             dtype={'Dpr': np.str, 'DIP': np.str})

#d = {
#    'A': '1',
#    'R': '2',
#    'N': '3',
#    'D': '4',
#    'B': '5',
#    'C': '6',
#    'E': '7',
#    'Q': '8',
#    'Z': '9',
#    'G': '10',
#    'H': '11',
#    'I': '12',
#    'L': '13',
#    'K': '14',
#    'M': '15',
#    'F': '16',
#    'P': '17',
#    'T': '18',
#    'W': '19',
#    'Y': '20',
#    'V': '21',
#    'S': '22',
#    'X': '23'
#    }

my_dict = dict(A='1',  R='2',  N='3',  D='4',  B='5',  C='6',  E='7',  Q='8',  Z='9',
               G='10', H='11', I='12', L='13', K='14', M='15', F='16', P='17', T='18',
               W='19', Y='20', V='21', S='22', X='23')

s = pd.read_table('top-10-MI-cognate-combined.txt', delim_whitespace=True, names=('Dpr', 'DIP'))

# unstack the pandas table into a list of all pairs
# and return the counts of each pair value
s_grp = (s.groupby(['Dpr', 'DIP']).size()
         .unstack(fill_value=0)
         .stack()
         .to_frame('size').reset_index())

# turn it into a float, scale for plotting
# currently normalized by the highest pair count in the protein
x = s_grp[['size']].values.astype(float)
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
s_normalized = pd.DataFrame(x_scaled)
#rounded = x.round(2)
#print(s_normalized)

s_grp = s_grp.assign(s_normalized=s_normalized.round(1).values)
#s_grp = s_grp.assign(rounded=rounded.values)

print(s_grp)

#df = pd.DataFrame(s_grp, columns=['Dpr','DIP','size'])
#print(df)

#heatmap_data_cognate = pd.pivot_table(s_grp, values='size', index=['Dpr'], columns='DIP')

heatmap_data_cognate = pd.pivot_table(s_grp, values='s_normalized', index=['Dpr'], columns='DIP')

mask = heatmap_data_cognate == 0
sns.heatmap(heatmap_data_cognate,cmap='coolwarm',annot=True,mask=mask)
plt.gca().invert_yaxis()  #put y-axis labels in reverse order
#plt.show()
plt.savefig("plotmat-cg-normalized.png", bbox_inches='tight')

#---------------------------scatter
#fig, ax = plt.subplots()
#sc = ax.scatter(df.Dpr, df.DIP, cmap='coolwarm')
#ax.set_aspect("equal")
#plt.show()
