import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv

#pre-processing of teh data and overall visualisation
name = 'spike_mutagensis_data.csv'
df = pd.read_csv(name)
df = df[: 2340]

original = df[df['Top_rep1'] == 0]['Substitution']

top_col = df[['Top_rep1','Top_rep2']].mean(axis=1).to_frame(name='top_means')
bottom_col = df[['Bottom_rep1','Bottom_rep2']].mean(axis=1).to_frame(name='bottom_means')
df = pd.concat([df, top_col, bottom_col], axis=1, sort=False)
seq = df['Residue #'].drop_duplicates()

#plotting mean top values against mean bottom value.
'''fig = plt.figure().gca(projection = '3d')
fig.scatter(df['Residue #'], df['top_means'],df['bottom_means'])
fig.set_xlabel('Residue Number')
fig.set_ylabel('Top_rep')
fig.set_zlabel('Bottom_rep')
plt.show()'''

#creating dataframes for hydrophilic/hydrophobic analysis
data = [['R', 3], ['D', 3], ['E', 3], ['K', 3], ['S', 0.3],
        ['N',0.2], ['Q', 0.2], ['G', 0],['P', 0],['T', -0.4],
        ['A', -0.5],['H', -0.5], ['C', -1], ['M', -1.3], ['V', -1.5],
        ['I', -1.8], ['L', -1.8], ['Y', -2.3], ['F', -2.5],['W', -3.4]]
hydro = pd.DataFrame(data, columns = ['original', 'Hydrophobicity_original'])
hydro2 = pd.DataFrame(data, columns = ['mutated', 'Hydrophobicity_mutated'])

#finding out the substitution that gives the greatest binding performance at each position
list= []
top_list = []
bottom_list = []
max_list = []
poss = []
aa = []
for i in seq:
    new_df = df[df['Residue #'] == i]

    top_max = new_df['top_means'].max()
    bottom_min = new_df['bottom_means'].min()

    top_aa = new_df[new_df['top_means'] == top_max].iloc[0]['Substitution']
    bottom_aa = new_df[new_df['bottom_means'] == bottom_min].iloc[0]['Substitution']
    likelihood = new_df[new_df['bottom_means'] == bottom_min].iloc[0]['top_means']

    top_list.append(top_aa)
    bottom_list.append(bottom_aa)
    max_list.append(top_max)
    aa.append(top_aa)
    if top_aa == bottom_aa:
        list.append(top_aa)
        poss.append(likelihood)
    else:
        list.append('NA')
        poss.append('NA')

result = pd.DataFrame(zip(seq,original,list,poss), columns = ('sequence', 'original', 'mutated', 'likelihood'))
region_join_df = pd.merge(result,hydro, on='original', how = 'left')
region_join_df = pd.merge(region_join_df,hydro2, on='mutated', how = 'left')


for i in range(len(region_join_df['sequence'])):
    region_join_df['hydrophobicity_difference'.format(i)] = region_join_df['Hydrophobicity_original'] - region_join_df['Hydrophobicity_mutated']
region_join_df.to_csv('result.csv')

#assigning hydrophobicity to the substitution yileding the greatest Binding
max = pd.DataFrame(zip(seq,original,aa), columns = ('sequence', 'original', 'mutated'))

join_df = pd.merge(max,hydro, on='original', how = 'left')
join_df = pd.merge(join_df,hydro2, on='mutated', how = 'left')

for i in range(len(join_df['sequence'])):
    join_df['hydrophobicity_difference'.format(i)] = join_df['Hydrophobicity_original'] - join_df['Hydrophobicity_mutated']

#finding out how both the top and the bottom value changes in terms of position
plt.plot(df['Residue #'].drop_duplicates(), max_list)
plt.plot(df['Residue #'].drop_duplicates(), join_df['hydrophobicity_difference'].to_list())
plt.xlabel('Residue Number')
plt.ylabel('Binding Index & hydrophobicity')
plt.legend(('mutation probability', 'hydrophobicity change'))
plt.show()
