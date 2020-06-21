import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import csv
import scipy.stats as stats

#loading data
name = 'covid_scRNA_metadata.csv'
df = pd.read_csv(name)

#calculate ratio between total_RNA_counts and total_number_genes
cons = df['Condition'].drop_duplicates()
celltypes = df['Celltype_tentative'].drop_duplicates()
patients = df['Patient'].drop_duplicates()
for i in range(len(df['Condition'])):
    df['Ratio'.format(i)] = df['total_RNA_counts']/df['total_number_genes']

#conduct one-way ANOVA to evaluate the significance of variabilities of each celltype across conditions
dict = {}
i = 0
dict2 = {}

# calculate one-way ANOVA significance
for celltype in celltypes:
    df1 = df[df['Celltype_tentative'] == celltype]
    F,p = stats.f_oneway(df1[df1['Condition'] == 'Mild COVID19']['Ratio'],
                        df1[df1['Condition'] == 'Severe COVID19']['Ratio'],
                        df1[df1['Condition'] == 'Healthy control']['Ratio'])

    F1,p1 = stats.f_oneway(df1[df1['Condition'] ==  'Mild COVID19']['percent_mitochondrial'],
                        df1[df1['Condition'] == 'Severe COVID19']['percent_mitochondrial'],
                        df1[df1['Condition'] == 'Healthy control']['percent_mitochondrial'])
    dict[celltype] = [p]
    dict[celltype].append(p1)

pd = pd.DataFrame(dict.items(),columns = ('Celltype',
                'significance level for RNA counts to gene number ratio and significance level for mitochondrial gene proportion'))

pd.to_csv('significance value.csv')



counts = df.groupby(['Patient','Condition','Celltype_tentative']).size().reset_index().rename(columns={0:'cell_number'})
count1 = df.groupby(['Patient','Condition','Celltype_tentative']).size()
countpercent = count1.groupby(level=0).apply(lambda x: 100*x /x.sum()).reset_index().rename(columns={0:'cell_percentage'})
countpercent = countpercent.groupby(['Condition', 'Celltype_tentative'])['cell_percentage'].mean().reset_index()
countpercent.to_csv('percentage_of_cell.csv')



#catplot of the ratio across celltypes and conditions
sns.set(style = 'ticks')
sns.set_context('paper')
g = sns.catplot(x ='Celltype_tentative', y = 'Ratio', hue='Condition',
                data = df, kind = 'violin',
                height = 10, aspect = 1.0)
plt.show()

#catplot of mitochondria gene proportion
sns.set(style = 'ticks')
sns.set_context('paper')
g = sns.catplot(x ='Celltype_tentative', y = 'percent_mitochondrial', hue='Condition',
                data = df, kind = 'violin',
                height = 10, aspect = 2.0)
plt.show()

#barplot percentage of cellfor different conditions
g = sns.catplot(x = 'Celltype_tentative', y = 'cell_percentage', hue = 'Condition', data = countpercent,
                kind = 'bar',height= 10, aspect =3.0)
plt.show()
