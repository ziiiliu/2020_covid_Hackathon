import pandas as pd
import csv
from scipy.stats import chisquare

input_filename  = "/Users/aldricgoh/COVID_Hackathon/Datasets/Patient_lung_biopsies/normalized_human_biopsy_data.csv"
output_filename_1 = "/Users/aldricgoh/human_biopsy_analysis_upreg.csv"
output_filename_2 = "/Users/aldricgoh/human_biopsy_analysis_downreg.csv"

#Read files
df = pd.read_csv(input_filename, usecols = ['gene_symbol','HealthyLungBiopsy_2','HealthyLungBiopsy_1', "COVID19Lung_2", 'COVID19Lung_1'])
#Sort file by gene_symbol
df.sort_values(by='gene_symbol', inplace=True) 

#Function to write csv file
def write_csv(filename, header, data):
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(data)

#Insert values into lists
genes = df['gene_symbol'].values
healthy_2 = df['HealthyLungBiopsy_2'].values
healthy_1 = df['HealthyLungBiopsy_1'].values
covid_2 = df['COVID19Lung_2'].values
covid_1 = df['COVID19Lung_1'].values

#Set up empty lists
sig_diff_down = []
sig_diff_up = []

for i in range(len(genes)):
    if healthy_2[i] == 0 and healthy_1[i] == 0 and covid_2[i] == 0 and covid_1[i] == 0:
        pass
    
    else:
        #Check is there is significant difference between the 2 healthy lungs genes themselves,
        #and between the 2 dammaged lung genes themselves
        #Only proceed when no significant differences
        #May encounter division by 0, but can be ignored
        #At a 0.1% level
        if chisquare([healthy_1[i], healthy_2[i]])[0] < 6.635 and chisquare([covid_1[i], covid_2[i]])[0] < 6.635:
            a = (healthy_1[i] + healthy_2[i])/2
            b = (covid_1[i] + covid_2[i])/2
            chi = chisquare([a, b])[0] #Degree of freedom of 1, and obtain value os chisquare

            #If difference is significant, only insert into list for monitoring
            #At a 0.5% level
            if chi > 7.879:
                if (healthy_2[i] + healthy_1[i])/2 < (covid_2[i] + covid_1[i])/2:
                    sig_diff_up.append([genes[i], healthy_2[i], healthy_1[i], covid_2[i], covid_1[i], chi])
                else:
                    sig_diff_down.append([genes[i], healthy_2[i], healthy_1[i], covid_2[i], covid_1[i], chi])
        
        else:
            pass

    
# sort the chisquare value
def take_chi(elem):
    return elem[5]
sig_diff_up = sorted(sig_diff_up, key=take_chi,reverse=True)
sig_diff_down = sorted(sig_diff_down, key=take_chi,reverse=True)

print(len(sig_diff_up))
print(len(sig_diff_down))

#Output data to csv file
write_csv(filename=output_filename_1, header=['gene_symbol', 'HealthyLungBiopsy_2','HealthyLungBiopsy_1', "COVID19Lung_2", 'COVID19Lung_1', 'Chisquare_Value'], data=sig_diff_up)
write_csv(filename=output_filename_2, header=['gene_symbol', 'HealthyLungBiopsy_2','HealthyLungBiopsy_1', "COVID19Lung_2", 'COVID19Lung_1', 'Chisquare_Value'], data=sig_diff_down)