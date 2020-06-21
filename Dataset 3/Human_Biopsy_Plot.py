#Import modules
import pandas as pd
import csv
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

#Import files
in_filename_1 = "/Users/aldricgoh/human_biopsy_analysis_upreg.csv"
in_filename_2 = "/Users/aldricgoh/human_biopsy_analysis_downreg.csv"

#Read files
df_up = pd.read_csv(in_filename_1, usecols = ['gene_symbol','HealthyLungBiopsy_2','HealthyLungBiopsy_1', "COVID19Lung_2", 'COVID19Lung_1', 'Chisquare_Value'])
df_down = pd.read_csv(in_filename_2, usecols = ['gene_symbol','HealthyLungBiopsy_2','HealthyLungBiopsy_1', "COVID19Lung_2", 'COVID19Lung_1', 'Chisquare_Value'])

genes_up = df_up['gene_symbol'].values
genes_down = df_down['gene_symbol'].values
chi_up = df_up['Chisquare_Value'].values
chi_down = df_down['Chisquare_Value'].values

genes_up = genes_up[:10]
genes_down = genes_down[:10]
chi_up = chi_up[:10]
chi_down = chi_down[:10]

with PdfPages('/Users/aldricgoh/COVID HACK/Human Biopsy Analysis (new).pdf') as pdf:
# Create an array with the position of each bar along the x-axis
    x_pos = np.arange(len(genes_up))

# Produce bar plot
    plt.bar(x_pos, chi_up, align='center', color='r');

# Replace the x ticks with the Tripos name, and rotate labels 30 degrees
    plt.xticks(x_pos, genes_up, rotation=90)

# Add axis labels 
    plt.xlabel('Genes')
    plt.ylabel('Chisquare Count')
    plt.title('Up regulated genes from Human Biopsy Data')
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    plt.bar(x_pos, chi_up, align='center', color='r')
    plt.semilogy(genes_up, chi_up)
    plt.xlabel('Genes')
    plt.ylabel('Chisquare Count')
    plt.xticks(x_pos, genes_up, rotation=90)
    plt.title('Up regulated genes from Human Biopsy Data (Log)')
    plt.yscale('log')
    plt.tight_layout()
    pdf.savefig()
    plt.close()


    x_pos = np.arange(len(genes_down))
    plt.bar(x_pos, chi_down, align='center', color='b');
    plt.xticks(x_pos, genes_down, rotation=90)
    plt.xlabel('Genes')
    plt.ylabel('Chisquare Count')
    plt.title('Down regulated genes from Human Biopsy Data')
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    plt.bar(x_pos, chi_down, align='center', color='b')
    plt.semilogy(genes_down, chi_down)
    plt.xlabel('Genes')
    plt.ylabel('Chisquare Count')
    plt.xticks(x_pos, genes_down, rotation=90)
    plt.title('Down regulated genes from Human Biopsy Data (Log)')
    plt.tight_layout()
    pdf.savefig()
    plt.close()