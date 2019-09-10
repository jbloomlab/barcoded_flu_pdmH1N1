"""
This script parses samples from the 'samples.csv' file, and creates a 'simple_samples.csv' file for passing information to cellranger mkfastq.

The format for the cellranger simple-samples.csv file is:
Lane,Sample,Index
1,test_sample,SI-GA-A3

For more information, see https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq

"""

import os
import pandas as pd

#Read in sample information
samples = pd.read_csv('samples.csv',comment="#")

#Only keep relevant columns
samples_clean = samples[['sample','index']]

#Add additional column info
samples_clean.loc[:, 'Lane'] = 1
samples_clean['index'] = 'SI-GA-' + samples_clean.loc[:, 'index'].astype(str)


#Rename columns
samples_clean.columns = ['Sample','Index','Lane']

#Reorder columns
samples_clean = samples_clean[['Lane', 'Sample', 'Index']]

#Write to file
samples_clean.to_csv('simple-samples.csv',index=False)
print('Sample data written to simple-samples.csv. Process done.')