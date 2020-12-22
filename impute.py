#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 11:29:10 2020

@author: dale
"""

import miceforest as mf
import pandas as pd
import numpy as np

raw = pd.read_csv('cul5_proteomics_new.csv')
raw = raw.replace('Filtered',np.nan)

#Extract copy number columns and divide into respective sample groups.
header_names = list(raw)
loc = [i for i, s in enumerate(header_names) if '.Quantity' in s]
raw_cleaned = raw.iloc[:,loc]
raw_cleaned = raw_cleaned.apply(pd.to_numeric)
raw_cleaned['index'] = raw['Genes'].str.upper()
raw_cleaned['genes'] = raw['Genes'].str.upper()
raw_cleaned = raw_cleaned.dropna(subset = ['genes'])
raw_cleaned = raw_cleaned.set_index('index')
KOc = raw_cleaned[['KOc_01.Quantity','KOc_02.Quantity','KOc_03.Quantity']]
KOd = raw_cleaned[['KOd_01.Quantity','KOd_02.Quantity','KOd_03.Quantity']]
WTc = raw_cleaned[['WTc_01.Quantity','WTc_02.Quantity','WTc_03.Quantity']]
WTd = raw_cleaned[['WTd_01.Quantity','WTd_02.Quantity','WTd_03.Quantity','genes']]

median_norm = True

KOd1_sum = KOd.iloc[:,0].median()
KOd2_sum = KOd.iloc[:,1].median()
KOd3_sum = KOd.iloc[:,2].median()
KOc1_sum = KOc.iloc[:,0].median()
KOc2_sum = KOc.iloc[:,1].median()
KOc3_sum = KOc.iloc[:,2].median()
WTd1_sum = WTd.iloc[:,0].median()
WTd2_sum = WTd.iloc[:,1].median()
WTd3_sum = WTd.iloc[:,2].median()
WTc1_sum = WTc.iloc[:,0].median()
WTc2_sum = WTc.iloc[:,1].median()
WTc3_sum = WTc.iloc[:,2].median()

if median_norm == True:
    KOd.iloc[:,0] = np.log2(KOd.iloc[:,0]/KOd1_sum)
    KOd.iloc[:,1] = np.log2(KOd.iloc[:,1]/KOd2_sum)
    KOd.iloc[:,2] = np.log2(KOd.iloc[:,2]/KOd3_sum)
    KOc.iloc[:,0] = np.log2(KOc.iloc[:,0]/KOc1_sum)
    KOc.iloc[:,1] = np.log2(KOc.iloc[:,1]/KOc2_sum)
    KOc.iloc[:,2] = np.log2(KOc.iloc[:,2]/KOc3_sum)
    WTd.iloc[:,0] = np.log2(WTd.iloc[:,0]/WTd1_sum)
    WTd.iloc[:,1] = np.log2(WTd.iloc[:,1]/WTd2_sum)
    WTd.iloc[:,2] = np.log2(WTd.iloc[:,2]/WTd3_sum)
    WTc.iloc[:,0] = np.log2(WTc.iloc[:,0]/WTc1_sum)
    WTc.iloc[:,1] = np.log2(WTc.iloc[:,1]/WTc2_sum)
    WTc.iloc[:,2] = np.log2(WTc.iloc[:,2]/WTc3_sum)
else:
    KOd.iloc[:,1] = np.log2(KOd.iloc[:,1])
    KOd.iloc[:,2] = np.log2(KOd.iloc[:,2])
    KOd.iloc[:,3] = np.log2(KOd.iloc[:,3])
    KOc.iloc[:,1] = np.log2(KOc.iloc[:,1])
    KOc.iloc[:,2] = np.log2(KOc.iloc[:,2])
    KOc.iloc[:,3] = np.log2(KOc.iloc[:,3])
    WTd.iloc[:,1] = np.log2(WTd.iloc[:,1])
    WTd.iloc[:,2] = np.log2(WTd.iloc[:,2])
    WTd.iloc[:,3] = np.log2(WTd.iloc[:,3])
    WTc.iloc[:,1] = np.log2(WTc.iloc[:,1])
    WTc.iloc[:,2] = np.log2(WTc.iloc[:,2])
    WTc.iloc[:,3] = np.log2(WTc.iloc[:,3])
    
print('Data loaded and cleaned.')
def target_append(sample):
    sample1 = sample + '_01.Quantity'
    sample2 = sample + '_02.Quantity'
    sample3 = sample + '_03.Quantity'
    targets = pd.DataFrame({'genes':['target'],sample1:[sample],
                            sample2:[sample],sample3:[sample]})
    return targets

targets = target_append('WTd')
WTd = WTd.append(targets)
WTd = WTd.drop_duplicates(subset=['genes'])
loc = []
for i in range(len(WTd)):
    temp = WTd.iloc[i,:]
    del temp['genes']
    temp = temp.dropna()
    if len(temp) > 1:
        loc = loc + [i]
        
WTd = WTd.iloc[loc,:]
WTd = WTd.set_index('genes')
WTd = WTd.transpose()

print('Creating kernel.')
# Create kernel. 
kds = mf.KernelDataSet(
  WTd,
  save_all_iterations=True,
)

print('Imputing.')
# Run the MICE algorithm for 3 iterations
kds.mice(3)

# Return the completed kernel data
completed_data = kds.complete_data()