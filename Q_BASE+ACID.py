#%%

import numpy as np
import scipy
import matplotlib.pyplot as plt
import pandas as pd
import math
import os
from sympy import symbols, Eq, solve, nsolve, re,S, Abs
from matplotlib import colors
colors.XKCD_COLORS
from matplotlib import patheffects
from scipy.optimize import curve_fit
from matplotlib import cm, ticker
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
from matplotlib import font_manager
from matplotlib import rcParams
from matplotlib.patches import FancyBboxPatch
from scipy.optimize import minimize



DIC = 2.5e-3       # mol/L
k_co2 = 10**(-6.4)
k_co3 = 10**(-10.3)

#original ocean concentration:
i = 8.1
H = 10 ** (-i)
OH = 10 ** (-14) / H
CO2_0 = DIC / (1 + k_co2 / H + k_co2*k_co3 / (H) ** 2)
HCO3_0 = DIC / (1 + H / k_co2 + k_co3 / H)
CO3_0 = DIC / (1 + H / k_co3 + (H) ** 2 / k_co2*k_co3)

F = 96485.3321 #s A / mol



# %%
# from scipy.optimize import fsolve
# df = pd.read_excel("base & acid_mol.xlsx")
# pH_values = []
# Q_base = []
# Q_acid = []
#%%
# for index, row in df.iterrows():
#     pH_final = row.iloc[0]  # pH is the first column
#     CO2 = row.iloc[1]  
#     HCO3 = row.iloc[2]
#     CO3 = row.iloc[3]
#     CaMgCO32 = row.iloc[6]
#     MgOH2 = row.iloc[7]  # MgOH2 is the 8th column
#     CaCO3 = row.iloc[8]
#     CaOH2 = row.iloc[9]
#     Q = (2 * MgOH2 + 2 * CaOH2 + 2*CaMgCO32 + CaCO3 + CO3- CO3_0 - (CO2-CO2_0) - 10**(8.1-14) + 10**(pH_final-14))*F
#     if Q < 0:
#         Q = None
#         Q1 = ((CO2-CO2_0) -CO3+ CO3_0 -2*CaMgCO32 + 10**(-pH_final) - 10**(-8.1))*F

#     else:
#         Q1 = None 
    
#     pH_values.append(pH_final)
#     Q_base.append(Q)
#     Q_acid.append(Q1)


# df_output = pd.DataFrame({
#     'pH': pH_values,
#     'Q_base': Q_base,
#     'Q_acid': Q_acid,
# })

# df_output.to_excel("output_Q.xlsx", index=False)

#%%
from scipy.optimize import fsolve

# Read the Excel file
df = pd.read_excel("base & acid_mol.xlsx")

# Initialize lists to store results
pH_values = []
Q_base = []
Q_acid = []
CO2_values = []
HCO3_values = []
CO3_values = []
CaMgCO32_values = []
MgOH2_values = []
CaCO3_values = []
CaOH2_values = []

# Function to extract relevant columns from a row
def extract_columns(row):
    pH_final = row.iloc[0]  # pH is the first column
    CO2 = row.iloc[1]  
    HCO3 = row.iloc[2]
    CO3 = row.iloc[3]
    CaMgCO32 = row.iloc[6]
    MgOH2 = row.iloc[7]  # MgOH2 is the 8th column
    CaCO3 = row.iloc[8]
    CaOH2 = row.iloc[9]
    return pH_final, CO2, HCO3, CO3, CaMgCO32, MgOH2, CaCO3, CaOH2

# Iterate over each row in the DataFrame
for index, row in df.iterrows():
    pH_final, CO2, HCO3, CO3, CaMgCO32, MgOH2, CaCO3, CaOH2 = extract_columns(row)
    
    # Calculate Q and Q1
    Q = (2 * MgOH2 + 2 * CaOH2 + 2 * CaMgCO32 + CaCO3 + CO3 - CO3_0 - (CO2 - CO2_0) - 10**(8.1 - 14) + 10**(pH_final - 14)) * F
    if Q < 0:
        Q = None
        Q1 = ((CO2 - CO2_0) - CO3 + CO3_0 - 2 * CaMgCO32 + 10**(-pH_final) - 10**(-8.1)) * F
    else:
        Q1 = None
    
    # Append results to lists
    pH_values.append(pH_final)
    Q_base.append(Q)
    Q_acid.append(Q1)
    CO2_values.append(CO2)
    HCO3_values.append(HCO3)
    CO3_values.append(CO3)
    CaMgCO32_values.append(CaMgCO32)
    MgOH2_values.append(MgOH2)
    CaCO3_values.append(CaCO3)
    CaOH2_values.append(CaOH2)

# Create a DataFrame from the results
df_output = pd.DataFrame({
    'pH': pH_values,
    'Q_base': Q_base,
    'Q_acid': Q_acid,
    'CO2': CO2_values,
    'HCO3': HCO3_values,
    'CO3': CO3_values,
    'CaMgCO32': CaMgCO32_values,
    'MgOH2': MgOH2_values,
    'CaCO3': CaCO3_values,
    'CaOH2': CaOH2_values,
})

# Write the results to an Excel file
df_output.to_excel("output_Q.xlsx", index=False)
#%%

#%%
#arrange Q-pH data, Q from small to large
df = pd.read_excel("output_Q.xlsx")

mask_base = df['Q_base'].notnull()
df_base = df.rename(columns = {'pH':'pH_base'})
df_base = df_base[mask_base][['Q_base', 'pH_base']]
df_base.to_excel("Q_base.xlsx", index=False)


mask_acid = df['Q_acid'].notnull()
df_acid = df.rename(columns = {'pH':'pH_acid'})
df_acid = df_acid[mask_acid][['Q_acid', 'pH_acid']]
df_acid = df_acid.sort_values(by='Q_acid')
df_acid.to_excel("Q_acid.xlsx", index=False)


df_base = pd.read_excel("Q_base.xlsx")
Q_base = pd.Series(df_base['Q_base'].values)
pH_base = pd.Series(df_base['pH_base'].values)

df_acid = pd.read_excel("Q_acid.xlsx")
Q_acid = pd.Series(df_acid['Q_acid'].values)
pH_acid = pd.Series(df_acid['pH_acid'].values)

df = pd.concat([Q_base, pH_base, Q_acid, pH_acid], axis=1)
df.columns = ['Q_base', 'pH_base', 'Q_acid', 'pH_acid']
df.to_excel("Q_pH.xlsx", index=False)



#%%
##find Q_base more near to which Q_acid
df_base = pd.read_excel("Q_pH.xlsx")
def find_nearest(row):
    q_base = row['Q_base']
    differences = np.abs(df_acid['Q_acid'] - q_base)
    nearest_index = differences.idxmin()
    nearest_q_acid = df_acid.loc[nearest_index, 'Q_acid']
    return nearest_q_acid

df_base['nearest_Q_acid'] = df_base.apply(find_nearest, axis=1)
df_base.rename(columns={'nearest_Q_acid':'Q'},inplace=True)
base_copy = df['pH_base'].copy()
acid_copy = df['pH_acid'].copy()
df_base['pH.base'] = base_copy
df_base['pH.acid'] = acid_copy

print(df_base.head())
df_base.to_excel("Q_pH.xlsx", index=False)

#%%
import pandas as pd

q_ph_df = pd.read_excel('Q_pH.xlsx')
output_q_df = pd.read_excel('output_Q.xlsx')
q_ph_data = q_ph_df.iloc[:, [4, 5, 6]]  # Selecting 5th, 6th, 7th columns

q_base_products = []
q_acid_products = []

for _, row in q_ph_data.iterrows():
    q_value = row[0]  # Q value (5th column in Q_pH.xlsx)
    base_ph_value = row[1]  # 6th column in Q_pH.xlsx
    acid_ph_value = row[2]  # 7th column in Q_pH.xlsx
    
    # Check for matching pH in output_Q.xlsx for base products
    base_match = output_q_df[output_q_df.iloc[:, 0] == base_ph_value]
    if not base_match.empty:
        base_match = base_match.copy()
        base_match['Q'] = q_value
        q_base_products.append(base_match)
    
    # Check for matching pH in output_Q.xlsx for acid products
    acid_match = output_q_df[output_q_df.iloc[:, 0] == acid_ph_value]
    if not acid_match.empty:
        acid_match = acid_match.copy()
        acid_match['Q'] = q_value
        q_acid_products.append(acid_match)

# Concatenate the results and sort by Q value
q_base_products_df = pd.concat(q_base_products).sort_values(by='Q')
q_acid_products_df = pd.concat(q_acid_products).sort_values(by='Q')

q_base_products_df.to_excel('Q_base_products.xlsx', index=False)
q_acid_products_df.to_excel('Q_acid_products.xlsx', index=False)
# Display the resulting dataframes




#%%
# Extract columns 4 to 10 as product columns from both acid and base dataframes
q_acid_merged = q_acid_products_df[['Q', 'pH'] + list(q_acid_products_df.columns[3:10])].copy()
q_base_merged = q_base_products_df[['Q', 'pH'] + list(q_base_products_df.columns[3:10])].copy()

# Rename the columns to differentiate between acid and base product data
q_acid_merged.columns = ['Q', 'pH_acid'] + list(q_acid_merged.columns[2:])
q_base_merged.columns = ['Q', 'pH_base'] + list(q_base_merged.columns[2:])

# Concatenate the acid and base product dataframes on Q
merged_df = pd.concat([q_acid_merged, q_base_merged], axis=1)

merged_df.to_excel('Q_merge.xlsx', index=False)



#%%

##Q-pH- products
df1 = pd.read_excel("Q_pH.xlsx")
df2 = pd.read_excel("base & acid_mol.xlsx")
df = pd.concat([df1, df2], axis=1)
df = df.drop(columns=['Q_base','pH_base', 'Q_acid','pH_acid'])
df.to_excel("Q_products.xlsx", index=False)
df = pd.read_excel("Q_products.xlsx")
##creat a new to save the products
##Q-pH base- products
df_result = pd.DataFrame()

for i, value in df['pH.base'].items():
    # find pH.base = pH
    df_same = df[df['pH'] == value]
    df_extracted = df_same.loc[:, 'CO2':'MgCO3']
    df_result = pd.concat([df_result, df_extracted])

df_result.reset_index(drop=True, inplace=True)
df.reset_index(drop=True, inplace=True)
df_result.insert(0, 'Q', df['Q'])
df_result.insert(1, 'pH.base', df['pH.base'])
print(df_result.head())
df_result.to_excel("Q_base_products.xlsx", index=False)



##Q-pH acid- products

df_result = pd.DataFrame()

for i, value in df['pH.acid'].items():
    # find pH.base = pH
    df_same = df[df['pH'] == value]
    df_extracted = df_same.loc[:, 'CO2':'MgCO3']
    df_result = pd.concat([df_result, df_extracted])

df_result.reset_index(drop=True, inplace=True)
df.reset_index(drop=True, inplace=True)
df_result.insert(0, 'Q', df['Q'])
df_result.insert(1, 'pH.acid', df['pH.acid'])
print(df_result.head())
df_result.to_excel("Q_acid_products.xlsx", index=False)
# %%
