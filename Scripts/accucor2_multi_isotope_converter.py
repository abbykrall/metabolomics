# %%
import pandas as pd
from tkinter.filedialog import askopenfilename
import numpy as np

# %% [markdown]
# Upload neg and pos csv files separately from Mzmine 2 with dual isotope labelling, and neg and pos metabolite info sheets separately from SQL DB.
# Gives pos and neg data files in excel format and charge info sheets in csv format ready to upload to Accucor 2 in R

# %%
file_name = 'all_mets_accucor2'
export_dir = "C:/Users/nmatulionis/Desktop/Python_Script_Dump/"

# %%
#Upload csv file after exporting from Mzmine 2 
print("Upload the metabolite csv file!")
path = askopenfilename()

# %%
#Rename met column and drop duplicates and NaNs
df_data = pd.read_csv(path)
df_data = df_data.rename(columns={'row identity (main ID)': 'name'})
df_data = df_data[~df_data['name'].isin(['row identity (main ID)'])]
df_data.drop_duplicates(subset="name", keep=False, inplace=True)
df_data.dropna(axis="columns",inplace=True)

# %%
#Insert columns to match the template of Accucor 2 data file
df_data.insert(loc=1, column='parent', value=0)
df_data.insert(loc=2, column='13C#', value=0)
df_data.insert(loc=3, column='15N#', value=0)
df_data.insert(loc=4, column='Expected?', value=1)

# %%
met_info_path = askopenfilename()
df_met_info = pd.read_csv(met_info_path)




# %%
df_mass_info = df_met_info[['name','mz']].copy()

# %%
#Match parent mass to compound isotopomer
df_data = df_data.merge(df_mass_info, on='name', how='inner')

df_data['parent'] = df_data['mz']
df_data = df_data.drop(columns='mz')


# %%
df_data['13C#'] = df_data['name'].astype(str).str.split(' C').str[1]
df_data['15N#'] = df_data['13C#'].astype(str).str.split('N').str[1]

# %%
df_data['name'] = df_data['name'].astype(str).str.split(' C').str[0]
df_data['13C#'] = df_data['13C#'].astype(str).str.split('N').str[0]

df_data.rename(columns = {'name': 'compoundId'})

# %%
#Charge file creation
df_charge_info = pd.DataFrame(columns=['compound', 'formula', 'charge'])

df_charge_info['compound'] = df_met_info['name']
df_charge_info['formula'] = df_met_info['formula']
df_charge_info['charge'] = 0

if 'pos' in met_info_path:
    df_charge_info['charge'] = 1
    charge_name = 'pos'
elif 'neg' in met_info_path:
    df_charge_info['charge'] = -1
    charge_name = 'neg'
else:
    print('Met info file must include pos or neg in name')
    charge_name = 'something went wrong'

# %%
df_charge_info['compound'] = df_charge_info['compound'].astype(str).str.split(' C').str[0]
df_charge_info.drop_duplicates(inplace=True)
df_charge_info.reset_index(drop=True)

# %%
df_data.to_excel(export_dir + file_name + '_' + charge_name + '.xlsx', index=False)
df_charge_info.to_excel(export_dir + 'charge_info' + '_' + charge_name + '.xlsx', index=False)
print("Files exported to: \n" + export_dir)