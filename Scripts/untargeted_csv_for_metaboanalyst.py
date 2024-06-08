import pandas as pd
import os
from tkinter import Tk
import re
from tkinter.filedialog import askopenfilename

# Hide the main tkinter window
Tk().withdraw()

def replace_p_with_dot(string):
    return re.sub(r'(\d)p(\d)', r'\1.\2', string)

# Function to modify the CSV file
def modify_csv(filepath):
    # Check if the file exists and is a CSV file
    if not os.path.isfile(filepath) or not filepath.lower().endswith('.csv'):
        raise ValueError("The file must exist and should be a CSV file.")

    df = pd.read_csv(filepath)
    
    # Process column names
    df.columns = [re.sub('___neg.*', '', col) for col in df.columns]
    df.columns = [re.sub('___pos.*', '', col) for col in df.columns]
    df.columns = [replace_p_with_dot(col) for col in df.columns]

    if 'row m/z' not in df.columns or 'row retention time' not in df.columns:
        raise ValueError("The CSV file must contain 'row m/z' and 'row retention time' columns.")

    # Combine 'row m/z' and 'row retention time' into a new 'Sample' column
    df['Sample'] = df['row m/z'].astype(str) + '__' + df['row retention time'].astype(str)
    
    # Remove the original 'row m/z' and 'row retention time' columns
    df.drop(['row m/z', 'row retention time'], axis=1, inplace=True)

    # Move 'Sample' column to the first position
    cols = df.columns.tolist()
    cols.insert(0, cols.pop(cols.index('Sample')))
    df = df.reindex(columns=cols)
    
    # Remove any unnamed columns
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    
    # Remove duplicate rows based on 'Sample' column keeping the one with the largest sum of the row
    df['row_sum'] = df.iloc[:, 1:].sum(axis=1)
    df = df.sort_values('row_sum', ascending=False).drop_duplicates('Sample').drop('row_sum', axis=1)

    # Insert the 'Group' row after the header
    group_values = ['Group'] + ['insert group'] * (len(df.columns) - 1)
    group_df = pd.DataFrame([group_values], columns=df.columns)
    df = pd.concat([group_df, df.iloc[1:]]).reset_index(drop=True)
    df.at[0, 'Sample'] = 'Group'  # Ensure the first cell is 'Group'

    # Create the new output path by appending '_metaboanalyst' before the file extension
    base, extension = os.path.splitext(filepath)
    output_path = base + '_metaboanalyst' + extension
    df.to_csv(output_path, index=False)

    return output_path

# Ask the user to select a file
filepath = askopenfilename()
if filepath:
    output_path = modify_csv(filepath)
    print(f"Modified file saved to: {output_path}")
else:
    print("No file selected.")