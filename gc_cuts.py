"""A python file to create a list of clusters from the harris catalog,
filtered by mass, dec, and distance

Created on Thur Jan 30 5:30 PM 2025

@author: Annika Deutsch
@date: 1/30/25
@title: gc_cuts.py
@CPU: MacBook Pro Apple M3
@Operating System: Sonoma 14.6.1
@Interpreter and version no.: python 3.12.2
"""

#--------------------imports and function definitions-----------------------

import astropy.table as at
from astropy.io import ascii
from astropy.table import Table
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Function to parse declination strings and convert to decimal degrees
def parse_declination(dec_str):
    """Parse declination strings and convert to decimal degrees

    Parameters:
    dec_str (str): Declination in the format of degrees, minutes, and seconds, 
                   separated by spaces (e.g., "-30 15 20.5" or "+10 45 00").

    Returns:
    float: Declination converted to decimal degrees. A negative value represents 
           southern declinations, and a positive value represents northern declinations.
    """
    sign = -1 if dec_str.startswith('-') else 1
    dec_str = dec_str.strip('+-')  # Remove leading '+' or '-' for splitting
    parts = dec_str.split()
    degrees = int(parts[0])
    minutes = int(parts[1])
    seconds = float(parts[2])
    return sign * (degrees + minutes / 60 + seconds / 3600)



#------------------------------Read in Harris file and parse-------------------------------

# Read the file and extract the relevant part
with open('/Users/adeutsch/Desktop/UVA/25 Spring/research-pulsar/harris_cat.txt', 'r') as f:
    lines = f.readlines()

# Find the start and end indices for Part I data---------------------
start_index_part1 = lines.index('   ID        Name           RA   (2000)   DEC         L       B     R_Sun  R_gc    X     Y     Z\n') - 1
end_index_part1 = lines.index('_________________________________________________________________________________________________\n')

# Read the data part into an Astropy Table with fixed-width columns
col_names_part1 = ['ID', 'Name', 'RA(2000)', 'DEC', 'L', 'B', 'R_Sun', 'R_gc', 'X', 'Y', 'Z']
col_starts_part1 = [0, 10, 25, 38, 51, 58, 66, 73, 80, 87, 93]
col_ends_part1 = [9, 24, 37, 50, 57, 65, 72, 79, 86, 92, None]

part1 = ascii.read(lines[start_index_part1:end_index_part1], format='fixed_width', names=col_names_part1, col_starts=col_starts_part1, col_ends=col_ends_part1)

# Find the start and end indices for Part II data--------------------
start_index_part2 = lines.index('   ID       [Fe/H] wt  E(B-V) V_HB (m-M)V V_t   M_V,t   U-B   B-V   V-R   V-I  spt   ellip\n') - 1
end_index_part2 = lines.index('            Part III:  Velocities and Structural Parameters\n')

# Read the data part into an Astropy Table with fixed-width columns
col_names_part2 = ['ID', '[Fe/H]', 'wt', 'E(B-V)', 'V_HB', '(m-M)V', 'V_t', 'M_V,t', 'U-B', 'B-V', 'V-R', 'V-I', 'spt', 'ellip']
col_starts_part2 = [0, 10, 18, 22, 28, 34, 40, 46, 53, 60, 68, 74, 78, 82] 
col_ends_part2 = [9, 17, 21, 27, 33, 39, 45, 52, 59, 65, 73, 79, 81, None] 

part2 = ascii.read(lines[start_index_part2:end_index_part2], format='fixed_width', names=col_names_part2, col_starts=col_starts_part2, col_ends=col_ends_part2)

# Find the start and end indices for Part III data--------------------
start_index_part3 = lines.index('    ID         v_r   +/-    v_LSR    sig_v  +/-    c        r_c   r_h    mu_V   rho_0 lg(tc) lg(th)\n') - 1
end_index_part3 = lines.index(' NGC 7492    -177.5   0.6  -176.2     1.2   1.0   0.72      0.86  1.15   20.68   1.27   9.60  9.44\n') + 1

# Read the data part into an Astropy Table with fixed-width columns
col_names_part3 = ['ID', 'v_r', '+/-', 'v_LSR', 'sig_v', 'p/m', 'c', 'c_c', 'r_c', 'r_h', 'mu_V', 'rho_0', 'lg(tc)', 'lg(th)']
col_starts_part3 = [0, 12, 20, 26, 33, 41, 47, 54, 60, 66, 73, 80, 87, 92]
col_ends_part3 = [11, 19, 25, 32, 39, 46, 53, 59, 65, 72, 79, 86, 91, None]

part3 = ascii.read(lines[start_index_part3:end_index_part3], format='fixed_width', names=col_names_part3, col_starts=col_starts_part3, col_ends=col_ends_part3)

# Now join all 3 parts into one big table

# Convert part2 and part3 to dictionaries with ID as the key
part2_dict = {row['ID']: row for row in part2}
part3_dict = {row['ID']: row for row in part3}

# Prepare data for the combined table
combined_data = []

for row in part1:
    combined_row = {**row}  # Start with part1 row data
    
    # Add data from part2 if available
    if row['ID'] in part2_dict:
        for col in part2.colnames:
            if col != 'ID':  # Skip the ID column to avoid duplication
                combined_row[col] = part2_dict[row['ID']][col]
    
    # Add data from part3 if available
    if row['ID'] in part3_dict:
        for col in part3.colnames:
            if col != 'ID':  # Skip the ID column to avoid duplication
                combined_row[col] = part3_dict[row['ID']][col]
    
    combined_data.append(combined_row)

# Define column names and create the combined table
combined_col_names = list(part1.colnames) + [col for col in part2.colnames if col != 'ID'] + [col for col in part3.colnames if col != 'ID']
combined_table = Table(rows=combined_data, names=combined_col_names)

# update the cluster ID's to match those in the Baumgardt masses
combined_table['ID'] = [s.replace(" ", "_") for s in combined_table['ID']]



#-------------------------------Read in the Baumgardt file, parse, and combine-------------------------------
# Read the file and extract the relevant part
with open('/Users/adeutsch/Desktop/UVA/25 Spring/research-pulsar/baumgardt_masses.txt', 'r') as f:
    lines = f.readlines()

# Find the start and end indices for the data section
start_index = 2  # Data starts after the second line (header lines)
end_index = len(lines)  # Read until the end of the file

# Read the data part into an Astropy Table with fixed-width columns
col_names = ['Cluster', 'RA', 'DEC', 'R_Sun', 'DRSun', 'R_GC', 'DRGC', 'N_RV', 'N_PM', 'Mass']
col_starts = [0, 14, 25, 36, 43, 52, 59, 65, 71, 78]
col_ends = [13, 24, 35, 42, 51, 58, 64, 70, 77, 90]

# Reading the specific data into the table
table = ascii.read(
    lines[start_index:end_index],
    format='fixed_width',
    names=col_names,
    col_starts=col_starts,
    col_ends=col_ends
)

# Keep only 'Cluster' and 'Mass' columns
baumgardt_masses = table['Cluster', 'Mass']

# Make the mass column be log mass
baumgardt_masses['Mass'] = np.log10(baumgardt_masses['Mass'].astype(float))

log_masses = np.zeros(len(combined_table))
combined_table['Log_Mtot'] = log_masses

# Add any masses from the one table to the other:
for i in range(0, len(combined_table)):
    for j in range(0, len(baumgardt_masses)):
        if combined_table[i]['ID'] == baumgardt_masses[j]['Cluster']:
            combined_table[i]['Log_Mtot'] = baumgardt_masses[j]['Mass']


#-----------------------Filter by dec, dist, mass, and save to excel file-----------------------
# Parse declination values and create a new column with parsed values
combined_table['DEC'] = [parse_declination(dec) for dec in combined_table['DEC']]

# Filter rows with user specified declination: 
user_dec_min_input = input("Enter the minimum (degrees) (press Enter to skip): ")
user_dec_max_input = input("Enter the maximum (degrees) (press Enter to skip): ")

# Convert inputs to floats or use default values
user_dec_min = float(user_dec_min_input) if user_dec_min_input else -np.inf
user_dec_max = float(user_dec_max_input) if user_dec_max_input else np.inf

filtered_dec = combined_table[(combined_table['DEC'] >= user_dec_min) & (combined_table['DEC'] <= user_dec_max)]

# Convert to Astropy Table
filtered_dec_table = at.Table(filtered_dec)

# Ask user for max distance:
user_dist_input = input("Enter the maximum distance (kpc) (press Enter to skip): ")
user_dist = float(user_dist_input) if user_dist_input else np.inf

# filter out clusters with a distance >= user input
filtered_dist_table = filtered_dec_table[filtered_dec_table['R_Sun'] <= user_dist]

# filter out GMRT proposed clusters 
# GMRT_clusters = ['NGC 5986', 'NGC 6624', 'NGC 6626', 'NGC 1851', 'NGC 6121', 'NGC 6266', 'NGC 6656']
# mask = ~np.isin(filtered_dist_table['ID'], GMRT_clusters)
# filtered_GMRT_table = filtered_dist_table[mask]

# to make easier for comparing tables, remove spaces from cluster IDs: 
filtered_dist_table['ID'] = [s.replace(" ", "_") for s in filtered_dist_table['ID']]
filtered_dist_table['ID'] = [s.replace("zan", "") for s in filtered_dist_table['ID']]

# remove clusters with known log_masses < 5, sort by declination, and only keep relevant columns (final_table_short)

# Ask user for minimum log mass:
user_mass_input = input("Enter the minimum log total mass: (press Enter to skip): ")
user_mass = float(user_mass_input) if user_mass_input else -np.inf

final_table = filtered_dist_table[filtered_dist_table['Log_Mtot'] >= user_mass]
final_table.sort('DEC')

# save as pandas df and excel:
excel_name_input = input("Enter the name to save the cluster excel sheet as: ")
excel_name = excel_name_input + ".xlsx"

df_final = final_table.to_pandas()
df_final.to_excel(excel_name, index=False)