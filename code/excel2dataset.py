import os
import os.path as op
import pandas as pd
import numpy as np
from itertools import islice
from openpyxl import load_workbook
import pickle
from nimare.dataset import Dataset


def get_worksheet_as_df(sheet_name):
    data = sheet_name.values
    cols = next(data)[0:]
    cols_rm_idx = []
    for i, col in enumerate(cols):
        if col == None:
            cols_rm_idx.append(i)
    cols = np.delete(cols, cols_rm_idx)
    data = list(data)
    data = np.delete(data, cols_rm_idx, axis=1)
    data = (islice(r, 0, None) for r in data)
    df = pd.DataFrame(data, columns=cols)
    df = df.replace('na', 'NaN')
    df = df.replace('ns', 'NaN')
    return df


project_directory = '/Users/miriedel/Desktop/GitHub/meta-analysis_implicit-learning'

wb = load_workbook(filename = op.join(project_directory, 'code', 'IL Meta-Analysis Coding for NiMARE 3-4-21.xlsx'))
contrasts_df = get_worksheet_as_df(wb['Papers and Contrasts'])
coordinates_df = get_worksheet_as_df(wb['Coordinates'])

#First, create a dictionary of the the coordinate spreadsheet
#The dictionary is called contrast_dict
#The dictionary is indexed by  paper (or study) ID, e.g.: contrast_dict[2] for paper/study ID 2
#Contrasts can be indexed by the paper (or study) ID and contrast ID, e.g.: contrast_dict[2]['Contrasts'][3] for paper/study ID 2,contrast 3
#The dictionary has the following structure:
#Paper ID
#   -> Author
#   -> Year
#   -> Contrasts
#       -> Contrast ID
#         -> Contrast Name
#         -> X
#         -> Y
#         -> Z
#         -> Cluster Extent (mm^3)
#         -> stat
#         -> statistic type
#         -> Cluster threshold (mm^3)

study_idx = np.where(coordinates_df['Study ID'].notna())[0]
contrast_dict = {}
for i in range(len(study_idx)):
    start_row = study_idx[i]
    if i == len(study_idx) - 1:
        stop_row = np.shape(coordinates_df)[0]
    else:
        stop_row = study_idx[i+1]
    study_info = coordinates_df.iloc[start_row]
    paper_id = study_info['Study ID']
    contrast_id = study_info['Contrast ID']

    if paper_id not in contrast_dict.keys():
        contrast_dict[paper_id] = {}
        contrast_dict[paper_id].update(study_info[['Author', 'Year']].to_dict())

    if 'contrasts' not in contrast_dict[paper_id].keys():
        contrast_dict[paper_id]['contrasts'] = {}

    contrast_dict[paper_id]['contrasts'][contrast_id] = {}
    contrast_info = coordinates_df[['X', 'Y', 'Z', 'Cluster Extent (mm^3)', 'stat', 'statistic type', 'Cluster threshold (mmm^3)']][start_row+1:stop_row]
    contrast_dict[paper_id]['contrasts'][contrast_id].update(contrast_info.to_dict())
    for col in contrast_info.columns:
        contrast_dict[paper_id]['contrasts'][contrast_id][col] = list(contrast_info[col].values)
    for a in ['X', 'Y', 'Z']:
        contrast_dict[paper_id]['contrasts'][contrast_id][a.lower()] = contrast_dict[paper_id]['contrasts'][contrast_id].pop(a)
    contrast_dict[paper_id]['contrasts'][contrast_id]['Contrast Name'] = study_info['Contrast']

#Second, parase the papers and contrasts spreadsheet
study_idx = np.where(contrasts_df['Paper ID'].notna())[0]
study_dict = {}
for i in range(len(study_idx)-1):
    start_row = study_idx[i]
    if i == len(study_idx) - 1:
        stop_row = np.shape(contrasts_df)[0]
    else:
        stop_row = study_idx[i+1]
    study_info = contrasts_df.iloc[start_row]
    paper_id = study_info['Paper ID']
    study_info = study_info.drop(labels=['Paper ID', 'Include Contrast? (0,1)', 'Contrast ID', 'Contrast Name', 'Contrast Type', 'Coordinate Space/Template', 'Sample Size (n)', 'Meta-Analysis Group'])
    study_dict[paper_id] = {}
    #study_dict[paper_id].update(study_info.to_dict())
    contrast_info = contrasts_df[['Include Contrast? (0,1)', 'Contrast ID', 'Contrast Name', 'Contrast Type', 'Coordinate Space/Template', 'Sample Size (n)', 'Meta-Analysis Group']][start_row+1:stop_row]
    study_dict[paper_id]['contrasts'] = {}
    for j, row in contrast_info.iterrows():
        study_dict[paper_id]['contrasts'][row['Contrast ID']] = {}
        study_dict[paper_id]['contrasts'][row['Contrast ID']]['metadata'] = {}
        study_dict[paper_id]['contrasts'][row['Contrast ID']]['labels'] = {}
        study_dict[paper_id]['contrasts'][row['Contrast ID']]['coords'] = {}
        study_dict[paper_id]['contrasts'][row['Contrast ID']]['coords']['space'] = {}

        study_dict[paper_id]['contrasts'][row['Contrast ID']]['metadata'].update(study_info[['Author', 'Year', 'Citation']].to_dict())
        study_dict[paper_id]['contrasts'][row['Contrast ID']]['metadata']['# Foci'] = np.shape(contrast_dict[paper_id]['contrasts'][row['Contrast ID']]['x'])[0]
        study_dict[paper_id]['contrasts'][row['Contrast ID']]['metadata']['Contrast Name'] = contrast_dict[paper_id]['contrasts'][row['Contrast ID']]['Contrast Name']
        study_dict[paper_id]['contrasts'][row['Contrast ID']]['metadata']['sample_sizes'] = row['Sample Size (n)']
        study_dict[paper_id]['contrasts'][row['Contrast ID']]['labels'].update(study_info.drop(labels=['Author', 'Year', 'Citation']).to_dict())
        study_dict[paper_id]['contrasts'][row['Contrast ID']]['labels'].update(row.drop(labels=['Contrast ID', 'Coordinate Space/Template', 'Sample Size (n)']).to_dict())
        study_dict[paper_id]['contrasts'][row['Contrast ID']]['coords']['space'] = row['Coordinate Space/Template']
        study_dict[paper_id]['contrasts'][row['Contrast ID']]['coords'].update(contrast_dict[paper_id]['contrasts'][row['Contrast ID']])

#save the dictionary
with open(op.join(project_directory, 'code', 'dataset_dictionary.pkl'), 'wb') as fo:
    pickle.dump(study_dict, fo)

#now convert to a nimare dataset and save
dset = Dataset(study_dict)
#save the dictionary
with open(op.join(project_directory, 'code', 'nimare_dataset.pkl'), 'wb') as fo:
    pickle.dump(dset, fo)
