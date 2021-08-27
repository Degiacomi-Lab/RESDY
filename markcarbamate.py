import pandas as pd

from numpy import true_divide

avg_df = pd.read_csv('Output/results_likely.csv')
avg_df = avg_df.drop(['Unnamed: 0'], axis=1)

carbam_df = pd.read_csv('Output/results_2608/results_carbamates_marked.csv')
carbam_df = carbam_df.drop(['Unnamed: 0'], axis=1)
carbam_df_true = carbam_df.where(carbam_df['carbamylated'] == True)
carbam_df_true = carbam_df_true[carbam_df_true['carbamylated'].notna()]
resid_list = carbam_df_true['resid'].to_list()
pdb_list = carbam_df_true['PDB Code'].to_list()
key_list = []
data = [carbam_df_true['PDB Code'], carbam_df_true['resid'], carbam_df_true['carbamylated']]
headers = ['PDB Code', 'resid', 'carbamylated']
map_df = pd.concat(data, axis=1, keys=headers)
print(map_df)

merge_df = map_df.merge(avg_df, how='outer', left_on=['PDB Code', 'resid'], right_on=['PDB Code', 'resid'])
merge_df['carbamylated'].fillna(False, inplace=True)
merge_df.to_csv('Output/results_2708_marked.csv')
