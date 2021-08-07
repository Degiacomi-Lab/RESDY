import inquirer
import urllib.request, urllib.parse, urllib.error
import re
from bs4 import BeautifulSoup
import csv
import pandas as pd
import os
import subprocess
import numpy as np
import biobox as bb
from biobox.measures.calculators import sasa

def calculate_sasa(code, df):
    try:
        if df.empty:
            print('No lysine residues with an epsilon amino group pKa < 9 were found in ' + code + '\n')

        else:
            path = 'assembled/' + code 
            print('The following lysine residues have epsilon amino group pKa value(s) < 9' + '\n')
            print(df)
            M = bb.Molecule()
            M.import_pdb(path)
            df_2 = M.data

            resid_list = df['resid'].tolist()
            chain_list = df['chain'].tolist()

            df.index = pd.RangeIndex(len(df.index))
            df.index = range(len(df.index))
            print('\n Obtaining sasa data for ' + code + '\n')
            list_of_sasa = list()

            for i in range(len(resid_list)):
                try:
                    percentage_through = (i/(len(resid_list)))*100
                    percentage_through = (str(percentage_through))[:4]
                    print(percentage_through + '%')


                    chain = chain_list[i]
                    ResidueID = resid_list[i]
                    ResidueID = [int(ResidueID)]
                    pts, indices = M.atomselect(chain,  ResidueID, ["CB", "CG", "CD", "CE", "NZ"],  use_resname=False, get_index=True)
                    x = sasa(M, targets=indices, probe=1.4, n_sphere_point=960, threshold=0.05)
                    acc_surf_area = str(x[0])
                    list_of_sasa.append(acc_surf_area)
                except Exception as e:
                    print("Error %s"%e)
                    print("Failed to calculate sasa for " + code)
                    df = df.drop(df.index[i])
                    continue




            df['sasa'] = list_of_sasa
            print(df)
                #s = bb.Structure(p = x[0])
                #s.write_pdb('teststructure.pdb')

    except Exception as e:
        print("ERROR: %s"%e)
        print("Failed to calculate sasa for " + code)

    return(df)