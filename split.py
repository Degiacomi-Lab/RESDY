import glob
import os
import pandas as pd
import subprocess
import numpy as np
import biobox as bb
import re


#This code checks to see if there are multiple conformations available (usually NMR structures have multiple conformations in one file).

def split_confs(pdb):
    
    print('Checking for alternate conformations')
    try:
        #Firstly it downloads the .pdb file.
        cwd = os.getcwd()
        print(cwd)

        os.chdir('curate_PDB/conformations')
        subprocess.check_call("wget https://files.rcsb.org/download/" + pdb + '.pdb', shell=True)

        number = 1
        name = pdb + '-alt-' + str(number) + '.pdb'
    except Exception as e:
        print('Error %s'%e)
        print('Print failure opening file.')
        return
#Next it checks if the phrase 'ENDMDL' is present (i.e. are there  multiple conformations).
    f = open(pdb + '.pdb')
    f_write = open(name, 'w')
    endmdls = list()
    for line in f:
        if re.search('ENDMDL', line):
            endmdls.append(line)
    f = open(pdb + '.pdb')

    if len(endmdls) == 0:
        for line in f:
            f_write.write(line)
        f_write.close()
        os.remove(pdb + '.pdb')

    else:
        for line in f:
            try:
                line = str(line)

                if re.search('MASTER', line):
                    os.remove(name)
                    os.remove(pdb + '.pdb')
                    break

                else:

                    if re.search('ENDMDL', line):
                        f_write.write(line)
                        f_write.close()
                        number = number + 1
                        name = pdb + '-alt-' + str(number) + '.pdb'
                        f_write = open(name, 'w')

                    else:
                        f_write.write(line)
            except:
                continue
    
    os.chdir(cwd)
    cwd = os.getcwd()
    print(cwd)
    print('done')
    return()








def split_confs_2():
    files = np.array(glob.glob("curate_PDB/conformations/*pdb"))
    for f in files:
        print(f)
        ABC_list = ['A', 'B', 'C', 'D']
        ABC_dict = {"A": 0, "B": 0, "C": 0}
        for i in range(len(ABC_list)):
            read = open(f)
            for line in read:
                if (line[:4] == 'ATOM') and (line[16] == ABC_list[i]):
                    ABC_dict[ABC_list[i]] = 1
        

        list_of_values = ABC_dict.values()
        if 1 not in list_of_values:
            continue
        else:
            print('FOUND ALTERNATE CONFORMATION')


        for i in range(len(ABC_list)):
            try:
                if ABC_dict.get(ABC_list[i]) == 1:
                    name = f[:-6] + f[-5] + ABC_list[i] + '.pdb'
                    print(name)
                    f_write = open(name, 'w')

                    target_letter = ABC_list[i]
                    non_target_letters = []
                    for letter in ABC_list:
                        if letter != target_letter:
                            non_target_letters.append(letter)

                    read = open(f)
                    for line in read:
                        
                        if (line[:4] == 'ATOM') and (line[16] == target_letter):
                            newline = line[:16] + ' ' + line[17:]
                            f_write.write(newline)
                            continue
                        if (line[:4] == 'ATOM') and (line[16] in non_target_letters):
                            continue
                        if (line[:4] == 'ATOM'):
                            f_write.write(line)
                        if (line[:6] == 'HETATM'):
                            f_write.write(line)
                        if (line[:3] == 'TER'):
                            f_write.write(line)
                    f_write.close()
            except Exception as e:
                print("Error %s"%e)
        os.remove(f)
    return





