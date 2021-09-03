from typing import ChainMap
import biobox as bb
import pandas as pd
import glob
import numpy as np


#This fixes any resid numbering issues that arise as a result of the autopatcher.

def correct_resid(pdb, chain):

#Firtly it opens each patched file in clean.

    cleanfiles = np.array(glob.glob("curate_PDB/clean/*pdb"))

    for cleanfile in cleanfiles:
        name = pdb + '_' + chain + '_patched'
        

        if (cleanfile[17:-4]) == name:
            print('*********FIXING PATCHED FILE RESID*************')
            try:
#Next it gets the first resid.

                M = bb.Molecule()
                M.import_pdb(cleanfile)
                cleandf = M.data
                cleanresid = cleandf.at[0, 'resid']
                print(cleanresid)

#Next it opens the corresponding raw file.

                rawfiles = np.array(glob.glob("curate_PDB/raw/*pdb"))

                for rawfile in rawfiles:
                    name2 = pdb + '_' + chain

#It opens it in biobox and gets the first resid

                    if rawfile[15:-4] == name2:
                        S = bb.Molecule()
                        S.import_pdb(rawfile)
                        rawdf = S.data
                        rawresid = rawdf.at[0, 'resid']
                        
#If they aren't teh same, it shifts every resid in the patched file so they match and writes a new pdb file.

                        if cleanresid != rawresid:
                            print('CHANGING')
                            cleandf['resid'] = cleandf['resid'] + rawresid - 1
                            M.write_pdb(cleanfile)
                        else:
                            return
                        
            except Exception as e:
                print('Error %s'%e)
                print('Issue correcting resid for' + pdb + chain)
    return
        

if __name__ == "__main__":
    try:
        correct_resid('4XBJ-alt1B', 'C')
    except Exception as e:
        print("ERROR: %s"%e)



