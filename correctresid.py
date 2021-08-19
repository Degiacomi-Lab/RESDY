from typing import ChainMap
import biobox as bb
import pandas as pd
import glob
import numpy as np



def correct_resid(pdb, chain):

    cleanfiles = np.array(glob.glob("curate_PDB/clean/*pdb"))
    print(cleanfiles)
    for cleanfile in cleanfiles:
        name = pdb + '_' + chain + '_patched'
        print(name)
        print(cleanfile[17:-4])
        if (cleanfile[17:-4]) == name:
            print('*********FIXING PATCHED FILE RESID*************')

            M = bb.Molecule()
            M.import_pdb(cleanfile)
            cleandf = M.data
            cleanresid = cleandf.at[0, 'resid']
            print(cleanresid)

            rawfiles = np.array(glob.glob("curate_PDB/raw/*pdb"))
            for rawfile in rawfiles:
                name2 = pdb + '_' + chain
                print(rawfile[15:-4])
                if rawfile[15:-4] == name2:
                    S = bb.Molecule()
                    S.import_pdb(rawfile)
                    rawdf = S.data
                    rawresid = rawdf.at[0, 'resid']
                    

                    if cleanresid != rawresid:
                        print('CHANGING')
                        cleandf['resid'] = cleandf['resid'] + rawresid - 1
                        M.write_pdb(cleanfile)
                    else:
                        return
    return
        

if __name__ == "__main__":
    try:
        correct_resid('4XBJ-alt1B', 'C')
    except Exception as e:
        print("ERROR: %s"%e)



