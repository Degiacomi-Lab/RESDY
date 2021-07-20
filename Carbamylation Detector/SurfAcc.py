import numpy as np
import biobox as bb
import pandas as pd
from biobox.measures.calculators import sasa

M = bb.Molecule()
M.import_pdb('1HPX.pdb')
df = M.data

only_lys_CA = df.where((df['resname'] == 'LYS') & (df['name'] == 'CA'))
only_lys_CA = only_lys_CA.dropna()
resid_list = only_lys_CA['resid'].tolist()

#idx = M.atomselect("*", "LYS", "NZ", use_resname=True, get_index=True)
#print(M)

#M = bb.Molecule("1UBQ.pdb")
#For loop here
for ResidueID in resid_list:
    ResidueID = str(ResidueID)
    pts, indices = M.atomselect("*",  ResidueID, ["CA"],  use_resname=False, get_index=True)
# indices will contain the position of atoms we care about
#print(indices)
    x = sasa(M, targets=indices, probe=0.5, n_sphere_point=960, threshold=0.05)
    print(x[0])
#s = bb.Structure(p=x[1])
#s.write_pdb('cloud.pdb')
