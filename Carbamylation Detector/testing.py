import numpy as np
import biobox as bb
import pandas as pd
from biobox.measures.calculators import sasa

M = bb.Molecule()
M.import_pdb('3bg3.pdb')
df = M.data

df['resid'] = pd.to_numeric(df['resid'])
x = (df['resid'])[0]

if x != 1:
    df['resid'] = df['resid'] - (x-1)
    print(df)