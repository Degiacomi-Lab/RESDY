import biobox as bb
import pandas as pd
def clean(pdb_code, chain):
    print('CLEANING ' + pdb_code + '_' + chain)
    path = 'curate_PDB/raw/' + pdb_code + '_' + chain + '.pdb'
    M = bb.Molecule()
    M.import_pdb(path)


    df = M.data
    name_column = df['name'].tolist()
    query_list= list()
    for name in name_column:
        names = name.split()
        if (len(names) == 1):
            if len(name) < 4:
                query_list.append(name)
            elif len(name) >= 4:
                if name[-1:] == 'B':
                    query_list.append(name)
                else:
                    pass
            
        if (len(names) != 1):
            if (names[1] == 'B'):
                query_list.append(name)
    pos, idx = M.query('name != ["H"] and name ==' + str(query_list), get_index=True)


    try:
        M.write_pdb("tmp", index=idx)
    except:
        M.write_pdb("tmp", index=idx, split_struc=False)
    return()