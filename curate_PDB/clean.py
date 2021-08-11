from os import remove
import biobox as bb
import pandas as pd
import fileinput
def clean(pdb_code, chain):
    print('CLEANING ' + pdb_code + '_' + chain)
    remove_kcx(pdb_code, chain)
    path = 'curate_PDB/raw/' + pdb_code + '_' + chain + '.pdb'
    M = bb.Molecule()
    M.import_pdb(path, include_hetatm=True)


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
    try:
        for line in fileinput.input('tmp', inplace = 1):
            print(line.replace('BLYS', ' LYS'))

    except Exception as e:
        print('Error %s'%e)

    return()


#This code removes the carbamate from any lysine that is already carbamylated in the .pdb structure.

def remove_kcx(pdb_code, chain):
    try:
    #Firstly it imports the pdb file from raw.
        prot_code = pdb_code + '_' + chain + '.pdb'
        path = 'curate_PDB/raw/' + prot_code
        M = bb.Molecule()
        M.import_pdb(path, include_hetatm=True)

        #Next it removes the carbon and two oxygens from the structure

        idxs = []
        pos ,idx = M.atomignore('*', 'KCX', 'CX', get_index=True, use_resname=True)
        A = M.get_subset(idxs=idx)
        pos, idx = A.atomignore('*', 'KCX', 'OQ1', get_index=True, use_resname=True)
        B = A.get_subset(idxs=idx)
        pos, idx = B.atomignore('*', 'KCX', 'OQ2', get_index=True, use_resname=True)
        C = B.get_subset(idxs=idx)

        #It then rewrites the structure without the carbamate

        C.write_pdb(path, split_struc=False)
        filename = path

        #Lastly it renames the residue as a lysine and relabels it as an atom instead of a hetereoatom.

        with fileinput.FileInput(filename, inplace = True) as f:
            for line in f:
                if("KCX" in line):
                    line = line.replace("HETATM","ATOM  ")
                    line = line.replace('KCX', 'LYS')
                    print(line, end ='') 
                else:
                    print(line, end ='') 
    except:
        return()

    return()


if __name__ == "__main__":
    try:
        print(clean('4E3T', 'A'))
    except Exception as e:
        print("ERROR: %s"%e)