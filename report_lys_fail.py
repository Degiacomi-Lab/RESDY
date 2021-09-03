import os
import math
import pandas as pd
import biobox as bb

def remove_problematic(code, propka_lys_fails, chain_resid_near_failed_chain, pka_sasa_res_df):

    def report_lys_fail(code, propka_lys_fails, chain_resid_clash, chain_resid_near_failed_chain):

        if not os.path.exists('Failed_Lysines'):
            os.mkdir('Failed_Lysines')

        code_no_pdb = code[:-4]
        path = 'Failed_Lysines/' + code_no_pdb + '_fails.txt'
        
        f = open(path, 'w')
        for entry in propka_lys_fails:
            f.write(entry)
            f.write('-Failed propka- check propka log for explanation')
            f.write('\n')

        for entry in chain_resid_clash:
            f.write(entry)
            f.write('-Clash error- likely due to autopatch adding a section which causes clashing')
            f.write('\n')

        for entry in chain_resid_near_failed_chain:
            f.write(entry)
            f.write('-Near chain error- likely due to autopatch adding a section which causes chains to overlap')
            f.write('\n')
        return




    def remove_failures(code, propka_lys_fails, chain_resid_near_failed_chain, pka_sasa_res_df):
    #The program then removes residues that have been selected as problematic along the line.
        print('Removing problematic residues...')
        print('\n')


    #Firstly, those which were selected as problematic in the propka report are removed.
        try:
            inverse_boolean_series = ~pka_sasa_res_df.Chain_Resid.isin(propka_lys_fails)
            pka_sasa_res_df = pka_sasa_res_df[inverse_boolean_series]
        except Exception as e:
            print('Error %s'%e)
            print('Error removing problematic propka residues')
            pka_sasa_res_df = pka_sasa_res_df[0:0]
            return(pka_sasa_res_df)


    #Next chain clash is investigated, if two chains clash the relevant residues are removed

        if code[:2] != 'AF':

            try:
                chain_resid_clash = check_clash(code)
                print(chain_resid_clash)
                inverse_boolean_series = ~pka_sasa_res_df.Chain_Resid.isin(chain_resid_clash)
                pka_sasa_res_df = pka_sasa_res_df[inverse_boolean_series]
                print(pka_sasa_res_df)

            except Exception as e:
                print('Error %s'%e)
                print('Error removing clashing atoms')
                pka_sasa_res_df = pka_sasa_res_df[0:0]
                return(pka_sasa_res_df)

    #Next residues are removed if they are exposed to any chain that failed the autopatch (as without the chain they are usually exposed to the results are not reliable).
            try:

                inverse_boolean_series = ~pka_sasa_res_df.Chain_Resid.isin(chain_resid_near_failed_chain)
                pka_sasa_res_df = pka_sasa_res_df[inverse_boolean_series]

            except Exception as e:
                print('Error %s'%e)
                print('Error removing files near failed chain')
                pka_sasa_res_df = pka_sasa_res_df[0:0]
                return(pka_sasa_res_df)
        else:
            chain_resid_clash = list()
        
        print('Residues removed:')
        print(propka_lys_fails)
        print(chain_resid_clash)
        print(chain_resid_near_failed_chain)



        return chain_resid_clash, pka_sasa_res_df



    #This module checks for chain clash.
    #Chain clash occurs when there is an issue during the assembly which means that two chains can partially overlap.
    #Lysines in clashing regions are ignored.

    def check_clash(pdb_code):
        list_of_chains = list()
        list_of_resid = list()

    #Firstly a df is constructed to store results in.

        columns_2 = ['PDB Code', 'Lysine Index', 'Score']
        score_df = pd.DataFrame(columns=columns_2)
        print('BREAKING UP MOLECULE')

    #Next the assembled file is opened in biobox and the data is stored in a df.
        try:
            M = bb.Molecule()
            path = 'assembled/' + pdb_code

            M.import_pdb(path, include_hetatm=True)
            df = M.data

    #Next the corrdinates and index for every lysine is found.

            lys_coords, lys_idx = M.atomselect('*','LYS', 'CA', use_resname=True, get_index=True)

    #Next the chain and resid for each lysine are put into lists.

            for entry in lys_idx:
                chain = df.at[entry, 'chain']
                list_of_chains.append(chain)

            for entry in lys_idx:
                resid = df.at[entry, 'resid']
                list_of_resid.append(resid)

    #Next the coordinates and index for every atom are found.

            all_coords, idx = M.atomselect('*','*','*', get_index=True)

        except Exception as e:
            print('Error %s'%e)
            print('Failure checking chain clash')
            return []

    #The program then cycles through each lysine and finds the distnace to every atom.

        for j in range(len(lys_coords)):
            try:
                lys_chain = list_of_chains[j]
                score = 0

                try:

                    
                    for i in range(len(all_coords)):
                        x_dist = (((lys_coords[j])[0] - (all_coords[i])[0])**2)
                        y_dist = (((lys_coords[j])[1] - (all_coords[i])[1])**2)
                        z_dist = (((lys_coords[j])[2] - (all_coords[i])[2])**2)
                        distance = math.sqrt(x_dist + y_dist + z_dist)
                        
                        #If the distance to an atom on another chain is found the be less than 4.5 an extra point is added to the score.

                        if distance < 4.5:
                            index_of_aa = idx[i]
                            test_chain = df.at[index_of_aa, 'chain']

                            if test_chain != lys_chain:
                                
                                score = score + 1



                except Exception as e:
                    print('Error %s'%e)
                
                #The score for each lysine is then added to a df.
                chain_resid = str(list_of_chains[j]) + str(list_of_resid[j])

                data = ({'PDB Code': pdb_code, 'Lysine Index': lys_idx[j], 'Resid': list_of_resid[j], 'Chain': list_of_chains[j], 'Chain_Resid': chain_resid,'Score': score})
                score_df = score_df.append(data, ignore_index=True)

            except Exception as e:
                print('Error %s'%e)
                print('Failure chekcing chain clash')
                score = 0
                data = ({'PDB Code': pdb_code, 'Lysine Index': lys_idx[j], 'Resid': list_of_resid[j], 'Chain': list_of_chains[j], 'Chain_Resid': chain_resid,'Score': score})
                score_df = score_df.append(data, ignore_index=True)
            
                continue


#If the score is greater than 1 the lysine is removed (actual removal happens in main) from the results.
            try:
                score_df = score_df.where(score_df['Score'] > 1)
                score_df = score_df.dropna()

                chain_resid_clash = score_df['Chain_Resid'].tolist()
            except Exception as e:
                print('Error %s'%e)
                print('Failure chekcing chain clash')
                return []

            return(chain_resid_clash)

    
    chain_resid_clash, pka_sasa_res_df = remove_failures(code, propka_lys_fails, chain_resid_near_failed_chain, pka_sasa_res_df)

    report_lys_fail(code, propka_lys_fails, chain_resid_clash, chain_resid_near_failed_chain)


    return(pka_sasa_res_df)
    

        




