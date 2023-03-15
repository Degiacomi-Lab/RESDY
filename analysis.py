import re
import urllib.request, urllib.parse, urllib.error
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import threading
import concurrent.futures
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

class Analysis(object):
    
    def __init__(self, df):
        self.df = df.dropna(subset=['pKa', 'sasa'])
        
        self.df_aggregated = pd.DataFrame(columns = ['Uniprot_Entry','Resid','Num','pKa mean','pKa std','pKa range', 
        'SASA mean','SASA std','SASA range'])
        
        self.df_sub = pd.DataFrame(columns = ['Uniprot_Entry','PDB Code','Method','Resolution','Chain','Resid','pKa','sasa']) 
        
        self.GO_dict = {} # code as the key
        
        self.name_to_code = {}
        self.code_to_name = None
    
    def get_data(self, uniprot_entry, resid):
        df_query = self.df[(self.df['Uniprot_Entry'] == uniprot_entry) & (self.df['Resid'] == resid)]
        return df_query
    
    def get_data_alphafold(self):
        df_query = self.df[self.df['Method'] == 'Predicted']
        return df_query
    

    def GO_search_term(self, df, code = '', name = ''):
        '''
        List the subset of UNIPROT codes associated with a GO Term
        '''
        
        if code == '' and name == '':
            return 'Insufficient input!'
        elif code == '' and name != '':
            code = self.name_to_code[name]
        elif code != '' and name != '':
            # check if they match
            if code != self.name_to_code[name]:
                return f'Unmatched GO term code and name, wrong input code {code}, should be {self.GO_decode_dict[name]}.'
        
        uni_list = self.GO_dict[code]
        df_out = pd.DataFrame()
        for uni in uni_list:
            cdf = df[df['Uniprot_Entry'] == uni]
            df_out = pd.concat([df_out, cdf], ignore_index=True)
            
        return df_out
    
    def GO_search_protein(self, uniprot_entry):
        '''
        List all the GO Terms associated with a UNIPROT code
        '''
        
        GO_list = list()
        for code, uni_list in self.GO_dict.items():
            if uniprot_entry in uni_list:
                GO_list.append(code)
        GO_list = [self.code_to_name[code] for code in GO_list]
        return GO_list
    
    def _GO_get_data(self, uniprot_code, lock, index, total):
        '''
        worker of the self.GO_get_data method
        '''
        
        print(f'Searching for {index}/{total} protein.')
        try:
            url_2 = 'https://www.uniprot.org/uniprot/' + uniprot_code + '.txt'
            html_2 = urllib.request.urlopen(url_2)
        except Exception as e:
            print('Failed to obtain UNIPROT data for %s. %s'%(uniprot_code, e))
        
        try:
            for line in html_2:
                line = str(line)
                find_go = re.findall('GO; GO', line)
                if len(find_go) > 0:
                    # GO; GO:0030089; C:phycobilisome; IEA:UniProtKB-KW.
                    # GO; GO:0102834; F:1-18:1-2-16:0-monogalactosyldiacylglycerol acyl-lipid omega-6 desaturase activity; IEA:UniProtKB-EC.
                    GO_code = line.split('; ')[1].split(':')[1]
                    idx = line.split('; ')[2].index(':') + 1
                    GO_word = line.split('; ')[2][idx:]
                    with lock:
                        if GO_code not in self.GO_dict.keys():
                            self.GO_dict[GO_code] = [uniprot_code]
                        else:
                            self.GO_dict[GO_code].append(uniprot_code)

                        # construct name2code dict
                        if GO_word not in self.name_to_code.keys():
                            self.name_to_code[GO_word] = GO_code
                            
        except Exception as e:
            print(uniprot_code)
            print('Error %s'%e)   
                    
    def GO_get_data(self):
        uniprot_codes = self.df['Uniprot_Entry'].unique()
        num = len(uniprot_codes)
        locks = [threading.Lock()]*num
        with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.map(self._GO_get_data, uniprot_codes, locks, range(num), [num]*num)

        self.code_to_name = {v: k for k, v in self.name_to_code.items()}
        
    def aggregate(self):
        df_temp = self.df.drop_duplicates(subset=['Uniprot_Entry','Resid'])
        for idx, row in df_temp.iterrows():
            entry = row['Uniprot_Entry']
            resid = row['Resid']
            df_query = self.df[(self.df['Uniprot_Entry'] == entry) & (self.df['Resid'] == resid)]
            num = len(df_query)
            pka_mean = round(df_query['pKa'].mean(),2)
            pka_std = round(df_query['pKa'].std(),2)
            pka_range = [round(df_query['pKa'].min(),2),round(df_query['pKa'].max(),2)]
            sasa_mean = round(df_query['sasa'].mean(),2)
            sasa_std = round(df_query['sasa'].std(),2)
            sasa_range = [round(df_query['sasa'].min(),2),round(df_query['sasa'].max(),2)]
           
            data = {'Uniprot_Entry': entry,
                        'Resid': resid,
                        'Num': num,
                        'pKa mean': pka_mean,
                        'pKa std': pka_std,
                        'pKa range': pka_range,
                        'SASA mean': sasa_mean,
                        'SASA std': sasa_std,
                        'SASA range': sasa_range}
            
            self.df_aggregated = pd.concat([self.df_aggregated, pd.DataFrame(data)], ignore_index=True)
                

    def plot_graph(self, plot_type, feature, uniprot_entry = False, resid = False):
        '''
        Basic plot, show either a histogram or a boxplot
        '''
        
        try:
            plt.clf()
        except:
            pass
        
                    
        if not uniprot_entry and not resid:
            try:
                x = self.df[feature]
            except:
                print(f'could not find feature {feature}')
                return
                
            if plot_type == 'histogram':
                sns.displot(x, kde=True)
        
            elif plot_type == 'boxplot':
                sns.boxplot(x=x)
            else:
                print('No Such Plot Available.')
                return
            
        else:    
            if uniprot_entry and resid:
                if uniprot_entry not in self.df['Uniprot_Entry'].unique():
                    print('Wrong Uniprot_Entry.')
                    return
                
                else:
                    if resid not in self.df[self.df['Uniprot_Entry'==uniprot_entry]]['Resid'].unique():
                        print('Wrong Resid.')
                        return

                df_query = self.df[(self.df['Uniprot_Entry'] == uniprot_entry) & (self.df['Resid'] == resid)]
                
                try:
                    x = df_query[feature]
                except:
                    print(f'could not find feature {feature}')
                    return

                if plot_type == 'histogram':
                    sns.displot(x, kde=True)
                elif plot_type == 'boxplot':
                    sns.boxplot(x=x)
                else:
                    print('No Such Plot Available.')
                    return
            else:
                print('Lack of Input Information.')
                return
                
        plt.show()
    
    def get_outliers(self, uniprot_entry, resid, feature, whis = 1.5):
        df_query = self.df[(self.df['Uniprot_Entry'] == uniprot_entry) & (self.df['Resid'] == resid)]
        x = df_query[feature]
        
        Q1 = x.quantile(0.25)
        Q3 = x.quantile(0.75)
        IQR = Q3 - Q1
        
        lower = Q1 - whis * IQR
        upper = Q3 + whis * IQR
        
        df_outlier = df_query[(df_query[feature] < lower) | (df_query[feature] > upper)]
        return df_outlier
    
    def get_extreme_values(self, feature, lower = 1, upper = 14):
        df_query = self.df[(self.df[feature] < lower) | (self.df[feature] > upper)]
        return df_query
    
    def remove_df(self, df_to_remove):
        L1 = len(self.df)
        remove_list = df_to_remove.index.tolist()
        self.df = self.df.drop(index = remove_list)
        L2 = len(self.df)
        print(f'Original num of rows: {L1}\nCurrent num of rows: {L2}\nNum of rows removed: {len(df_to_remove)}')
        
    def subset(self, df, weight = 0.5, method = 'average'):
        
        if method != 'average' and method != 'south_east':
            raise ValueError('Wrong input method, try average or south_east.')
        
        # function to evaluate trade-off between pka and sasa
        def low_pka_large_sasa(df, weight = 0.5):
            if len(df) > 1:
                df = df.reset_index(drop = True)

                # situation 1: same pka and same sasa
                if (len(df['pKa'].unique()) == 1) and (len(df['sasa'].unique()) == 1):
                    return df.iloc[0,:]

                # situation 2: same pka but different sasa
                elif (len(df['pKa'].unique()) == 1) and (len(df['sasa'].unique()) != 1):
                    index = df['sasa'].idxmax()
                    return df.iloc[index,:]

                # situation 3: different pka but same sasa
                elif (len(df['pKa'].unique()) != 1) and (len(df['sasa'].unique()) == 1):
                    index = df['pKa'].idxmin()
                    return df.iloc[index,:]

                # situation 4: different pka and different sasa
                else:
                    pka_max = df['pKa'].max()
                    pka_list = [(pka_max-i) for i in df['pKa'].tolist()]
                    sasa_list = df['sasa'].tolist()

                    # compute mean and std
                    pka_m, pka_std = np.mean(pka_list), np.std(pka_list)
                    sasa_m, sasa_std = np.mean(sasa_list), np.std(sasa_list)

                    # standardize two lists
                    pka_list = (pka_list - pka_m) / pka_std
                    sasa_list = (sasa_list - sasa_m) / sasa_std

                    w_pka = weight
                    w_sasa = 1 - weight
                    index = -1
                    base = 0
                    for i in range(len(pka_list)):
                        weighted_sum = w_pka * round(pka_list[i],2) + w_sasa * round(sasa_list[i],2)
                        if weighted_sum >= base:
                            index = i
                            base = weighted_sum
                    return df.iloc[index,:]
            else:
                return df.iloc[0,:]

        df_temp = df.drop_duplicates(subset=['Uniprot_Entry','Resid'])
        for idx, row in df_temp.iterrows():
            entry = row['Uniprot_Entry']
            resid = row['Resid']
            df_query = df[(df['Uniprot_Entry'] == entry) & (df['Resid'] == resid)]
            
            if method == 'south_east':
                row_to_append = low_pka_large_sasa(df_query, weight = weight)
                self.df_sub = self.df_sub.append(row_to_append, ignore_index = True)
            else:
                pka_mean = df_query['pKa'].mean()
                sasa_mean = df_query['sasa'].mean()
                data = {'Uniprot_Entry':entry, 'PDB Code':np.nan, 'Method':np.nan, 'Resolution':np.nan, 'Chain':np.nan,
                           'Resid':resid, 'pKa':pka_mean, 'sasa':sasa_mean}
                self.df_sub = pd.concat([self.df_sub, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)


    def get_contingency_table(self, GO_code, my_list, reference):
        '''
        Compute contingency table given a GO Term, the list of interest, and a reference list
        '''
        
        BP_list = 0
        for uni in my_list:
            if uni in self.GO_dict[GO_code]:
                BP_list += 1
            
        BP_not_list = len(self.GO_dict[GO_code]) - BP_list
            
        not_BP_list = len(my_list) - BP_list
            
        not_BP_not_list = 0
        for uni in reference:
            if (uni not in my_list) and (uni not in self.GO_dict[GO_code]):
                not_BP_not_list += 1
            
        table = [[BP_list, BP_not_list], [not_BP_list, not_BP_not_list]]
        
        return table
    
    
    def enrichment_analysis(self, pka_range, sasa_range, uniprot_cnt_cutoff=1):
        '''
        Analyse prevalence of GO-terms in sub-regions of the SASA vs pKa graph
        '''
        
        # get all the uniprot codes inside the range and the reference uniprot code list
        pka_l, pka_u = pka_range[0], pka_range[1]
        sasa_l, sasa_u = sasa_range[0], sasa_range[1]
        selected_df = self.df_sub[(self.df_sub['pKa'] >= pka_l) & (self.df_sub['pKa'] <= pka_u)]
        selected_df = selected_df[(selected_df['sasa'] >= sasa_l) & (selected_df['sasa'] <= sasa_u)]
        my_list = selected_df['Uniprot_Entry'].unique()
        reference = self.df_sub['Uniprot_Entry'].unique()
        
        # get all the GO Terms in the background (that are associated with more than uniprot_cnt_cutoff)
        GO_bacgou = [code for code in list(self.GO_dict.keys()) if len(self.GO_dict[code]) >= uniprot_cnt_cutoff]
        
        p_val_dict = {}
        
        # for each GO Term, compute a contingency table
        for GO_code in GO_bacgou:
            # compute contigency table
            table = self.get_contingency_table(GO_code, my_list, reference)
            # compute p values and store them into the dictionary
            oddsratio, pvalue = fisher_exact(table, alternative='greater')
            p_val_dict[GO_code] = pvalue
        
        # sort the dictionary based on p-values
        p_val_dict = dict(sorted(p_val_dict.items(), key = lambda x: x[1]))
        p_val_list = [x[1] for x in p_val_dict.items()]
        
        # adjust p-values using the BH method
        y=multipletests(pvals=p_val_list, alpha=0.05, method="fdr_bh")
        out_dict = {}
        
        if sum(y[0]) != 0: # if there is enrichment
            for i in range(len(y[0])):
                if y[0][i]: # get p-values below 0.05
                    go_code = list(p_val_dict.items())[i][0]
                    raw_p = p_val_list[i]
                    adj_p = y[1][i]
                    out_dict[go_code] = (raw_p, adj_p, self.get_contingency_table(go_code, my_list, reference))
        else:
            for i in range(len(y[0])): # if no enrichment at all, still output the p values
                go_code = list(p_val_dict.items())[i][0]
                raw_p = p_val_list[i]
                adj_p = y[1][i]
                out_dict[go_code] = (raw_p, adj_p, self.get_contingency_table(go_code, my_list, reference))
        
        df = pd.DataFrame(columns = ['GO ID', 'GO Term', 'raw p value', 'FDR',
                             'num in the region', 'num in the bkgd'])
        for data in out_dict.items():
            ID = data[0]
            Term = self.code_to_name[ID]
    
            r_p = round(data[1][0],6)
            fdr = round(data[1][1],6)
    
            n1 = sum(data[1][2][0])
            d1 = sum(data[1][2][0]) + sum(data[1][2][1])                    
            bkgd = str(n1) + '/' + str(d1)
            n2 = data[1][2][0][0]
            d2 = data[1][2][1][0] + data[1][2][0][0]
            reg = str(n2) + '/' + str(d2)
    
            dt = {'GO ID':ID, 'GO Term':Term, 'raw p value':r_p, 'FDR':fdr,
                             'num in the bkgd':bkgd, 'num in the region':reg}
            df_dictionary = pd.DataFrame([dt])
            df = pd.concat([df, df_dictionary], ignore_index=True)
        
        return df
