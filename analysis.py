#!/usr/bin/env python
# coding: utf-8

# In[13]:

import re
import urllib.request, urllib.parse, urllib.error
from bs4 import BeautifulSoup
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# In[12]:

class Analysis(object):
    def __init__(self, df):
        self.df = df.dropna(subset=['pKa', 'sasa'])
        
        self.aggregated_df = pd.DataFrame(columns = ['Uniprot Entry','Resid','Num','pKa mean','pKa std','pKa range', 
        'SASA mean','SASA std','SASA range'])
        
        self.df_concise = pd.DataFrame(columns = ['Uniprot Entry','PDB Code','Method','Resolution','Chain','Resid','pKa','sasa']) 
        
        self.GO_dict = {}
        self.GO_decode_dict = {}
        self.GO_decode_dict_reverse = {}
        
    def get_unique_uniprot_entry(self):
        return self.df['Uniprot Entry'].unique()
    
    def get_unique_resid(self, uniprot_entry):
        current_df = self.df[self.df['Uniprot Entry'] == uniprot_entry]
        return current_df['Resid'].unique()
    
    def get_data(self, uniprot_entry, resid):
        current_df = self.df[(self.df['Uniprot Entry'] == uniprot_entry) & (self.df['Resid'] == resid)]
        return current_df
    
    def get_data_alphafold(self):
        current_df = self.df[self.df['Method'] == 'Predicted']
        return current_df
    
    # list all the uniprot codes associated with a GO Term
    def GO_Search_Term(self, df, code = '', name = ''):
        
        if code == '' and name == '':
            return 'Insufficient input!'
        elif code == '' and name != '':
            code = self.GO_decode_dict[name]
        elif code != '' and name != '':
            # check if they match
            if code != self.GO_decode_dict[name]:
                return f'Unmatched GO Term code and name; wrong input code: {code}; Should be: {self.GO_decode_dict[name]}'
        
        uni_list = self.GO_dict[code]
        df_out = pd.DataFrame()
        for uni in uni_list:
            cdf = df[df['Uniprot Entry'] == uni]
            df_out = pd.concat([df_out, cdf], ignore_index=True)
        return df_out
    
    # list all the GO Terms associated with a uniprot code
    def GO_Search_Protein(self, uniprot_entry):
        
        GO_list = list()
        for code, uni_list in self.GO_dict.items():
            if uniprot_entry in uni_list:
                GO_list.append(code)
        GO_list = [self.GO_decode_dict_reverse[code] for code in GO_list]
        return GO_list
      
    def GO_Get_Data(self):
        
        uniprot_codes = self.df['Uniprot Entry'].unique()
        for uniprot_code in uniprot_codes:
            try:
                url_2 = 'https://www.uniprot.org/uniprot/' + uniprot_code + '.txt'
                html_2 = urllib.request.urlopen(url_2)
            except Exception as e:
                raise Exception('Failed to obtain UNIPROT data. %s'%e)
            
            try:
                for line in html_2:
                    line = str(line)
                    GO_entry = re.findall('GO;', line)

                    if len(GO_entry) > 0:
                        # GO; GO:0030089; C:phycobilisome; IEA:UniProtKB-KW.
                        # GO; GO:0102834; F:1-18:1-2-16:0-monogalactosyldiacylglycerol acyl-lipid omega-6 desaturase activity; IEA:UniProtKB-EC.
                        GO_code = line.split('; ')[1].split(':')[1]
                        GO_word_idx = line.split('; ')[2].index(':') + 1
                        GO_word = line.split('; ')[2][GO_word_idx:]

                        # put it into self.GO_dict
                        if GO_code not in self.GO_dict.keys():
                            self.GO_dict[GO_code] = [uniprot_code]
                        else:
                            self.GO_dict[GO_code].append(uniprot_code)

                        # put it into self.GO_decode_dict
                        if GO_word not in self.GO_decode_dict.keys():
                            self.GO_decode_dict[GO_word] = GO_code
            
            except Exception as e:
                print('Error %s'%e)
                continue
        
        self.GO_decode_dict_reverse = {v: k for k, v in self.GO_decode_dict.items()}
        
        
    def aggregate(self):

        # get unique uniport entry
        unique_uniport_entries = self.get_unique_uniprot_entry()
        
        index = 0
        for entry in unique_uniport_entries:
            
            # get unique resid
            unique_resid = self.get_unique_resid(entry)
            
            for resid in unique_resid:
                current_df = self.df[(self.df['Uniprot Entry'] == entry) & (self.df['Resid'] == resid)]
                num_of_instances = len(current_df)

                # append a row to new_df
                pka_mean = round(current_df['pKa'].mean(),2)
                pka_std = round(current_df['pKa'].std(),2)
                pka_range = [round(current_df['pKa'].min(),2),round(current_df['pKa'].max(),2)]
                sasa_mean = round(current_df['sasa'].mean(),2)
                sasa_std = round(current_df['sasa'].std(),2)
                sasa_range = [round(current_df['sasa'].min(),2),round(current_df['sasa'].max(),2)]
                self.aggregated_df.loc[index] = [entry, resid, num_of_instances, pka_mean, pka_std, pka_range,
                                                 sasa_mean, sasa_std, sasa_range]
                index += 1
        
        print(f'In total {index} pairs of Uniport Entry and Resid identified')
    
    # plot either a histogram with kernal density estimation or a boxplot
    def plot_graph(self, plot_type, feature, uniprot_entry = False, resid = False):
        
        if not uniprot_entry and not resid:
            x = self.df[feature]
            if plot_type == 'histogram':
                sns.displot(x, kde=True)
        
            elif plot_type == 'boxplot':
                sns.boxplot(x=x)

            else:
                print('Sorry No Such Plot Available!')
            
        else:    
            if uniprot_entry and resid:
                if uniprot_entry not in self.get_unique_uniprot_entry():
                    return 'Wrong Uniprot Entry!'
                else:
                    if resid not in self.get_unique_resid(uniprot_entry):
                        return 'Wrong Resid!'

                current_df = self.df[(self.df['Uniprot Entry'] == uniprot_entry) & (self.df['Resid'] == resid)]
                x = current_df[feature]

                if plot_type == 'histogram':
                    sns.displot(x, kde=True)

                elif plot_type == 'boxplot':
                    sns.boxplot(x=x)

                else:
                    print('Sorry No Such Plot Available!')
            else:
                print('Lack of Information!')
    
    
    def get_outliers(self, uniprot_entry, resid, feature, whis = 1.5):
        current_df = self.df[(self.df['Uniprot Entry'] == uniprot_entry) & (self.df['Resid'] == resid)]
        x = current_df[feature]
        
        Q1 = x.quantile(0.25)
        Q3 = x.quantile(0.75)
        IQR = Q3 - Q1
        
        lower = Q1 - whis * IQR
        upper = Q3 + whis * IQR
        
        outlier_df = current_df[(current_df[feature] < lower) | (current_df[feature] > upper)]
        return outlier_df
    
    def get_extreme_values(self, feature, lower = 1, upper = 14):
        current_df = self.df[(self.df[feature] < lower) | (self.df[feature] > upper)]
        return current_df
    
    
    def remove_df(self, df_to_remove):
        L1 = len(self.df)
        remove_list = df_to_remove.index.tolist()
        self.df = self.df.drop(index = remove_list)
        L2 = len(self.df)
        print(f'Original length: {L1}\nCurrent length: {L2}\nNum of rows removed: {len(df_to_remove)}')
    
    def concise(self, df, weight = 0.5, method = 'average'):
        
        if method != 'average' and method != 'south_east':
            return 'Wrong input method, try average or south_east.'
        
        self.df_concise = pd.DataFrame(columns = ['Uniprot Entry','PDB Code','Method','Resolution','Chain','Resid','pKa','sasa'])
        
        # create the function to quantify the tradeoff between pka and sasa
        def low_pka_large_sasa(df, weight = 0.5):
    
            if len(df) > 1:

                df = df.reset_index(drop = True)

                # situation 1: same pka and same sasa
                if (len(df['pKa'].unique()) == 1) and (len(df['sasa'].unique()) == 1):
                    return df.iloc[0,:]

                # situation 2: same pka but different sasa maybe unlikely
                elif (len(df['pKa'].unique()) == 1) and (len(df['sasa'].unique()) != 1):
                    index = df['sasa'].idxmax()
                    return df.iloc[index,:]

                # situation 3: different pka but same sasa sometimes  
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

        unique_uni_entries = df['Uniprot Entry'].unique()
        for uni_entry in unique_uni_entries:
            current_df = df[df['Uniprot Entry'] == uni_entry]

            unique_resids = current_df['Resid'].unique()
            for resid in unique_resids:
                current_df = df[(df['Uniprot Entry'] == uni_entry) & (df['Resid'] == resid)]
                
                if method == 'south_east':
                    row_to_append = low_pka_large_sasa(current_df, weight = weight)
                    self.df_concise = self.df_concise.append(row_to_append, ignore_index = True)
                
                else:
                    pka_mean = current_df['pKa'].mean()
                    sasa_mean = current_df['sasa'].mean()
                    data = {'Uniprot Entry':uni_entry, 'PDB Code':np.nan, 'Method':np.nan, 'Resolution':np.nan, 'Chain':np.nan,
                           'Resid':resid, 'pKa':pka_mean, 'sasa':sasa_mean}
                    df_dictionary = pd.DataFrame([data])
                    self.df_concise = pd.concat([self.df_concise, df_dictionary], ignore_index=True)
          
        return self.df_concise
