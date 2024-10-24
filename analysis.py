from cgi import test
import os
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

#### TODO Section #### - for general todos in this file, may be more further down
# TODO GW 13.09.24 - most of the GO term analysis currently only works for propka not pkaani, look into adding this in
# TODO GW 13.09.24 - look into the GO term functions and see if these stilll actually work with all the extra stuff added in





class Analysis(object):
    
    def __init__(self, df,  outdir="result"):
        self.df = df.dropna(subset=['propka', 'sasa'])
        
        self.df_aggregated = pd.DataFrame(columns = ['Uniprot_Entry','Resid','Num','propka mean','propka std','propka range', 
        'SASA mean','SASA std','SASA range', 'Depth mean', 'Depth std', 'Depth range'])
        
        #self.df_sub = pd.DataFrame(columns = ['Uniprot_Entry','Resid','propka','sasa', 'depth']) 
        self.df_sub = pd.DataFrame(columns=['Uniprot_Entry'])
        
        self.GO_dict = {} # code as the key
        
        self.name_to_code = {}
        self.code_to_name = None

        self.outdir = outdir
    
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
        #df_temp is just used to get the table of uniprot codes and associated resids for repeating over, the aggregation code uses the full dataset
        df_temp = self.df.drop_duplicates(subset=['Uniprot_Entry','Resid'])
        for idx, row in df_temp.iterrows():
            # extract the uniprot and resid of interest
            entry = row['Uniprot_Entry']
            resid = row['Resid']
            # create a subset of the dataframe of measurements where the uniprot and resid match
            df_query = self.df[(self.df['Uniprot_Entry'] == entry) & (self.df['Resid'] == resid)]
            # calculate all the values of interest from the subset dataframe (df_query)
            num = len(df_query)
            propka_mean = round(df_query['propka'].mean(),2)
            propka_std = round(df_query['propka'].std(),2)
            propka_range_values = [round(df_query['propka'].min(),2),round(df_query['propka'].max(),2)]
            propka_range = propka_range_values[1] - propka_range_values[0]
            # TODO GW 13.09.24 - Do we need to split up analysis into Propka and pkaANI, at the moment it just takes forward the propka data
            sasa_mean = round(df_query['sasa'].mean(),2)
            sasa_std = round(df_query['sasa'].std(),2)
            sasa_range_values = [round(df_query['sasa'].min(),2),round(df_query['sasa'].max(),2)]
            sasa_range = sasa_range_values[1] - sasa_range_values[0]
            depth_mean = round(df_query['depth'].mean(), 2)
            depth_std = round(df_query['depth'].std(), 2)
            depth_range_values = [round(df_query['depth'].min(), 2), round(df_query['depth'].max(), 2)]
            depth_range = depth_range_values[1] - depth_range_values[0]
           
            data = {'Uniprot_Entry': entry,
                        'Resid': resid,
                        'Num': num,
                        'propka mean': propka_mean,
                        'propka std': propka_std,
                        'propka range': propka_range,
                        'SASA mean': sasa_mean,
                        'SASA std': sasa_std,
                        'SASA range': sasa_range,
                        'Depth mean': depth_mean,
                        'Depth std': depth_std,
                        'Depth range': depth_range}
            
            # add in the data to the aggregated dataframe, index provided due to only scalar values being used before being ignored when it is added in. 
            self.df_aggregated = pd.concat([self.df_aggregated, pd.DataFrame(data, index=[0])], ignore_index=True)
                

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

    # function added by GW 09.11.23 to remove all measures that were done on residues that aren't in a set of data
    def remove_not_important_residues(self, req_resid_table):
        print('>> Removing unrequired residues')
        # duplicate the req_resid_table to allow to delete rows with testing
        test_table = req_resid_table
        initial_data_one = len(self.df)
        self.df = self.df.drop_duplicates()
        duplicate_rows_removed = initial_data_one - len(self.df)
        print(f'Removed {duplicate_rows_removed} rows of duplicates')
        initial_full_data_rows = len(self.df)
        # remove rows which have a UNIPROT code which isnt required
        uniprot_codes = test_table['Uniprot_Entry'].drop_duplicates().tolist()
        entries_to_remove = []
        for i, r in self.df.iterrows():
            if r['Uniprot_Entry'] not in uniprot_codes:
                entries_to_remove.append(i)
        self.df = self.df.drop(index=entries_to_remove)
        uniprot_rows_removed = initial_full_data_rows - len(self.df)
        print(f'Removed {uniprot_rows_removed} rows of Uniprot codes which were not mentioned in the required residues file')
        # iterate over each set of residues of a protein
        while len(test_table) > 0:
            # prints the number of rows left in the hits sheet updating how far through you are
            print("Number of rows left: " + str(len(test_table)))
            # read in the uniprot code at the top of the hits sheet
            test_uniprot = test_table["Uniprot_Entry"][test_table.first_valid_index()]
            # print out which one you are finding, mainly just for checking
            print("test_uniprot: " + str(test_uniprot))
            # find all the desired residues from the particular uniprot code and put into a list
            # automatically removes duplicates from this (doesn't retain order)
            desired_residues = list(set(test_table[test_table["Uniprot_Entry"].str.contains(test_uniprot)]["Resid"].tolist()))
            # search the measures spreadsheet for all rows containing the desired uniprot code
            search_uniprot = self.df[self.df["Uniprot_Entry"].str.contains(test_uniprot.strip())][["Uniprot_Entry", "Resid"]]
            all_search_rows = search_uniprot.index.tolist()
            # go over each row of search_uniprot, see if the residue matches one of the desired ones
            wanted_rows = search_uniprot[search_uniprot["Resid"].isin(desired_residues)].index.tolist()
            not_wanted_rows = [x for x in all_search_rows if x not in wanted_rows]
            print("Rows removed: " + str(len(not_wanted_rows)))
            # remove the rows which aren't wanted from the main data set
            self.df = self.df.drop(index = not_wanted_rows)
            # remove rows which contain the uniprot code that has been searched from test_table
            test_table = test_table.drop(index = test_table[test_table["Uniprot_Entry"] == test_uniprot].index.tolist())

        final_full_data_rows = len(self.df)
        diff_rows = initial_full_data_rows - final_full_data_rows
        print(f'Original num of rows: {initial_full_data_rows}\nCurrent num of rows: {final_full_data_rows}\nNum of rows removed: {diff_rows}')
        self.df.to_csv(os.path.join(self.outdir, "measures_cut.csv"), index_label=False, index=False) 

    def subset(self, df, weight = 0.5, method = 'average', metrics=['propka', 'sasa', 'depth']):
        
        if method != 'average' and method != 'south_east':
            raise ValueError('Wrong input method, try average or south_east.')
        
        # catch errors with bad input weights and transfer into list format if not there already
        try:
            if type(weight) == float or type(weight) == int:
                if weight == 0.5:
                    print('Weight remains as default and equal for all metrics')
                #transfer weight to list to make sure 
                weight = [weight]
            elif type(weight) == list:
                if weight[0] == 0.5:
                    print('Weight remains as default and equal for all metrics')
       
        except:
            print('Failed reading weight')

                
        # reset the dataframe incase it is rerun with the other option, stops the dataframe getting bigger and bigger
        self.df_sub = pd.DataFrame(columns=['Uniprot_Entry'])


        # function to evaluate trade-off between pka and sasa - OLD function
        '''
        def low_pka_large_sasa(df, weight = 0.5):
            if len(df) > 1:
                df = df.reset_index(drop = True)

                # situation 1: both pKa and sasa only have 1 unique value each -> take the first row as all the same
                if (len(df['propka'].unique()) == 1) and (len(df['sasa'].unique()) == 1):
                    return df.iloc[[0],:]

                # situation 2: only 1 unique pKa value but more than 1 unique sasa value
                elif (len(df['propka'].unique()) == 1) and (len(df['sasa'].unique()) != 1):
                    # as all pKa the same, just find the max value for sasa and return the row which has this
                    index = df['sasa'].idxmax()
                    return df.iloc[[index],:]

                # situation 3: more than 1 unique pKa value, only 1 unique sasa value
                elif (len(df['propka'].unique()) != 1) and (len(df['sasa'].unique()) == 1):
                    # all sasa values the same, so just find the lowest pKa value and return the row that this is on
                    index = df['propka'].idxmin()
                    return df.iloc[[index],:]

                # situation 4: more than 1 unique value for both pKa and sasa
                else:
                    pka_max = df['propka'].max()
                    pka_list = [(pka_max-i) for i in df['propka'].tolist()]  # this effectively inverts the values, eg a lower pKa now had a higher value, allows normalised comparison later
                    sasa_list = df['sasa'].tolist()

                    # compute mean and std
                    pka_m, pka_std = np.mean(pka_list), np.std(pka_list)
                    sasa_m, sasa_std = np.mean(sasa_list), np.std(sasa_list)

                    # standardize two lists
                    pka_list = (pka_list - pka_m) / pka_std
                    #print(pka_list)
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
                    return df.iloc[[index],:]
            else:
                return df.iloc[[0],:]
        '''
        
        # general function for investigating the trade off between 2 metrics
        def relative_best_2D(df, weight_list, metrics):
            # sort the weights into a list of length 2
            try:
                if len(weight_list) == 1:
                    # only 1 weight provided for comparison - assume equal weighting for each metric
                    new_weights = [0.5, 0.5]
                elif len(weight_list) == 2:
                    # 2 values provided for weights, check they sum to 1 otherwise normalise
                    if weight_list[0] + weight_list[1] != 1:
                        new_weights = [round((i/sum(weight_list)), 2) for i in weight_list]
                elif len(weight_list) > 2:
                    print('More weights provided than needed, taking the first two through.')
                    new_weights = weight_list[0,1]
            except:
                print('Error sorting the weights for the analysis')

            # if data for the lysine is only 1 row, don't need to run analysis on it
            if len(df) > 1:
                df = df.reset_index(drop = True)

                # situation 1: both metric1 and metric2 only have 1 unique value each -> take the first row as all the same
                if (len(df[metrics[0]].unique()) == 1) and (len(df[metrics[1]].unique()) == 1):
                    return df.iloc[[0],:]

                # situation 2: only 1 unique metric1 value but more than 1 unique metric2 value
                elif (len(df[metrics[0]].unique()) == 1) and (len(df[metrics[1]].unique()) != 1):
                    # as all metric1 values the same, just find the optimal for metric2 and return the row which has this
                    match metrics[1]:
                        case 'propka':
                            index = df[metrics[1]].idxmin()
                        case 'sasa':
                            index = df[metrics[1]].idxmax()
                        case 'depth':
                            index = df[metrics[1]].idxmin()
                    return df.iloc[[index],:]

                # situation 3: more than 1 unique metric1 value, only 1 unique metric1 value
                elif (len(df[metrics[0]].unique()) != 1) and (len(df[metrics[1]].unique()) == 1):
                    # all metric2 values the same, so just find the lowest pKa value and return the row that this is on
                    match metrics[0]:
                        case 'propka':
                            index = df[metrics[0]].idxmin()
                        case 'sasa':
                            index = df[metrics[0]].idxmax()
                        case 'depth':
                            index = df[metrics[0]].idxmin()
                    return df.iloc[[index],:]

                # situation 4: more than 1 unique value for both pKa and sasa
                else:
                    # link together the metrics with the weights
                    try:
                        metric_and_weights = zip(metrics, new_weights)
                    except:
                        print('Failed to link metrics with weights, were the metrics entered correctly?')
                    
                    # setup temporary list to house the metrics list after they have been calculated
                    metric_calculated_values_list_temp = []

                    try:
                        for metric, met_weight in metric_and_weights:
                            match metric:
                                case 'propka':
                                    try:
                                        # preference for lower pKa -> invert list
                                        pka_max = df['propka'].max()
                                        pka_list = [(pka_max-i) for i in df['propka'].tolist()]
                                        # compute mean and std
                                        pka_avg, pka_std = np.mean(pka_list), np.std(pka_list)
                                        # standardise list
                                        pka_list_standardised = (pka_list - pka_avg) / pka_std
                                        # weight the list
                                        pka_list_weighted = pka_list_standardised * met_weight
                                        # append the list to the temporary list
                                        metric_calculated_values_list_temp.append(pka_list_weighted)
                                    except:
                                        print('Failed to load the data for the metric: propka')
                                case 'sasa':
                                    try:
                                        # preference for highest sasa -> just take list
                                        sasa_list = df['sasa'].tolist()
                                        # compute mean and std
                                        sasa_avg, sasa_std = np.mean(sasa_list), np.std(sasa_list)
                                        # standardise list
                                        sasa_list_standardised = (sasa_list - sasa_avg) / sasa_std
                                        # weight the list
                                        sasa_list_weighted = sasa_list_standardised * met_weight
                                        # append the list to the temporary list
                                        metric_calculated_values_list_temp.append(sasa_list_weighted)
                                    except:
                                        print('Failed to load the data for the metric: sasa')
                                case 'depth':
                                    try:
                                        # preference for lower depth -> invert list
                                        depth_max = df['depth'].max()
                                        depth_list = [(depth_max-i) for i in df['depth'].tolist()]
                                        # compute mean and std
                                        depth_avg, depth_std = np.mean(depth_list), np.std(depth_list)
                                        # standardise list
                                        depth_list_standardised = (depth_list - depth_avg) / depth_std
                                        # weight the list
                                        depth_list_weighted = depth_list_standardised * met_weight
                                        # append the list to the temporary list
                                        metric_calculated_values_list_temp.append(depth_list_weighted)
                                    except:
                                        print('Failed to load the data for the metric: depth')
                                case _:
                                    print('Make sure metrics entered are correct: accepted metrics are currently propka, sasa and depth')
                    except:
                        print('Error loading data for the metrics provided')

                    try:
                        index = -1
                        base = 0
                        for i in range(len(metric_calculated_values_list_temp[0])):
                            weighted_sum = metric_calculated_values_list_temp[0][i] + metric_calculated_values_list_temp[1][i]
                            if weighted_sum >= base:
                                index = i
                                base = weighted_sum
                        return df.iloc[[index],:]
                    except:
                        print('Failed to find the best row for the desired trade off between the two metrics')
                        return df.iloc[[0],:]
            else:
                return df.iloc[[0],:]

        # This function gives the best row of values for that lysine from the dataframe based on relative trade off between metrics
        def relative_best_3D(df, weight_list, metrics):
            uniprot_temp = str(df['Uniprot_Entry'].values[0])
            lysine_temp = str(df['Resid'].values[0])

            # TODO GW-16.09.24 - as the 2D function is used within the 3D function, this affects the error messages for the weightings, 

            # sort the weights into a list of length 2
            try:
                if len(weight_list) == 1:
                    # only 1 weight provided for comparison -> assume equal weighting for each metric
                    new_weights = [0.33, 0.33, 0.33]
                elif len(weight_list) == 2:
                    # only 2 weights provided for comparison -> assume equal weighting for each metric
                    print('Only 2 weightings given for 3 metrics, assuming equal weighting for all')
                    new_weights = [0.33, 0.33, 0.33]
                elif len(weight_list) == 3:
                    # 3 values provided for weights -> check they sum to 1 otherwise normalise
                    if sum(weight_list) != 1:
                        new_weights = [round((i/sum(weight_list)), 2) for i in weight_list]
                elif len(weight_list) > 3:
                    print('More weights provided than needed, taking the first three through.')
                    new_weights = weight_list[0,1,2]
            except:
                print('Error sorting the weights for the analysis')

            # if data for the lysine is only 1 row, don't need to run analysis on it
            if len(df) > 1:
                df = df.reset_index(drop = True)

                # situation 1: all 3 metrics only have 1 unique value each -> take the first row as all the same
                if (len(df[metrics[0]].unique()) == 1) and (len(df[metrics[1]].unique()) == 1) and (len(df[metrics[2]].unique()) == 1):
                    try:
                        return df.iloc[[0],:]
                    except:
                        print(f'Failed 3 metric analysis for situation 1 on protein: {uniprot_temp}, Lysine: {lysine_temp}')

                # situation 2.1: only 1 unique value for metric1 and metric2 but more than 1 unique metric3 value
                elif (len(df[metrics[0]].unique()) == 1) and (len(df[metrics[1]].unique()) == 1) and (len(df[metrics[2]].unique()) != 1):
                    try:
                        # as all metric1 values the same, just find the max value for metric2 and return the row which has this
                        match metrics[2]:
                            case 'propka':
                                index = df[metrics[2]].idxmin()
                            case 'sasa':
                                index = df[metrics[2]].idxmax()
                            case 'depth':
                                index = df[metrics[2]].idxmin()
                        return df.iloc[[index],:]
                    except:
                        print(f'Failed 3 metric analysis for situation 2.1 on protein: {uniprot_temp}, Lysine: {lysine_temp}')

                # situation 2.2: only 1 unique value for metric1 and metric3 but more than 1 unique metric2 value
                elif (len(df[metrics[0]].unique()) == 1) and (len(df[metrics[1]].unique()) != 1) and (len(df[metrics[2]].unique()) == 1):
                    try:
                        # all sasa values the same, so just find the lowest pKa value and return the row that this is on
                        match metrics[1]:
                            case 'propka':
                                index = df[metrics[1]].idxmin()
                            case 'sasa':
                                index = df[metrics[1]].idxmax()
                            case 'depth':
                                index = df[metrics[1]].idxmin()
                        return df.iloc[[index],:]
                    except:
                        print(f'Failed 3 metric analysis for situation 2.2 on protein: {uniprot_temp}, Lysine: {lysine_temp}')

                # situation 2.3: only 1 unique value for metric2 and metric3 but more than 1 unique metric1 value
                elif (len(df[metrics[0]].unique()) != 1) and (len(df[metrics[1]].unique()) == 1) and (len(df[metrics[2]].unique()) == 1):
                    try:
                        # all sasa values the same, so just find the lowest pKa value and return the row that this is on
                        match metrics[0]:
                            case 'propka':
                                index = df[metrics[0]].idxmin()
                            case 'sasa':
                                index = df[metrics[0]].idxmax()
                            case 'depth':
                                index = df[metrics[0]].idxmin()
                        return df.iloc[[index],:]
                    except:
                        print(f'Failed 3 metric analysis for situation 2.3 on protein: {uniprot_temp}, Lysine: {lysine_temp}')

                # situation 3.1: 1 metric has only unique values, 2 metrics have different values -> go to the function for 2D analysis
                elif (len(df[metrics[0]].unique()) != 1) and (len(df[metrics[1]].unique()) != 1) and (len(df[metrics[2]].unique()) == 1):
                    try:
                        new_weight_list = [metrics[0], metrics[1]]
                        new_metric_list = [new_weights[0], new_weights[1]]
                        return relative_best_2D(df, metrics=new_metric_list, weight_list=new_weight_list)
                    except:
                        print(f'Failed 3 metric analysis for situation 3.1 on protein: {uniprot_temp}, Lysine: {lysine_temp}')
                
                # situation 3.2: 1 metric has only unique values, 2 metrics have different values -> go to the function for 2D analysis
                elif (len(df[metrics[0]].unique()) != 1) and (len(df[metrics[1]].unique()) == 1) and (len(df[metrics[2]].unique()) != 1):
                    try:
                        new_weight_list = [metrics[0], metrics[2]]
                        new_metric_list = [new_weights[0], new_weights[2]]
                        return relative_best_2D(df, metrics=new_metric_list, weight_list=new_weight_list)
                    except:
                        print(f'Failed 3 metric analysis for situation 3.2 on protein: {uniprot_temp}, Lysine: {lysine_temp}')
                
                # situation 3.3: 1 metric has only unique values, 2 metrics have different values -> go to the function for 2D analysis
                elif (len(df[metrics[0]].unique()) == 1) and (len(df[metrics[1]].unique()) != 1) and (len(df[metrics[2]].unique()) != 1):
                    try:
                        new_weight_list = [metrics[1], metrics[2]]
                        new_metric_list = [new_weights[1], new_weights[2]]
                        return relative_best_2D(df, metrics=new_metric_list, weight_list=new_weight_list)
                    except:
                        print(f'Failed 3 metric analysis for situation 3.3 on protein: {uniprot_temp}, Lysine: {lysine_temp}')
                
                # situation 4: more than 1 unique value for both all 3 metrics
                else:
                    # link together the metrics with the weights
                    try:
                        metric_and_weights = zip(metrics, new_weights)
                    except:
                        print('Failed to link metrics with weights, were the metrics entered correctly?')
                    
                    # setup temporary list to house the metrics list after they have been calculated
                    metric_calculated_values_list_temp = []
                    try:
                        for metric, met_weight in metric_and_weights:
                            match metric:
                                case 'propka':
                                    try:
                                        # preference for lower pKa -> invert list
                                        pka_max = df['propka'].max()
                                        pka_list = [(pka_max-i) for i in df['propka'].tolist()]
                                        # compute mean and std
                                        pka_avg, pka_std = np.mean(pka_list), np.std(pka_list)
                                        # standardise list
                                        pka_list_standardised = (pka_list - pka_avg) / pka_std
                                        # weight the list
                                        pka_list_weighted = pka_list_standardised * met_weight
                                        # append the list to the temporary list
                                        metric_calculated_values_list_temp.append(pka_list_weighted)
                                    except:
                                        print('Failed to load the data for the metric: propka')
                                case 'sasa':
                                    try:
                                        # preference for highest sasa -> just take list
                                        sasa_list = df['sasa'].tolist()
                                        # compute mean and std
                                        sasa_avg, sasa_std = np.mean(sasa_list), np.std(sasa_list)
                                        # standardise list
                                        sasa_list_standardised = (sasa_list - sasa_avg) / sasa_std
                                        # weight the list
                                        sasa_list_weighted = sasa_list_standardised * met_weight
                                        # append the list to the temporary list
                                        metric_calculated_values_list_temp.append(sasa_list_weighted)
                                    except:
                                        print('Failed to load the data for the metric: sasa')
                                case 'depth':
                                    try:
                                        # preference for lower depth -> invert list
                                        depth_max = df['depth'].max()
                                        depth_list = [(depth_max-i) for i in df['depth'].tolist()]
                                        # compute mean and std
                                        depth_avg, depth_std = np.mean(depth_list), np.std(depth_list)
                                        # standardise list
                                        depth_list_standardised = (depth_list - depth_avg) / depth_std
                                        # weight the list
                                        depth_list_weighted = depth_list_standardised * met_weight
                                        # append the list to the temporary list
                                        metric_calculated_values_list_temp.append(depth_list_weighted)
                                    except:
                                        print('Failed to load the data for the metric: depth')
                                case _:
                                    print('Make sure metrics entered are correct: accepted metrics are currently propka, sasa and depth')
                    except:
                        print('Error loading data for the metrics provided')

                    try:
                        index = -1
                        base = 0
                        for i in range(len(metric_calculated_values_list_temp[0])):
                            weighted_sum = 0
                            for weighted_value_index in range(len(metric_calculated_values_list_temp)):
                                weighted_sum += metric_calculated_values_list_temp[weighted_value_index][i]
                            if weighted_sum >= base:
                                index = i
                                base = weighted_sum
                        return df.iloc[[index],:]
                    except:
                        print('Failed to find the best row for the desired trade off between the three metrics')
                        return df.iloc[[0],:]
            else:
                return df.iloc[[0],:]

        # TODO GW-16.09.24 -if more metrics are added, the above 2 functions can be combined into 1 where user can request any number of metrics to be evaluated against each other
        #                   code is almost there for this, just needs sorting of the weights incase differnet number to the number of metrics


        df_temp = df.drop_duplicates(subset=['Uniprot_Entry','Resid'])
        for idx, row in df_temp.iterrows():
            entry = row['Uniprot_Entry']
            resid = row['Resid']
            df_query = df[(df['Uniprot_Entry'] == entry) & (df['Resid'] == resid)]
            
            if method == 'south_east':
                if len(metrics) == 3:
                    row_to_append = relative_best_3D(df_query, weight_list = weight, metrics=metrics)
                elif len(metrics) == 2:
                    row_to_append = relative_best_2D(df_query, weight_list = weight, metrics=metrics)
                else:
                    print('Please enter at least 2 metrics to compare the trade off between')
                # TODO GW-16.09.24 - is it worth adding a function here where a user could take through a row which has the best value for just 1 metric
                #self.df_sub = self.df_sub.append(row_to_append, ignore_index = True)  # old line which doesnt work anymore, replaced by new one
                self.df_sub = pd.concat([self.df_sub, row_to_append], axis=0, ignore_index=True)
                self.df_sub['Resid'] = self.df_sub['Resid'].astype(int)
            else:
                propka_mean = df_query['propka'].mean()
                sasa_mean = df_query['sasa'].mean()
                depth_mean = df_query['depth'].mean()
                # GW: Have changed the data entry from the following line to the one after to; not worth including the NaN values in this dataframe when they dont add anything to it
                #data = {'Uniprot_Entry':entry, 'PDB_Code':np.nan, 'Method':np.nan, 'Resolution':np.nan, 'Chain':np.nan, 'Resid':resid, 'propka':propka_mean, 'sasa':sasa_mean, 'depth':depth_mean}
                data = {'Uniprot_Entry':entry, 'Resid':resid, 'propka mean':propka_mean, 'sasa mean':sasa_mean, 'depth mean':depth_mean}
                self.df_sub = pd.concat([self.df_sub, pd.DataFrame(data, index=[0])], ignore_index=True)



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
        selected_df = self.df_sub[(self.df_sub['propka'] >= pka_l) & (self.df_sub['propka'] <= pka_u)]
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
