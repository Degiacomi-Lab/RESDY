import os
import re
import urllib.request
import urllib.parse
import urllib.error
import threading
import concurrent.futures
from ast import literal_eval
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

#### TODO Section #### - for general todos in this file, may be more further down
# TODO GW 13.09.24 - most of the GO term analysis currently only works for propka
#                    not pkaani, look into adding this in
# TODO GW 13.09.24 - look into the GO term functions and see if these still actually
#                    work with all the extra stuff added in
# TODO GW 16.04.25 - removed dropna function on init, need to add in function which cleans the dataframe at the start instead. Dont want to blanket remove all null rows incase only null for some measurements and these arent being used


class Analysis(object):

    def __init__(self, df, outdir="result", features_to_analyse = []):
        '''
        Initialise the Analysis class which allows you to create graph and go
        over other metrics such as GO terms

        Parameters
        ----------
        df -> dataframe
            Dataframe of measurements to go over the analysis for
        outdir -> string
            Name of the directory to write to
        features_to_analyse -> list
            List of features which should be analysed over
        '''
        if features_to_analyse == []:
            self.df = df
        else:
            self.df = df.dropna(subset=features_to_analyse)

        self.df_aggregated = pd.DataFrame(columns = ['Uniprot_Entry','Resid','Num'])

        self.df_sub = pd.DataFrame(columns=['Uniprot_Entry'])

        self.GO_dict = {} # code as the key

        self.name_to_code = {}
        self.code_to_name = None

        self.outdir = outdir


    # GW 05.12.24 function potentially unused - remove?
    def get_data(self, uniprot_entry, resid):
        df_query = self.df[(self.df['Uniprot_Entry'] == uniprot_entry) & (self.df['Resid'] == resid)]
        return df_query

    # GW 05.12.24 function potentially unused - remove?
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
            url_2 = f'https://www.uniprot.org/uniprot/{uniprot_code}.txt'
            with urllib.request.urlopen(url_2, timeout=10) as response:
                html_2 = response.read()
        except Exception as e:
            print(f'Failed to obtain UNIPROT data for {uniprot_code}. {e}')

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

                        if GO_word not in self.name_to_code.keys():
                            self.name_to_code[GO_word] = GO_code

        except Exception as e:
            print(f'Error analysis GO data for Uniprot: {uniprot_code} with error: {e}')

    def GO_get_data(self):
        uniprot_codes = self.df['Uniprot_Entry'].unique()
        num = len(uniprot_codes)
        locks = [threading.Lock()]*num
        with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.map(self._GO_get_data, uniprot_codes, locks, range(num), [num]*num)

        self.code_to_name = {v: k for k, v in self.name_to_code.items()}


    def plot_graph(self, plot_type, feature, uniprot_entry = False, resid = False):
        '''
        Basic plot, show either a histogram or a boxplot
        '''
        try:
            plt.clf()
        except Exception:
            pass

        if not uniprot_entry and not resid:
            try:
                x = self.df[feature]
            except Exception:
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
                    if resid not in self.df[self.df['Uniprot_Entry'] == uniprot_entry]['Resid'].unique():
                        print('Wrong Resid.')
                        return

                df_query = self.df[(self.df['Uniprot_Entry'] == uniprot_entry) & (self.df['Resid'] == resid)]

                try:
                    x = df_query[feature]
                except Exception:
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
        '''
        Extract outliers from a dataset using quartiles based on a given feature
        and specific uniprot code and resid number.

        Parameters
        ----------
        uniprot_entry -> string
            Uniprot code of interest
        resid -> string
            Residue number of the uniprot code of interest
        feature -> string
            Feature of interest to extract outliers over
        whis -> float
        '''
        df_query = self.df[(self.df['Uniprot_Entry'] == uniprot_entry) & (self.df['Resid'] == resid)]
        x = df_query[feature]

        q_one = x.quantile(0.25)
        q_three = x.quantile(0.75)
        iqr = q_three - q_one

        lower = q_one - whis * iqr
        upper = q_three + whis * iqr

        df_outlier = df_query[(df_query[feature] < lower) | (df_query[feature] > upper)]
        return df_outlier

    def get_extreme_values(self, feature, lower = 1, upper = 14):
        df_query = self.df[(self.df[feature] < lower) | (self.df[feature] > upper)]
        return df_query

    def remove_df(self, df_to_remove):
        '''
        
        '''
        len_one = len(self.df)
        remove_list = df_to_remove.index.tolist()
        self.df = self.df.drop(index = remove_list)
        len_two = len(self.df)
        print(f'Original num of rows: {len_one}')
        print(f'Current num of rows: {len_two}')
        print(f'Num of rows removed: {len(df_to_remove)}')


    def add_extra_measures(self, extra_measures_filename, write_new_file = False, out_filename='measures_new.csv'):
        '''
        Function to add in extra measurements to the measures frame that has been
        autoloaded into the analysis class on defining this. This will match up the
        measurements in each case and hold in for the analysis. A new measures file
        will be written with the new filename that has been passed into the function.

        Parameters
        ----------
        extra_measures_filename -> string
            The name of the new measures file written of the combination of both
            measures dataframe.
        write_new_file -> bool
            True/False option for writing a new measures.csv file when the new data
            has been added in. Auto set to False. 
        out_filename -> string
            The name of the new measures.csv file that you want to be produced. Auto
            set to be measures_new.csv
        '''

        # Step 1: read in new dataframe, extract column names, check for overlap and
        #         deal if is, otherwise add new column in
        try:
            new_df = pd.read_csv(extra_measures_filename)
            base_columns = ['Uniprot_Entry', 'PDB_Code', 'Method', 'Resolution', 'Chain', 'Resid']
            orig_measures_columns = [a for a in self.df.columns if a not in base_columns]
            new_measures_columns = [a for a in new_df.columns if a not in base_columns]
            overlap_columns = [a for a in new_measures_columns if a in orig_measures_columns]
            new_nonoverlap_columns = [a for a in new_measures_columns if a not in orig_measures_columns]
            for column in overlap_columns:
                proper_answer = False
                print(f'Measurement {column} already present in the loaded dataframe, do you want to replace it?')

                while not proper_answer:
                    col_to_keep = input('Enter old or new for data to keep: ')
                    match col_to_keep.lower():
                        case 'new':
                            print(f'Keeping new measurements for {column}')
                            for i, r in self.df.iterrows():

                                protein_code = r['PDB_Code']
                                chain_value = r["Chain"]
                                resid_value = r["Resid"]

                                idx = np.where((new_df["PDB_Code"] == protein_code) & (new_df["Chain"] == chain_value) & (new_df["Resid"].astype(int) == resid_value))
                                if len(idx[0]) == 0:
                                    continue

                                # account for measurements that have special cases
                                # aevs - add the list of aevs in one column to the overall dataframe
                                if column == 'aev':
                                    self.df['aev'] = self.df['aev'].astype('object')
                                elif column == 'aev_legolas':
                                    self.df['aev_legolas'] = self.df['aev_legolas'].astype('object')
                                else:
                                    self.df.at[i, column] = new_df.loc[idx[0][0], column]
                            proper_answer = True
                        case 'old':
                            print(f'Keeping old measurements for {column}')
                            proper_answer = True
                        case _:
                            print(f'{col_to_keep} was not recognised')

            num_lines_selfdf = len(self.df)

            for column in new_nonoverlap_columns:
                print(f'\nAdding new measurement {column} to the dataframe \n')
                for i, r in self.df.iterrows():

                    protein_code = r['PDB_Code']
                    chain_value = r["Chain"]
                    resid_value = r["Resid"]

                    idx = np.where((new_df["PDB_Code"] == protein_code) & (new_df["Chain"] == chain_value) & (new_df["Resid"].astype(int) == resid_value))
                    if len(idx[0]) == 0:
                        continue

                    # account for measurements that have special cases
                    # aevs - add the list of aevs in one column to the overall dataframe
                    if column == 'aev':
                        new_df['aev'] = new_df['aev'].astype('object')
                    elif column == 'aev_legolas':
                        new_df['aev_legolas'] = new_df['aev_legolas'].astype('object')
                    self.df.at[i, column] = new_df.loc[idx[0][0], column]
                    print(f'Progress adding {column} data: {round(((i+1)/num_lines_selfdf)*100, 2)} %\r', end='', flush=True)

            print(f'\nFinished adding extra measures data to dataframe')

        except Exception as e:
            print(f'Failed to load in the new measures dataframe (name: {extra_measures_filename}), error: {e}')

        # Step 2: write new measures.csv file if required
        if write_new_file:
            self.df.to_csv(out_filename)
            print(f'New measures file written with name: {out_filename}')


    def remove_not_important_residues(self, req_resid_table, outname='measures_cut.csv'):
        '''
        Function to take the input file documenting which residues are required to
        keep due to being of interest and remove anything from the dataframe that
        isnt in this list. This is required due to the codebase calculating data
        for every possible resid in the structure.

        Method
        ------
        Extract the list of residues and taking data for these. Goes over the dataframe
        and extracts any residues which are not present within the required residues.
        Removes these from the dataframe and then writes a new dataframe with the
        updated data.

        Parameters
        ----------
        req_resid_table -> dataframe
            Dataframe containing all the measured data inputted into the analysis class
        
        outname -> string
            The name of the file to give in output for the new updated measures file.
            Auto set to measures_cut.csv
        '''
        print('>> Removing unrequired residues')
        # Remove duplicated data from the measurements
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

        while len(test_table) > 0:
            print("Number of rows left: " + str(len(test_table)))
            test_uniprot = test_table["Uniprot_Entry"][test_table.first_valid_index()]
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
        print(f'Original num of rows: {initial_full_data_rows}')
        print(f'Current num of rows: {final_full_data_rows}')
        print(f'Num of rows removed: {diff_rows}')
        self.df.to_csv(os.path.join(self.outdir, outname), index_label=False, index=False)


    def relative_best(self, df, weights, features=['depth']):
        '''
        Function to extract the best relative list of features for all combinations of uniprot
        entry and residues based on a given list of metrics to do the calculation over and
        the desired weightings for each of those features.

        Parameters
        ----------
        df -> dataframe
            Dataframe of measurements to do the analysis over
        weights -> list
            List of floats which sum to 1 of the weights for each of the given features. Length
            should match the list of features given
        features -> list
            List of features that the analyis should extract the relative best values for. 
        '''
        # go over the weights to make sure the values are good and then matches the number of features
        try:
            if isinstance(weights, float) or isinstance(weights, int):
                if weights == 0.5:
                    print('>> Weight remains as default and equal for all metrics')
                    if len(features) == 1: weights = [1]
                    else: weights = [1/len(features)] * len(features)
            elif isinstance(weights, list):
                if len(set(weights)) == 1:
                    print('>> Weight remains as default and equal for all metrics')
                    weights = [1/len(features)] * len(features)
                elif len(weights) < len(features):
                    print('>> Number of weights given is lower than the number of features. Adjusting weights to include')
                    while len(weights) != len(features): weights.append(1/len(features))
                    weights = [round((i/sum(weights)), 2) for i in weights]
                elif len(weights) > len(features):
                    print(f'>> Number of weights given is greater than the number of features. Taking first {len(features)} through')
                    weights = weights[:len(features)]
                    weights = [round((i/sum(weights)), 2) for i in weights]
                elif sum(weights) != 1:
                    print('>> Weights do not sum to 1. Adjusting weights as required.')
                    weights = [round((i/sum(weights)), 2) for i in weights]
                else: print('>> Weights fit, taking through as given')
        except Exception as e:
                print(f'Error sorting the weights for the analysis: {e}')

        df = df.dropna(subset=features)
        self.df_sub = pd.DataFrame(columns=['Uniprot_Entry'])

        # if features not known which trend is best, ask user which is required
        min_feats = ['propka', 'pkaani', 'legolas', 'depth']
        max_feats = ['sasa', 'das', 'seqcharge', 'frustration']
        for feat in features:
            if feat in ['aev', 'aev_legolas']:
                print(f'>> Feature {feat} is not supported with this analysis. Dropping feature from ')
            if feat not in min_feats + max_feats:
                pref = 'tmp'
                while pref not in ['min', 'max']:
                    pref = input('Bias towards min or max for curvature (enter "min" or "max"): ')
                if pref == 'min': min_feats.append(feat)
                elif pref == 'max': max_feats.append(feat)

        df_temp = df.drop_duplicates(subset=['Uniprot_Entry', 'Resid'])
        for i, r in df_temp.iterrows():
            uniprot_tmp, resid_tmp = r['Uniprot_Entry'], r['Resid']
            df_query = df[(df['Uniprot_Entry'] == uniprot_tmp) & (df['Resid'] == resid_tmp)]

            if len(df_query) > 1:
                df_query = df_query.reset_index(drop=True)

                try:
                    features_and_weights = zip(features, weights)
                except Exception as e:
                    print(f'Failed to link metrics with weights, were the metrics entered correctly? {e}')
                
                metric_calculated_values_list_temp = []
                try:
                    for feature, weight in features_and_weights:
                        try:
                            match feature:
                                case x if x in min_feats:
                                    # features need min value, invert values in list
                                    feat_max = df_query[feature].max()
                                    feat_list = [(feat_max-i) for i in list(df_query[feature])]
                                case x if x in max_feats:
                                    # features need max value
                                    feat_list = list(df_query[feature])
                                case _:
                                    print(f'>> Feature ({feature}) unknown, taking minimum')
                                    feat_list = list(df_query[feature])
                            feat_avg, feat_std = np.mean(feat_list), np.std(feat_list)
                            feat_list_standardised = (feat_list - feat_avg) / feat_std
                            feat_list_weighted = feat_list_standardised * weight
                            metric_calculated_values_list_temp.append(feat_list_weighted)
                        except Exception as e:
                            print(f'>> Failed to load the data for the metric: {feature} with error: {e}')
                except Exception as e:
                    print(f'Error loading data for the metrics provided: {e}')
                
                try:
                    index = -1
                    base = 0
                    for i, val in enumerate(metric_calculated_values_list_temp[0]):
                        weighted_sum = 0
                        for weighted_value_index, metric_value in enumerate(metric_calculated_values_list_temp):
                            weighted_sum += metric_value[i]
                        if weighted_sum >= base:
                            index = i
                            base = weighted_sum
                    row_to_append = df_query.iloc[[index],:]
                except Exception as e:
                    print(f'Failed to find the best row for the desired trade off between the three metrics: {e}')
                    row_to_append = df_query.iloc[[0],:]
            else:
                row_to_append = df_query.iloc[[0],:]

            self.df_sub = pd.concat([self.df_sub, row_to_append], axis=0, ignore_index=True)
            self.df_sub['Resid'] = self.df_sub['Resid'].astype(int)


    def get_contingency_table(self, GO_code, my_list, reference):
        '''
        Compute contingency table given a GO Term, the list of interest, and a reference list
        '''

        bp_list = 0
        for uni in my_list:
            if uni in self.GO_dict[GO_code]:
                bp_list += 1

        bp_not_list = len([p for p in reference if GO_code in self.GO_dict[p]]) - bp_list
        not_bp_list = len(my_list) - bp_list

        not_bp_not_list = 0
        for uni in reference:
            if (uni not in my_list) and (uni not in self.GO_dict[GO_code]):
                not_bp_not_list += 1

        table = [[bp_list, bp_not_list], [not_bp_list, not_bp_not_list]]

        return table


    def enrichment_analysis(self, feature_one = ['propka', 7, 11], feature_two = ['sasa', 0, 10], uniprot_cnt_cutoff=1):
        '''
        Analyse prevalence of GO-terms in sub-regions of the SASA vs pKa graph
        '''

        # get all the uniprot codes inside the range and the reference uniprot code list
        feat_one, feat_one_low, feat_one_upper = str(feature_one[0]), float(feature_one[1]), float(feature_one[2])
        feat_two, feat_two_low, feat_two_upper = str(feature_two[0]), float(feature_two[1]), float(feature_two[2])
        selected_df = self.df_sub[(self.df_sub[feat_one] >= feat_one_low) & (self.df_sub[feat_one] <= feat_one_upper)]
        selected_df = selected_df[(selected_df[feat_two] >= feat_two_low) & (selected_df[feat_two] <= feat_two_upper)]
        uniprot_selected = selected_df['Uniprot_Entry'].unique()
        uniprot_reference = self.df_sub['Uniprot_Entry'].unique()

        # get all the GO Terms in the background (that are associated with more than uniprot_cnt_cutoff)
        GO_bacgou = [code for code in list(self.GO_dict.keys()) if len(self.GO_dict[code]) >= uniprot_cnt_cutoff]

        p_val_dict = {}

        # for each GO Term, compute a contingency table
        for GO_code in GO_bacgou:
            table = self.get_contingency_table(GO_code, uniprot_selected, uniprot_reference)
            # compute p values and store them into the dictionary
            oddsratio, pvalue = fisher_exact(table, alternative='greater')
            p_val_dict[GO_code] = pvalue

        # sort the dictionary based on p-values
        p_val_dict = dict(sorted(p_val_dict.items(), key = lambda x: x[1]))
        p_val_list = [x[1] for x in p_val_dict.items()]

        # adjust p-values using the BH method
        y=multipletests(pvals=p_val_list, alpha=0.05, method="fdr_bh")
        out_dict = {}

        # no matter enrichment, output p values to the dictionary
        for i in range(len(y[0])):
            go_code = list(p_val_dict.items())[i][0]
            raw_p = p_val_list[i]
            adj_p = y[1][i]
            out_dict[go_code] = (raw_p, adj_p, self.get_contingency_table(go_code, uniprot_selected, uniprot_reference))

        df = pd.DataFrame(columns = ['GO ID', 'GO Term', 'raw p value', 'FDR',
                                    'num in the region', 'num in the bkgd'])
        for data in out_dict.items():
            id = data[0]
            term = self.code_to_name[id]

            r_p = round(data[1][0],6)
            fdr = round(data[1][1],6)

            n1 = sum(data[1][2][0])
            d1 = sum(data[1][2][0]) + sum(data[1][2][1])
            bkgd = str(n1) + '/' + str(d1)
            n2 = data[1][2][0][0]
            d2 = data[1][2][1][0] + data[1][2][0][0]
            reg = str(n2) + '/' + str(d2)

            dt = {'GO ID':id, 'GO Term':term, 'raw p value':r_p, 'FDR':fdr,
                'num in the bkgd':bkgd, 'num in the region':reg}
            df_dictionary = pd.DataFrame([dt])
            df = pd.concat([df, df_dictionary], ignore_index=True)

        return df


if __name__ == '__main__':
    df_test = pd.read_csv(f'data{os.sep}measures_CannData_all_12.05.25.csv')
    analysis = Analysis(df_test, outdir='result')
    analysis.relative_best(analysis.df, weights=0.5, features=['depth', 'sasa', 'phi', 'das'])
    print(analysis.df_sub)
