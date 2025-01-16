import os
import random
import statistics
from ast import literal_eval
from copy import deepcopy as dc
#data handling
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn
from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score

cwd = os.getcwd()
pd.set_option('display.max_rows', 500)
pd.set_option('display.width', 1000)
#pd.set_option('expand_frame_repr', False)
pd.set_option('display.min_rows', 50)

class Model(object):
    '''
    Class for the random forest model which we are producing on the dataset to see how good the data is
    '''

    def __init__(self, pos_measures_files,
                 neg_measures_files,
                 features_to_include=['aev'],
                 aggregation_method='avg',
                 subtract_avg_aev=True,
                 num_aev_features_req = 100):
        
        self.pos_measures_files = pos_measures_files
        self.neg_measures_files = neg_measures_files
        self.X_all = pd.DataFrame()
        self.X_agg = pd.DataFrame()
        self.X_final = pd.DataFrame()
        self.y = []

        self.features_to_include = features_to_include
        self.aggregation_method = aggregation_method
        self.subtract_avg_aev = subtract_avg_aev
        self.top_n_features = 0
        self.num_aev_features_req = num_aev_features_req


    def _prepare_aevs(self):
        '''
        Function to prepare aevs correctly from the overall set for the aevs
        Takes the column of AEVs, splits the AEV on commas and changes it into seperate columns
        Removes the brackets left over from splitting the string of the list

        Parameters
        ----------
        df : dataframe
            measures dataframe read in from file

        Returns
        -------
        df : dataframe
            a dataframe containing the split aevs 

        Example
        -------
        >> model.prepare_aevs(df)
        '''
        for i, r in self.X_all.iterrows():
            new_aev = literal_eval(r['aev'])
            self.X_all.at[i, 'aev'] = new_aev

        if self.subtract_avg_aev:
            with open('aev_data_randomforestmodel/average_aevs_overall_noH.csv', 'r') as avg_aev_file:
                average_lys_aev = avg_aev_file.readline().split(',')
            #average_lys_aev = list(pd.read_csv('aev_data_randomforestmodel/average_aevs_overall_noH.csv', header=None))
            average_lys_aev = [float(a) for a in average_lys_aev]
            print(average_lys_aev)
            for i, r in self.X_all.iterrows():
                old_aev = r['aev']
                new_aev = [(a - b) for a, b in zip(old_aev, average_lys_aev)]
                self.X_all.at[i, 'aev'] = new_aev



    def prepare_dataset(self):
        '''
        Function to prepare the dataset from a list of measures.csv files
        For pos and neg list of measures.csv files

        Returns
        -------
        X : dataframe
            datafame containing the prepared dataset ready for using within the RF model

        Example
        -------
        >> model.prepare(df)
        '''
        # read in data
        pos_data = pd.DataFrame()
        neg_data = pd.DataFrame()
        for file in self.pos_measures_files:
            pos_data = pd.concat([pos_data, pd.read_csv(file)], ignore_index=True)
        for file in self.neg_measures_files:
            neg_data = pd.concat([neg_data, pd.read_csv(file)], ignore_index=True)

        pos_data['class'] = 1
        neg_data['class'] = 0
        print([len(pos_data), len(neg_data)])
        self.X_all = pd.concat([pos_data, neg_data], ignore_index=True)

        #remove any rows with rubbish depth measurements
        df_suspicious = self.get_extreme_values(feature='depth', lower=0, upper=20)
        self.remove_df(df_suspicious)

        original_cols = self.X_all.columns.values
        cols_required = ['Uniprot_Entry', 'PDB_Code', 'Resid', 'class']
        for feature in self.features_to_include:
            if feature in original_cols:
                cols_required.append(feature)
            else:
                print(f'Feature given as input is not available in all input files, will not be included: {feature}')
        self.X_all = self.X_all[cols_required]

        # prepare aevs
        if 'aev' in self.features_to_include:
            self._prepare_aevs()
            #self._reduce_aevs()
        # aggregate based on aggregation method - will do both pos and neg data automatically
        self._aggregate()

        # create X
        # prepare the random numbers to get a random subset of the negative dataset to balance out the data
        num_pos_data = len(self.X_agg[self.X_agg['class'] == 1])
        print([num_pos_data, (len(self.X_agg) - num_pos_data)])
        indices_random_neg_data = random.sample(range(num_pos_data, len(self.X_agg)), num_pos_data)

        # prepare the X dataset
        self.X_final = pd.concat([self.X_agg[self.X_agg['class'] == 1], self.X_agg.iloc[indices_random_neg_data]])

        print([len(self.X_final[self.X_agg['class'] == 1]), len(self.X_final[self.X_agg['class'] == 0])])
        self.X_final = self.X_final.drop(['Uniprot_Entry', 'Resid', 'class'], axis=1)
        print(self.X_final.head())
        self._reduce_aevs_after()
        self._cut_columns()

        # create Y set - the classification set, 0 is negative, 1 is positive, positives go into the set first
        for i in range(num_pos_data):
            self.y.append(1)
        while len(self.y) < len(self.X_final):
            self.y.append(0)


    def _reduce_aevs_before(self):
        '''
        Model performs worse when more of the features of the AEV are taken through to training.
        Reduce the AEVs down to the required number of features based on one of the methods
        chosen below. The first method is to use the standard deviation of the individual features
        of the AEV to workout which features show variation and will be likely to be good choices
        to take through to the model. This is currently setup to find the top 100 from the base
        aevs fed in.
        
        Method
        ------
        Take the dataframe and perform a standard deviation over the AEVs, take the top
        required number of structures in terms of standard deviation as the new input
        dataframe going forward. 
        
        Example
        -------
        >> self.reduce_aevs()
        '''
        # TODO work on X_all - will eventually change all these functions so that they are general rather than single use
        temp_aevs_df = self.X_all['aev'].tolist()
        X_all_std = np.std(temp_aevs_df, axis=0)
        df_std = pd.DataFrame({'aev_std': X_all_std})
        df_std = df_std.sort_values(by=['aev_std'], ascending=False)
        top_n_features = df_std.index.values[:self.num_aev_features_req]
        self.top_n_features = top_n_features
        print('Top 100 aev features: ', self.top_n_features)


    def _reduce_aevs_after(self):
        '''
        Model performs worse when more of the features of the AEV are taken through to training.
        Reduce the AEVs down to the required number of features based on one of the methods chosen
        below. The first method is to use the standard deviation of the individual features of the
        AEV to workout which features show variation and will be likely to be good choices to take
        through to the model. This method does it after aggregation to ensure that the maximal
        variance after aggregation is taken through into the model.
        
        Method
        ------
        Take the dataframe and perform a standard deviation over the AEVs, take the top
        required number of structures in terms of standard deviation as the new input
        dataframe going forward. 
        
        Example
        -------
        >> self.reduce_aevs()
        '''
        # TODO work on X_all - will eventually change all these functions so that they are general rather than single use
        aev_stds = {}
        for column in self.X_final.columns:
            if 'AEV' in column:
                aev_stds[column] = self.X_final[column].std(ddof=0)

        df_std = pd.DataFrame({'aev_std': aev_stds})
        df_std = df_std.sort_values(by=['aev_std'], ascending=False)
        top_n_features = list(df_std.index.values[:self.num_aev_features_req])
        self.top_n_features = top_n_features
        print('Top 100 aev features: ', self.top_n_features)


    def _aggregate(self):
        '''
        Match aggregation type up to the relevant aggregation function
        
        Method
        ------
        Uses the input parameter of aggregation_method and calls the relevant function.
        If no cases match, assumes average and prints to terminal to state this.
        
        Example
        -------
        >> self._aggregate()
        '''
        match self.aggregation_method:
            case 'avg':
                self._aggregate_avg()
            case 'random':
                self._aggregate_random()
            case 'max':
                self._aggregate_max()
            case 'min':
                self._aggregate_min()
            case 'average subtract aev':
                self._aggregate_avg_less_avgaev()
            case 'mixmatch':
                self._aggregate_mixmatch()
            case _:
                print('Aggregation method not recognised; using average values')
                self._aggregate_avg()


    def add_average_sd(self, data_set):
        average = statistics.fmean(data_set)
        sd = statistics.stdev(data_set)
        data_set = np.append(data_set, average)
        data_set = np.append(data_set, sd)
        return data_set


    def _cut_columns(self):
        # first remove all the null columns
        num_cols_to_cut = sum((self.X_final != 0).any(axis=0))
        print(f'{num_cols_to_cut} columns kept from the AEV input data')
        self.X_final = self.X_final.loc[:, (self.X_final != 0).any(axis=0)]

        # then remove the least important columns according to the method given
        # drop all the columns that were not required through the std deviation
        print(f'Initial number of columns: {len(self.X_final.columns.tolist())}')
        if 'aev' in self.features_to_include and len(self.top_n_features) != 0:
            #self.top_n_features = [f'AEV_{x}' for x in self.top_n_features]
            current_feature_columns = self.X_final.columns.tolist()
            for feature in current_feature_columns:
                if 'AEV' in feature and feature not in self.top_n_features:
                    self.X_final = self.X_final.drop(feature, axis=1)
        print(f'Final number of columns: {len(self.X_final.columns.tolist())}')

    def get_extreme_values(self, feature, lower = 1, upper = 14):
        df_query = self.X_all[(self.X_all[feature] < lower) | (self.X_all[feature] > upper)]
        return df_query

    def remove_df(self, df_to_remove):
        len_one = len(self.X_all)
        remove_list = df_to_remove.index.tolist()
        self.X_all = self.X_all.drop(index = remove_list)
        len_two = len(self.X_all)
        print(f'Original num of rows: {len_one}')
        print(f'Current num of rows: {len_two}')
        print(f'Num of rows removed: {len(df_to_remove)}')



    # ------------------------------------------------------------------------------------
    #  Aggregation methods
    def _aggregate_avg(self):
        seperate_lys = self.X_all.drop_duplicates(subset=['Uniprot_Entry', 'Resid', 'class'])
        df_data_agg = pd.DataFrame(columns=['Uniprot_Entry', 'Resid', 'class'])
        for idx, row in seperate_lys.iterrows():
            entry = row['Uniprot_Entry']
            resid = row['Resid']
            class_val = row['class']
            # create a subset of the dataframe of measurements where the uniprot and resid match
            df_query = self.X_all[(self.X_all['Uniprot_Entry'] == entry) & (self.X_all['Resid'] == resid)]
            # calculate all the values of interest from the subset dataframe (df_query)
            data = {'Uniprot_Entry': entry,
                    'Resid' : resid,
                    'class' : class_val}
            for feature in self.features_to_include:
                if feature == 'aev':
                    list_aevs = df_query['aev'].tolist()
                    aev_avg = np.average(list_aevs, axis=0)
                    data['aev'] = aev_avg
                else:
                    temp_avg = round(df_query[feature].mean(),2)
                    data[feature] = temp_avg

            df_data_agg = pd.concat([df_data_agg, pd.DataFrame([data])], ignore_index=True)

        if 'aev' in self.features_to_include:
            df_out = pd.DataFrame(df_data_agg['aev'].to_list())
            df_out = df_out.add_prefix('AEV_')
            self.X_agg = pd.concat([df_data_agg, df_out], axis=1)
            self.X_agg = self.X_agg.drop('aev', axis=1)
            #aev_kept_cols = [f'AEV_{pos}' for pos in aev_kept_cols]
            #df_out = df_out.set_axis(aev_kept_cols, axis=1)
        else:
            self.X_agg = df_data_agg

    def _aggregate_max(self):
        seperate_lys = self.X_all.drop_duplicates(subset=['Uniprot_Entry', 'Resid', 'class'])
        df_data_agg = pd.DataFrame(columns=['Uniprot_Entry', 'Resid', 'class'])
        for idx, row in seperate_lys.iterrows():
            entry = row['Uniprot_Entry']
            resid = row['Resid']
            class_val = row['class']
            # create a subset of the dataframe of measurements where the uniprot and resid match
            df_query = self.X_all[(self.X_all['Uniprot_Entry'] == entry) & (self.X_all['Resid'] == resid)]
            # calculate all the values of interest from the subset dataframe (df_query)
            data = {'Uniprot_Entry': entry,
                    'Resid' : resid,
                    'class' : class_val}
            for feature in self.features_to_include:
                if feature == 'aev':
                    list_aevs = df_query['aev'].tolist()
                    aev_avg = np.max(list_aevs, axis=0)
                    data['aev'] = aev_avg
                else:
                    temp_max = round(df_query[feature].max(),2)
                    data[feature] = temp_max

            df_data_agg = pd.concat([df_data_agg, pd.DataFrame([data])], ignore_index=True)

        if 'aev' in self.features_to_include:
            df_out = pd.DataFrame(df_data_agg['aev'].to_list())
            df_out = df_out.add_prefix('AEV_')
            self.X_agg = pd.concat([df_data_agg, df_out], axis=1)
            self.X_agg = self.X_agg.drop('aev', axis=1)
            #aev_kept_cols = [f'AEV_{pos}' for pos in aev_kept_cols]
            #df_out = df_out.set_axis(aev_kept_cols, axis=1)
        else:
            self.X_agg = df_data_agg

    def _aggregate_min(self):
        seperate_lys = self.X_all.drop_duplicates(subset=['Uniprot_Entry', 'Resid', 'class'])
        df_data_agg = pd.DataFrame(columns=['Uniprot_Entry', 'Resid', 'class'])
        for idx, row in seperate_lys.iterrows():
            entry = row['Uniprot_Entry']
            resid = row['Resid']
            class_val = row['class']
            # create a subset of the dataframe of measurements where the uniprot and resid match
            df_query = self.X_all[(self.X_all['Uniprot_Entry'] == entry) & (self.X_all['Resid'] == resid)]
            # calculate all the values of interest from the subset dataframe (df_query)
            data = {'Uniprot_Entry': entry,
                    'Resid' : resid,
                    'class' : class_val}
            for feature in self.features_to_include:
                if feature == 'aev':
                    list_aevs = df_query['aev'].tolist()
                    aev_avg = np.min(list_aevs, axis=0)
                    data['aev'] = aev_avg
                else:
                    temp_min = round(df_query[feature].min(),2)
                    data[feature] = temp_min

            df_data_agg = pd.concat([df_data_agg, pd.DataFrame([data])], ignore_index=True)

        if 'aev' in self.features_to_include:
            df_out = pd.DataFrame(df_data_agg['aev'].to_list())
            df_out = df_out.add_prefix('AEV_')
            self.X_agg = pd.concat([df_data_agg, df_out], axis=1)
            self.X_agg = self.X_agg.drop('aev', axis=1)
            #aev_kept_cols = [f'AEV_{pos}' for pos in aev_kept_cols]
            #df_out = df_out.set_axis(aev_kept_cols, axis=1)
        else:
            self.X_agg = df_data_agg

    def _aggregate_avg_less_avgaev(self):
        seperate_lys = self.X_all.drop_duplicates(subset=['Uniprot_Entry', 'Resid', 'class'])
        df_data_agg = pd.DataFrame(columns=['Uniprot_Entry', 'Resid', 'class'])
        for idx, row in seperate_lys.iterrows():
            entry = row['Uniprot_Entry']
            resid = row['Resid']
            class_val = row['class']
            # create a subset of the dataframe of measurements where the uniprot and resid match
            df_query = self.X_all[(self.X_all['Uniprot_Entry'] == entry) & (self.X_all['Resid'] == resid)]
            # calculate all the values of interest from the subset dataframe (df_query)
            data = {'Uniprot_Entry': entry,
                    'Resid' : resid,
                    'class' : class_val}
            for feature in self.features_to_include:
                if feature == 'aev':
                    list_aevs = df_query['aev'].tolist()
                    aev_avg = np.average(list_aevs, axis=0)
                    data['aev'] = aev_avg
                else:
                    temp_avg = round(df_query[feature].mean(),2)
                    data[feature] = temp_avg

            df_data_agg = pd.concat([df_data_agg, pd.DataFrame([data])], ignore_index=True)

        if 'aev' in self.features_to_include:
            df_out = pd.DataFrame(df_data_agg['aev'].to_list())
            df_out = df_out.add_prefix('AEV_')
            self.X_agg = pd.concat([df_data_agg, df_out], axis=1)
            self.X_agg = self.X_agg.drop('aev', axis=1)
            #aev_kept_cols = [f'AEV_{pos}' for pos in aev_kept_cols]
            #df_out = df_out.set_axis(aev_kept_cols, axis=1)
        else:
            self.X_agg = df_data_agg

    def _aggregate_mixmatch(self):
        '''
        This one will do the best case for each individual feature
        Max sasa, das
        Min propka, pkaANI, depth
        Average AEV
        '''
        seperate_lys = self.X_all.drop_duplicates(subset=['Uniprot_Entry', 'Resid', 'class'])
        df_data_agg = pd.DataFrame(columns=['Uniprot_Entry', 'Resid', 'class'])
        max_features = ['sasa', 'das']
        min_features = ['propka', 'pkaANI', 'depth']
        for idx, row in seperate_lys.iterrows():
            entry = row['Uniprot_Entry']
            resid = row['Resid']
            class_val = row['class']
            # create a subset of the dataframe of measurements where the uniprot and resid match
            df_query = self.X_all[(self.X_all['Uniprot_Entry'] == entry) & (self.X_all['Resid'] == resid)]
            # calculate all the values of interest from the subset dataframe (df_query)
            data = {'Uniprot_Entry': entry,
                    'Resid' : resid,
                    'class' : class_val}
            for feature in self.features_to_include:
                if feature == 'aev':
                    list_aevs = df_query['aev'].tolist()
                    aev_avg = np.average(list_aevs, axis=0)
                    data['aev'] = aev_avg
                else:
                    if feature in max_features:
                        temp_agg = round(df_query[feature].max(),2)
                    elif feature in min_features:
                        temp_agg = round(df_query[feature].min(),2)
                    data[feature] = temp_agg

            df_data_agg = pd.concat([df_data_agg, pd.DataFrame([data])], ignore_index=True)

        if 'aev' in self.features_to_include:
            df_out = pd.DataFrame(df_data_agg['aev'].to_list())
            df_out = df_out.add_prefix('AEV_')
            self.X_agg = pd.concat([df_data_agg, df_out], axis=1)
            self.X_agg = self.X_agg.drop('aev', axis=1)
            #aev_kept_cols = [f'AEV_{pos}' for pos in aev_kept_cols]
            #df_out = df_out.set_axis(aev_kept_cols, axis=1)
        else:
            self.X_agg = df_data_agg

    def _aggregate_random(self):
        seperate_lys = self.X_all.drop_duplicates(subset=['Uniprot_Entry', 'Resid', 'class'])
        df_data_agg = pd.DataFrame(columns=['Uniprot_Entry', 'Resid', 'class'])
        for idx, row in seperate_lys.iterrows():
            entry = row['Uniprot_Entry']
            resid = row['Resid']
            class_val = row['class']
            # create a subset of the dataframe of measurements where the uniprot and resid match
            df_query = self.X_all[(self.X_all['Uniprot_Entry'] == entry) & (self.X_all['Resid'] == resid)]
            # calculate all the values of interest from the subset dataframe (df_query)
            data = {'Uniprot_Entry': entry,
                    'Resid' : resid,
                    'class' : class_val}
            idx_row_chosen = 0

            if len(df_query) >= 1:
                idx_row_chosen = random.sample(range(0, len(df_query)), 1)[0]

            row_to_append = df_query[idx_row_chosen:idx_row_chosen+1]

            for feature in self.features_to_include:
                data[feature] = row_to_append[feature].iloc[0]

            df_data_agg = pd.concat([df_data_agg, pd.DataFrame([data])], ignore_index=True)

        if 'aev' in self.features_to_include:
            df_out = pd.DataFrame(df_data_agg['aev'].to_list())
            df_out = df_out.add_prefix('AEV_')
            self.X_agg = pd.concat([df_data_agg, df_out], axis=1)
            self.X_agg = self.X_agg.drop('aev', axis=1)
            #aev_kept_cols = [f'AEV_{pos}' for pos in aev_kept_cols]
            #df_out = df_out.set_axis(aev_kept_cols, axis=1)
        else:
            self.X_agg = df_data_agg



    # ------------------------------------------------------------------------------------
    # model testing
    def rf_five_fold(self):
        '''
        Function to perform a 5 fold test on the data that has been presented to the
        model.

        Method
        ------
        Takes the data that has been prepared using the earlier functions and uses the
        SciKitLearn RandomForestClassifier to perform predictions on this. This uses the
        standard setup that it comes with without any optimisation. Metrics are calculated
        for each of the models and reported both individually and as an average after
        all have been completed. Functions for graphings the data are called to display
        this after.

        Example
        -------
        model.rf_five_fold()
        '''
        scores = ['accuracy', 'f1', 'precision']
        kf = KFold(n_splits=5, shuffle=True, random_state=False)

        clf = RandomForestClassifier(n_estimators=100)

        def confusion_matrix_scorer(model, X, y):
            y_pred = model.predict(X)
            cm = metrics.confusion_matrix(y, y_pred)
            return {'tn': cm[0, 0], 'fp': cm[0, 1],
                'fn': cm[1, 0], 'tp': cm[1, 1]}

        model_output = sklearn.model_selection.cross_validate(clf, self.X_final, self.y, cv=kf, scoring=confusion_matrix_scorer, return_estimator = True)
        accuracy_output = cross_val_score(clf, self.X_final, self.y, cv=kf, scoring='accuracy')
        f1_output = cross_val_score(clf, self.X_final, self.y, cv=kf, scoring='f1')
        precision = cross_val_score(clf, self.X_final, self.y, cv=kf, scoring='precision')

        df_input_importance = pd.DataFrame()
        for idx,estimator in enumerate(model_output['estimator']):
            feature_importances = pd.DataFrame(estimator.feature_importances_,
                                            index = self.X_final.columns.values,
                                            columns=[f'Run {idx + 1}'])
            df_input_importance = pd.concat([df_input_importance, feature_importances], axis=1)

        df_input_importance['Total'] = df_input_importance.sum(axis=1)


        test_tn = self.add_average_sd(model_output['test_tn'].astype(float))
        test_fp = self.add_average_sd(model_output['test_fp'].astype(float))
        test_fn = self.add_average_sd(model_output['test_fn'].astype(float))
        test_tp = self.add_average_sd(model_output['test_tp'].astype(float))

        # "Accuracy", "F1_Score", "Precision", "tn", "fp", "fn", "tp", "sensitivity", "specificity"
        output = pd.DataFrame(["1", "2", "3", "4", "5", "Average", "s.d."], columns=["Run"])
        accuracy_output = self.add_average_sd(accuracy_output)
        f1_output = self.add_average_sd(f1_output)
        precision = self.add_average_sd(precision)
        output['Accuracy'] = accuracy_output
        output['F1_Score'] = f1_output
        output['Precision'] = precision
        output['TN'] = test_tn
        output['FP'] = test_fp
        output['FN'] = test_fn
        output['TP'] = test_tp

        output['Sensitivity'] = ''
        output['Specificity'] = ''

        for i, r in output.iterrows():
            if i < 5:
                sensitivity = r['TP'] / (r['TP'] + r['FN'])
                specificity = r['TN'] / (r['TN'] + r['FP'])
                output.at[i, 'Sensitivity'] = sensitivity
                output.at[i, 'Specificity'] = specificity

        output.at[5, 'Sensitivity'] = statistics.fmean(output['Sensitivity'][:5])
        output.at[6, 'Sensitivity'] = statistics.stdev(output['Sensitivity'][:5])
        output.at[5, 'Specificity'] = statistics.fmean(output['Specificity'][:5])
        output.at[6, 'Specificity'] = statistics.stdev(output['Specificity'][:5])


        print(output)

        if len(self.features_to_include) == 1 and 'aev' in self.features_to_include:
            self.graph_aev_importance(df_input_importance)
        else:
            self.graph_top_ten_importance(df_input_importance)
        self.find_feature_importance_pos(df_input_importance)


    def rf_five_fold_optimised(self):
        '''
        Function to perform a 5 fold test on the data that has been presented to the
        model. A grid_search is performed before to optimise the RF model to the
        data it will be trained on.

        Method
        ------
        Takes the data that has been prepared using the earlier functions and uses a
        grid_search from SciKitLearn to find the optimised RF model from a seiries of
        trial runs on this. The SciKitLearn RandomForestClassifier is then used to get
        an accurate reading of how the model has performed. This uses the standard setup
        that it comes with without any optimisation. Metrics are calculated for each of
        the models and reported both individually and as an average after all have been
        completed. Functions for graphings the data are called to display this after.

        Example
        -------
        model.rf_five_fold_optimised()
        '''
        #define classifier
        rf = RandomForestClassifier()

        # define parameters and values to test
        params = {
            'max_depth': [2, 3, 5, 10, 20],
            'min_samples_leaf': [5, 10, 20, 50, 100, 200],
            'n_estimators': [10, 25, 30, 50, 100, 200]
        }

        # grid search over parameters space
        grid_search = GridSearchCV(estimator=rf,
                                param_grid=params,
                                cv = 4,
                                n_jobs=-1, verbose=1, scoring="accuracy")

        grid_search.fit(self.X_final, self.y)

        # best estimator found with grid search
        RF_best = grid_search.best_estimator_


        scores = ['accuracy', 'f1', 'precision']
        kf = KFold(n_splits=5, shuffle=True, random_state=False)

        def confusion_matrix_scorer(model, X, y):
            y_pred = model.predict(X)
            cm = metrics.confusion_matrix(y, y_pred)
            return {'tn': cm[0, 0], 'fp': cm[0, 1],
                'fn': cm[1, 0], 'tp': cm[1, 1]}

        model_output = sklearn.model_selection.cross_validate(RF_best, self.X_final, self.y, cv=kf, scoring=confusion_matrix_scorer, return_estimator = True)
        accuracy_output = cross_val_score(RF_best, self.X_final, self.y, cv=kf, scoring='accuracy')
        f1_output = cross_val_score(RF_best, self.X_final, self.y, cv=kf, scoring='f1')
        precision = cross_val_score(RF_best, self.X_final, self.y, cv=kf, scoring='precision')

        df_input_importance = pd.DataFrame()
        for idx,estimator in enumerate(model_output['estimator']):
            feature_importances = pd.DataFrame(estimator.feature_importances_,
                                            index = self.X_final.columns.values,
                                            columns=[f'Run {idx + 1}'])
            df_input_importance = pd.concat([df_input_importance, feature_importances], axis=1)

        df_input_importance['Total'] = df_input_importance.sum(axis=1)


        test_tn = self.add_average_sd(model_output['test_tn'].astype(float))
        test_fp = self.add_average_sd(model_output['test_fp'].astype(float))
        test_fn = self.add_average_sd(model_output['test_fn'].astype(float))
        test_tp = self.add_average_sd(model_output['test_tp'].astype(float))

        # "Accuracy", "F1_Score", "Precision", "tn", "fp", "fn", "tp", "sensitivity", "specificity"
        output = pd.DataFrame(["1", "2", "3", "4", "5", "Average", "s.d."], columns=["Run"])
        accuracy_output = self.add_average_sd(accuracy_output)
        f1_output = self.add_average_sd(f1_output)
        precision = self.add_average_sd(precision)
        output['Accuracy'] = accuracy_output
        output['F1_Score'] = f1_output
        output['Precision'] = precision
        output['TN'] = test_tn
        output['FP'] = test_fp
        output['FN'] = test_fn
        output['TP'] = test_tp

        output['Sensitivity'] = ''
        output['Specificity'] = ''

        for i, r in output.iterrows():
            if i < 5:
                sensitivity = r['TP'] / (r['TP'] + r['FN'])
                specificity = r['TN'] / (r['TN'] + r['FP'])
                output.at[i, 'Sensitivity'] = sensitivity
                output.at[i, 'Specificity'] = specificity

        output.at[5, 'Sensitivity'] = statistics.fmean(output['Sensitivity'][:5])
        output.at[6, 'Sensitivity'] = statistics.stdev(output['Sensitivity'][:5])
        output.at[5, 'Specificity'] = statistics.fmean(output['Specificity'][:5])
        output.at[6, 'Specificity'] = statistics.stdev(output['Specificity'][:5])


        print(output)

        if len(self.features_to_include) == 1 and 'aev' in self.features_to_include:
            self.graph_aev_importance(df_input_importance)
        else:
            self.graph_top_ten_importance(df_input_importance)
        self.find_feature_importance_pos(df_input_importance)


    def rf_ubq_test(self):
        '''
        Find the importance ranking of the non AEV features within the features used in the model


        Returns
        ----------
        df_ubq_analysis : dataframe
            Table containing a comparison between the true modfications for ubiquitin
            against the predicted ones for easy comparison.

        Example
        -------
        >> model.rf_ubq_test()
        '''
        # read in data about the ubiquitin data
        ubq_data = pd.read_csv('aev_data_randomforestmodel/measures_ubq.csv')
        y_test = [1, 0, 0, 0, 1, 1, 1]
        ubq_data['class'] = y_test
        original_cols = ubq_data.columns.values
        cols_required = ['Uniprot_Entry', 'PDB_Code', 'Resid', 'class']
        for feature in self.features_to_include:
            if feature in original_cols:
                cols_required.append(feature)
            else:
                print(f'Feature given as input is not available in all input files, will not be included: {feature}')
        ubq_data = ubq_data[cols_required]


        if 'aev' in self.features_to_include:
            for i, r in ubq_data.iterrows():
                new_aev = literal_eval(r['aev'])
                ubq_data.at[i, 'aev'] = new_aev

            df_ubq_analysis = pd.DataFrame(ubq_data['aev'].to_list())
            df_ubq_analysis = df_ubq_analysis.add_prefix('AEV_')
            ubq_data = pd.concat([ubq_data, df_ubq_analysis], axis=1)
            ubq_data = ubq_data.drop('aev', axis=1)

        ubq_data = ubq_data.drop(['Uniprot_Entry', 'PDB_Code', 'Resid', 'class'], axis=1)

        cols_to_keep = self.X_final.columns.values
        cols_to_remove = []
        for col_name in ubq_data.columns.values:
            if col_name not in cols_to_keep:
                cols_to_remove.append(col_name)
        ubq_data = ubq_data.drop(cols_to_remove, axis=1)

        # create the model and test on ubq data

        clf = RandomForestClassifier(n_estimators=100)
        clf.fit(self.X_final, self.y)
        y_pred = clf.predict(ubq_data)

        df_ubq_analysis = pd.DataFrame({'Lysines': [6, 11, 27, 29, 33, 48, 63],
                               'True Modifications': y_test,
                               'Predicted Modifications': y_pred})

        print(df_ubq_analysis)
        return df_ubq_analysis


    # --------------------------------------------------------------------------------------------
    # Importance graphing functions
    def graph_aev_importance(self, df_input_importance):
        '''
        Creates a graph over the entire AEV showing the importance rather than individual AEV
        values itself. Any parameters of the AEV which haven't been taken forward are set at -0.002
        in order to differentiate from the features that were fed into the model, but came out with
        importances of 0. The graph created shows 5 AEVs,
        one for each of the runs in the 5 fold test.

        Parameters
        ----------
        df_input_importance : dataframe
            dataframe containing all the features and the corresponding importances from the model

        Example
        -------
        >> self.find_feature_importance_pos(df_input_importance)
                                 Metric  Value
        0           das_importance_rank      1
        1  Num zero importance features    102
        '''
        # setup the grid for the graphs
        total_num_plots = 5
        fig, axs = plt.subplots(total_num_plots, sharex=True)
        plot_no = 0
        x = [idx_val.split('_')[-1] for idx_val in df_input_importance.index.tolist()]
        for i in range(total_num_plots):
            y = df_input_importance[f'Run {i + 1}'].tolist()
            plot_y = []
            plot_x = []
            temp_y_idx = 0
            for idcol in range(1008):
                if str(idcol) in x:
                    plot_y.append(y[temp_y_idx])
                    temp_y_idx += 1
                else:
                    plot_y.append(-0.002)
                plot_x.append(idcol)
            #print((plot_x, plot_y))
            axs[plot_no].bar(x=plot_x, height=plot_y, label=f'Run {i + 1}', color=(0.408, 0.141, 0.427))
            #axs[plotNo].set_ylabel(f'Time Progression (TimeStep : {frameStep}ns)')
            #axs[plotNo].set_title(f'Lysine No: {lys_res_nos[i]}', y=1, pad=-20)
            #axs[plot_no].set_ylim(top=0.5)
            axs[plot_no].set_xlabel('AEV')
            axs[plot_no].legend(loc='center right')


            plot_no += 1

            #pdbCode = testStructure.split('.')[0]
            #plt.savefig(f'testAEVimshow_{pdbCode}.png')
        fig.supylabel(f'Importance values of inputs from AEVs', ha='right')
        plt.subplots_adjust(wspace=0, hspace=0)
        plt.show()

    def graph_top_ten_importance(self, df_input_importance):
        '''
        Create a graph for the top ten most important features from the model used. 

        Parameters
        ----------
        df_input_importance : dataframe
            dataframe containing all the features and the corresponding importances from the model

        Example
        -------
        >> self.graph_top_ten_importance(df_input_importance)
        '''
        sorted_df = df_input_importance.sort_values('Total', ascending=False)
        top_ten = sorted_df.iloc[:10]

        fig, ax = plt.subplots(layout='constrained')
        ax = top_ten.plot(kind='bar')
        #fig = ax.get_figure()
        ax.set_xlabel('Top Ten Features')
        ax.set_ylabel('Feature Importance')
        plt.show()

    def find_feature_importance_pos(self, df_input_importance):
        '''
        Find the importance ranking of the non AEV features within the features used in the model

        Method
        ------
        Sort the total importance column to give an order.
        Loop over the features that have been included (excluding the aev)
        and find the position of each feature.
        Also, do a search over the dataframe to identify all the features with 0 importance.

        Parameters
        ----------
        df_input_importance : dataframe
            dataframe containing all the features and the corresponding importances from the model

        Example
        -------
        >> self.find_feature_importance_pos(df_input_importance)
                                 Metric  Value
        0           das_importance_rank      1
        1  Num zero importance features    102
        '''
        sorted_df = df_input_importance.sort_values('Total', ascending=False)
        index_positions = sorted_df.index.tolist()
        feature_pos = {}
        temp_features = dc(self.features_to_include)
        temp_features.remove('aev')
        for feature in temp_features:
            pos = index_positions.index(feature)
            feature_pos[f'{feature}_importance_rank'] = pos + 1
        zero_total_df = sorted_df[sorted_df['Total'] == 0]
        feature_pos['Num zero importance features'] = len(zero_total_df)
        feature_pos_df = pd.DataFrame(list(feature_pos.items()), columns=['Metric', 'Value'])
        print(feature_pos_df)


# --------------------------------------------------------------------------------------------------
# Testing section


if __name__ == "__main__":

    pos_measure_files = ['data/measures_cut_CannPositiveData_all.csv',
                        'data/measures_cut_KingHighConfData_all.csv',
                        'data/measures_cut_Ecoli(hCit)_all_3.csv',
                        'data/measures_cut_Synecho(hCit)_all_2.csv']
    #pos_measure_files = ['data/measures_cut_Synecho(hCit)_all_2.csv']
    neg_measure_files = ['data/measures_cut_KingAllNegative_all.csv']
    features_to_include = ['aev']
    model = Model(pos_measures_files=pos_measure_files,
                neg_measures_files=neg_measure_files,
                features_to_include=features_to_include,
                aggregation_method='min',
                subtract_avg_aev=False,
                num_aev_features_req=100)
    model.prepare_dataset()
    #print(model.X_final)

    model.rf_ubq_test()
    model.rf_five_fold()
    model.rf_five_fold_optimised()
