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
from sklearn.decomposition import PCA
from preprocessing import Preprocessing
from aggregation import Aggregation

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
                 subtract_avg_aev=False,
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
        self.X_all.rename(columns={'aev': 'aev_old'}, inplace=True)
        self.X_all['aev'] = 1
        self.X_all['aev'] = self.X_all['aev'].astype('object')

        for i, r in self.X_all.iterrows():
            #new_aev = literal_eval(r['aev'].values[0])
            print(literal_eval(r['aev_old'].values[0]))
            self.X_all.at[i, 'aev'] = literal_eval(r['aev_old'].values[0])
            #self.X_all.loc[[i], 'aev'] = pd.Series([new_aev], index=self.X_all.index[[i]])

        self.X_all.drop(['aev_old'], axis=1)

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
        pos_data = pd.DataFrame()
        neg_data = pd.DataFrame()
        for file in self.pos_measures_files:
            pos_data = pd.concat([pos_data, pd.read_csv(file)], ignore_index=True)
        for file in self.neg_measures_files:
            neg_data = pd.concat([neg_data, pd.read_csv(file)], ignore_index=True)

        pos_data['class'] = 1
        neg_data['class'] = 0
        print(f'Len positve data: {len(pos_data)}, Len negative data: {len(neg_data)}')
        self.X_all = pd.concat([pos_data, neg_data], ignore_index=True)

        if 'aev_legolas' in self.X_all.columns:
            self.X_all.drop(columns=['aev'], inplace=True)
            self.X_all.rename(columns={'aev_legolas': 'aev'}, inplace=True)
        self.X_all.dropna(subset=self.features_to_include, inplace=True)

        self._check_for_overlap()

        #remove any rows with rubbish depth measurements - can be expanded into other features if needed too
        if 'depth' in self.features_to_include:
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
        self.y_all = list(self.X_all['class'])

        agg = Aggregation(self.X_all, aggregation_method='minmax', features_to_include=['all'])
        self.X_agg = agg.aggregate_data()

        # create X
        # prepare the random numbers to get a random subset of the negative dataset to balance out the data
        print(self.X_agg)
        #self.X_agg.to_csv('X_agg.csv')
        num_pos_data = len(self.X_agg[self.X_agg['class'] == 1])
        print([num_pos_data, (len(self.X_agg) - num_pos_data)])
        indices_random_neg_data = random.sample(range(num_pos_data, len(self.X_agg)), num_pos_data)

        # prepare the X dataset
        self.X_final = pd.concat([self.X_agg[self.X_agg['class'] == 1], self.X_agg.iloc[indices_random_neg_data]])

        print([len(self.X_final[self.X_agg['class'] == 1]), len(self.X_final[self.X_agg['class'] == 0])])
        self.X_final = self.X_final.drop(['Uniprot_Entry', 'Resid', 'class'], axis=1)
        print(self.X_final.head())


        # create Y set - the classification set, 0 is negative, 1 is positive, positives go into the set first
        for i in range(num_pos_data):
            self.y.append(1)
        while len(self.y) < len(self.X_final):
            self.y.append(0)


    def _check_for_overlap(self):
        '''
        Function to check for overlap between the positive and negative parts of the dataset
        If there are any overlaps, remove the data from the negative dataset.
        '''

        all_overlaps = pd.DataFrame(columns=['Uniprot_Entry', 'PDB_Code', 'Resid'])

        for i, r in self.X_all[self.X_all['class'] == 1].iterrows():
            pdb = r['PDB_Code']
            resid = r['Resid']
            # create a subset of the dataframe of measurements where the uniprot and resid match
            df_query = self.X_all[(self.X_all['PDB_Code'] == pdb) & (self.X_all['Resid'] == resid)]
            if len(df_query) > 1:
                #print(f'Overlap found: {pdb} {resid}')
                all_overlaps = pd.concat([all_overlaps, df_query[['Uniprot_Entry', 'PDB_Code', 'Resid', 'class']]], ignore_index=True)
                self.X_all = self.X_all.drop(index=df_query.index[1:])

        print(f'Number of overlaps: {len(all_overlaps)}')
        #print(all_overlaps)


    def add_average_sd(self, data_set):
        average = statistics.fmean(data_set)
        sd = statistics.stdev(data_set)
        data_set = np.append(data_set, average)
        data_set = np.append(data_set, sd)
        return data_set


    def _cut_columns(self):
        # first remove all the null columns
        num_cols_to_cut = sum((self.X_final != 0).any(axis=0))
        #print([i for i, a in enumerate(list((self.X_final != 0).any(axis=0))) if a is False])
        print(f'{num_cols_to_cut} columns kept from the AEV input data that are non zero')
        self.X_final = self.X_final.loc[:, (self.X_final != 0).any(axis=0)]

        # then remove the least important columns according to the method given
        # drop all the columns that were not required through the std deviation
        print(f'Initial number of columns: {len(self.X_final.columns.tolist())}')
        if 'aev' in self.features_to_include and len(self.top_n_features) != 0:
            #self.top_n_features = [f'AEV_{x}' for x in self.top_n_features]
            current_feature_columns = self.X_final.columns.tolist()
            for feature in current_feature_columns:
                #if 'AEV' in feature and int(feature.split('_')[-1]) not in self.top_n_features:
                if 'AEV' in feature and feature not in self.top_n_features:
                    self.X_final = self.X_final.drop(feature, axis=1)
        print(f'Final number of columns after reduction: {len(self.X_final.columns.tolist())}')

    def get_extreme_values(self, feature, lower = 1, upper = 14):
        '''
        Obtain any values within the dataframe that are deemed to be values that are likely wrong.

        Parameters
        ----------
        feature -> string
            The feature of interest to get the extreme values from
        lower -> float
            The lowest acceptable value for the feature of interest
        upper -> float
            The highest acceptable value for the feature of interest
        
        Returns
        -------
        df_query -> dataframe
            A subset of the dataframe which are the extreme values within the dataframe X_all
        '''
        df_query = self.X_all[(self.X_all[feature] < lower) | (self.X_all[feature] > upper)]
        return df_query

    def remove_df(self, df_to_remove):
        '''
        Function for removing values from a dataframe give a dataframe of the items to remove

        Parameters
        ----------
        df_to_remove -> dataframe
            dataframe of items from the overall dataset you want to remove from X_all
        '''
        len_one = len(self.X_all)
        remove_list = df_to_remove.index.tolist()
        self.X_all = self.X_all.drop(index = remove_list)
        len_two = len(self.X_all)
        print(f'Original num of rows: {len_one}')
        print(f'Current num of rows: {len_two}')
        print(f'Num of rows removed: {len(df_to_remove)}')


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
        self.avg_output = output.iloc[5].to_dict()

        '''
        if len(self.features_to_include) == 1 and 'aev' in self.features_to_include:
            self.graph_aev_importance(df_input_importance)
        else:
            self.graph_top_ten_importance(df_input_importance)
        '''
        self.find_feature_importance_pos(df_input_importance)


    def rf_five_fold_optimised_expanded(self):
        '''
        Function to perform a 5 fold test on the data that has been presented to the
        model. A grid_search is performed before to optimise the RF model to the
        data it will be trained on.
        This one also takes the datasets that are produced throught the rf model and
        expands them such that the training and test data include all the measurements
        that are available without including individual lysines in both training and test

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

        # sort the data out for this
        x_all_agg = self.X_all.drop_duplicates(subset=['Uniprot_Entry', 'Resid', 'class'])
        print(x_all_agg[x_all_agg.duplicated(subset=['Uniprot_Entry', 'Resid'], keep=False)].sort_values(by=['Uniprot_Entry', 'Resid']))
        x_all_agg = x_all_agg.drop_duplicates(subset=['Uniprot_Entry', 'Resid'], keep='first', inplace=False).dropna()
        print(f'len of each; pos:{len(x_all_agg[x_all_agg["class"] == 1])}, neg: {len(x_all_agg[x_all_agg["class"] == 0])}')

        num_pos_data = len(x_all_agg[x_all_agg['class'] == 1])
        print([num_pos_data, (len(x_all_agg) - num_pos_data)])
        indices_random_neg_data = random.sample(range(num_pos_data, len(x_all_agg)), num_pos_data)
        x_all_agg = pd.concat([x_all_agg[x_all_agg['class'] == 1], x_all_agg.iloc[indices_random_neg_data]])
        y_all_agg = x_all_agg['class']

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

        grid_search.fit(x_all_agg.drop(['Uniprot_Entry', 'PDB_Code', 'Resid', 'class'], axis=1), y_all_agg)

        # best estimator found with grid search
        RF_best = grid_search.best_estimator_


        def confusion_matrix_scorer(model, X, y):
            y_pred = model.predict(X)
            cm = metrics.confusion_matrix(y, y_pred)
            return {'tn': cm[0, 0], 'fp': cm[0, 1],
                    'fn': cm[1, 0], 'tp': cm[1, 1]}

        scores = ['accuracy', 'f1', 'precision']
        df_input_importance = pd.DataFrame()
        kf = KFold(n_splits=5, shuffle=True, random_state=False)
        for train, test in kf.split(x_all_agg, y_all_agg):
            print(f'train: {len(train)}, test: {len(test)}, total: {len(train) + len(test)}')

            # extend the indices into the full dataset
            train_extended = []
            test_extended = []
            for i in train:
                entry = x_all_agg.iloc[i]['Uniprot_Entry']
                resid = x_all_agg.iloc[i]['Resid']
                df_query = self.X_all[(self.X_all['Uniprot_Entry'] == entry) & (self.X_all['Resid'] == resid)]
                train_extended += list(df_query.index)
            for i in test:
                entry = x_all_agg.iloc[i]['Uniprot_Entry']
                resid = x_all_agg.iloc[i]['Resid']
                df_query = self.X_all[(self.X_all['Uniprot_Entry'] == entry) & (self.X_all['Resid'] == resid)]
                test_extended += list(df_query.index)

            print(f'train: {len(train_extended)}, test: {len(test_extended)}, total: {len(train_extended) + len(test_extended)}')            
            train_data_X = self.X_all.loc[train_extended].drop(['Uniprot_Entry', 'PDB_Code', 'Resid', 'class'], axis=1)
            self.X_all['class'] = self.y_all
            train_data_y = list(self.X_all.loc[train_extended]['class'])
            print(f'Balance of training data: 1: {len([a for a in train_data_y if a == 1])}, 0: {len([a for a in train_data_y if a == 0])}')

            #print(train_data_X)

            RF_best.fit(train_data_X, train_data_y)

            test_data_X = self.X_all.loc[test_extended].drop(['Uniprot_Entry', 'PDB_Code', 'Resid', 'class'], axis=1)
            test_data_y = list(self.X_all.loc[test_extended]['class'])

            test_y_pred = RF_best.predict(test_data_X)
            acc = metrics.accuracy_score(test_data_y, test_y_pred)
            print(f'Accuracy: {acc}')

        df_input_importance['Total'] = df_input_importance.sum(axis=1)


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
        ubq_data = pd.read_csv('data/measures_ubq.csv')
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
        #plt.show()
        plt.close()

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
        #plt.show()
        plt.close()

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
        self.sorted_df_importances = sorted_df['Total'].to_dict()
        index_positions = sorted_df.index.tolist()
        feature_pos = {}
        temp_features = dc(self.features_to_include)
        if 'aev' in temp_features:
            temp_features.remove('aev')
        
        '''
        if self.aggregation_method == 'minmax':
            temp_features = [f'{feature}_min' for feature in temp_features] + [f'{feature}_max' for feature in temp_features]
        '''
            
        for feature in temp_features:
            pos = index_positions.index(feature)
            feature_pos[f'{feature}_importance_rank'] = pos + 1
        zero_total_df = sorted_df[sorted_df['Total'] == 0]
        feature_pos['Num zero importance features'] = len(zero_total_df)
        feature_pos_df = pd.DataFrame(list(feature_pos.items()), columns=['Metric', 'Value'])
        print(feature_pos_df.sort_values('Value', ascending=True))


    def multi_run_test(self, num_runs=1):
        '''
        Run the complete setup model multiple times to get an average of the results
        including average importance rankings

        Parameters
        ----------
        num_runs : int
            The number of times to run the model

        Example
        -------
        >> self.multi_run_test(num_runs=10)
        '''
        df_output = pd.DataFrame(columns=['Run', 'Accuracy', 'F1_Score', 'Precision', 'TN', 'FP', 'FN', 'TP', 'Sensitivity', 'Specificity'])
        df_importance = pd.DataFrame()
        for i in range(num_runs):
            self.rf_five_fold_optimised()
            df_output = pd.concat([df_output, pd.DataFrame([self.avg_output])], ignore_index=True)
            df_importance = pd.concat([df_importance, pd.DataFrame([self.sorted_df_importances])], ignore_index=True)

        df_output = df_output.drop(columns=['Run'])
        df_output.loc['avg'] = df_output.mean()
        df_importance.loc['avg'] = df_importance.mean()
        print(df_output)
        print(df_importance.T.sort_values('avg', ascending=False))
        return df_output, df_importance.T.sort_values('avg', ascending=False)


# --------------------------------------------------------------------------------------------------
# Testing section


def create_comparison_between_datasets(num_runs=5):
    '''
    Function to create a comparison between the different datasets that are available
    for the model. This is to see the difference in performance between the different
    datasets.

    Example
    -------
    >> create_comparison_between_datasets()
    '''

    features_to_include = ['propka', 'sasa', 'das', 'seqcharge', 'curvature', 'writhing', 'torsion', 'arc_length', 'phi', 'psi', 'legolas', 'aev']
    pos_measure_files = [['data/measures_cut_CannData_all_12.05.25.csv',
                        'data/measures_cut_KingHighConf_all_01.05.25_joined.csv',
                        'data/measures_cut_Ecoli(hCit)_all_01.05.25_joined.csv',
                        'data/measures_cut_Synecho(hCit)_all_01.05.25_joined.csv'],
                        ['data/measures_cut_CannData_all_12.05.25.csv'],
                        ['data/measures_cut_KingHighConf_all_01.05.25_joined.csv',
                        'data/measures_cut_Ecoli(hCit)_all_01.05.25_joined.csv',
                        'data/measures_cut_Synecho(hCit)_all_01.05.25_joined.csv']]
    neg_measure_files = ['data/measures_cut_KingAllNegative_all_01.05.25_joined.csv']
    rankings = []
    scores_accuracy = pd.DataFrame()
    scores_f1score = pd.DataFrame()
    scores_sensitivity = pd.DataFrame()
    scores_specificity = pd.DataFrame()

    set_names = ['Cann + King Data', 'Cann Data', 'King Data']
    for a, file in enumerate(pos_measure_files):
        tmp_scores_accuracy = []
        tmp_scores_f1score = []
        tmp_scores_sensitivity = []
        tmp_scores_specificity = []
        tmp_rankings = []
        for i in range(num_runs):
            model = Model(pos_measures_files=file,
                        neg_measures_files=neg_measure_files,
                        features_to_include=features_to_include,
                        aggregation_method='minmax',
                        subtract_avg_aev=False,
                        num_aev_features_req=100)
            model.prepare_dataset()
            multi_run_output, multi_run_importance = model.multi_run_test(num_runs=1)
            tmp_scores_accuracy.append(multi_run_output['Accuracy'].loc['avg'])
            tmp_scores_f1score.append(multi_run_output['F1_Score'].loc['avg'])
            tmp_scores_sensitivity.append(multi_run_output['Sensitivity'].loc['avg'])
            tmp_scores_specificity.append(multi_run_output['Specificity'].loc['avg'])
            tmp_rankings.append(multi_run_importance)
        print(tmp_scores_accuracy)
        print(tmp_scores_f1score)
        print(tmp_scores_sensitivity)
        print(tmp_scores_specificity)
        #new_setname = file.split('/')[-1].split('_')[2]
        new_setname = set_names[a]
        scores_accuracy = pd.concat([scores_accuracy, pd.DataFrame({new_setname: tmp_scores_accuracy})], axis=1)
        scores_f1score = pd.concat([scores_f1score, pd.DataFrame({new_setname: tmp_scores_f1score})], axis=1)
        scores_specificity = pd.concat([scores_specificity, pd.DataFrame({new_setname: tmp_scores_specificity})], axis=1)
        scores_sensitivity = pd.concat([scores_sensitivity, pd.DataFrame({new_setname: tmp_scores_sensitivity})], axis=1)
        rankings.append(tmp_rankings)
        set_names.append(new_setname)

    print(scores_accuracy)
    print(scores_f1score)
    print(scores_sensitivity)
    print(scores_specificity)
    '''
    #set_names = ['All', 'CannData']
    for i, df in enumerate(scores_accuracy):
        sets_acc = pd.concat([sets_acc, df['Accuracy']], axis=1)
        sets_acc.rename(columns={'Accuracy': set_names[i]}, inplace=True)
    
    print(sets_acc)
    '''

    def plot_comparison(set, metric):
        plt.close()
        fig, ax = plt.subplots()
        set.boxplot(grid=False)
        plt.xlabel('Positive Datasets')
        plt.ylabel(metric)
        #column=sets_acc.columns, ax=ax
        #plt.show()
        plt.savefig(f'dataset_comparison_3way_2_{num_runs}_{metric}_repreparedata.png')
        plt.savefig(f'dataset_comparison_3way_2_{num_runs}_{metric}_repreparedata.svg')
    
    datasets = [scores_accuracy, scores_f1score, scores_sensitivity, scores_specificity]
    metrics = ['Accuracy', 'F1_Score', 'Sensitivity', 'Specificity']
    for dataset, metric in zip(datasets, metrics):
        plot_comparison(dataset, metric)

if __name__ == "__main__":

    pos_measure_files = ['data/measures_cut_CannData_all_12.05.25.csv',
                        'data/measures_cut_KingHighConf_all_01.05.25_joined.csv',
                        'data/measures_cut_Ecoli(hCit)_all_01.05.25_joined.csv',
                        'data/measures_cut_Synecho(hCit)_all_01.05.25_joined.csv']
    #pos_measure_files = ['data/measures_cut_Synecho(hCit)_all_2.csv']
    neg_measure_files = ['data/measures_cut_KingAllNegative_all_01.05.25_joined.csv']
    #pos_measure_files = ['data/measures_cut_Synecho(hCit)_all_01.05.25.csv']
    #neg_measure_files = ['April25_negONLY.csv']
    features_to_include = ['propka', 'sasa', 'das', 'seqcharge', 'curvature', 'writhing', 'torsion', 'arc_length', 'phi', 'psi', 'legolas', 'aev']
    model = Model(pos_measures_files=pos_measure_files,
                neg_measures_files=neg_measure_files,
                features_to_include=features_to_include,
                aggregation_method='minmax',
                subtract_avg_aev=False,
                num_aev_features_req=100)
    model.prepare_dataset()
    #print(model.X_final)

    #model.rf_ubq_test()
    #model.rf_five_fold()
    model.rf_five_fold_optimised()
    #model.multi_run_test(num_runs=1)
    #model.rf_five_fold_optimised_expanded()

    #create_comparison_between_datasets(3)
