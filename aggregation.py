import os
import random
import statistics
from ast import literal_eval
from copy import deepcopy as dc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.decomposition import PCA
from preprocessing import Preprocessing


class Aggregation:
    '''
    Class to handle the different aggregation methods for sorting over the
    different measurements that have been calculated. This will take a
    dataframe of measurements, and return a dataframe with only one set
    of measurements per lysine residue according to the aggregation method.
    '''

    def __init__(self, df_measurements, aggregation_method='minmax',
                 features_to_include=['all'], aev_red_method='pca',
                 num_sd_aev_features=100):
        '''
        Initialisation of the Aggregation class.

        Parameters
        ----------
        df_measurements : DataFrame
            The dataframe of measurements which need to be aggregated for use
        aggregation_method : str
            The aggregation method chosen to reduce the measurements into a usable format
            for training on. Options are:
            - 'avg': Take the average of all measurements for each lysine
            - 'random': Take a random measurement out of all measurements for the lysine
            - 'max': Take the maximum value of each feature for the lysine
            - 'min': Take the minimum value of each feature for the lysine
            - 'average subtract aev': Take the average of all aev features except the aevs WHAT IS USED HERE
            - 'mixmatch': Takes the predicted metric which will work best for each feature. The
                          max is used for features where a high value is likely to be important
                          and min for features where less of it is required.
            - 'minmax': Takes both the min and max values of each feature which duplicates the
                        the feature space allowing the model to use both as needed
            - 'minmaxavg': Takes the min, max and average of each feature, tripling feature space
        features_to_include : list, optional
            The list of features that are to be included in the aggregation. The default
            for this is taken to be all of them.
        aev_red_method : string, optional
            The dimensionality reduction method for reducing the size of the AEVs, reducing
            clouding. PCA is taken as default for this. Options:
            - 'PCA': Create a PCA which represents 99% of the variance of the data and use
                     these new vectors to represent the AEVs instead.
            - 'sd': Take the standard deviation of all the AEV columns and work out which the
                    top n are taken through for use
            - 'vif': Uses variance inflation factors to calculate the decorrelation between
                     the different columns of the AEV. Only takes through the features which
                     are shown to not be correlated.
            - 'autoencoder': Uses an autoencoder to reduce the dimensions, better for non-linear
                             data.
            - 'null': Removes all the null columns from the AEV
        num_sd_aev_features : int, optional
            The number of features to keep from the aevs when the standard deviation method is
            used. Defualt is set to 100.
        '''
        self.df_measurements = df_measurements
        self.aggregation_method = aggregation_method
        self.features_to_include = features_to_include
        self.aev_red_method = aev_red_method
        self.num_sd_aev_features = num_sd_aev_features

        if self.features_to_include == ['all']:
            self.features_to_include = [a for a in self.df_measurements.columns if a not in ['Uniprot_Entry', 'PDB_Code', 'Resid']]


    def aggregate_data(self):
        '''
        Match aggregation type up to the relevant aggregation function
        
        Method
        ------
        Uses the input parameter of aggregation_method and calls the relevant function.
        If no cases match, assumes average and prints to terminal to state this.
        
        Example
        -------
        >> self.aggregate_data()
        '''
        match self.aggregation_method:
            case 'avg':
                return self._aggregate_avg()
            case 'random':
                return self._aggregate_random()
            case 'max':
                return self._aggregate_max()
            case 'min':
                return self._aggregate_min()
            case 'average subtract aev':
                return self._aggregate_avg_less_avgaev()
            case 'mixmatch':
                return self._aggregate_mixmatch()
            case 'minmax':
                return self._aggregate_minmax()
            case 'minmaxavg':
                return self._aggregate_minmaxavg()
            case _:
                print('Aggregation method not recognised; using minmax values')
                return self._aggregate_minmax()


    def _reduce_aev_dimensions(self):
        '''
        Due to curse of dimenstionality the model performs worse when the AEVs are clouding the
        data as it cannot work out which features are actually important. Therefore, this function
        will call the required method to reduce the AEVs down to a specified number of features
        depending on the method chosen. Method options:
        - 'pca': Principal Component Analysis taking 99% of the variance within the data
        - 'sd': Standard Deviation of the AEV
        - 'vif': Variance Inflation Factor Correlation analysis to remove features which are correlated
        - 'autoencoder': Autoencoder dimension reduction method, try and capture any non-linearity
        - 'null': Remove all the columns within the AEVs which are always zero
        '''
        if 'aev' in self.features_to_include:
            df_aevs = pd.DataFrame(list([literal_eval(aev) for aev in self.df_measurements['aev']]))
            df_aevs = df_aevs.add_prefix('AEV_')
            self.df_measurements = pd.concat([self.df_measurements.reset_index(), df_aevs.reset_index()], axis=1)
            self.df_measurements = self.df_measurements.drop(['index', 'aev'], axis=1)
            match self.aev_red_method:
                case 'pca':
                    self._prepare_pca()
                case 'sd':
                    self._prepare_aev_sd()
                case 'vif':
                    self._prepare_vif_aev()
                case 'autoencoder':
                    i = 1
                    # TODO: implement the autoencoder method here 21.08.25
                case 'null':
                    self._cut_null_aev_columns()
                case _:
                    print(f'>> AEV dimensionality reduction method: {self.aev_red_method}, was not recognised, using PCA method.')
                    self._prepare_pca()


    def _prepare_aev_sd(self):
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
        print('>> Reducing AEV dimensions with standard deviation...')
        aev_stds = {}
        aev_cols = [a for a in self.df_measurements.columns if 'AEV_' in a]
        for col in aev_cols:
            aev_stds[col] = self.df_measurements[col].std(ddof=0)

        df_std = pd.DataFrame({'aev_std': aev_stds})
        df_std = df_std.sort_values(by=['aev_std'], ascending=False)
        top_n_features = list(df_std.index.values[:self.num_sd_aev_features])
        aev_cols = [a for a in self.df_measurements.columns if 'AEV_' in a]
        cols_to_remove = [a for a in aev_cols if a not in top_n_features]
        self.df_measurements.drop(cols_to_remove, axis=1, inplace=True)
        print(f'>> Removed {len(cols_to_remove)} from AEVs, {len(top_n_features)} kept instead of the AEVs.')



    def _cut_null_aev_columns(self):
        '''
        Function takes all the AEV columns and removes any that are always null which would
        add nothing to the model except noise.
        '''
        print('>> Removing null AEV columns...')
        initial_num_cols = len(self.df_measurements.columns)
        num_cols_to_keep = sum((self.df_measurements != 0).any(axis=0))
        self.df_measurements = self.df_measurements.loc[:, (self.df_measurements != 0).any(axis=0)]
        print(f'>> Removed {initial_num_cols - len(self.df_measurements.columns)} null AEV columns, {num_cols_to_keep} AEV colums left.')


    def _prepare_vif_aev(self):
        '''
        Use the preprocessing module to calculate the variance inflation factor values
        for each of the columns within the AEVs and remove the columns which are highly
        correlated together such that it is the minimum number of columns without
        correlation. Non-correlation was taken to be a VIF value of less than 5.
        
        Method
        ------
        Take the dataframe of columns of the AEV and plug this into the preprocessing
        module to calculate the VIF values for each of these. This will return a list
        of decorrelated columns which can be used to cut down the dataframe.

        Example
        -------
        >> self._prepare_vif_aev()
        '''
        # The columns_to_keep that is commented out is the oriignal set calculated by Phong
        # over the negative dataset
        '''
        # this columns_to_keep is the original set calculated with the original AEVs
        columns_to_keep = [12,120,135,16,17,182,19,197,20,211,22,228,23,231,26,27,28,29,
                            30,31,339,34,344,35,351,365,366,367,37,370,372,38,387,39,396,
                            399,40,407,41,415,42,425,428,43,431,44,442,443,45,46,463,47,
                            50,51,52,53,54,543,544,555,556,557,558,563,57,573,579,58,580,
                            583,588,59,590,591,60,61,62,622,63,689,696,699,70,704,705,706,
                            709,711,716,719,73,74,745,75,751,76,77,78,79,9]
        '''
        print('>> Reducing AEV dimensions using VIF analysis...')
        # use the preprocesing module to come up with exact columns to keep via VIF analysis
        non_aev_cols = [a for a in self.df_measurements.columns if 'AEV_' not in a]
        P = Preprocessing(self.df_measurements, non_aev_cols)
        columns_to_keep = P.calculate_diff_features(self.df_measurements)

        top_aev_feature_names = ['AEV_' + str(a) for a in columns_to_keep]
        cols_to_remove = [a for a in self.df_measurements.columns if a not in top_aev_feature_names]
        self.df_measurements.drop(cols_to_remove, axis=1, inplace=True)
        print(f'>> VIF analysis reduced AEV dimensions from {len(non_aev_cols)} to {len(columns_to_keep)}.')


    def _prepare_pca(self):
        '''
        Short function to transform the AEV data using PCA to reduce the dimensions.
        Changes the data in self.X_all ready for the aggregation to actually reduce
        the dataset for training on
        '''
        print('>> Calculating PCA on AEV data...')
        n_components = 0.99
        pca = PCA(n_components=n_components)
        aev_col_names = [a for a in self.df_measurements.columns if 'AEV_' in a]
        aev_data = self.df_measurements[aev_col_names].values.tolist()
        pca_data = pca.fit_transform(aev_data)

        #classification_values = self.X_final['class'].values
        variance_data = pca.explained_variance_ratio_

        self.df_measurements.drop(aev_col_names, axis=1, inplace=True)
        df_out = pd.DataFrame(pca_data)
        df_out = df_out.add_prefix('AEV_')
        self.df_measurements = pd.concat([self.df_measurements, df_out], axis=1)
        print(f'>> PCA used {pca.n_components_} components to explain {n_components} variance. AEVs are now in reduced dimension format.')


    def _aggregate_avg(self):
        seperate_lys = self.df_measurements.drop_duplicates(subset=['Uniprot_Entry', 'Resid', 'class'])
        df_data_agg = pd.DataFrame(columns=['Uniprot_Entry', 'Resid', 'class'])

        self._reduce_aev_dimensions()

        for idx, row in seperate_lys.iterrows():
            entry = row['Uniprot_Entry']
            resid = row['Resid']
            class_val = row['class']
            # create a subset of the dataframe of measurements where the uniprot and resid match
            df_query = self.df_measurements[(self.df_measurements['Uniprot_Entry'] == entry) & (self.df_measurements['Resid'] == resid)]
            # calculate all the values of interest from the subset dataframe (df_query)
            data = {'Uniprot_Entry': entry,
                    'Resid' : resid,
                    'class' : class_val}
            df_query = df_query.copy()
            for feature in self.features_to_include:
                if feature == 'aev':
                    list_aevs = df_query['aev'].tolist()
                    aev_avg = np.average(list_aevs, axis=0)
                    data['aev'] = aev_avg
                else:
                    temp_avg = round(df_query[feature].mean(),2)
                    data[feature] = temp_avg

            df_data_agg = pd.concat([df_data_agg, pd.DataFrame([data])], ignore_index=True)

        return df_data_agg

    def _aggregate_max(self):
        seperate_lys = self.df_measurements.drop_duplicates(subset=['Uniprot_Entry', 'Resid', 'class'])
        df_data_agg = pd.DataFrame(columns=['Uniprot_Entry', 'Resid', 'class'])

        self._reduce_aev_dimensions()

        for idx, row in seperate_lys.iterrows():
            entry = row['Uniprot_Entry']
            resid = row['Resid']
            class_val = row['class']
            # create a subset of the dataframe of measurements where the uniprot and resid match
            df_query = self.df_measurements[(self.df_measurements['Uniprot_Entry'] == entry) & (self.df_measurements['Resid'] == resid)]
            # calculate all the values of interest from the subset dataframe (df_query)
            data = {'Uniprot_Entry': entry,
                    'Resid' : resid,
                    'class' : class_val}
            df_query = df_query.copy()
            for feature in self.features_to_include:
                if feature == 'aev':
                    list_aevs = df_query['aev'].tolist()
                    aev_avg = np.max(list_aevs, axis=0)
                    data['aev'] = aev_avg
                else:
                    temp_max = round(df_query[feature].max(),2)
                    data[feature] = temp_max

            df_data_agg = pd.concat([df_data_agg, pd.DataFrame([data])], ignore_index=True)

        return df_data_agg

    def _aggregate_min(self):
        seperate_lys = self.df_measurements.drop_duplicates(subset=['Uniprot_Entry', 'Resid', 'class'])
        df_data_agg = pd.DataFrame(columns=['Uniprot_Entry', 'Resid', 'class'])

        self._reduce_aev_dimensions()

        for idx, row in seperate_lys.iterrows():
            entry = row['Uniprot_Entry']
            resid = row['Resid']
            class_val = row['class']
            # create a subset of the dataframe of measurements where the uniprot and resid match
            df_query = self.df_measurements[(self.df_measurements['Uniprot_Entry'] == entry) & (self.df_measurements['Resid'] == resid)]
            # calculate all the values of interest from the subset dataframe (df_query)
            data = {'Uniprot_Entry': entry,
                    'Resid' : resid,
                    'class' : class_val}
            df_query = df_query.copy()
            for feature in self.features_to_include:
                if feature == 'aev':
                    list_aevs = df_query['aev'].tolist()
                    aev_avg = np.min(list_aevs, axis=0)
                    data['aev'] = aev_avg
                else:
                    temp_min = round(df_query[feature].min(),2)
                    data[feature] = temp_min

            df_data_agg = pd.concat([df_data_agg, pd.DataFrame([data])], ignore_index=True)

        return df_data_agg

    def _aggregate_avg_less_avgaev(self):
        seperate_lys = self.df_measurements.drop_duplicates(subset=['Uniprot_Entry', 'Resid', 'class'])
        df_data_agg = pd.DataFrame(columns=['Uniprot_Entry', 'Resid', 'class'])
        for idx, row in seperate_lys.iterrows():
            entry = row['Uniprot_Entry']
            resid = row['Resid']
            class_val = row['class']
            # create a subset of the dataframe of measurements where the uniprot and resid match
            df_query = self.df_measurements[(self.df_measurements['Uniprot_Entry'] == entry) & (self.df_measurements['Resid'] == resid)]
            # calculate all the values of interest from the subset dataframe (df_query)
            data = {'Uniprot_Entry': entry,
                    'Resid' : resid,
                    'class' : class_val}
            df_query = df_query.copy()
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
            return df_data_agg

    def _aggregate_mixmatch(self):
        '''
        This one will do the best case for each individual feature
        Max sasa, das
        Min propka, pkaANI, depth
        Average AEV
        '''
        seperate_lys = self.df_measurements.drop_duplicates(subset=['Uniprot_Entry', 'Resid', 'class'])
        df_data_agg = pd.DataFrame(columns=['Uniprot_Entry', 'Resid', 'class'])
        max_features = ['sasa', 'das', 'frustration', 'seqcharge', ]
        min_features = ['propka', 'pkaANI', 'depth', 'density']

        self._reduce_aev_dimensions()

        for idx, row in seperate_lys.iterrows():
            entry = row['Uniprot_Entry']
            resid = row['Resid']
            class_val = row['class']
            # create a subset of the dataframe of measurements where the uniprot and resid match
            df_query = self.df_measurements[(self.df_measurements['Uniprot_Entry'] == entry) & (self.df_measurements['Resid'] == resid)]
            # calculate all the values of interest from the subset dataframe (df_query)
            data = {'Uniprot_Entry': entry,
                    'Resid' : resid,
                    'class' : class_val}
            df_query = df_query.copy()
            features = [a for a in self.df_measurements.columns if a not in ['Uniprot_Entry', 'PDB_Code', 'Resid', 'class']]
            for feature in features:
                if 'AEV_' in feature:
                    temp_avg = round(df_query[feature].mean(),2)
                    data[feature + '_avg'] = temp_avg
                else:
                    if feature in max_features:
                        temp_max = round(df_query[feature].max(),2)
                        data[feature + '_max'] = temp_max
                    elif feature in min_features:
                        temp_min = round(df_query[feature].min(),2)
                        data[feature + '_min'] = temp_min
                    else:
                        temp_max = round(df_query[feature].max(),2)
                        temp_min = round(df_query[feature].min(),2)
                        data[feature + '_max'] = temp_max
                        data[feature + '_min'] = temp_min

            df_data_agg = pd.concat([df_data_agg, pd.DataFrame([data])], ignore_index=True)

        self.features_to_include = [a for a in list(self.X_agg.columns) if a not in ['Uniprot_Entry', 'Resid', 'class']]
        return df_data_agg

    def _aggregate_random(self):
        seperate_lys = self.df_measurements.drop_duplicates(subset=['Uniprot_Entry', 'Resid', 'class'])
        df_data_agg = pd.DataFrame(columns=['Uniprot_Entry', 'Resid', 'class'])
        for idx, row in seperate_lys.iterrows():
            entry = row['Uniprot_Entry']
            resid = row['Resid']
            class_val = row['class']
            # create a subset of the dataframe of measurements where the uniprot and resid match
            df_query = self.df_measurements[(self.df_measurements['Uniprot_Entry'] == entry) & (self.df_measurements['Resid'] == resid)]
            # calculate all the values of interest from the subset dataframe (df_query)
            data = {'Uniprot_Entry': entry,
                    'Resid' : resid,
                    'class' : class_val}
            df_query = df_query.copy()
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
            return df_data_agg


    def _aggregate_minmax(self):
        '''
        Aggregation method which will take the max and min values of the features and duplicate
        the features to allow for this without having problems with the data overlapping and
        biasing the output.
        '''
        seperate_lys = self.df_measurements.drop_duplicates(subset=['Uniprot_Entry', 'Resid', 'class'])
        df_data_agg = pd.DataFrame(columns=['Uniprot_Entry', 'Resid', 'class'])

        self._reduce_aev_dimensions()

        for idx, row in seperate_lys.iterrows():
            entry = row['Uniprot_Entry']
            resid = row['Resid']
            class_val = row['class']
            # create a subset of the dataframe of measurements where the uniprot and resid match
            df_query = self.df_measurements[(self.df_measurements['Uniprot_Entry'] == entry) & (self.df_measurements['Resid'] == resid)]
            # calculate all the values of interest from the subset dataframe (df_query)
            data = {'Uniprot_Entry': entry,
                    'Resid' : resid,
                    'class' : class_val}
            df_query = df_query.copy()
            features = [a for a in self.features_to_include if a not in ['Uniprot_Entry', 'PDB_Code', 'Resid', 'class']]
            for feature in features:
                if feature == 'aev':
                    df_query['sumaev'] = df_query[[a for a in df_query.columns if 'AEV_' in a]].sum(axis=1)
                    min_row = df_query['sumaev'].idxmin()
                    max_row = df_query['sumaev'].idxmax()
                    for feat in [a for a in df_query.columns if 'AEV_' in a]:
                        data[feat + '_min'] = round(df_query[feat].loc[min_row], 2)
                        data[feat + '_max'] = round(df_query[feat].loc[max_row], 2)

                else:
                    temp_min = round(df_query[feature].min(),2)
                    data[feature + '_min'] = temp_min
                    temp_max = round(df_query[feature].max(),2)
                    data[feature + '_max'] = temp_max

            df_data_agg = pd.concat([df_data_agg, pd.DataFrame([data])], ignore_index=True)


        self.features_to_include = [a for a in list(df_data_agg.columns) if a not in ['Uniprot_Entry', 'Resid', 'class']]
        return df_data_agg


    def _aggregate_minmaxavg(self):
        '''
        Aggregation method which will take the max and min values of the features and duplicate
        the features to allow for this without having problems with the data overlapping and
        biasing the output.
        '''
        print('MinMax aggregation type only takes the min for aev, works for all other features: GW FIX ME!')
        seperate_lys = self.df_measurements.drop_duplicates(subset=['Uniprot_Entry', 'Resid', 'class'])
        df_data_agg = pd.DataFrame(columns=['Uniprot_Entry', 'Resid', 'class'])

        self._reduce_aev_dimensions()

        for idx, row in seperate_lys.iterrows():
            entry = row['Uniprot_Entry']
            resid = row['Resid']
            class_val = row['class']
            # create a subset of the dataframe of measurements where the uniprot and resid match
            df_query = self.df_measurements[(self.df_measurements['Uniprot_Entry'] == entry) & (self.df_measurements['Resid'] == resid)]
            # calculate all the values of interest from the subset dataframe (df_query)
            data = {'Uniprot_Entry': entry,
                    'Resid' : resid,
                    'class' : class_val}
            df_query = df_query.copy()
            features = [a for a in self.df_measurements.columns if a not in ['Uniprot_Entry', 'PDB_Code', 'Resid', 'class']]
            for feature in features:
                temp_min = round(df_query[feature].min(),2)
                temp_max = round(df_query[feature].max(),2)
                temp_avg = round(df_query[feature].mean(),2)
                data[feature + '_min'] = temp_min
                data[feature + '_max'] = temp_max
                data[feature + '_avg'] = temp_avg

            df_data_agg = pd.concat([df_data_agg, pd.DataFrame([data])], ignore_index=True)

        self.features_to_include = [a for a in list(self.X_agg.columns) if a not in ['Uniprot_Entry', 'Resid', 'class']]
        return df_data_agg
