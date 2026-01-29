import os
import random
import statistics
import random
from ast import literal_eval
from copy import deepcopy as dc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
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
                        the feature space allowing the model to use both as needed (default)
            - 'minmaxavg': Takes the min, max and average of each feature, tripling feature space
            - 'all': Takes all potential statistical features that have been coded to be calculated
            - 'choose': Allows the user to choose which statistic for the feature they want
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
        # Note: current preference for using aev_legolas as easier to obtain - change here if necessary
        self.df_measurements = df_measurements
        self.aggregation_method = aggregation_method
        self.features_to_include = features_to_include
        self.aev_red_method = aev_red_method
        self.num_sd_aev_features = num_sd_aev_features

        if 'aev_legolas' in self.df_measurements.columns:
            if 'aev' in self.df_measurements.columns:
                self.df_measurements.drop(columns=['aev'], inplace=True)
            self.df_measurements.rename(columns={'aev_legolas': 'aev'}, inplace=True)

        if self.features_to_include == ['all']:
            self.df_measurements = self.df_measurements.loc[:, ~self.df_measurements.columns.str.contains('^Unnamed')]
            self.features_to_include = [a for a in self.df_measurements.columns if a not in ['Uniprot_Entry', 'PDB_Code', 'Method', 'Resolution', 'Chain', 'Resid', 'class']]

        data_cols_entered = self.df_measurements.columns.values
        cols_required = ['Uniprot_Entry', 'PDB_Code', 'Resid', 'class']
        for feat in self.features_to_include:
            if feat in data_cols_entered:
                cols_required.append(feat)
            else:
                print(f'Feature given as input not available in all input files, will not be included: {feat}')
        self.df_measurements = self.df_measurements[cols_required]

        print(self.df_measurements.columns)
        if 'method' in self.df_measurements.columns: self.df_measurements = self.df_measurements.drop(columns='Method')
        if 'Resolution' in self.df_measurements.columns: self.df_measurements = self.df_measurements.drop(columns='Resolution')
        print(len(self.df_measurements))
        self.df_measurements = self.df_measurements.dropna(subset=self.features_to_include)  # TODO GW 29.01.26 - add somethign to let you know how mnay lines have been removed and if many of them are from one specific feature
        print(self.df_measurements.columns)
        print(len(self.df_measurements))

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
        self._data_tidying()
        df_stats = self._calculate_statistics()
        match self.aggregation_method:
            case 'avg':
                return df_stats.drop(columns=[a for a in [b for b in df_stats.columns if 'avg' not in b] if a not in ['Uniprot_Entry', 'Resid', 'class']])
            case 'random':
                return df_stats.drop(columns=[a for a in [b for b in df_stats.columns if 'rand' not in b] if a not in ['Uniprot_Entry', 'Resid', 'class']])
            case 'max':
                return df_stats.drop(columns=[a for a in [b for b in df_stats.columns if 'max' not in b] if a not in ['Uniprot_Entry', 'Resid', 'class']])
            case 'min':
                return df_stats.drop(columns=[a for a in [b for b in df_stats.columns if 'min' not in b] if a not in ['Uniprot_Entry', 'Resid', 'class']])
            case 'average subtract aev':
                return self._aggregate_avg_less_avgaev()
            case 'mixmatch':
                max_features = ['sasa', 'das', 'frustration', 'seqcharge']
                min_features = ['propka', 'pkaANI', 'depth', 'density', 'legolas']
                max_feat_cols = [a for a in self.features_to_include if any(b in a for b in max_features) and 'max' in a]
                min_feat_cols = [a for a in self.features_to_include if any(b in a for b in min_features) and 'min' in a]
                avg_features = [a for a in [b for b in self.features_to_include if not any (c in b for c in max_features) and not any (c in b for c in min_features)] if 'avg' in a]
                return df_stats.drop(columns=[a for a in self.features_to_include if a not in max_feat_cols + min_feat_cols + avg_features])
            case 'minmax':
                return df_stats.drop(columns=[a for a in [b for b in df_stats.columns if not any(c in b for c in ['min', 'max'])] if a not in ['Uniprot_Entry', 'Resid', 'class']])
            case 'minmaxavg':
                return df_stats.drop(columns=[a for a in [b for b in df_stats.columns if not any(c in b for c in ['min', 'max', 'avg'])] if a not in ['Uniprot_Entry', 'Resid', 'class']])
            case 'all':
                return df_stats
            case 'choose':
                return self._aggregate_choose(df_stats)
            case _:
                print('Aggregation method not recognised; using minmax values')
                return df_stats.drop(a for a in df_stats.columns if a not in ['Uniprot_Entry', 'Resid', 'class'] or any(b in a for b in ['min', 'max']))

    def _data_tidying(self):
        '''
        Go over the provided data and remove any poor data from this. 
        '''
        if 'depth' in self.df_measurements.columns:
            self._remove_bad_data('depth', 0, 20)

    def _remove_bad_data(self, feature, lower, upper):
        '''
        Easy function for removing bad rows of data from the aggregated data

        Parameters
        ----------
        feature -> string
            The feature to investigate bad values for
        lower -> float
            The lower bound of bad values to accept
        upper -> float
            The upper bound of bad values to accept
        '''
        df_suspicious = self.df_measurements[(self.df_measurements[feature] < lower) | (self.df_measurements[feature] > upper)]
        len_one = len(self.df_measurements)
        remove_list = df_suspicious.index.tolist()
        self.df_measurements = self.df_measurements.drop(index = remove_list)
        print(f'>> Removed {len(df_suspicious)} rows from the measurements data, new length {len(self.df_measurements)} (old length: {len_one} rows)')


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


    def _calculate_statistics(self):
        '''
        Function for creating a dataframe which includes all the potential statistics which
        could then be used for aggregation later on. This can then be shortened as desired
        based on which method of aggregation is required for this.
        '''
        print('>> Calculating statistics for measurements data provided...')
        if 'class' not in self.df_measurements.columns:
            self.df_measurements['class'] = -1  # set to -1 as unsure if pos or neg

        seperate_lys = self.df_measurements.drop_duplicates(subset=['Uniprot_Entry', 'Resid', 'class'])
        df_stats = pd.DataFrame(columns=['Uniprot_Entry', 'Resid', 'class'])

        self._reduce_aev_dimensions()

        for idx, row in seperate_lys.iterrows():
            entry = row['Uniprot_Entry']
            resid = row['Resid']
            class_val = row['class']
            df_query = self.df_measurements[(self.df_measurements['Uniprot_Entry'] == entry) & (self.df_measurements['Resid'] == resid)].copy().reset_index()
            data = {'Uniprot_Entry': entry,
                    'Resid' : resid,
                    'class' : class_val}
            features = [a for a in self.features_to_include if a not in ['Uniprot_Entry', 'PDB_Code', 'Resid']]
            for feature in features:
                if feature == 'aev' or feature == 'aev_legolas':
                    df_query['sumaev'] = df_query[[a for a in df_query.columns if 'AEV_' in a]].sum(axis=1)
                    min_row = df_query['sumaev'].idxmin()
                    max_row = df_query['sumaev'].idxmax()
                    rand_row = random.randrange(0, len(df_query))
                    for feat in [a for a in df_query.columns if 'AEV_' in a]:

                        data[feat + '_min'] = round(float(df_query[feat].loc[min_row]), 2)  # min value at this position in min AEV
                        data[feat + '_max'] = round(float(df_query[feat].loc[max_row]), 2)  # max value at this position in max AEV
                        data[feat + '_rand'] = round(df_query[feat].loc[rand_row], 2)

                        data[feat + '_absmin'] = round(df_query[feat].min(), 2)  # min value of any AEV at this position in the AEV
                        data[feat + '_absmax'] = round(df_query[feat].max(), 2)  # max value of any AEV at this position in the AEV
                        data[feat + '_avg'] = round(df_query[feat].mean(), 2)
                        data[feat + '_sd'] = round(df_query[feat].std(), 2)
                        data[feat + '_range'] = round(df_query[feat].max(), 2) - round(df_query[feat].min(), 2)

                else:
                    data[feature + '_min'] = round(float(df_query[feature].min()), 2)
                    data[feature + '_max'] = round(float(df_query[feature].max()), 2)
                    data[feature + '_avg'] = round(df_query[feature].mean(),2)
                    data[feature + '_sd'] = round(df_query[feature].std(),2)
                    data[feature + '_range'] = data[feature + '_max'] - data[feature + '_min']
                    data[feature + '_rand'] = df_query[feature].iloc[random.randrange(0, len(df_query))]

            df_stats = pd.concat([df_stats, pd.DataFrame([data])], ignore_index=True)

        self.features_to_include = [a for a in list(df_stats.columns) if a not in ['Uniprot_Entry', 'Resid', 'class']]
        return df_stats

    def _aggregate_choose(self, df_stats):
        '''
        Choose method of aggregation allowing the user to choose the statistics they require 
        '''
        max_features = []; min_features = []; avg_features = []; sd_features = []; range_features = []; rand_features = []
        print('>> For the following list of features, please choosen which statistics you would like to include in the aggregation')
        print('>> The following statistics can be chosen (enter word or number seperated by ;): all (1), max (2), min (3), avg (4), sd (5), range (6), random (7)')
        for feat in list(set([a.split('_')[0] for a in self.features_to_include if 'max' in a])):
            tmp_stats = input(f'>> Choose statistics for feature: {feat}')
            try:
                stats = [a.strip() for a in tmp_stats.split(';')]
                for stat in stats:
                    match stat:
                        case 'all' | '1':
                            if feat not in max_features: max_features.append(feat)
                            if feat not in min_features: min_features.append(feat)
                            if feat not in avg_features: avg_features.append(feat)
                            if feat not in sd_features: sd_features.append(feat)
                            if feat not in range_features: range_features.append(feat)
                            if feat not in rand_features: rand_features.append(feat)
                        case 'max' | '2':
                            if feat not in max_features: max_features.append(feat)
                        case 'min' | '3':
                            if feat not in min_features: min_features.append(feat)
                        case 'avg' | '4':
                            if feat not in avg_features: avg_features.append(feat)
                        case 'sd' | '5':
                            if feat not in sd_features: sd_features.append(feat)
                        case 'range' | '6':
                            if feat not in range_features: range_features.append(feat)
                        case 'random' | '7':
                            if feat not in rand_features: rand_features.append(feat)
                        case _:
                            print(f'>> Could not evaluate input ({tmp_stats}) for feature: {feat}, using all metrics for {feat} instead')
                            if feat not in max_features: max_features.append(feat)
                            if feat not in min_features: min_features.append(feat)
                            if feat not in avg_features: avg_features.append(feat)
                            if feat not in sd_features: sd_features.append(feat)
                            if feat not in range_features: range_features.append(feat)
                            if feat not in rand_features: rand_features.append(feat)
            except Exception as e:
                print(f'>> Could not evaluate input ({tmp_stats}) for feature: {feat} with error {e}, using all metrics for {feat} instead')
        max_feat_cols = [a for a in self.features_to_include if any(b in a for b in max_features) and 'max' in a]
        min_feat_cols = [a for a in self.features_to_include if any(b in a for b in min_features) and 'min' in a]
        avg_feat_cols = [a for a in self.features_to_include if any(b in a for b in avg_features) and 'avg' in a]
        sd_feat_cols = [a for a in self.features_to_include if any(b in a for b in sd_features) and 'sd' in a]
        range_feat_cols = [a for a in self.features_to_include if any(b in a for b in range_features) and 'range' in a]
        rand_feat_cols = [a for a in self.features_to_include if any(b in a for b in rand_features) and 'rand' in a]
        return df_stats.drop(columns=[a for a in self.features_to_include if a not in max_feat_cols + min_feat_cols + avg_feat_cols + sd_feat_cols + range_feat_cols + rand_feat_cols])

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



if __name__ == "__main__":
    
    test_dataframe_name = 'data/measures_cut_Ecoli(hCit)_all_01.05.25_joined.csv'
    test_measures_dataframe = pd.read_csv(test_dataframe_name)
    #print(test_measures_dataframe.head())
    print(len(test_measures_dataframe.dropna().drop_duplicates(subset=['Uniprot_Entry', 'Resid'])))
    agg = Aggregation(test_measures_dataframe, aggregation_method='minmax', features_to_include=['all'], aev_red_method='pca')
    test_agg_df = agg.aggregate_data()
    print(test_agg_df)
    #test_agg_df.to_csv('AllNegative_aggregated_minmaxavg.csv')