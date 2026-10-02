import os
import numbers
import random
from ast import literal_eval
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from .preprocessing import Preprocessing


class Aggregation:
    '''
    Class to handle the different aggregation methods for sorting over the different measurements
    that have been calculated. This will take a dataframe of measurements, and return a dataframe
    with only one set of measurements per lysine residue according to the aggregation method.
    '''

    def __init__(self, df_measurements, outdir='result', aggregation_method='minmax',
                 features_to_include=['all'], aev_red_method='pca',
                 num_sd_aev_features=100, include_chain=False,
                 max_feature_nan_fraction=0.5, get_nan_df=False,
                 aev_pca_variance=0.99, vif_threshold=5.0):
        '''
        Initialisation of the Aggregation class.

        :param df_measurements: The dataframe of measurements which need to be aggregated for use
        :type df_measurements: pandas.DataFrame
        :param outdir: Name of the directory to write any data to
        :type outdir: str
        :param aggregation_method: The aggregation method chosen to reduce the measurements into a
            usable format for training on. Options are:

            - 'avg': Take the average of all measurements for each lysine
            - 'random': Take a random measurement out of all measurements for the lysine
            - 'max': Take the maximum value of each feature for the lysine
            - 'min': Take the minimum value of each feature for the lysine
            - 'median': Take the median value of each feature for the lysine
            - 'mixmatch': Takes the predicted metric which will work best for each feature. The
              max is used for features where a high value is likely to be important and min for
              features where less of it is required.
            - 'minmax': Takes both the min and max values of each feature which duplicates the
              feature space allowing the model to use both as needed (default)
            - 'minmaxavg': Takes the min, max and average of each feature, tripling feature
              space
            - 'all': Takes all potential statistical features that have been coded to be
              calculated
            - custom: Allows the user to choose which statistic for the feature they want. You
                will be prompted on the command line to specify options for each feature passed.
        :type aggregation_method: str
        :param features_to_include: The list of features that are to be included in the aggregation.
            The default for this is taken to be all of them.
        :type features_to_include: list, optional
        :param aev_red_method: The dimensionality reduction method for reducing the size of the
            AEVs, reducing clouding. PCA is taken as default for this. Options:

            - 'PCA': Create a PCA which represents a fraction aev_pca_variance of the variance
              of the data (99% by default) and use these new vectors to represent the AEVs
              instead.
            - 'sd': Take the standard deviation of all the AEV columns and work out which the
              top n are taken through for use
            - 'vif': Uses variance inflation factors to calculate the decorrelation between the
              different columns of the AEV. Only takes through the features whose variance
              inflation factor falls below vif_threshold.
            - 'autoencoder': Uses an autoencoder to reduce the dimensions, better for non-linear
              data.
            - 'null': Removes all the null columns from the AEV
        :type aev_red_method: str, optional
        :param num_sd_aev_features: The number of features to keep from the aevs when the standard
            deviation method is used. Default is set to 100.
        :type num_sd_aev_features: int, optional
        :param include_chain: Option to include chain in the aggregation key, if True it will
            aggregate on 'Uniprot_Entry', 'Chain', 'Resid' else will aggregate on 'Uniprot_Entry',
            'Resid' (default)
        :type include_chain: bool, optional
        :param max_feature_nan_fraction: The max fraction of values in a feature column that would
            allow a feature to remain in the measures dataframe for aggregation. If set to 1.0 then
            the feature will be kept and all rows containing NaN will be removed before aggregation.
        :type max_feature_nan_fraction: float
        :param get_nan_df: Option to create a dataframe (saved as csv) which contains all the rows
            that are being removed when aggregating, this allows curation of the data being removed
            for investigations into potential problems.
        :type get_nan_df: bool, optional
        :param aev_pca_variance: Used when aev_red_method is 'pca'. A fraction between 0 and 1
            keeps the fewest principal components that explain that fraction of the variance
            of the AEVs; an integer keeps that number of components. Default 0.99.
        :type aev_pca_variance: float or int, optional
        :param vif_threshold: Used when aev_red_method is 'vif'. AEV columns are removed until
            every remaining one has a variance inflation factor below this value. Default 5.
        :type vif_threshold: float, optional
        '''
        if isinstance(aev_pca_variance, bool) or not (
                (isinstance(aev_pca_variance, numbers.Integral) and aev_pca_variance >= 1)
                or (isinstance(aev_pca_variance, numbers.Real)
                    and not isinstance(aev_pca_variance, numbers.Integral)
                    and 0.0 < aev_pca_variance < 1.0)):
            raise ValueError(f'aev_pca_variance must be a fraction between 0 and 1, or a '
                             f'number of components of at least 1, got {aev_pca_variance}')
        if not vif_threshold > 1:
            raise ValueError(f'vif_threshold must be greater than 1, the smallest possible '
                             f'variance inflation factor, got {vif_threshold}')
        self.aev_pca_variance = aev_pca_variance
        self.vif_threshold = vif_threshold
        if isinstance(df_measurements, str):
            self.df_measurements = pd.read_csv(df_measurements)
        else:
            self.df_measurements = df_measurements

        if isinstance(features_to_include, str):
            self.features_to_include = [features_to_include]
        else:
            self.features_to_include = features_to_include

        self.outdir = outdir
        os.makedirs(outdir, exist_ok=True)
        self.aggregation_method = aggregation_method
        self.aev_red_method = aev_red_method
        self.num_sd_aev_features = num_sd_aev_features
        self.include_chain = include_chain
        self.get_nan_df = get_nan_df
        self.df_agg = pd.DataFrame()
        self.max_feature_nan_fraction = max_feature_nan_fraction

        self.non_feature_cols = ['Uniprot_Entry', 'PDB_Code', 'Chain', 'Modified', 'Method', 'Source',
                                 'Resolution', 'Resid', 'class', 'PLDDT', 'Largest_Gap']

        if self.features_to_include == ['all']:
            self.df_measurements = self.df_measurements.loc[:, ~self.df_measurements.columns.str.contains('^Unnamed')]
            self.features_to_include = [a for a in self.df_measurements.columns if a not in self.non_feature_cols]

        if 'Uniprot_Entry' in self.df_measurements.columns:
            self.lys_key = ['Uniprot_Entry', 'Resid']
            if self.include_chain:
                self.lys_key.append('Chain')
        elif 'PDB_Code' in self.df_measurements.columns:
            self.lys_key = ['PDB_Code', 'Resid']
            if self.include_chain:
                self.lys_key.append('Chain')
        else:
            return KeyError('>> Measurements dataframe must contain either a Uniprot_Entry '
                            'or PDB_Code column in measurements file')

        data_cols_entered = self.df_measurements.columns.values
        cols_required = [a for a in self.df_measurements.columns if a in self.lys_key + ['PDB_Code', 'class']]
        cols_removed = [f for f in self.features_to_include if f not in data_cols_entered]
        self.features_to_include = [f for f in self.features_to_include if f in data_cols_entered]
        if cols_removed:
            print(f'The following features given as input not available in all '
                  f'input files, will not be included: {cols_removed}')
        self.df_measurements = self.df_measurements[cols_required + self.features_to_include]

        self.meta_data_cols = [a for a in ['Method', 'Resolution'] if a in self.df_measurements.columns]
        self.features_to_include = [a for a in self.features_to_include if a not in self.meta_data_cols]

        print('>> Finding numbers of na values present in each feature column in the dataframe')
        always_na_cols = []
        for col in self.features_to_include:
            tmp_na_vals = self.df_measurements[col].isna().sum()
            print(f'Number of na values for feature: {col}: num na: {tmp_na_vals}')
            if len(self.df_measurements) > 0 and tmp_na_vals == len(self.df_measurements):
                always_na_cols.append(col)
        if always_na_cols:
            print(f'>> Columns {", ".join(always_na_cols)} are always na and provide '
                  f'no information, removing from dataframe')
            self.df_measurements = self.df_measurements.drop(columns=always_na_cols, axis=1)
            self.features_to_include = [a for a in self.features_to_include if a not in always_na_cols]

        len_before_df = len(self.df_measurements)

        if self.get_nan_df:
            df_nan = self.df_measurements[self.df_measurements[self.features_to_include].isna().any(axis=1)]
            df_nan.to_csv(f'{self.outdir}{os.sep}measures_nan_feature_data_removed.csv', index=False)
            print(f'>> Wrote {len(df_nan)} rows carrying a NaN of data which have been removed from the '
                  f'aggregation data to \'measures_nan_feature_data_removed.csv\'')

        na_per_feature = {c: int(self.df_measurements[c].isna().sum()) for c in self.features_to_include}
        n_meas = len(self.df_measurements)
        sparse_features = [c for c, k in na_per_feature.items() if n_meas and
                                k / n_meas > self.max_feature_nan_fraction]
        if sparse_features:
            print(f'>> Features {", ".join(sparse_features)} are missing values on more than '
                  f'{self.max_feature_nan_fraction:.0%} of rows; dropping the features instead '
                  f'of the rows. If you would like to keep the feature and get rid of the rows '
                  f'pass 1.0 as the value for max_feature_nan_fraction.')
            self.features_to_include = [feat for feat in self.features_to_include
                                        if feat not in sparse_features]

        keep = self.df_measurements.dropna(subset=self.features_to_include)
        del_lysines = (self.df_measurements[self.lys_key].drop_duplicates().shape[0] - keep[self.lys_key].drop_duplicates().shape[0])
        if len_before_df and keep.empty:
            raise ValueError('>> Every row was removed by the NaN filter. The features with the most missing '
                             'values are: ' + ', '.join(f'{c} ({n})' for c, n in
                            sorted(((c, self.df_measurements[c].isna().sum()) for c in self.features_to_include), key=lambda t: -t[1])[:5]))
        else:
            print(f'>> Removed {len_before_df - len(keep)} rows of {len_before_df} possible measurement '
                  f'rows containing nan values; which cost {del_lysines} distinct residues. The worst '
                  f'features for this are: ' + ', '.join(f'{c} ({n})' for c, n in sorted(na_per_feature.items(), key=lambda t: -t[1])[:3]))
            self.df_measurements = keep


    def aggregate_data(self):
        '''
        Match aggregation type up to the relevant aggregation function. Calls
        _calculate_statistics() to get the dataframe of statistics data on each residue key.
        Uses the input parameter of aggregation_method and matches this to a case, this case
        the drops any columns in the dataframe which aren't related to the aggregation
        method.

        :returns: Dataframe with aggregated data accoriding to the aggregation method
        :rtype: pandas.DataFrame
        '''

        bad_feature_sets = [('depth', 0, 20)]  # add to as more confinements on features needed
        for bad_feat, feat_low, feat_up in bad_feature_sets:
            if bad_feat in self.df_measurements.columns:
                self._remove_bad_data(feature=bad_feat, lower=feat_low, upper=feat_up)
        df_stats = self._calculate_statistics()
        not_cols = self.lys_key + ['class']
        match self.aggregation_method:
            case 'avg':
                self.df_agg = df_stats.drop(columns=[a for a in [b for b in df_stats.columns if 'avg' not in b] if a not in not_cols])
            case 'random':
                self.df_agg = df_stats.drop(columns=[a for a in [b for b in df_stats.columns if 'rand' not in b] if a not in not_cols])
            case 'max':
                self.df_agg = df_stats.drop(columns=[a for a in [b for b in df_stats.columns if 'max' not in b] if a not in not_cols])
            case 'min':
                self.df_agg = df_stats.drop(columns=[a for a in [b for b in df_stats.columns if 'min' not in b] if a not in not_cols])
            case 'median':
                self.df_agg = df_stats.drop(columns=[a for a in [b for b in df_stats.columns if 'med' not in b] if a not in not_cols])
            case 'mixmatch':
                max_features = ['sasa', 'das', 'frustration', 'seqcharge']
                min_features = ['propka', 'pkaANI', 'depth', 'density', 'legolas']
                max_feat_cols = [a for a in self.features_to_include if any(b in a for b in max_features) and 'max' in a]
                min_feat_cols = [a for a in self.features_to_include if any(b in a for b in min_features) and 'min' in a]
                avg_features = [a for a in [b for b in self.features_to_include if not any (c in b for c in max_features) and not any (c in b for c in min_features)] if 'avg' in a]
                self.df_agg = df_stats.drop(columns=[a for a in self.features_to_include if a not in max_feat_cols + min_feat_cols + avg_features])
            case 'minmax':
                self.df_agg = df_stats.drop(columns=[a for a in [b for b in df_stats.columns if not any(c in b for c in ['min', 'max'])] if a not in not_cols])
            case 'minmaxavg':
                self.df_agg = df_stats.drop(columns=[a for a in [b for b in df_stats.columns if not any(c in b for c in ['min', 'max', 'avg'])] if a not in not_cols])
            case 'all':
                self.df_agg = df_stats
            case 'choose':
                self.df_agg = self._aggregate_choose(df_stats)
            case _:
                print('Aggregation method not recognised; using minmax values')
                self.aggregation_method = 'minmax'
                self.df_agg = df_stats.drop(columns=[a for a in [b for b in df_stats.columns if not any(c in b for c in ['min', 'max'])] if a not in not_cols])
        return self.df_agg


    def _remove_bad_data(self, feature, lower, upper):
        '''
        Easy function for removing bad rows of data from the aggregated data

        :param feature: The feature to investigate bad values for
        :type feature: str
        :param lower: The lower bound of bad values to accept
        :type lower: float
        :param upper: The upper bound of bad values to accept
        :type upper: float
        '''
        df_suspicious = self.df_measurements[(self.df_measurements[feature] < lower) | (self.df_measurements[feature] > upper)]
        #self.df_measurements[feature] = self.df_measurements[feature].where(self.df_measurements[feature].between(lower, upper))  # sets out of bounds entries to NaN
        self.df_measurements = self.df_measurements[self.df_measurements[feature].between(lower, upper)]
        #print(f'>> Changed {len(df_suspicious)} values for {feature} from the measurements data which did not fall inside the bounds, replaced with NaN')
        print(f'>> Removed {len(df_suspicious)} rows which contain out of bounds entries for {feature} from the measurements data.')


    def _reduce_aev_dimensions(self):
        '''
        Due to curse of dimensionality the model performs worse when the AEVs are clouding the data
        as it cannot work out which features are actually important. Therefore, this function will
        call the required method to reduce the AEVs down to a specified number of features depending
        on the method chosen. Method options:

        - 'pca': Principal Component Analysis taking 99% of the variance within the data
        - 'sd': Standard Deviation of the AEV
        - 'vif': Variance Inflation Factor Correlation analysis to remove features which are
          correlated
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
                case 'null':
                    self._cut_null_aev_columns()
                case _:
                    print(f'>> AEV dimensionality reduction method: {self.aev_red_method}, '
                          f'was not recognised, using PCA method.')
                    self._prepare_pca()
        if 'evolution' in self.features_to_include:
            df_evolution = pd.DataFrame(list([literal_eval(evol) for evol in self.df_measurements['evolution']]))
            df_evolution = df_evolution.add_prefix('EVL_')
            self.df_measurements = pd.concat([self.df_measurements.reset_index(), df_evolution.reset_index()], axis=1)
            self.df_measurements = self.df_measurements.drop(['index', 'evolution'], axis=1)



    def _prepare_aev_sd(self):
        '''
        Model performs worse when more of the features of the AEV are taken through to training.
        Reduce the AEVs down to the required number of features based on one of the methods chosen
        below. The first method is to use the standard deviation of the individual features of the
        AEV to work out which features show variation and will be likely to be good choices to take
        through to the model. This is currently setup to find the top 100 from the base aevs fed in.
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
        Function takes all the AEV columns and removes any that are always null which would add
        nothing to the model except noise.
        '''
        print('>> Removing null AEV columns...')
        aev_col_names = [a for a in list(self.df_measurements.columns) if 'AEV_' in a]
        cols_to_remove = [a for a in aev_col_names if (self.df_measurements[a] == 0).all()]
        self.df_measurements = self.df_measurements.drop(columns=cols_to_remove, axis=0)
        print(f'>> Removed {len(cols_to_remove)} null AEV columns, '
              f'{len(aev_col_names) - len(cols_to_remove)} AEV columns left.')


    def _prepare_vif_aev(self):
        '''
        Use the preprocessing module to calculate the variance inflation factor values for each of
        the columns within the AEVs and remove the columns which are highly correlated together such
        that it is the minimum number of columns without correlation. Non-correlation was taken to
        be a VIF value of less than vif_threshold (5 by default).
        '''
        print('>> Reducing AEV dimensions using VIF analysis...')
        aev_cols = [a for a in self.df_measurements.columns if 'AEV_' in a]
        P = Preprocessing(self.df_measurements, aev_cols)
        columns_to_keep = P.calculate_diff_features(self.df_measurements,
                                                    vif_threshold=self.vif_threshold)

        cols_to_remove = [a for a in self.df_measurements.columns if ('AEV_' in a) and (a not in columns_to_keep)]
        self.df_measurements.drop(cols_to_remove, axis=1, inplace=True)
        print(f'>> VIF analysis reduced AEV dimensions from {len(aev_cols)} to {len(columns_to_keep)}')


    def _prepare_pca(self):
        '''
        Short function to transform the AEV data using PCA to reduce the dimensions. Changes the
        data in self.X_all ready for the aggregation to actually reduce the dataset for training on
        '''
        print('>> Calculating PCA on AEV data...')
        n_components = self.aev_pca_variance
        pca = PCA(n_components=n_components)
        aev_col_names = [a for a in self.df_measurements.columns if 'AEV_' in a]
        aev_data = self.df_measurements[aev_col_names].values.tolist()
        pca_data = pca.fit_transform(aev_data)

        variance_data = pca.explained_variance_ratio_

        self.df_measurements.drop(aev_col_names, axis=1, inplace=True)
        df_out = pd.DataFrame(pca_data)
        df_out = df_out.add_prefix('AEV_')
        self.df_measurements = pd.concat([self.df_measurements, df_out], axis=1)
        print(f'>> PCA used {pca.n_components_} components to explain '
              f'{variance_data.sum():.4f} of the variance. AEVs are now in reduced dimension '
              f'format.')


    def _calculate_statistics(self):
        '''
        Function for creating a dataframe which includes all the potential statistics which could
        then be used for aggregation later on. This can then be shortened as desired based on which
        method of aggregation is required for this.
        AEVs following the handling that was given in the class creation. Evolution vectors are
        treated as almost identical due to being a sequence effect, therefore take AF file which
        uses canonical sequence as the one for the aggregated dataframe.
        '''
        print('>> Calculating statistics for measurements data provided...')
        if 'class' not in self.df_measurements.columns:
            self.df_measurements['class'] = -1  # set to -1 as unsure if pos or neg

        cols_key = self.lys_key + ['class']
        df_stats = pd.DataFrame(columns=cols_key)

        if 'aev' in self.features_to_include or 'evolution' in self.features_to_include:
            self._reduce_aev_dimensions()

        for row_key, df_query in self.df_measurements.groupby(self.lys_key):
            if not isinstance(row_key, tuple):
                row_key = (row_key,)
            df_query = df_query.reset_index(drop=True)
            if len(set(df_query['class'])) != 1:
                print(f'>> Not all instances of resid assigned to same class (instances: '
                      f'{list(set(df_query["class"]))}), using class -1 instead')
                class_val = -1
            else:
                class_val = df_query['class'].iloc[0]  # take first value of class as overall class for resid
            data = dict(zip(self.lys_key, row_key))
            data['class'] = class_val
            features = [a for a in self.features_to_include if a not in self.non_feature_cols]
            rand_row = random.randrange(0, len(df_query))
            for feature in features:
                if feature == 'aev':
                    df_query['sumaev'] = df_query[[a for a in df_query.columns if 'AEV_' in a]].sum(axis=1)
                    min_row = df_query['sumaev'].idxmin()
                    max_row = df_query['sumaev'].idxmax()
                    for feat in [a for a in df_query.columns if 'AEV_' in a]:

                        data[feat + '_min'] = round(float(df_query[feat].loc[min_row]), 2)  # min value at this position in min AEV
                        data[feat + '_max'] = round(float(df_query[feat].loc[max_row]), 2)  # max value at this position in max AEV
                        data[feat + '_rand'] = round(df_query[feat].loc[rand_row], 2)

                        data[feat + '_absmin'] = round(df_query[feat].min(), 2)  # min value of any AEV at this position in the AEV
                        data[feat + '_absmax'] = round(df_query[feat].max(), 2)  # max value of any AEV at this position in the AEV
                        data[feat + '_med'] = round(df_query[feat].median(), 2)
                        data[feat + '_avg'] = round(df_query[feat].mean(), 2)
                        data[feat + '_sd'] = round(df_query[feat].std(ddof=0), 2)
                        data[feat + '_range'] = round(df_query[feat].max(), 2) - round(df_query[feat].min(), 2)
                elif feature == 'evolution':
                    pred_query = df_query[df_query.get('Method', pd.Series(dtype=object)) == 'Predicted']
                    if not pred_query.empty:
                        for feat in [a for a in df_query.columns if 'EVL_' in a]:
                            data[feat] = pred_query[feat].iloc[0]
                    else:
                        for feat in [a for a in df_query.columns if 'EVL_' in a]:
                            data[feat] = df_query[feat].iloc[0]

                else:
                    data[feature + '_min'] = round(float(df_query[feature].min()), 2)
                    data[feature + '_max'] = round(float(df_query[feature].max()), 2)
                    data[feature + '_med'] = round(df_query[feature].median(), 2)
                    data[feature + '_avg'] = round(df_query[feature].mean(),2)
                    data[feature + '_sd'] = round(df_query[feature].std(ddof=0),2)
                    data[feature + '_range'] = data[feature + '_max'] - data[feature + '_min']
                    data[feature + '_rand'] = df_query[feature].iloc[rand_row]

            df_stats = pd.concat([df_stats, pd.DataFrame([data])], ignore_index=True)

        self.features_to_include = [a for a in list(df_stats.columns) if a not in cols_key]
        return df_stats

    def _aggregate_choose(self, df_stats):
        '''
        Choose method of aggregation allowing the user to choose the statistics they require
        '''
        max_features = []; min_features = []; med_features = []; avg_features = []; sd_features = []; range_features = []; rand_features = []
        print('>> For the following list of features, please choosen which statistics you would like to include in the aggregation')
        print('>> The following statistics can be chosen (enter word or number seperated by ;): all (1), max (2), min (3), med (4), avg (5), sd (6), range (7), random (8)')
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
                        case 'med' | '4':
                            if feat not in med_features: med_features.append(feat)
                        case 'avg' | '5':
                            if feat not in avg_features: avg_features.append(feat)
                        case 'sd' | '6':
                            if feat not in sd_features: sd_features.append(feat)
                        case 'range' | '7':
                            if feat not in range_features: range_features.append(feat)
                        case 'random' | '8':
                            if feat not in rand_features: rand_features.append(feat)
                        case _:
                            print(f'>> Could not evaluate input ({tmp_stats}) for feature: {feat}, using all metrics for {feat} instead')
                            if feat not in max_features: max_features.append(feat)
                            if feat not in min_features: min_features.append(feat)
                            if feat not in med_features: med_features.append(feat)
                            if feat not in avg_features: avg_features.append(feat)
                            if feat not in sd_features: sd_features.append(feat)
                            if feat not in range_features: range_features.append(feat)
                            if feat not in rand_features: rand_features.append(feat)
            except Exception as e:
                print(f'>> Could not evaluate input ({tmp_stats}) for feature: {feat} with error {e}, using all metrics for {feat} instead')
        max_feat_cols = [a for a in self.features_to_include if any(b in a for b in max_features) and 'max' in a]
        min_feat_cols = [a for a in self.features_to_include if any(b in a for b in min_features) and 'min' in a]
        med_feat_cols = [a for a in self.features_to_include if any(b in a for b in med_features) and 'med' in a]
        avg_feat_cols = [a for a in self.features_to_include if any(b in a for b in avg_features) and 'avg' in a]
        sd_feat_cols = [a for a in self.features_to_include if any(b in a for b in sd_features) and 'sd' in a]
        range_feat_cols = [a for a in self.features_to_include if any(b in a for b in range_features) and 'range' in a]
        rand_feat_cols = [a for a in self.features_to_include if any(b in a for b in rand_features) and 'rand' in a]
        return df_stats.drop(columns=[a for a in self.features_to_include if a not in max_feat_cols + min_feat_cols + med_feat_cols + avg_feat_cols + sd_feat_cols + range_feat_cols + rand_feat_cols])


    def save_state(self, outname="measures_aggregated.csv"):
        '''
        Function saves a copy of the aggregated dataframe to a csv

        :param outname: the name of the csv file that the output is written to
        :type outname: str
        '''
        self.df_agg.to_csv(os.path.join(self.outdir, outname), index_label=False, index=False)


if __name__ == "__main__":
    test_dataframe_name = 'data/measures_cut_Ecoli(hCit)_all_01.05.25_joined.csv'
    test_measures_dataframe = pd.read_csv(test_dataframe_name)
    agg = Aggregation(test_measures_dataframe, aggregation_method='median', features_to_include=['all'], aev_red_method='pca')
    agg.aggregate_data()
    print(agg.df_agg)
    #agg.save_state()
