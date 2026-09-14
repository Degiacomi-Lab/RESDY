'''
General notes on the work still outstanding in this file. There may be more further down.

.. todo::

   Look into the GO term functions and see if these still actually work with all the extra
   material that has been added in (GW, 13.09.24).

.. todo::

   The dropna function was removed on init, a function that cleans the dataframe at the start is
   needed instead. All null rows should not be removed indiscriminately, in case they are only
   null for some measurements and those measurements are not being used (GW, 16.04.25).
'''

import os
import re
import urllib.request
import threading
import concurrent.futures
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from .aggregation import Aggregation
from scipy.stats import fisher_exact

try:
    from statsmodels.stats.multitest import multipletests
    statsmodel_available = True
except Exception as e:
    statsmodel_available = False
    print(f'>> Packages available for calculating enrichment analysis '
          f'(statsmodel.multipletests) are not available. Won\'t be '
          f'able to calculate this. Error: {e}')

class Analysis(object):
    '''
    Class to handle nalysis of a measurements dataframe. GO term analysis, 
    '''

    def __init__(self, df, outdir="result", features_to_analyse = []):
        '''
        Initialise the Analysis class which allows you to create graph and go over other metrics
        such as GO terms

        :param df: Dataframe of measurements to go over the analysis for
        :type df: pandas.DataFrame
        :param outdir: Name of the directory to write to
        :type outdir: str
        :param features_to_analyse: List of features which should be analysed over. If ['all'] is
            passed, all features available in the table will be used.
        :type features_to_analyse: list
        '''
        if isinstance(df, str):
            self.df = pd.read_csv(df)
        else:
            self.df = df

        if features_to_analyse == ['all']:
            self.df = self.df
        else:
            self.df = self.df.dropna(subset=features_to_analyse)

        standard_cols = ['Uniprot_Entry', 'Chain', 'Resid', 'Resolution', 'Method', 'Class',
                        'Modified', 'PLDDT']
        self.features_to_analyse = [a for a in self.df.columns if a not in standard_cols]

        self.df_aggregated = pd.DataFrame(columns = ['Uniprot_Entry','Resid','Num'])

        self.df_sub = pd.DataFrame(columns=['Uniprot_Entry'])

        self.GO_dict = {} # code as the key

        self.name_to_code = {}
        self.code_to_name = None

        self.outdir = outdir
        self.lys_key = ['Uniprot_Entry', 'Chain', 'Resid']


    def get_data(self, uniprot_entry, resid):
        '''
        Retrieve data from the given dataframe which match up to a given uniprot
        entry and residue.

        :param uniprot_entry: Code for the uniprot entry to subset
        :type  uniprot_entry: str
        :param resid: Residue number to subset
        :type resid: str

        :returns: Dataframe containing only information from an alphafold structure
        :rtype: pandas.DataFrame
        '''
        df_query = self.df[(self.df['Uniprot_Entry'] == uniprot_entry) & (self.df['Resid'] == resid)]
        return df_query


    def get_data_alphafold(self):
        '''
        Retrieve data from the given dataframe to the class which were produced by alphafold.

        :returns: Dataframe containing only information from an alphafold structure
        :rtype: pandas.DataFrame
        '''
        df_query = self.df[self.df['Method'] == 'Predicted']
        return df_query


    def GO_search_term(self, df, code = '', name = ''):
        '''
        List the subset of UNIPROT codes associated with a GO Term

        :param df: Dataframe to get matches to GO term from
        :type  df: pandas.DataFrame
        :param code: Code of the Go term to analyse; format: 
        :type code: str
        :param name: Name corresponding to the GO term
        :type name: str
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
                html_2 = response.read().decode('utf-8')
        except Exception as e:
            print(f'Failed to obtain UNIPROT data for {uniprot_code}. {e}')

        try:
            for line in html_2.splitlines():
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


    def plot_feature_histogram(self, plot_type='all', feature='', agg_type='', uniprot='',
                               chain='', resid='', save_name=''):
        '''
        Create a basic plot showing the distribution of values for the specified feature
        across the set of measurements passed into the analysis class. An aggregation
        method will be employed to reduce biases, default here is avg unless other
        method passed.

        :param plot_type: Asks for while type of plot you want; options are 'all' which
            will include graphs for all possible features, 'single' which will just plot
            a single graph for the selected feature or 'agg' which will plot a collation
            of histograms of distribution of data after different aggregation methods
            for a single feature.
        :type plot_type: str
        :param feature: Name of the feature to produce histogram of data for
        :type feature: str
        :param agg_type: Name of aggregation method to use when plotting the histogram.
            For more information, see documentation of aggregation class.
        :type agg_type: str
        :param uniprot: Uniprot code given if a specific analysis of feature data for a
            uniprot code is required. Default is left as '' and will take full dataframe
            unless code given.
        :type uniprot: str, optional
        :param chain: Chain given if a specific analysis of feature data for a chain is 
            required. Default is left as '' and will take full dataframe unless code given.
        :type chain: str, optional
        :param resid: Resid given if a specific analysis of feature data for a resid is 
            required. Default is left as '' and will take full dataframe unless code given.
        :type resid: str, optional
        :param save_name: File name to save the histogram to.
        :type save_name: str, optional
        '''
        try:
            plt.clf()
        except Exception:
            pass

        if plot_type == 'single' and feature == '':
            print('Selected single feature for plot type but no feature given as input. '
                  'Please enter a feature when calling the function.')
            return
        elif plot_type == 'single' and feature not in self.df.columns:
            print(f'Selected single feature for plot type but feature given as input (input: '
                  f'{feature}) not in the dataframe given to the class. Please enter a feature '
                  f'when calling the function.')
            return
        elif plot_type == 'agg' and feature == '':
            print('Selected agg for plot type but no feature given as input. Please '
                  'enter a feature when calling the function.')
            return

        if uniprot != '' and chain != '' and resid != '':
            df_plot = self.df[(self.df['Uniprot_Entry'] == uniprot) & (self.df['Chain'] == chain) & (self.df['Resid'] == resid)]
        elif uniprot != '' and chain != '' and resid == '':
            df_plot = self.df[(self.df['Uniprot_Entry'] == uniprot) & (self.df['Chain'] == chain)]
        elif uniprot != '' and chain == '' and resid != '':
            df_plot = self.df[(self.df['Uniprot_Entry'] == uniprot) & (self.df['Resid'] == resid)]
        elif uniprot == '' and chain != '' and resid != '':
            df_plot = self.df[(self.df['Chain'] == chain) & (self.df['Resid'] == resid)]
        elif uniprot != '' and chain == '' and resid == '':
            df_plot = self.df[(self.df['Uniprot_Entry'] == uniprot)]
        elif uniprot == '' and chain != '' and resid == '':
            df_plot = self.df[(self.df['Chain'] == chain)]
        elif uniprot == '' and chain == '' and resid != '':
            df_plot = self.df[(self.df['Resid'] == resid)]
        else:
            df_plot = self.df

        feature_labels = {'depth': 'Depth (Å)',
                        'sasa': 'SASA (Å\u00b2)',
                        'propka': 'pKa (propka)',
                        'flexibility': 'Flexibility (B-factor) (Å\u00b2)',
                        'curvature': 'Curvature',
                        'arc_length': 'Arc Length',
                        'das': 'Dynamically Accessible Surface',
                        'legolas': '15N NMR Backbone Shift (ppm) (LEGOLAS)',
                        'phi': 'Phi',
                        'psi': 'Psi',
                        'seqcharge': 'Sequence Charge',
                        'torsion': 'Torsion Angle',
                        'writhing': 'Writhing',
                        'frustration': 'Frustration',
                        'rmsf': 'RMSF'}

        agg_labels = {'avg': 'Average',
                      'med': 'Median',
                      'min': 'Minimum',
                      'max': 'Maximum',
                      'rand': 'Random',
                      'sd': 'Standard Deviation',
                      'range': 'Range'}

        feature_histwidths = {'depth': 0.1,
                            'sasa': 2,
                            'propka': 0.2,
                            'flexibility': 0.1,
                            'curvature': 0.1,
                            'arc_length': 0.1,
                            'das': 2,
                            'legolas': 2,
                            'phi': 2,
                            'psi': 2,
                            'seqcharge': 0.5,
                            'torsion': 0.2,
                            'writhing': 0.02,
                            'frustration': 0.5,
                            'rmsf': 0.2}

        match plot_type:
            case 'single':
                if agg_type == '':
                    print('No aggregation type given, using avg')
                    agg_type = 'avg'

                agg_feature = f'{feature}_{agg_type}'

                agg = Aggregation(df_measurements=df_plot,
                                aggregation_method=agg_type,
                                features_to_include=[feature])
                df_plot = agg.aggregate_data()

                fig, ax = plt.subplots()
                fig.set_figheight(8)
                fig.set_figwidth(8)

                x_left = df_plot[agg_feature].min()
                x_right = df_plot[agg_feature].max()
                step=feature_histwidths[feature]

                palatinate_colour = '#682860'
                ax.hist(df_plot[agg_feature], rwidth=1, density=True, histtype='bar',
                        bins=np.arange(x_left, x_right, step), color=palatinate_colour,
                        alpha=0.5, label='')
                sns.kdeplot(df_plot[agg_feature], color=palatinate_colour, clip=(x_left, x_right), ax=ax)
                ax.set_xlabel(feature_labels[feature])
                ax.set_xlim(x_left, x_right)

                ax.set_title(f'Histogram Feature Analysis: {feature_labels[feature]}')

            case 'all':
                non_feat_cols = ['Uniprot_Entry', 'PDB_Code', 'Chain', 'Resid', 'PLDDT', 'Method',
                                 'Resolution', 'Modified', 'class']
                feat_cols = [a.split('_')[0] for a in df_plot.columns if a not in non_feat_cols]
                if 'arc' in feat_cols:
                    feat_cols.remove('arc')
                    feat_cols.append('arc_length')
                if 'aev' in feat_cols:
                    print(f'> Feature \'aev\' in dataframe, ignoring for scalar feature histogram creation')
                    feat_cols.remove('aev')
                if 'evolution' in feat_cols:
                    print(f'> Feature \'evolution\' in dataframe, ignoring for scalar feature histogram creation')
                    feat_cols.remove('evolution')

                if agg_type == '':
                    print('No aggregation type given, using avg')
                    agg_type = 'avg'

                agg = Aggregation(df_measurements=df_plot,
                                aggregation_method=agg_type,
                                features_to_include=feat_cols)
                df_plot = agg.aggregate_data()

                if len(feat_cols) <= 9: col_len = 3
                elif len(feat_cols) <= 16: col_len = 4
                else: col_len = 5

                row_len = math.ceil(len(feat_cols) / col_len)

                fig, axs = plt.subplots(row_len, col_len)
                fig.set_figheight(12)
                fig.set_figwidth(12)
                fig.suptitle(f'Histogram All Features Analysis: Aggregation Type={agg_type}')
                fig.subplots_adjust(left=0.07, right=0.98, top=0.95, bottom=0.05, wspace=0.27, hspace=0.25)

                used_plots = []
                for i, feat in enumerate(feat_cols):
                    agg_feature = f'{feat}_{agg_type}'
                    row = math.floor(i / col_len)
                    col = i % col_len
                    used_plots.append([row, col])

                    x_left = df_plot[agg_feature].min()
                    x_right = df_plot[agg_feature].max()
                    step=feature_histwidths[feat]

                    palatinate_colour = '#682860'
                    axs[row,col].hist(df_plot[agg_feature], rwidth=1, density=True, histtype='bar',
                            bins=np.arange(x_left, x_right, step), color=palatinate_colour,
                            alpha=0.5, label='')
                    sns.kdeplot(df_plot[agg_feature], color=palatinate_colour, clip=(x_left, x_right), ax=axs[row,col])
                    axs[row,col].set_xlabel(feature_labels[feat])
                    axs[row,col].set_xlim(x_left, x_right)

                for i in range(row_len * col_len):
                    row = math.floor(i / col_len)
                    col = i % col_len
                    if [row, col] not in used_plots:
                        fig.delaxes(axs[row][col])

            case 'agg':
                agg_types = ['avg', 'med', 'sd', 'range', 'rand', 'max', 'min']
                if agg_type == '':
                    print('No aggregation type given, using all for aggregation histogram plot')
                    agg_type = 'all'

                agg = Aggregation(df_measurements=df_plot,
                                aggregation_method=agg_type,
                                features_to_include=[feature])
                df_plot = agg.aggregate_data()

                non_feat_cols = ['Uniprot_Entry', 'PDB_Code', 'Chain', 'Resid', 'PLDDT', 'Method',
                                'Resolution', 'Modified', 'class']
                agg_types = [a.split('_')[-1] for a in df_plot.columns if a not in non_feat_cols]

                if len(agg_types) <= 4: col_len = 2
                else: col_len = 3

                row_len = math.ceil(len(agg_types) / col_len)

                fig, axs = plt.subplots(row_len, col_len)
                fig.set_figheight(12)
                fig.set_figwidth(12)
                fig.suptitle(f'Histogram All Aggregation Analysis; Feature: {feature_labels[feature]}')
                fig.subplots_adjust(left=0.05, right=0.98, top=0.92, bottom=0.05, wspace=0.27, hspace=0.27)

                used_plots = []
                for i, agg in enumerate(agg_types):
                    agg_feature = f'{feature}_{agg}'
                    row = math.floor(i / col_len)
                    col = i % col_len
                    used_plots.append([row, col])

                    x_left = df_plot[agg_feature].min()
                    x_right = df_plot[agg_feature].max()
                    step=feature_histwidths[feature]

                    palatinate_colour = '#682860'
                    axs[row,col].hist(df_plot[agg_feature], rwidth=1, density=True, histtype='bar',
                            bins=np.arange(x_left, x_right, step), color=palatinate_colour,
                            alpha=0.5, label='')
                    sns.kdeplot(df_plot[agg_feature], color=palatinate_colour, clip=(x_left, x_right), ax=axs[row,col])
                    axs[row,col].set_xlabel(feature_labels[feature])
                    axs[row,col].set_xlim(x_left, x_right)
                    axs[row,col].set_title(agg_labels[agg])

                for i in range(row_len * col_len):
                    row = math.floor(i / col_len)
                    col = i % col_len
                    if [row, col] not in used_plots:
                        fig.delaxes(axs[row][col])
            case _:
                print(f'>> Plot type give ({plot_type}) not recognised, please choose either all, single, agg')

        if save_name != '':
            plt.savefig(save_name)

        plt.show()


    def plot_feature_violins(self, features='', agg_type='', uniprot='',
                               chain='', resid='', save_name=''):
        '''
        Create a basic plot showing the distribution of values for the specified feature
        across the set of measurements passed into the analysis class. An aggregation
        method will be employed to reduce biases, default here is avg unless other
        method passed.

        :param features: Name of the feature or list of features to produce violins of
            data for
        :type features: str or list
        :param agg_type: Name of aggregation method to use when plotting the histogram.
            For more information, see documentation of aggregation class.
        :type agg_type: str
        :param uniprot: Uniprot code given if a specific analysis of feature data for a
            uniprot code is required. Default is left as '' and will take full dataframe
            unless code given.
        :type uniprot: str, optional
        :param chain: Chain given if a specific analysis of feature data for a chain is 
            required. Default is left as '' and will take full dataframe unless code given.
        :type chain: str, optional
        :param resid: Resid given if a specific analysis of feature data for a resid is 
            required. Default is left as '' and will take full dataframe unless code given.
        :type resid: str, optional
        :param save_name: File name to save the histogram to.
        :type save_name: str, optional
        '''
        try:
            plt.clf()
        except Exception:
            pass

        if agg_type == '':
            print('No aggregation type given as input; will include all options in the graphs')
            agg_type = 'all'

        if uniprot != '' and chain != '' and resid != '':
            df_plot = self.df[(self.df['Uniprot_Entry'] == uniprot) & (self.df['Chain'] == chain) & (self.df['Resid'] == resid)]
        elif uniprot != '' and chain != '' and resid == '':
            df_plot = self.df[(self.df['Uniprot_Entry'] == uniprot) & (self.df['Chain'] == chain)]
        elif uniprot != '' and chain == '' and resid != '':
            df_plot = self.df[(self.df['Uniprot_Entry'] == uniprot) & (self.df['Resid'] == resid)]
        elif uniprot == '' and chain != '' and resid != '':
            df_plot = self.df[(self.df['Chain'] == chain) & (self.df['Resid'] == resid)]
        elif uniprot != '' and chain == '' and resid == '':
            df_plot = self.df[(self.df['Uniprot_Entry'] == uniprot)]
        elif uniprot == '' and chain != '' and resid == '':
            df_plot = self.df[(self.df['Chain'] == chain)]
        elif uniprot == '' and chain == '' and resid != '':
            df_plot = self.df[(self.df['Resid'] == resid)]
        else:
            df_plot = self.df

        feature_labels = {'depth': 'Depth (Å)',
                        'sasa': 'SASA (Å\u00b2)',
                        'propka': 'pKa (propka)',
                        'flexibility': 'Flexibility (B-factor) (Å\u00b2)',
                        'curvature': 'Curvature',
                        'arc_length': 'Arc Length',
                        'das': 'Dynamically Accessible Surface',
                        'legolas': '15N NMR Backbone Shift (ppm) (LEGOLAS)',
                        'phi': 'Phi',
                        'psi': 'Psi',
                        'seqcharge': 'Sequence Charge',
                        'torsion': 'Torsion Angle',
                        'writhing': 'Writhing',
                        'frustration': 'Frustration',
                        'rmsf': 'RMSF'}

        agg_labels = {'avg': 'Average',
                      'med': 'Median',
                      'min': 'Minimum',
                      'max': 'Maximum',
                      'rand': 'Random',
                      'sd': 'Standard Deviation',
                      'range': 'Range'}

        non_feat_cols = ['Uniprot_Entry', 'PDB_Code', 'Chain', 'Resid', 'PLDDT', 'Method',
                        'Resolution', 'Modified', 'class']

        if features == '' or features == []:
            print('No features given as input, using all possible scalar features')
            features = [a for a in df_plot.columns if a not in non_feat_cols]

        if isinstance(features, str):
            features = [features]

        features = [a for a in features if a not in ['aev', 'evolution']]
        num_feats = len(features)

        if num_feats == 1:
            fig, ax = plt.subplots()
            fig.set_figheight(6)
            fig.set_figwidth(8)

            for i, feat in enumerate(features):
                agg = Aggregation(df_measurements=df_plot,
                                aggregation_method=agg_type,
                                features_to_include=feat)
                df_plot_feat = agg.aggregate_data()

                agg_non_feat_cols = [a for a in df_plot_feat.columns if a in non_feat_cols]
                df_plot_feat = df_plot_feat.drop(columns=agg_non_feat_cols)
                df_plot_feat = df_plot_feat.rename(columns={k: v for k, v in zip(list(df_plot_feat.columns),
                                        [agg_labels[a.split('_')[-1]] for a in list(df_plot_feat.columns)])})

                sns.violinplot(data=df_plot_feat,
                            inner="quart",
                            fill=False,
                            color='#682860',
                            ax=ax)

                if i == num_feats - 1:
                    ax.set(xlabel='Aggregation Types')
                ax.set(ylabel=feature_labels[feat])
                ax.set(title=f'Feature Aggregation Analysis: {feature_labels[features[0]]}')

            plt.xticks(rotation=20)

        else:
            fig, axs = plt.subplots(num_feats)
            fig.set_figheight(3*num_feats)
            fig.set_figwidth(8)
            fig.subplots_adjust(left=0.07, right=0.98, top=0.95, bottom=0.05, wspace=0, hspace=0)

            for i, feat in enumerate(features):
                agg = Aggregation(df_measurements=df_plot,
                                aggregation_method=agg_type,
                                features_to_include=feat)
                df_plot_feat = agg.aggregate_data()

                agg_non_feat_cols = [a for a in df_plot_feat.columns if a in non_feat_cols]
                df_plot_feat = df_plot_feat.drop(columns=agg_non_feat_cols)
                df_plot_feat = df_plot_feat.rename(columns={k: v for k, v in zip(list(df_plot_feat.columns),
                                        [agg_labels[a.split('_')[-1]] for a in list(df_plot_feat.columns)])})

                sns.violinplot(data=df_plot_feat,
                            inner="quart",
                            fill=False,
                            color='#682860',
                            ax=axs[i])

                if i == num_feats - 1:
                    axs[i].set(xlabel='Aggregation Types')
                else:
                    axs[i].set(xticks=[])
                axs[i].set(ylabel=feature_labels[feat])

            fig.suptitle('Violin Feature Aggregation Analysis')
            plt.xticks(rotation=20)


        if save_name != '':
            plt.savefig(save_name)

        plt.show()


    def get_outliers(self, uniprot_entry, resid, feature, whis = 1.5):
        '''
        Extract outliers from a dataset using quartiles based on a given feature and specific
        uniprot code and resid number.

        :param uniprot_entry: Uniprot code of interest
        :type uniprot_entry: str
        :param resid: Residue number of the uniprot code of interest
        :type resid: str
        :param feature: Feature of interest to extract outliers over
        :type feature: str
        :param whis:
        :type whis: float
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
        '''
        Retrieve data from the given dataframe which are outside the lower and upper
        boundaries for that feature of interest.

        :param feature: Feature of interest to examine
        :type feature: str
        :param lower: Lower boundary for a value for the feature to not be considered extreme
        :type lower: int
        :param upper: Upper boundary for a value for the feature to not be considered extreme
        :type upper: int

        :returns: Dataframe containing only information from an alphafold structure
        :rtype: pandas.DataFrame
        '''
        df_query = self.df[(self.df[feature] < int(lower)) | (self.df[feature] > int(upper))]
        return df_query


    def add_extra_measures(self, extra_measures_filename, write_new_file = False,
                           out_filename='measures_new.csv'):
        '''
        Function to add in extra measurements to the measures frame that has been autoloaded into
        the analysis class on defining this. This will match up the measurements in each case and
        hold in for the analysis. A new measures file will be written with the new filename that has
        been passed into the function.

        :param extra_measures_filename: The name of the new measures file written of the combination
            of both measures dataframe.
        :type extra_measures_filename: str
        :param write_new_file: True/False option for writing a new measures.csv file when the new
            data has been added in. Auto set to False.
        :type write_new_file: bool
        :param out_filename: The name of the new measures.csv file that you want to be produced.
            Auto set to be measures_new.csv
        :type out_filename: str
        '''

        # Step 1: read in new dataframe, extract column names, check for overlap and
        #         deal if is, otherwise add new column in
        try:
            if isinstance(extra_measures_filename, str):
                new_df = pd.read_csv(extra_measures_filename)
            if 'Unnamed: 0' in new_df.columns: new_df = new_df.drop(columns='Unnamed: 0')
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
                            if column in ['aev', 'aev_legolas']:
                                self.df[column] = self.df[column].astype('object')
                            for i, r in self.df.iterrows():

                                protein_code = r['PDB_Code']
                                chain_value = r["Chain"]
                                resid_value = r["Resid"]

                                idx = np.where((new_df["PDB_Code"] == protein_code) &
                                               (new_df["Chain"] == chain_value) &
                                               (new_df["Resid"].astype(int) == resid_value))
                                if len(idx[0]) == 0:
                                    continue

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

                    idx = np.where((new_df["PDB_Code"] == protein_code) &
                                   (new_df["Chain"] == chain_value) &
                                   (new_df["Resid"].astype(int) == resid_value))
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
        Function to take the input file documenting which residues are required to keep due to being
        of interest and remove anything from the dataframe that is not in this list. This is required
        due to the codebase calculating data for every possible resid in the structure.

        .. rubric:: Method

        Extract the list of residues and taking data for these. Goes over the dataframe and extracts
        any residues which are not present within the required residues. Removes these from the
        dataframe and then writes a new dataframe with the updated data.

        :param req_resid_table: Dataframe containing all the measured data inputted into the
            analysis class
        :type req_resid_table: pandas.DataFrame
        :param outname: The name of the file to give in output for the new updated measures file.
            Auto set to measures_cut.csv
        :type outname: str
        '''
        print('>> Removing unrequired residues')
        # Remove duplicated data from the measurements
        if isinstance(req_resid_table, str):
            df_req_res = pd.read_csv(req_resid_table)
        else:
            df_req_res = req_resid_table
        initial_data_one = len(self.df)
        self.df = self.df.drop_duplicates()
        duplicate_rows_removed = initial_data_one - len(self.df)
        print(f'Removed {duplicate_rows_removed} rows of duplicates')
        initial_full_data_rows = len(self.df)

        # remove rows which have a UNIPROT code which isnt required
        uniprot_codes = df_req_res['Uniprot_Entry'].drop_duplicates().tolist()
        entries_to_remove = []
        for i, r in self.df.iterrows():
            if r['Uniprot_Entry'] not in uniprot_codes:
                entries_to_remove.append(i)
        self.df = self.df.drop(index=entries_to_remove)
        uniprot_rows_removed = initial_full_data_rows - len(self.df)
        print(f'Removed {uniprot_rows_removed} rows of Uniprot codes which were not mentioned in the required residues file')

        while len(df_req_res) > 0:
            print("Number of rows left: " + str(len(df_req_res)))
            test_uniprot = df_req_res["Uniprot_Entry"][df_req_res.first_valid_index()]
            print("test_uniprot: " + str(test_uniprot))

            # find all the desired residues from the particular uniprot code and put into a list
            # automatically removes duplicates from this (doesn't retain order)
            desired_residues = list(set(df_req_res[df_req_res["Uniprot_Entry"] == test_uniprot]["Resid"].tolist()))

            # search the measures spreadsheet for all rows containing the desired uniprot code
            search_uniprot = self.df[self.df["Uniprot_Entry"] == test_uniprot.strip()][["Uniprot_Entry", "Resid"]]
            all_search_rows = search_uniprot.index.tolist()

            # go over each row of search_uniprot, see if the residue matches one of the desired ones
            wanted_rows = search_uniprot[search_uniprot["Resid"].isin(desired_residues)].index.tolist()
            not_wanted_rows = [x for x in all_search_rows if x not in wanted_rows]
            print("Rows removed: " + str(len(not_wanted_rows)))

            # remove the rows which aren't wanted from the main data set
            self.df = self.df.drop(index = not_wanted_rows)
            # remove rows which contain the uniprot code that has been searched from test_table
            df_req_res = df_req_res.drop(index = df_req_res[df_req_res["Uniprot_Entry"] == test_uniprot].index.tolist())

        final_full_data_rows = len(self.df)
        diff_rows = initial_full_data_rows - final_full_data_rows
        print(f'Original num of rows: {initial_full_data_rows}')
        print(f'Current num of rows: {final_full_data_rows}')
        print(f'Num of rows removed: {diff_rows}')
        self.df.to_csv(os.path.join(self.outdir, outname), index_label=False, index=False)
        return self.df


    def relative_best(self, df, weights, features=['depth']):
        '''
        Function to extract the best relative list of features for all combinations of uniprot entry
        and residues based on a given list of metrics to do the calculation over and the desired
        weightings for each of those features.

        :param df: Dataframe of measurements to do the analysis over
        :type df: pandas.DataFrame
        :param weights: List of floats which sum to 1 of the weights for each of the given features.
            Length should match the list of features given
        :type weights: list
        :param features: List of features that the analysis should extract the relative best values
            for.
        :type features: list
        '''
        # go over the weights to make sure the values are good and then matches the number of features
        try:
            if isinstance(weights, float) or isinstance(weights, int):
                if weights == 0.5:
                    print('>> Weight remains as default and equal for all metrics')
                    if len(features) == 1: weights = [1]
                    else: weights = [1/len(features)] * len(features)
                elif len(features) == 2:
                    print(f'>> 1 weight given, 2 features, assigning given weight ({weights}) to first feature ({features[0]}), second feature ({features[1]}) will take {1-weights} for weighting.')
                    weights = [weights, 1-weights]
                else:
                    print('>> Not enough weights entered to work out weightings for features. Assuming equal weight for all metrics')
                    weights = [1/len(features)] * len(features)
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
        unsure_feats = ['phi', 'psi', 'curvture', 'writhing', 'arc_length']
        min_feats = ['propka', 'pkaANI', 'legolas', 'depth']
        max_feats = ['sasa', 'das', 'seqcharge', 'frustration'] + unsure_feats  # unsure feats added to max feats for now
        for feat in features:
            if feat in ['aev', 'aev_legolas']:
                print(f'>> Feature {feat} is not supported with this analysis. Dropping feature from ')
            if feat not in min_feats + max_feats:
                pref = 'tmp'
                while pref not in ['min', 'max']:
                    pref = input(f'Enter bias towards min or max for {feat} (enter "min" or "max"): ')
                if pref == 'min': min_feats.append(feat)
                elif pref == 'max': max_feats.append(feat)

        for (uniprot_tmp, chain_tmp, resid_tmp), df_query in df.groupby(self.lys_key):
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
                            if feat_std == 0:
                                print(f'Standard deviation of feature: {feature} is 0 as is always constant, '
                                      f'{feature} therefore contributes 0 to relative best.')
                                feat_list_standardised = np.zeros(len(feat_list))
                            else:
                                feat_list_standardised = (feat_list - feat_avg) / feat_std
                            feat_list_weighted = feat_list_standardised * weight
                            metric_calculated_values_list_temp.append(feat_list_weighted)
                        except Exception as e:
                            print(f'>> Failed to load the data for the metric: {feature} with error: {e}')
                except Exception as e:
                    print(f'Error loading data for the metrics provided: {e}')

                try:
                    index = 0
                    base = -np.inf
                    for i, val in enumerate(metric_calculated_values_list_temp[0]):
                        weighted_sum = 0
                        for weighted_value_index, metric_value in enumerate(metric_calculated_values_list_temp):
                            weighted_sum += metric_value[i]
                        if np.isfinite(weighted_sum) and weighted_sum >= base:
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

        bp_not_list = len([p for p in reference if p in self.GO_dict[GO_code]]) - bp_list
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

        if not statsmodel_available:
            raise ImportError(f'>> GO enrichment analysis requires statsmodel.multipletests which isn\' available')

        # get all the uniprot codes inside the range and the reference uniprot code list
        feat_one, feat_one_low, feat_one_upper = str(feature_one[0]), float(feature_one[1]), float(feature_one[2])
        feat_two, feat_two_low, feat_two_upper = str(feature_two[0]), float(feature_two[1]), float(feature_two[2])
        selected_df = self.df[(self.df[feat_one] >= feat_one_low) & (self.df[feat_one] <= feat_one_upper)]
        selected_df = selected_df[(selected_df[feat_two] >= feat_two_low) & (selected_df[feat_two] <= feat_two_upper)]
        uniprot_selected = selected_df['Uniprot_Entry'].unique()
        uniprot_reference = self.df['Uniprot_Entry'].unique()

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
