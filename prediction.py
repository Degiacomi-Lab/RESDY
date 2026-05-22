# File to house the generalised useage for calling a model for prediction on the aggregated data

import random
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.lines import Line2D

from aggregation import Aggregation
from models.pulearning import PUlearning

class Prediction:
    '''
    Class for handling the overall functions required for running prediction models on the data
    produced within the pipeline. Individual model scripts will be kept in the directory models
    as subclasses. 
    '''

    def __init__(self,
                 pos_aggregated_data=None,
                 neg_aggregated_data=None,
                 aggregated_data=None):
        '''
        Initialise the prediction class and complete some of the early data tasks before passing
        onto the desired model.
        '''
        self.pos_aggregated_data = pos_aggregated_data
        self.neg_aggregated_data = neg_aggregated_data
        self.aggregated_data = aggregated_data
        self.df_agg = pd.DataFrame()
        self.df_prepped = pd.DataFrame()
        self.y = []

        self._collate_data()
        if not self._check_aggregation_status():
            quit


    def _collate_data(self):
        '''
        Go over the inputs for the aggregated data and sort this into the overall dataframe for
        model input
        '''
        if self.pos_aggregated_data is not None:
            tmp_pos = self._type_checking_work(self.pos_aggregated_data)
            if 'class' not in tmp_pos.columns: tmp_pos['class'] = 1
            self.df_agg = pd.concat([self.df_agg, tmp_pos], ignore_index=True)
        if self.neg_aggregated_data is not None:
            tmp_neg = self._type_checking_work(self.neg_aggregated_data)
            if 'class' not in tmp_neg.columns: tmp_neg['class'] = 0
            self.df_agg = pd.concat([self.df_agg, tmp_neg], ignore_index=True)
        if self.aggregated_data is not None:
            self.df_agg = pd.concat([self.df_agg, self._type_checking_work(self.aggregated_data)], ignore_index=True)

    def _type_checking_work(self, data_to_sort):
        '''
        Case handling function for working out what is needed to do with processing the input
        of data to the prediction class

        Parameters
        ----------
        data_to_sort -> dataframe, string, list
            Set of data input which needs transforming into dataframes to use further
        
        Returns
        -------
        df_out -> dataframe
            Dataframe from all the input to use further
        '''
        df_out = pd.DataFrame()
        if isinstance(data_to_sort, str):
            df_out = pd.read_csv(data_to_sort)
        elif isinstance(data_to_sort, pd.DataFrame):
            df_out = data_to_sort
        elif isinstance(data_to_sort, list):
            for dataset in data_to_sort:
                if isinstance(dataset, str):
                    df_out = pd.concat([df_out, pd.read_csv(dataset)], ignore_index=True)
                elif isinstance(dataset, pd.DataFrame):
                    df_out = pd.concat([df_out, dataset], ignore_index=True)
        return df_out

    def _check_aggregation_status(self):
        '''
        Go over the dataset that has been provided to the prediction class and check if there are
        any duplicates left within this. Break at this point to avoid false results production.
        '''
        df_checking = self.df_agg.duplicated(subset=['Uniprot_Entry', 'Resid', 'class'], keep=False)
        list_checking = list(set(list(df_checking)))
        if all(a is False for a in list_checking):
            print('No duplicates in dataset, data has been properly aggregated')
            return True
        else:
            df_duplicated = self.df_agg[df_checking]
            print('>> Duplicates present in the dataset. Sort this out before proceeding.')
            print('>> The following is the table of duplicates:')
            print(df_duplicated)
            return False


    def predict_carbamates(self, prediction_model):
        '''
        Central function for managing the prediction methods

        Parameters
        ----------
        prediction_model -> string
            desired model to use for the prediction output, choose from:

        '''
        match prediction_model:
            case 'pulearning':
                PUlearn = PUlearning(self.df_agg, self.y)
                PUlearn.basic_pulearn()
    
    def produce_pca_2d(self):
        '''
        Create a 2D PCA based on the aggregated data
        Note this should be moved to viewer or analysis later on when those are sorted out

        '''
        X = self.df_agg.drop(['Uniprot_Entry', 'Resid', 'class'], axis=1)
        y = self.df_agg['class']
        print(y)
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        pca = PCA(n_components=2)
        X_pca = pca.fit_transform(X_scaled)
        variance_one = round(pca.explained_variance_[0], 2); variance_two = round(pca.explained_variance_[1], 2)
        cust_cmap = matplotlib.colors.ListedColormap([(1, 0.498, 0.055), (0.408, 0.141, 0.427)])
        plt.figure()
        plt.scatter(X_pca[:, 0], X_pca[:, 1], c=y, cmap=cust_cmap, edgecolors='Black')
        plt.xlabel(f'PC 1 ({variance_one}%)')
        plt.ylabel(f'PC 2 ({variance_two}%)')
        plt.title("PCA Aggregated Data")
        legend_elements = [Line2D([0], [0], color=(0.408, 0.141, 0.427), lw=4, label='Carbamylated'),
                           Line2D([0], [0], color=(1, 0.498, 0.055), lw=4, label='Not Carbamylated')]
        plt.legend(handles=legend_elements, loc='best', frameon=True)
        plt.savefig('pca_allfeatures_2d.svg')
        plt.show()

    def _prepare_dataset(self):
        '''
        Prepare the dataset from the aggregated dataset to a general form for model input.
        Remove any overlap between the datasets to avoid any problems with prediction.
        Remove part of the negative dataset to set to match the number of positives within
        the dataset. Remove the unnecessary columns.
        '''
        all_overlaps = pd.DataFrame(columns=['Uniprot_Entry', 'PDB_Code', 'Resid'])
        for i, r in self.df_agg[self.df_agg['class'] == 1].iterrows():
            pdb = r['PDB_Code']
            resid = r['Resid']
            df_query = self.df_agg[(self.df_agg['PDB_Code'] == pdb) & (self.df_agg['Resid'] == resid)]
            if len(df_query) > 1:
                all_overlaps = pd.concat([all_overlaps, df_query[['Uniprot_Entry', 'PDB_Code', 'Resid', 'class']]], ignore_index=True)
                self.df_agg = self.df_agg.drop(index=df_query.index[1:])  # drop all duplicates after the first occurance
        print(f'Number of overlaps: {len(all_overlaps)}')
        num_pos_data = len(self.df_agg[self.df_agg['class'] == 1])
        indices_random_neg_data = random.sample(range(num_pos_data, len(self.df_agg)), num_pos_data)
        self.df_prepped = pd.concat([self.df_agg[self.df_agg['class'] == 1], self.df_agg.iloc[indices_random_neg_data]])
        self.y = self.df_prepped['class']
        self.df_prepped = self.df_prepped.drop(['Uniprot_Entry', 'Resid', 'class'], axis=1)
    
    def produce_pca_3d(self):
        '''
        Create a 3D PCA based on the aggregated data
        Note this should be moved to viewer or analysis later on when those are sorted out

        '''
        X = self.df_agg.drop(['Uniprot_Entry', 'Resid', 'class'], axis=1)
        y = self.df_agg['class']
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        pca = PCA(n_components=3)
        X_pca = pca.fit_transform(X_scaled)
        variance_one = round(pca.explained_variance_[0], 2); variance_two = round(pca.explained_variance_[1], 2); variance_three = round(pca.explained_variance_[2], 2)
        cust_cmap = matplotlib.colors.ListedColormap([(1, 0.498, 0.055), (0.408, 0.141, 0.427)])
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(X_pca[:, 0], X_pca[:, 1], X_pca[:, 2], c=y, cmap=cust_cmap, edgecolors='Black')
        ax.set_xlabel(f'PC 1 ({variance_one}%)')
        ax.set_ylabel(f'PC 2 ({variance_two}%)')
        ax.set_zlabel(f'PC 3 ({variance_three}%)')
        plt.title("PCA Aggregated Data")
        legend_elements = [Line2D([0], [0], color=(0.408, 0.141, 0.427), lw=4, label='Carbamylated'),
                           Line2D([0], [0], color=(1, 0.498, 0.055), lw=4, label='Not Carbamylated')]
        plt.legend(handles=legend_elements, loc='best', frameon=True)
        plt.savefig('pca_allfeatures.svg')
        plt.show()


if __name__ == '__main__':
    pos_measure_files = ['data/measures_cut_CannData_all_12.05.25.csv',
                        'data/measures_cut_KingHighConf_all_01.05.25_joined.csv',
                        'data/measures_cut_Ecoli(hCit)_all_01.05.25_joined.csv',
                        'data/measures_cut_Synecho(hCit)_all_01.05.25_joined.csv']
    neg_measure_files = ['data/measures_cut_KingAllNegative_all_01.05.25_joined.csv']
    features_to_include = ['propka', 'sasa', 'das', 'seqcharge', 'curvature', 'writhing', 'torsion', 'arc_length', 'phi', 'psi', 'legolas', 'aev']
    df_pos = pd.DataFrame(); df_neg = pd.DataFrame()
    for file in pos_measure_files:
        df_pos = pd.concat([df_pos, pd.read_csv(file)], ignore_index=True)
    for file in neg_measure_files:
        df_neg = pd.concat([df_neg, pd.read_csv(file)], ignore_index=True)
    df_pos['class'] = 1; df_neg['class'] = 0
    print(f'Len positve data: {len(df_pos)}, Len negative data: {len(df_neg)}')
    df_collated = pd.concat([df_neg, df_pos], ignore_index=True)
    agg_col = Aggregation(df_collated, aggregation_method='minmax', features_to_include=features_to_include, aev_red_method='pca')
    df_col_agg = agg_col.aggregate_data()
    num_pos_data = len(df_col_agg[df_col_agg['class'] == 1])
    print([num_pos_data, (len(df_col_agg) - num_pos_data)])
    prediction = Prediction(aggregated_data=df_col_agg)
    #print(prediction.df_agg.isna().any(axis=1))
    prediction.produce_pca_3d()
    #prediction.predict_carbamates(prediction_model='pulearning')
