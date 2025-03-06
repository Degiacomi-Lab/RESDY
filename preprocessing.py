from copy import deepcopy
from ast import literal_eval
import numpy as np
import pandas as pd
from statsmodels.stats.outliers_influence import variance_inflation_factor as VIF
from statsmodels.tools.tools import add_constant
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, DBSCAN
from scipy.spatial.distance import euclidean
import matplotlib.pyplot as plt
import seaborn as sns


scaler = StandardScaler()

class Preprocessing:
    '''
    Class encompasing methods used for preprocesing the data before it is passed through
    to a model to train on the data. Functions include correlation analysis of the features
    using VIF and undersampling to be used to reduce the size of a dataset without losing
    key information.
    '''

    def __init__(self, df, features=[]):
        self.df = df
        self.n_obs = len(self.df)
        self.features = features

        self.df_cleaned, self.n_obs_cleaned = self.clean(df, features)
        self.data = self.df_cleaned[features]

        self.aev = np.zeros((len(self.df_cleaned), 1008))
        for i in range(self.n_obs_cleaned):
            self.aev[i] = np.array(self.df_cleaned['aev'][i].strip("[]").split(","), dtype=float)
        self.aev = pd.DataFrame(self.aev)
        self.aev = self.aev.loc[:, (self.aev != 0).any(axis=0)]

        self.data_full = pd.concat([self.data, self.aev], axis=1)

        # VIF variables
        self.vif = pd.DataFrame()

        # undersampling variables
        self.outliers = pd.DataFrame()
        self.common = pd.DataFrame()
        self.undersampled_data = pd.DataFrame()
        self.undersampled_data_common = pd.DataFrame()
        self.undersampled_data_outlier = pd.DataFrame()
        self.undersampled_data_normalised = pd.DataFrame()


    def clean(self, df, features):
        """
        Clean the input dataframe by removing rows containing NaN and duplicated rows.
        
        Parameters
        ----------
        df : Pandas DataFrame
        
        features : list
            List of features (columns) to be checked for NaN entries and duplication

        Returns
        -------
        df_cleaned : Pandas DataFrame
        n_obs_cleaned : int
            Number of rows in the cleaned dataframe.

        """
        # remove any features from measurements that
        print(f"Original number of observations: {len(df)}")
        df_cleaned = df.dropna(subset=features)
        print(f"Number of observations after removing rows containing NaN: {len(df_cleaned)}")
        df_cleaned = df_cleaned.drop_duplicates(subset=features)
        df_cleaned = df_cleaned.reset_index(drop=True)
        n_obs_cleaned = len(df_cleaned)
        print(f"Number of observations after removing duplicates: {n_obs_cleaned}")

        return df_cleaned, n_obs_cleaned



    def calculate_vif(self, data, features, multi_vif=False):
        """
        Calculate Variance Inflation Factors (VIFs) of the features selected.
        This allows for the n most decorrelated features to be selected later to take
        forward into the model, reducing the degrees of complexity.

        Parameters
        ----------
        data : Pandas DataFrame
            The overall measures dataframe
        features : list
            List of features to be considered for VIF calculations.

        """
        #X = scaler.fit_transform(data[features])
        if not multi_vif:
            data = data[features]
            data = add_constant(data)
            if 'aev' in data.columns:
                for i, r in data.iterrows():
                    new_aev = literal_eval(r['aev'])
                    data.at[i, 'aev'] = new_aev
                df_out = pd.DataFrame(data['aev'].to_list())
                df_out = df_out.add_prefix('AEV_')
                for col in df_out.columns:
                    col_vals = list(df_out[col])
                    if col_vals.count(0) == len(col_vals):
                        df_out = df_out.drop(col, axis=1)
                data = pd.concat([data, df_out], axis=1)
                data = data.drop('aev', axis=1)

        #print(data)

        data_correlation_matrix = data.corr()
        plt.figure(figsize=(10, 8))
        sns.heatmap(data_correlation_matrix, annot=True, cmap='coolwarm', vmin=-1, vmax=1)
        plt.title('Pairwise Correlation Matrix')
        plt.show()

        print(f'Number columns = {len(data.columns)}')

        #data = data.assign(const=1)
        vif_values = [VIF(data.values, i) for i in range(data.shape[1])]
        
        #vif_series = pd.Series([VIF(data.values, i) for i in range(data.shape[1])], index=data.columns)

        if self.vif.empty:
            self.vif['Parameter'] = data.columns
            self.vif.set_index('Parameter')
            new_col_header = 0
        elif len(list(self.vif.columns)) != 0:
            new_col_header = int(list(self.vif.columns)[-1]) + 1
        else:
            new_col_header = 0

        tmp_vif_df = pd.DataFrame({'Parameter': data.columns,
                                   new_col_header: vif_values},
                                   columns=['Parameter', new_col_header])

        self.vif = pd.merge(self.vif, tmp_vif_df, on='Parameter', how='left')



    def calculate_diff_features(self, data):
        '''
        Function to calculate the most decorrelated features from the measurements through
        VIF analysis. This will continuously call the VIF calculation until all the values
        returned are less than 5 (the commonly used value for decorrelation)
        '''

        all_decorrelated = False

        data = data[self.features]
        data = add_constant(data)
        print(data)
        if 'aev' in data.columns:
            for i, r in data.iterrows():
                new_aev = literal_eval(r['aev'])
                data.at[i, 'aev'] = new_aev
            df_out = pd.DataFrame(data['aev'].to_list())
            df_out = df_out.add_prefix('AEV_')
            for col in df_out.columns:
                col_vals = list(df_out[col])
                if col_vals.count(0) == len(col_vals):
                    df_out = df_out.drop(col, axis=1)
            data = pd.concat([data, df_out], axis=1)
            data = data.drop('aev', axis=1)

        while not all_decorrelated:
            # check for correlation and then change the data
            print(data)
            if not self.vif.empty:
                latest_vals = list(self.vif[list(self.vif.columns)[-1]])
                if all(x < 5 for x in latest_vals):
                    all_decorrelated = True
                
                # remove the column with the highest vif
                max_val_idx = self.vif[list(self.vif.columns)[-1]].idxmax()
                max_col = self.vif['Parameter'].iloc[max_val_idx]
                print(f'max_col: {max_col}, num cols:{len(list(data.columns))}')
                data = data.drop(max_col, axis='columns')
                print(f'num cols:{len(list(data.columns))}')

            # calculate new set of vif values
            self.calculate_vif(data, self.features, multi_vif=True)
            print(self.vif)

        print(self.vif)




    def _normalise(self, data_input="data", features=None):
        """
        Normalise all specified feature columns to mean zero, standard deviation 1.

        Parameters
        ----------
        data_input : Pandas DataFrame, optional
            If no dataframe specified, normalise the cleaned dataframe and 
            return the result to class attribute self.data_normalised
        features : list, optional
            List of features to be normalised. If none specified, normalise all 
            columns in the given dataframe.

        Returns (if a dataframe is given)
        -------
        data : Pandas DataFrame
            Normalised dataframe.

        """

        _case = 1
        if isinstance(data_input, pd.core.frame.DataFrame):
            data = deepcopy(self.data)
            _case = 0
        else:
            data = deepcopy(data_input)

        if features is None:
            features = data.columns
        for feature in features:
            data[feature] = (data[feature] - np.mean(data[feature])) / np.std(data[feature])

        if _case == 0:
            self.data_normalised = data
        else:
            return data



    def undersampling(self, features, n_cluster, n_init=100, max_iter=500,
                      iqr_reject_range=1.5, outlier_cluster_radius=0.6):
        '''
        Function to select a sample of points from a dataset which is representative
        of the entire dataset that has been fed, finding the most different points.
        
        Method
        ------
        First separate the dataset into two parts, the 'common' and the 
        'outliers', based on interquartile range (IQR). The sensitivity can be 
        changed by fixing the argument iqr_reject_range: all data points having 
        at least one feature lie in the rejection region are considered as 
        outliers; this region is defined by 
        [25th percentile - iqr_reject_range*IQR, 75th percentile + iqr_reject_range*IQR].
        
        +) Common part: use K-means to group data points into clusters, then 
        from each cluster choose the point closest to cluster centre.
        
        +) Outlier part: use DBSCAN algorithm to find small clusters in the set 
        of outliers, then from each cluster choose the point closest to cluster 
        centre. The sensitivity can be tuned by changing the argument 
        outlier_cluster_radius.
        

        Parameters
        ----------
        features : list
            List of features to be considered in the selection process.
        n_cluster : int
            Number of data points to be chosen from the common part.
        n_init : int, optional
            Number of runs for the K-means. The best run is chosen as the final 
            result. The default is 100.
        max_iter : int, optional
            Number of iterations in each run. The default is 500.
        iqr_reject_range : float, optional
            Determines the boundary between 'common' and 'outliers'. The 
            default is 1.5.
        outlier_cluster_radius : float, optional
            Determines the sensitivity of DBSCAN. The default is 0.6.


        Returns
        -------
        Selections are stored in the attribute self.undersampled_data

        '''


        data = self.data[features]
        outlier_indices = np.array([])
        for feature in features:
            quartiles = np.percentile(data[feature], [25,75])
            iqr = quartiles[1] - quartiles[0]
            new_outliers = np.where((data[feature] < quartiles[0]-iqr_reject_range*iqr)|(data[feature] > quartiles[1]+iqr_reject_range*iqr))[0]
            outlier_indices = np.unique(np.concatenate((outlier_indices, new_outliers)))

        data_outliers = data.iloc[outlier_indices]
        data_outliers = data_outliers.reset_index(drop=True)
        self.outliers = data_outliers
        data_outliers_normalised = self._normalise(data_outliers, features)

        data_common = data.drop(outlier_indices)
        data_common = data_common.reset_index(drop=True)
        self.common = data_common
        data_common_normalised = self._normalise(data_common, features)

        print("K-means clustering...")
        kmeans = KMeans(n_clusters=n_cluster, n_init=n_init, max_iter=max_iter).fit(data_common_normalised)
        avgdistance = np.sqrt(kmeans.inertia_/len(data_common_normalised))
        r = outlier_cluster_radius * avgdistance


        common_pts_idx = []
        for clust in range(kmeans.n_clusters):
            cluster_pts_indices = np.where(kmeans.labels_ == clust)[0]
            cluster_cen = kmeans.cluster_centers_[clust]
            min_idx = np.argmin([euclidean(data_common_normalised.iloc[idx], cluster_cen) for idx in cluster_pts_indices])
            #print(np.array([euclidean(data_common_normalised.iloc[idx], cluster_cen) for idx in cluster_pts_indices])/r)

            common_pts_idx.append(cluster_pts_indices[min_idx])


        print("Outliers clustering...")
        dbscan = DBSCAN(eps=r, min_samples=1, metric='euclidean').fit(data_outliers_normalised)
        oclusters = np.unique(dbscan.labels_)
        outlier_pts_idx = []
        for oclust in oclusters:
            ocluster_pts_indices = np.where(dbscan.labels_ == oclust)[0]
            if len(ocluster_pts_indices) == 1:
                outlier_pts_idx.append(ocluster_pts_indices[0])
                continue
            ocluster_cen = np.mean(data_outliers_normalised.iloc[ocluster_pts_indices], axis=0)
            min_idx = np.argmin([euclidean(data_outliers_normalised.iloc[idx], ocluster_cen) for idx in ocluster_pts_indices])

            outlier_pts_idx.append(ocluster_pts_indices[min_idx])

        common_final_list = data_common.iloc[common_pts_idx]
        outlier_final_list = data_outliers.iloc[outlier_pts_idx]
        common_normalised_final = data_common_normalised.iloc[common_pts_idx]
        outlier_normalised_final = data_outliers_normalised.iloc[outlier_pts_idx]

        self.undersampled_data = pd.concat([common_final_list, outlier_final_list], ignore_index=True)
        self.undersampled_data_common = common_final_list
        self.undersampled_data_outlier = outlier_final_list
        self.undersampled_data_normalised = pd.concat([common_normalised_final, outlier_normalised_final], ignore_index=True)
