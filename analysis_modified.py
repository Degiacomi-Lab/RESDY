#!/usr/bin/env python
# coding: utf-8

# In[13]:


import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go


# In[9]:

# this is a class that will be used to store a small subset of the data
class data(object):
    def __init__(self, df, uniprot_entry, resid):
        self.df = df
        self.uniprot_entry = uniprot_entry
        self.resid = resid
    
    def __str__(self):
        num_rows = self.df.shape[0]
        return f'{num_rows} instances under [{self.uniprot_entry}, {self.resid}]'


# In[12]:


class Analysis(object):
    def __init__(self, df):
        self.df = df.dropna() # this is the dataframe that will be aggregated
        self.df_withna = df # this is the dataframe that will be plotted in the interactive plot
        self.aggregated_df = pd.DataFrame(columns = ['Uniprot Entry','Resid','Object','pKa mean','pKa std','pKa range', 
                                                     'SASA mean','SASA std','SASA range'])    
    def get_unique_uniprot_entry(self):
        return self.df['Uniprot Entry'].unique()
    
    def get_unique_resid(self, uniprot_entry):
        current_df = self.df[self.df['Uniprot Entry'] == uniprot_entry]
        return current_df['Resid'].unique()
    
    # convenient way to get a subset of the data
    def get_data(self, entry, resid):
        current_df = self.df[(self.df['Uniprot Entry'] == entry) & (self.df['Resid'] == resid)]
        return current_df
    
    def aggregate(self):

        # get unique uniport entry
        unique_uniport_entries = self.get_unique_uniprot_entry()
        
        index = 0
        for entry in unique_uniport_entries:
            
            # get unique resid
            unique_resid = self.get_unique_resid(entry)
            
            for resid in unique_resid:
                current_df = self.df[(self.df['Uniprot Entry'] == entry) & (self.df['Resid'] == resid)]
                data_obj = data(current_df, entry, resid)

                # append a row to new_df
                pka_mean = round(current_df['pKa'].mean(),2)
                pka_std = round(current_df['pKa'].std(),2)
                pka_range = [round(current_df['pKa'].min(),2),round(current_df['pKa'].max(),2)]
                sasa_mean = round(current_df['sasa'].mean(),2)
                sasa_std = round(current_df['sasa'].std(),2)
                sasa_range = [round(current_df['sasa'].min(),2),round(current_df['sasa'].max(),2)]
                self.aggregated_df.loc[index] = [entry, resid, data_obj, pka_mean, pka_std, pka_range,
                                                 sasa_mean, sasa_std, sasa_range]
                index += 1
        
        print(f'In total {index} Uniport Entry and Resid combinations have been aggregated')
        return self.aggregated_df
     
    ########################################################################################
    
    # plot either a histogram with kde or a boxplot
    def plot_graph(self, plot_type, uniprot_entry, resid, feature):
        
        if uniprot_entry not in self.get_unique_uniprot_entry():
            return 'Wrong Uniprot Entry!'
        else:
            if resid not in self.get_unique_resid(uniprot_entry):
                return 'Wrong Resid!'
        
        current_df = self.df[(self.df['Uniprot Entry'] == uniprot_entry) & (self.df['Resid'] == resid)]
        x = current_df[feature]
        
        if plot_type == 'histogram':
            sns.displot(x, kde=True)
        
        elif plot_type == 'boxplot':
            sns.boxplot(x=x)
        
        else:
            print('Sorry No Such Plot Available!')
    
    # whis is a number controlling the extremeness of the outliers, conventionally it is 1.5
    def get_outliers(self, uniprot_entry, resid, feature, whis):
        current_df = self.df[(self.df['Uniprot Entry'] == uniprot_entry) & (self.df['Resid'] == resid)]
        x = current_df[feature]
        
        Q1 = x.quantile(0.25)
        Q3 = x.quantile(0.75)
        IQR = Q3 - Q1
        
        lower = Q1 - whis * IQR
        upper = Q3 + whis * IQR
        
        outlier_df = current_df[(current_df[feature] < lower) | (current_df[feature] > upper)]
        return outlier_df
    
    ##########################################################################################
    def interactive_plot(self):
        
        df_data = self.df_withna
        
        def update_point(trace, points, selector):
          
            # only work with gray points
            if len(points.point_inds) == 0:
                return  

            # gather uniprot and chain ID from full dataset clicked point
            idx = points.point_inds[0]
            my_uniprot = df_data.loc[idx, "Uniprot Entry"]
            my_resid = df_data.loc[idx, "Resid"]
            my_label = "%s(%i)"%(my_uniprot, my_resid)
            df_query = df_data[(df_data['Uniprot Entry'] == my_uniprot) & (df_data['Resid'] == my_resid)]

            
            #plot data associated with selection
            #labels = ["UNIPROT: %s<br>resid: %i"%(df_query["Uniprot Entry"].values[i], df_query["resid"].values[i]) for i in range(len(df_query))]
            labels = ["PDB: %s"%(df_query["PDB Code"].values[i]) for i in range(len(df_query))]

            # refresh the previous scatter plot
            if len(f.data)>1:
                f.data = [f.data[0]]


            f.add_scatter(x=df_query["sasa"], y=df_query["pKa"],
                          mode='markers', showlegend=False, name=my_label,
                          text = labels, hovertemplate='%{text}<br>SASA: %{x:.2f}<br>pKa: %{y:.2f}')

            
        # create scatter plot
        labels = ["UNIPROT: %s<br>resid: %i"%(df_data["Uniprot Entry"].values[i], df_data ["Resid"].values[i]) for i in range(len(df_data))]
        f = go.FigureWidget([go.Scatter(x=df_data["sasa"], y=df_data["pKa"],
                                        mode='markers', name="aggregate", showlegend=False, opacity=0.75,
                                        text = labels, hovertemplate='%{text}<br>SASA: %{x:.2f}<br>pKa: %{y:.2f}'
                                       )
                            ])

        # add labels
        f.update_layout(
            xaxis_title="SASA (A2)",
            yaxis_title="pKa")

        #set axes properties
        f.update_xaxes(range=[0, 130])
        f.update_yaxes(range=[-2, 20])
        f.update_xaxes(showspikes=True)
        f.update_yaxes(showspikes=True)

        #parameterize colours and triggers
        scatter = f.data[0]
        colors = ['#bae2be'] * len(df_data)
        scatter.marker.color = colors
        f.layout.hovermode = 'closest'
        scatter.on_click(update_point)
        
        return f


# In[ ]:




