from ipywidgets import HBox, VBox
from ipywidgets import widgets
import numpy as np
import pandas as pd
import nglview as nv
import plotly.express as px
import plotly.graph_objects as go
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import os
import webbrowser

class CoolPlots(object):
    
    def __init__(self, analysis, outdir = 'result'):
        self.analysis = analysis
        self.df = analysis.df
        self.df_concise = analysis.df_sub
        self.GO_dict = analysis.GO_dict
        
        # from Term to code
        self.GO_decode_dict = analysis.name_to_code
        
        # from code to Term
        self.GO_decode_dict_reverse = analysis.code_to_name
        
        # some temp dfs
        self.temp_df = pd.DataFrame()
        self.temp_df_2 = pd.DataFrame()
        
        # the path to store the regional data
        self.outdir = outdir
        self.export_path = ''
        
        # store the uniprot code clicked for the pop-up uniprot page
        self.uni_clicked = ''
        
        # plot
        self.plot = 'You have not called advanced_plot method.'
        
        # pka slider
        self.p = widgets.FloatRangeSlider(
                    value=[1, 14],
                    min=1,
                    max=14.0,
                    step=0.1,
                    description='pKa:',
                    disabled=False,
                    continuous_update=False,
                    orientation='horizontal',
                    readout=True,
                    readout_format='.1f',)
        
        # sasa slider
        self.s = widgets.FloatRangeSlider(
                    value=[0, 100],
                    min=0,
                    max=100.0,
                    step=0.1,
                    description='SASA:',
                    disabled=False,
                    continuous_update=False,
                    orientation='horizontal',
                    readout=True,
                    readout_format='.1f',)
        
        # GO Terms dropdown
        options = [f'{code}: {self.GO_decode_dict_reverse[code]}' for code in list(self.GO_dict.keys())]
        self.GO = widgets.Dropdown(
                    options=(['Welcome'] + sorted(options, key = lambda x: int(x.split(':')[0]))),
                    value='Welcome',
                    description='GO Terms: ',
                    disabled=False,)
        
        # PDB file text box
        self.PDB_box = widgets.Text(
                            value='Welcome',
                            placeholder='Type something',
                            description='PDB file:',
                            disabled=True)
        
        
        # call back function for the export buttom
        def call_back_buttom_export(b_export):
            if self.temp_df.empty:
                return
            df = self.temp_df
            columns = ['Uniprot Entry', 'Resid', 'Num', 'pKa_mean', 'sasa_mean', 'GO_Terms']
            df_out = pd.DataFrame(columns = columns)
            for idx, row in df.iterrows():
                my_uniprot = row['Uniprot Entry']
                my_resid = row['Resid']
                df_query = self.df[(self.df['Uniprot Entry'] == my_uniprot) & (self.df['Resid'] == my_resid)]
                GO_Terms = self.analysis.GO_Search_Protein(my_uniprot)
                data = ({'Uniprot Entry':my_uniprot, 'Resid':my_resid, 'Num':len(df_query), 'pKa_mean': round(df_query['pKa'].mean(),2), 'sasa_mean':round(df_query['sasa'].mean(),2), 'GO_Terms': GO_Terms})
                df_dictionary = pd.DataFrame([data])
                df_out = pd.concat([df_out, df_dictionary], ignore_index=True)
            
            df_out.to_csv(self.export_path, index = False)
        
        # export buttom for exporting the data within a region including GO Terms
        self.b_export = widgets.Button(
                    description='EXPORT TO CSV',
                    disabled=False,
                    button_style='info', # 'success', 'info', 'warning', 'danger' or ''
                    tooltip='Click me',
                    icon='check')
        self.b_export.on_click(call_back_buttom_export)
        
        # clear buttom for clearing the region
        def call_back_buttom(b):
            self.p.value, self.s.value = [1, 14], [0, 100]
            self.GO.options = (['Welcome'] + list(self.GO_dict.keys()))
            self.GO.value = 'Welcome'
            self.uni_clicked = ''
    
        self.b = widgets.Button(
                    description='RESET',
                    disabled=False,
                    button_style='info', # 'success', 'info', 'warning', 'danger' or ''
                    tooltip='Click me',
                    icon='check')
        
        self.b.on_click(call_back_buttom)
        
        # clear buttom for clearing the 3D visualisation
        self.b_2 = widgets.Button(
                        description='CLEAR',
                        disabled=False,
                        button_style='info', # 'success', 'info', 'warning', 'danger' or ''
                        tooltip='Click me',
                        icon='check')
        
        # buttom that once clicked will pop up the uniprot webpage
        def call_back_buttom_open_url(b_open_url):
            if self.uni_clicked == '':
                return
            url = 'https://www.uniprot.org/uniprot/' + self.uni_clicked
            try:
                webbrowser.open(url)
            except:
                print(f'access failed for {url}.')
                
        self.b_open_url = widgets.Button(
                        description='Go to Uniprot',
                        disabled=False,
                        button_style='info', # 'success', 'info', 'warning', 'danger' or ''
                        tooltip='Click me',
                        icon='check')
        
        self.b_open_url.on_click(call_back_buttom_open_url)
        
        # the very fundamental plot
        labels = ["UNIPROT: %s<br>resid: %i"%(self.df_concise["Uniprot Entry"].values[i], self.df_concise["Resid"].values[i]) for i in range(len(self.df_concise))]
        self.f = go.FigureWidget([go.Scatter(x=self.df_concise["sasa"], y=self.df_concise["pKa"],
                                mode='markers', name="aggregate", showlegend=False, opacity=0.75,
                                text = labels, hovertemplate='%{text}<br>SASA: %{x:.2f}<br>pKa: %{y:.2f}')])
        
        self.f.update_layout(
        xaxis_title="SASA (A2)",
        yaxis_title="pKa")
        
        self.f.update_xaxes(range=[0, 100])
        self.f.update_yaxes(range=[(int(self.df['pKa'].min())-1), (int(self.df['pKa'].max())+1)])
        self.f.update_xaxes(showspikes=True)
        self.f.update_yaxes(showspikes=True)
        
        scatter = self.f.data[0]
        colors = ['#7f7f7f'] * len(self.df_concise)
        scatter.marker.color = colors
        self.f.layout.hovermode = 'closest'
        
        # bar chart for displaying p values in enrichment analysis
        self.bar = go.FigureWidget()
        self.bar.update_layout(barmode='overlay',
                              xaxis_title='-log10(p value)',
                              yaxis_title='GO codes')
        
        
    # compute contingency table given a GO Term, the list of interest, and a reference list
    def get_contingency_table(self, GO_code, my_list, reference):
        
        BP_list = 0
        for uni in my_list:
            if uni in self.GO_dict[GO_code]:
                BP_list += 1
            
        BP_not_list = len(self.GO_dict[GO_code]) - BP_list
            
        not_BP_list = len(my_list) - BP_list
            
        not_BP_not_list = 0
        for uni in reference:
            if (uni not in my_list) and (uni not in self.GO_dict[GO_code]):
                not_BP_not_list += 1
            
        table = [[BP_list, BP_not_list], [not_BP_list, not_BP_not_list]]
        
        return table
    
    # enrichment analysis
    def enrichment_analysis(self, pka_range, sasa_range):
        
        # get all the uniprot codes inside the range and the reference uniprot code list
        pka_l, pka_u = pka_range[0], pka_range[1]
        sasa_l, sasa_u = sasa_range[0], sasa_range[1]
        selected_df = self.df_concise[(self.df_concise['pKa'] >= pka_l) & (self.df_concise['pKa'] <= pka_u)]
        selected_df = selected_df[(selected_df['sasa'] >= sasa_l) & (selected_df['sasa'] <= sasa_u)]
        my_list = selected_df['Uniprot Entry'].unique()
        reference = self.df_concise['Uniprot Entry'].unique()
        
        # get all the GO Terms in the background (that are associated with more than 4 uniprot codes)
        GO_bacgou = [code for code in list(self.GO_dict.keys()) if len(self.GO_dict[code]) >= 5]
        
        p_val_dict = {}
        
        # for each GO Term, compute a contingency table
        for GO_code in GO_bacgou:
            # compute contigency table
            table = self.get_contingency_table(GO_code, my_list, reference)
            # compute p values and store them into the dictionary
            oddsratio, pvalue = fisher_exact(table, alternative='greater')
            p_val_dict[GO_code] = pvalue
        
        # sort the dictionary based on p-values
        p_val_dict = dict(sorted(p_val_dict.items(), key = lambda x: x[1]))
        p_val_list = [x[1] for x in p_val_dict.items()]
        
        # adjust p-values using the BH method
        y=multipletests(pvals=p_val_list, alpha=0.05, method="fdr_bh")
        out_dict = {}
        
        if sum(y[0]) != 0: # if there is enrichment
            for i in range(len(y[0])):
                if y[0][i]: # get p-values below 0.05
                    go_code = list(p_val_dict.items())[i][0]
                    raw_p = p_val_list[i]
                    adj_p = y[1][i]
                    out_dict[go_code] = (raw_p, adj_p, self.get_contingency_table(go_code, my_list, reference))
        else:
            for i in range(len(y[0])): # if no enrichment at all, still output the p values
                go_code = list(p_val_dict.items())[i][0]
                raw_p = p_val_list[i]
                adj_p = y[1][i]
                out_dict[go_code] = (raw_p, adj_p, self.get_contingency_table(go_code, my_list, reference))
        
        df = pd.DataFrame(columns = ['GO ID', 'GO Term', 'raw p value', 'FDR',
                             'num in the region', 'num in the bkgd'])
        for data in out_dict.items():
            ID = data[0]
            Term = self.GO_decode_dict_reverse[ID]
    
            r_p = round(data[1][0],6)
            fdr = round(data[1][1],6)
    
            n1 = sum(data[1][2][0])
            d1 = sum(data[1][2][0]) + sum(data[1][2][1])                    
            bkgd = str(n1) + '/' + str(d1)
            n2 = data[1][2][0][0]
            d2 = data[1][2][1][0] + data[1][2][0][0]
            reg = str(n2) + '/' + str(d2)
    
            dt = {'GO ID':ID, 'GO Term':Term, 'raw p value':r_p, 'FDR':fdr,
                             'num in the bkgd':bkgd, 'num in the region':reg}
            df_dictionary = pd.DataFrame([dt])
            df = pd.concat([df, df_dictionary], ignore_index=True)
        
        return df
    
    # here we define a function that can easily add more data to the main data frame
    #def more_data(self, analysis_2):
        
    
    def advanced_plot(self, export_path = 'Regional_Data.csv'):
        
        self.export_path = os.path.join(self.outdir, export_path)
        
        def interact_slides(p,s):
            if len(self.f.data)>1:
                self.f.data = [self.f.data[0]]
    
            pka_l, pka_u = p[0], p[1]
            sasa_l, sasa_u = s[0], s[1]
            
            # dont do anything if we are at the initial state
            if ((pka_l, pka_u) == (1, 14.0)) and ((sasa_l, sasa_u) == (0, 100.0)):
                for i in range(len(self.bar.data)):
                    self.bar.data[i].visible = False # clear the bar chart
                return
            
            selected_df = self.df_concise[(self.df_concise['pKa'] >= pka_l) & (self.df_concise['pKa'] <= pka_u)]
            selected_df = selected_df[(selected_df['sasa'] >= sasa_l) & (selected_df['sasa'] <= sasa_u)]
            self.temp_df = selected_df
            
            # update GO Term dropdown options
            selected_uni_codes = selected_df['Uniprot Entry'].unique() # get unique uniprot codes
            new_options = list()
            for code, uni_list in self.GO_dict.items():
                for uni in selected_uni_codes:
                    if uni in uni_list:
                        new_options.append(code)
                        break
            new_options = [f'{code}: {self.GO_decode_dict_reverse[code]}' for code in new_options]
            self.GO.options = ['Welcome'] + sorted(new_options, key = lambda x: int(x.split(':')[0]))
            
            ###################################################################################################
            # enrichment analysis and update the barchart
            for i in range(len(self.bar.data)):
                self.bar.data[i].visible = False
            
            df_e = self.enrichment_analysis(p, s)
            
            labels = ['%s'%(self.GO_decode_dict_reverse[df_e['GO ID'].values[i]]) for i in range(len(df_e))]
            
            self.bar.add_trace(go.Bar(
                        y=df_e['GO ID'],
                        x=-np.log10(df_e['raw p value']),
                        name='raw p value',
                        orientation='h',
                        marker=dict(
                                color='rgba(246, 78, 139, 0.6)',
                        line=dict(color='rgba(246, 78, 139, 1.0)', width=3)),
                        text = labels,
                        hovertemplate = '%{x}<br>%{text}'))
            self.bar.add_trace(go.Bar(y=df_e['GO ID'],
                        x=-np.log10(df_e['FDR']),
                        name='FDR',
                        orientation='h',
                        marker=dict(
                                    color='rgba(58, 71, 80, 0.6)',
                                    line=dict(color='rgba(58, 71, 80, 1.0)', width=3)),
                                     text = labels,
                                     hovertemplate = '%{x}<br>%{text}'))
            self.bar.update_yaxes(autorange="reversed")
            ###################################################################################################
            
            label = 'pKa: %.1f-%.1f | SASA: %.1f-%.1f'%(pka_l, pka_u, sasa_l, sasa_u)
            labels = ["UNIPROT: %s<br>resid: %i"%(selected_df["Uniprot Entry"].values[i], selected_df["Resid"].values[i]) for i in range(len(selected_df))]
            self.f.add_scatter(x=selected_df["sasa"], y=selected_df["pKa"],
                    mode='markers', showlegend=False, name=label,
                    text = labels, hovertemplate='%{text}<br>SASA: %{x:.2f}<br>pKa: %{y:.2f}')
    
            scatter_1 = self.f.data[-1]
            scatter_1.marker.color = ['#1f77b4'] * len(selected_df)
            scatter_1.marker.size = [10] * len(selected_df)
        
        def interact_dropdown(GO,p,s):
            pka_l, pka_u = p[0], p[1]
            sasa_l, sasa_u = s[0], s[1]
    
            if GO == 'Welcome':
                # first case: having a region selected and have chosen a GO Term
                if len(self.f.data) == 3 and (((pka_l, pka_u) != (1, 14.0)) or ((sasa_l, sasa_u) != (0, 100.0))):
                    self.f.data = self.f.data[:-1]
                    return
        
                # second case: no region selected and have chosen a GO Term code
                elif ((len(self.f.data) == 2) or (len(self.f.data) == 3)) and ((pka_l, pka_u) == (1, 14.0)) and ((sasa_l, sasa_u) == (0, 100.0)):
                    self.f.data = [self.f.data[0]]
                    return
        
                # third case: a region and an option both selected and a point clicked
                elif len(self.f.data) == 4:
                    self.f.data = self.f.data[:2]
                    return
                
                # everything else    
                else:
                    return
            
            # no regional selection, but want to search for a GO Term code
            if (len(self.f.data) == 2) and ((pka_l, pka_u) == (1, 14.0)) and ((sasa_l, sasa_u) == (0, 100.0)):
                self.f.data = self.f.data[:-1]
    
            # there is a region selected and want to search for a GO Term code inside the region
            if len(self.f.data) == 3:
                self.f.data = self.f.data[:-1]
    
            # a region, an option, and a point all selected and want to update using the call_back function
            if len(self.f.data) == 4:
                self.f.data = self.f.data[:2]
    
            GO_code = GO.split(': ')[0]
            GO_name = self.GO_decode_dict_reverse[GO_code]
            
            selected_df = self.analysis.GO_Search_Term(df = self.df_concise, code = GO_code)
            
            # only want the points inside the region
            selected_df = selected_df[(selected_df['pKa'] >= pka_l) & (selected_df['pKa'] <= pka_u)]
            selected_df = selected_df[(selected_df['sasa'] >= sasa_l) & (selected_df['sasa'] <= sasa_u)]
            selected_df = selected_df.reset_index(drop=True) # this step is needed when attatching a call_back
            
            # store the dataframe for the call back function
            self.temp_df_2 = selected_df
            
            label = '%s: %s | pKa: %.1f-%.1f | SASA: %.1f-%.1f'%(GO_code, GO_name, pka_l, pka_u, sasa_l, sasa_u)
            labels = ["Uniprot Entry: %s"%(selected_df["Uniprot Entry"].values[i]) for i in range(len(selected_df))]
    
            self.f.add_scatter(x=selected_df["sasa"], y=selected_df["pKa"],
                    mode='markers', showlegend=False, name=label,
                    text = labels, hovertemplate='%{text}<br>SASA: %{x:.2f}<br>pKa: %{y:.2f}')
    
            scatter_2 = self.f.data[-1]
            scatter_2.marker.color = ['#d62728'] * len(selected_df)
            scatter_2.marker.size = [10] * len(selected_df)
            scatter_2.on_click(update_point_2)
        
        def interact_3D(PDB_box):
            if PDB_box == 'Welcome':
                return
            elif PDB_box == 'Cleared':
                return
            else:
                view = nv.show_file(PDB_box)
                view.clear_representations()
                view.add_representation('cartoon')
                view.add_licorice('LYS')
                display(view)
        
        def update_point_2(trace, points, selector):
            if len(points.point_inds) == 0:
                return
    
            idx = points.point_inds[0]
    
            #0000122: negative regulation of transcription by RNA polymerase II | pKa: 1.0-8.6 | SASA: 37.6-100.0
            pka_range = points.trace_name.split(' | ')[-2].split(': ')[-1].split('-')
            sasa_range = points.trace_name.split(' | ')[-1].split(': ')[-1].split('-')
    
            pka_l, pka_u = float(pka_range[0]), float(pka_range[1])
            sasa_l, sasa_u = float(sasa_range[0]), float(sasa_range[1])

            my_uniprot = self.temp_df_2.loc[idx, "Uniprot Entry"] # use the temp df defined in the data structure
            self.uni_clicked = my_uniprot
            my_resid = self.temp_df_2.loc[idx, "Resid"]
            my_label = "%s(%i)"%(my_uniprot, my_resid)
            df_query = self.df[(self.df['Uniprot Entry'] == my_uniprot) & (self.df['Resid'] == my_resid)]
            labels = ["PDB: %s"%(df_query["PDB Code"].values[i]) for i in range(len(df_query))]
            
            # a region, an option, and a point all selected
            if (len(self.f.data) == 4):
                self.f.data = self.f.data[:-1]
            
            # only an option and a point selected
            if (len(self.f.data) == 3) and ((pka_l, pka_u) == (1, 14.0)) and ((sasa_l, sasa_u) == (0, 100.0)):
                self.f.data = self.f.data[:-1]
    
            self.f.add_scatter(x=df_query["sasa"], y=df_query["pKa"],
                                    mode='markers', showlegend=False, name=my_label,
                                    text = labels, hovertemplate='%{text}<br>SASA: %{x:.2f}<br>pKa: %{y:.2f}')
            
            scatter_3 = self.f.data[-1]
            scatter_3.marker.color = ['#2ca02c'] * len(df_query)
            scatter_3.marker.size = [12] * len(df_query)
            scatter_3.on_click(update_point_3)
        
        def update_point_3(trace, points, selector):
            if len(points.point_inds) == 0:
                return
    
            idx = points.point_inds[0]
            # e.g. P19338(398)
            my_uniprot = points.trace_name.split('(')[0]
            my_resid = int(points.trace_name.split('(')[1].strip(')'))
            df_query = self.df[(self.df['Uniprot Entry'] == my_uniprot) & (self.df['Resid'] == my_resid)].reset_index(drop=True)
            pdb_file_path = df_query.loc[idx, 'PDB Code'] + '.pdb'
            self.PDB_box.value = pdb_file_path
        
        out_1 = widgets.interactive_output(interact_slides, {'p': self.p, 's':self.s})
        out_2 = widgets.interactive_output(interact_dropdown, {'GO': self.GO, 'p':self.p, 's':self.s})
        out_3 = widgets.interactive_output(interact_3D, {'PDB_box': self.PDB_box})
        
        # we have to define the callback for the second buttom here
        def call_back_buttom_2(b_2):
            out_3.clear_output()
            self.PDB_box.value = 'Cleared'
        self.b_2.on_click(call_back_buttom_2)
       
        # output format
        block0 = widgets.VBox([self.f, self.b_open_url, self.bar, out_3])
        block1 = widgets.HBox([self.b_export, widgets.VBox([self.p,self.s,out_1,out_2])])
        block2 = widgets.HBox([self.b,self.GO])
        block3 = widgets.HBox([self.b_2,self.PDB_box])
        
        self.plot = widgets.VBox([block0,block1,block2,block3])
        
        return widgets.VBox([block0,block1,block2,block3])                      
                       
