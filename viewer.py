from ipywidgets import widgets
import numpy as np
import pandas as pd
import nglview as nv
import plotly.graph_objects as go
import os
import webbrowser

class Viewer(object):
    
    def __init__(self, analysis, outdir = 'result'):
        self.analysis = analysis
        self.df = analysis.df
  
        # some temp dataframes
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
        options = [f'{code}: {self.analysis.code_to_name[code]}' for code in list(self.analysis.GO_dict.keys())]
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
            columns = ['Uniprot_Entry', 'Resid', 'Num', 'pKa_mean', 'sasa_mean', 'GO_Terms']
            df_out = pd.DataFrame(columns = columns)
            for idx, row in df.iterrows():
                my_uniprot = row['Uniprot_Entry']
                my_resid = row['Resid']
                df_query = self.df[(self.df['Uniprot_Entry'] == my_uniprot) & (self.df['Resid'] == my_resid)]
                GO_Terms = self.analysis.GO_search_protein(my_uniprot)
                data = ({'Uniprot_Entry':my_uniprot, 'Resid':my_resid, 'Num':len(df_query), 'pKa_mean': round(df_query['pKa'].mean(),2), 'sasa_mean':round(df_query['sasa'].mean(),2), 'GO_Terms': GO_Terms})
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
            self.GO.options = (['Welcome'] + list(self.analysis.GO_dict.keys()))
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
        labels = ["UNIPROT: %s<br>resid: %i"%(self.analysis.df_sub["Uniprot_Entry"].values[i], self.analysis.df_sub["Resid"].values[i]) for i in range(len(self.analysis.df_sub))]
        self.f = go.FigureWidget([go.Scatter(x=self.analysis.df_sub["sasa"], y=self.analysis.df_sub["pKa"],
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
        colors = ['#7f7f7f'] * len(self.analysis.df_sub)
        scatter.marker.color = colors
        self.f.layout.hovermode = 'closest'
        
        # bar chart for displaying p values in enrichment analysis
        self.bar = go.FigureWidget()
        self.bar.update_layout(barmode='overlay',
                              xaxis_title='-log10(p value)',
                              yaxis_title='GO codes')
        
        

    
    def advanced_plot(self, export_path = 'Regional_Data.csv', cutoff=0):
        
        self.export_path = os.path.join(self.outdir, export_path)
        
        def interact_slides(p, s):
            if len(self.f.data)>1:
                self.f.data = [self.f.data[0]]
    
            pka_l, pka_u = p[0], p[1]
            sasa_l, sasa_u = s[0], s[1]
            
            # dont do anything if we are at the initial state
            if ((pka_l, pka_u) == (1, 14.0)) and ((sasa_l, sasa_u) == (0, 100.0)):
                for i in range(len(self.bar.data)):
                    self.bar.data[i].visible = False # clear the bar chart
                return
            
            selected_df = self.analysis.df_sub[(self.analysis.df_sub['pKa'] >= pka_l) & (self.analysis.df_sub['pKa'] <= pka_u)]
            selected_df = selected_df[(selected_df['sasa'] >= sasa_l) & (selected_df['sasa'] <= sasa_u)]
            self.temp_df = selected_df
            
            # update GO Term dropdown options
            selected_uni_codes = selected_df['Uniprot_Entry'].unique() # get unique uniprot codes
            new_options = list()
            for code, uni_list in self.analysis.GO_dict.items():
                for uni in selected_uni_codes:
                    if uni in uni_list:
                        new_options.append(code)
                        break
            new_options = [f'{code}: {self.analysis.code_to_name[code]}' for code in new_options]
            self.GO.options = ['Welcome'] + sorted(new_options, key = lambda x: int(x.split(':')[0]))
            
            ###################################################################################################
            # enrichment analysis and update the barchart
            for i in range(len(self.bar.data)):
                self.bar.data[i].visible = False
            
            df_e = self.analysis.enrichment_analysis(p, s)
            
            # remove lines containing 0 examples in the region
            test = [int(v.split("/")[0])>cutoff for v in df_e["num in the region"].values]
            df_e = df_e[test]
            
            labels = ['%s'%(self.analysis.code_to_name[df_e['GO ID'].values[i]]) for i in range(len(df_e))]
            
            self.bar.add_trace(go.Bar(
                        y=df_e['GO ID'],
                        x=-np.log10(df_e['raw p value']),
                        name='raw p-value',
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
            labels = ["UNIPROT: %s<br>resid: %i"%(selected_df["Uniprot_Entry"].values[i], selected_df["Resid"].values[i]) for i in range(len(selected_df))]
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
            GO_name = self.analysis.code_to_name[GO_code]
            
            selected_df = self.analysis.GO_search_term(df = self.analysis.df_sub, code = GO_code)
            
            # only want the points inside the region
            selected_df = selected_df[(selected_df['pKa'] >= pka_l) & (selected_df['pKa'] <= pka_u)]
            selected_df = selected_df[(selected_df['sasa'] >= sasa_l) & (selected_df['sasa'] <= sasa_u)]
            selected_df = selected_df.reset_index(drop=True) # this step is needed when attatching a call_back
            
            # store the dataframe for the call back function
            self.temp_df_2 = selected_df
            
            label = '%s: %s | pKa: %.1f-%.1f | SASA: %.1f-%.1f'%(GO_code, GO_name, pka_l, pka_u, sasa_l, sasa_u)
            labels = ["Uniprot_Entry: %s"%(selected_df["Uniprot_Entry"].values[i]) for i in range(len(selected_df))]
    
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
                
                try:
                    display(view)
                except:
                    pass
        
        
        def update_point_2(trace, points, selector):
            if len(points.point_inds) == 0:
                return
    
            idx = points.point_inds[0]
    
            #0000122: negative regulation of transcription by RNA polymerase II | pKa: 1.0-8.6 | SASA: 37.6-100.0
            pka_range = points.trace_name.split(' | ')[-2].split(': ')[-1].split('-')
            sasa_range = points.trace_name.split(' | ')[-1].split(': ')[-1].split('-')
    
            pka_l, pka_u = float(pka_range[0]), float(pka_range[1])
            sasa_l, sasa_u = float(sasa_range[0]), float(sasa_range[1])

            my_uniprot = self.temp_df_2.loc[idx, "Uniprot_Entry"] # use the temp df defined in the data structure
            self.uni_clicked = my_uniprot
            my_resid = self.temp_df_2.loc[idx, "Resid"]
            my_label = "%s(%i)"%(my_uniprot, my_resid)
            df_query = self.df[(self.df['Uniprot_Entry'] == my_uniprot) & (self.df['Resid'] == my_resid)]
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
            df_query = self.df[(self.df['Uniprot_Entry'] == my_uniprot) & (self.df['Resid'] == my_resid)].reset_index(drop=True)
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
                       
