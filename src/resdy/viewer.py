import numpy as np
import pandas as pd
import nglview as nv
import plotly.express as px
import os
import webbrowser
from dash import Dash, dcc, html, Input, Output, callback
from .aggregation import Aggregation


class Viewer(object):
    '''
    Class to investigate the measurements data produced from measures using the plotly
    dash functionality
    '''

    def __init__(self, df_measures, outdir = 'result'):
        if isinstance(df_measures, str):
            self.df_measures = pd.read_csv(df_measures)
        else:
            self.df_measures = df_measures

        self.df_measures.columns.name = 'Feature'
        self.outdir = outdir

        self.non_feat_cols = ['Uniprot_Entry', 'PDB_Code', 'Chain', 'Resid',
                         'Method', 'Resolution', 'PLDDT', 'class']
        self.features = [a for a in self.df_measures.columns if a not in self.non_feat_cols]

        self.vector_features = ['aev', 'evolution']
        self.df_scalar_measures = self.df_measures.copy()
        for feat in self.vector_features:
            if feat in self.df_scalar_measures.columns:
                self.df_scalar_measures = self.df_scalar_measures.drop(columns=feat)
                print(f'>> This analysis can only investigate scalar features, removing {feat} from dataframe')
        self.features = [a for a in self.features if a not in self.vector_features]

        self.aggregation_types = {'Average': 'avg',
                                  'Standard Deviation': 'sd',
                                  'Median': 'med',
                                  'Minimum': 'min',
                                  'Maximum': 'max',
                                  'Range': 'range',
                                  'Random': 'rand'}

        self.feature_labels = {'pKa (PROPKA3)': 'propka',
                               'pKa (pKaANI)': 'pkaani',
                               'Depth': 'depth',
                               'Solvent Accessible Surface Area': 'sasa',
                               'Dynamically Accessible Area': 'das',
                               'Sequence Charge': 'seqcharge',
                               'Curvature': 'curvature',
                               'Writhing': 'writhing',
                               'Phi': 'phi',
                               'Psi': 'psi',
                               'Torsion': 'torsion',
                               'Arc Length': 'arc_length',
                               'Flexibility': 'flexibility',
                               'Root Mean Square Fluctuation': 'rmsf',
                               'Atomic Environment Vector': 'aev',
                               '15N NMR Shift (LEGOLAS)': 'legolas',
                               'Frustration': 'frustration',
                               'Density': 'density'}

        self.feature_labels_reverse = {v: k for k, v in self.feature_labels.items()}
        self.feature_labels_list = [self.feature_labels_reverse[a] for a in self.features]

        # create the aggregated table for investigations
        self.A = Aggregation(df_measurements=self.df_scalar_measures,
                             aggregation_method='all',
                             features_to_include='all',
                             aev_red_method='pca')
        self.df_agg = self.A.aggregate_data()
        self.df_agg.columns.name = 'Feature'

        self.df_meas_stack = self._stack_data(self.df_scalar_measures, ['depth', 0, 40], ['seqcharge', -5.0, 5.0])
        print(self.df_meas_stack)
        self.df_agg_stack = self._stack_data(self.df_agg, ['depth', 0, 40], ['seqcharge', -5.0, 5.0])
        print(self.df_agg_stack)

        self.data_viewer_app = Dash(__name__)
        self._setup_html()


    def _stack_data(self, df, feat_a=['depth', 0, 40], feat_b=['seqcharge', -5.0, 5.0]):
        '''
        Transform the data such that it allows for easy access into the viewer plotting sets
        
        Parameters
        ----------
        :param df: Dataframe containing all the measurements data for analysis
        :type df: pandas Dataframe
        :param feat_a: List containing 3 elements, first the feature column name, then the
            floats of the lower and upper bounds for the feature from the sliders
        :type feat_a: list
        :param feat_b: List containing 3 elements, first the feature column name, then the
            floats of the lower and upper bounds for the feature from the sliders
        :type feat_b: list
        '''
        feat_cols_df = [a for a in df.columns if a not in self.non_feat_cols]
        non_feat_cols_df = [a for a in df.columns if a in self.non_feat_cols]

        feat_one, feat_one_low, feat_one_upper = str(feat_a[0]), float(feat_a[1]), float(feat_a[2])
        feat_two, feat_two_low, feat_two_upper = str(feat_b[0]), float(feat_b[1]), float(feat_b[2])

        def _subset_data(df, feat, low, upper):
            matching_feats = [a for a in df.columns if feat in a]
            for feat in matching_feats:
                df = df[(df[feat] >= low) & (df[feat] <= upper)]
            return df

        selected_df = _subset_data(df, feat_one, feat_one_low, feat_one_upper)
        selected_df = _subset_data(df, feat_two, feat_two_low, feat_two_upper)

        selected_df = selected_df.set_index(non_feat_cols_df)
        df_stack = selected_df.stack()
        df_stack.name = 'Value'
        df_stack = df_stack.reset_index()

        return df_stack


    def launch_viewer(self):
        '''
        Launch a plotly web browser window which can then be used 
        '''
        self.data_viewer_app.run(debug=True)


    def _setup_html(self):
        '''
        Currently house all the setup of the plotly work here, may split up if possible to do so too.
        '''
        self.data_viewer_app.layout = html.Div([
            html.H1('RESDY Measurements Analysis', style={'text-align': 'center'}),
            dcc.Tabs(id="tabs_analysis", value='tab_scalar_analysis', children=[
                dcc.Tab(label='Scalar Measurements Analysis', value='tab_scalar_analysis'),
                dcc.Tab(label='Feature Histogram Analysis', value='tab_histogram_analysis'),
                dcc.Tab(label='GO Ontology Analysis', value='tab_GOterm_analysis'),
            ]),
            html.Div(id='tabs_content_analysis')
        ])

        @callback(Output('tabs_content_analysis', 'children'),
                  Input('tabs_analysis', 'value'))
        def render_content(tab):
            if tab == 'tab_scalar_analysis':
                # investigate 2D and 3D plots of aggregated analysis
                return html.Div([
                    html.Div([
                        html.Div([
                            html.H4('X-Axis Feature', style={'text-align': 'center'}),
                            dcc.Dropdown(
                                self.feature_labels_list,
                                self.feature_labels_list[0],
                                id='crossfilter-xaxis-column',
                            ),
                            dcc.RadioItems(
                                ['Linear', 'Log'],
                                'Linear',
                                id='crossfilter-xaxis-type',
                                labelStyle={'display': 'inline-block', 'marginTop': '5px'}
                            )
                        ],
                        style={'width': '32%', 'display': 'inline-block'}),

                        html.Div([
                            html.H4('Y-Axis Feature', style={'text-align': 'center'}),
                            dcc.Dropdown(
                                self.feature_labels_list,
                                self.feature_labels_list[1],
                                id='crossfilter-yaxis-column'
                            ),
                            dcc.RadioItems(
                                ['Linear', 'Log'],
                                'Linear',
                                id='crossfilter-yaxis-type',
                                labelStyle={'display': 'inline-block', 'marginTop': '5px'}
                            )
                        ],
                        style={'width': '32%', 'float': 'center', 'display': 'inline-block'}),

                        html.Div([
                            html.H4('Measurements Aggregation Type', style={'text-align': 'center'}),
                            dcc.Dropdown(
                                list(self.aggregation_types.keys()),
                                'Average',
                                id='crossfilter-agg-col'
                            ),
                            dcc.RadioItems(
                                ['On', 'Off'],
                                'On',
                                id='crossfilter-agg-type',
                                labelStyle={'display': 'inline-block', 'marginTop': '5px'}
                            )
                        ],
                        style={'width': '32%', 'float': 'right', 'display': 'inline-block'})
                    ], style={
                        'padding': '10px 5px'
                    }),

                    html.Div([
                        dcc.Graph(
                            id='crossfilter-indicator-scatter',
                            hoverData={'points': [{'customdata': self.df_measures['Uniprot_Entry'].iloc[0]}]}
                        )
                    ], style={'width': '49%', 'display': 'inline-block', 'padding': '0 20'}),
                    html.Div([
                        dcc.Graph(id='x-feat-hist'),
                        dcc.Graph(id='y-feat-hist'),
                    ], style={'display': 'inline-block', 'width': '49%'}),

                    html.Div([
                        html.H4('SASA Slider'),
                        dcc.Slider(
                        self.df_scalar_measures['sasa'].min(),
                        self.df_scalar_measures['sasa'].max(),
                        step=None,
                        id='crossfilter-sasa--slider',
                        value=self.df_scalar_measures['sasa'].max(),
                        marks=list(range(int(self.df_scalar_measures['sasa'].min()), int(self.df_scalar_measures['sasa'].max()), int(((self.df_scalar_measures['sasa'].max() - self.df_scalar_measures['sasa'].min())/ 5))))),
                    ], style={'width': '49%', 'padding': '0px 20px 20px 20px'})
                ])
            elif tab == 'tab_histogram_analysis':
                # investigate feature specific histogram and violins of different aggregation types
                return html.Div([
                    html.Div([
                        html.Div([
                            html.H3('Feature', style={'text-align': 'center'}),
                            dcc.Dropdown(
                                self.feature_labels_list,
                                self.feature_labels_list[0],
                                id='crossfilter-feature-hist',
                            ),
                        ],
                        style={'width': '49%', 'display': 'inline-block'}),

                        html.Div([
                            html.H3('Feature Aggregation Type', style={'text-align': 'center'}),
                            dcc.Dropdown(
                                ['None'] + list(self.aggregation_types.keys()),
                                'None',
                                id='crossfilter-agg-type-hist'
                            ),
                        ],
                        style={'width': '49%', 'horizontal-align': 'right', 'display': 'inline-block'}),
                    ]),

                    html.Div([
                        html.Div([
                            html.H4('Uniprot Entry', style={'text-align': 'center'}),
                            dcc.Dropdown(
                                ['All'] + list(self.df_measures['Uniprot_Entry'].unique()),
                                'All',
                                id='crossfilter-uniprot-hist',
                            ),
                        ],
                        style={'width': '32%', 'display': 'inline-block', 'vertical-align': 'top'}),

                        html.Div([
                            html.H4('Chain', style={'text-align': 'center'}),
                            dcc.Dropdown(
                                ['All'] + list(self.df_measures['Chain'].unique()),
                                'All',
                                id='crossfilter-chain-hist'
                            ),
                            dcc.RadioItems(
                                ['Together', 'Seperate'],
                                'Together',
                                id='crossfilter-chain-split-type',
                                labelStyle={'display': 'inline-block', 'marginTop': '5px'}
                            )
                        ],
                        style={'width': '32%', 'horizontal-align': 'center', 'display': 'inline-block'}),

                        html.Div([
                            html.H4('Resid', style={'text-align': 'center'}),
                            dcc.Dropdown(
                                ['All'] + list(self.df_measures['Resid'].unique()),
                                'All',
                                id='crossfilter-resid-hist'
                            ),
                            dcc.RadioItems(
                                ['Together', 'Seperate'],
                                'Together',
                                id='crossfilter-resid-split-type',
                                labelStyle={'display': 'inline-block', 'marginTop': '5px'}
                            )
                        ],
                        style={'width': '32%', 'float': 'right', 'display': 'inline-block'})
                    ], style={
                        'padding': '10px 5px'
                    }),
                    html.Div([
                        html.H4('Any potential warnings for plotting non aggregated feature data')
                    ], style={'width': '100%', 'display': 'inline-block', 'padding': '0 20'}),
                    html.Div([
                        dcc.Graph(
                            id='histogram-main'
                        )
                    ], style={'width': '100%', 'display': 'inline-block', 'padding': '0 20'}),
                ])
            
            elif tab == 'tab_GOterm_analysis':
                return html.Div([
                    html.H2('GO Term Analysis'),
                ])

        @callback(
            Output('crossfilter-indicator-scatter', 'figure'),
            Output('crossfilter-xaxis-column', 'value'),
            Output('crossfilter-yaxis-column', 'value'),
            Output('crossfilter-agg-col', 'value'),
            Input('crossfilter-xaxis-column', 'value'),
            Input('crossfilter-yaxis-column', 'value'),
            Input('crossfilter-xaxis-type', 'value'),
            Input('crossfilter-yaxis-type', 'value'),
            Input('crossfilter-agg-col', 'value'),
            Input('crossfilter-agg-type', 'value'),
            #Input('crossfilter-sasa--slider', 'value')
            )
        def update_graph(xaxis_column_name, yaxis_column_name,
                        xaxis_type, yaxis_type, agg_type, agg_on):
            if xaxis_column_name is None:
                xaxis_column_name = self.feature_labels_list[0]
            if yaxis_column_name is None:
                yaxis_column_name = self.feature_labels_list[1]
            if agg_type is None:
                agg_type = 'Average'
            #dff = self.df_agg[self.df_agg['sasa_avg'] == sasa_value]
            dff = self.df_agg

            xaxis = f'{self.feature_labels[xaxis_column_name]}_{self.aggregation_types[agg_type]}'
            yaxis = f'{self.feature_labels[yaxis_column_name]}_{self.aggregation_types[agg_type]}'

            fig = px.scatter(dff,
                             x=xaxis,
                    y=yaxis
                    )
            
            '''
            dff = self.df_agg_stack

            xaxis = f'{self.feature_labels[xaxis_column_name]}_{self.aggregation_types[agg_type]}'
            yaxis = f'{self.feature_labels[yaxis_column_name]}_{self.aggregation_types[agg_type]}'
            print(xaxis, yaxis)

            fig = px.scatter(x=dff.loc[dff['Feature'] == xaxis, 'Value'],
                    y=dff.loc[dff['Feature'] == yaxis, 'Value'],
                    hover_name=dff[dff['Feature'] == yaxis][['Uniprot_Entry', 'Resid']]
                    )
            '''
            

            fig.update_traces(customdata=dff['Uniprot_Entry'])

            fig.update_xaxes(title=xaxis_column_name, type='linear' if xaxis_type == 'Linear' else 'log')

            fig.update_yaxes(title=yaxis_column_name, type='linear' if yaxis_type == 'Linear' else 'log')

            fig.update_layout(margin={'l': 40, 'b': 40, 't': 10, 'r': 0}, hovermode='closest')

            return fig, xaxis_column_name, yaxis_column_name, agg_type


        def create_feature_hist(dff, feature, title, axis_type):

            fig = px.histogram(dff, x=feature)

            #fig.update_traces(mode='lines+markers')

            fig.update_xaxes(showgrid=False)

            fig.update_yaxes(type='linear' if axis_type == 'Linear' else 'log')

            fig.add_annotation(x=0, y=0.85, xanchor='left', yanchor='bottom',
                            xref='paper', yref='paper', showarrow=False, align='left',
                            text=title)

            fig.update_layout(height=225, margin={'l': 20, 'b': 30, 'r': 10, 't': 10},
                              xaxis_title_text=self.feature_labels_reverse[feature], yaxis_title_text='Count')

            return fig


        @callback(
            Output('x-feat-hist', 'figure'),
            Input('crossfilter-indicator-scatter', 'hoverData'),
            Input('crossfilter-xaxis-column', 'value'),
            Input('crossfilter-xaxis-type', 'value'))
        def update_x_hist(hoverData, xaxis_column_name, axis_type):
            if xaxis_column_name is None:
                xaxis_column_name = self.feature_labels_list[0]
            xaxis = self.feature_labels[xaxis_column_name]
            uniprot_entry = hoverData['points'][0]['customdata']
            dff = self.df_scalar_measures[self.df_scalar_measures['Uniprot_Entry'] == uniprot_entry]
            dff = dff[['Uniprot_Entry', 'Resid'] + [xaxis]]
            title = '<b>{}</b><br>{}'.format(uniprot_entry, xaxis_column_name)
            return create_feature_hist(dff, xaxis, title, axis_type)


        @callback(
            Output('y-feat-hist', 'figure'),
            Input('crossfilter-indicator-scatter', 'hoverData'),
            Input('crossfilter-yaxis-column', 'value'),
            Input('crossfilter-yaxis-type', 'value'))
        def update_y_hist(hoverData, yaxis_column_name, axis_type):
            if yaxis_column_name is None:
                yaxis_column_name = self.feature_labels_list[1]
            yaxis = self.feature_labels[yaxis_column_name]
            uniprot_entry = hoverData['points'][0]['customdata']
            dff = self.df_scalar_measures[self.df_scalar_measures['Uniprot_Entry'] == uniprot_entry]
            dff = dff[['Uniprot_Entry', 'Resid'] + [yaxis]]
            title = '<b>{}</b><br>{}'.format(uniprot_entry, yaxis_column_name)
            return create_feature_hist(dff, yaxis, title, axis_type)


        # tab 2 functions
        def create_main_hist(dff, feature, title):

            fig = px.histogram(dff, x=feature, marginal='rug')

            #fig.update_traces(mode='lines+markers')

            fig.update_xaxes(showgrid=False, )
            fig.update_yaxes(showgrid=False)

            '''
            fig.add_annotation(x=0, y=0.85, xanchor='left', yanchor='bottom',
                            xref='paper', yref='paper', showarrow=False, align='left',
                            text=title)
            '''

            fig.update_layout(height=225, margin={'l': 20, 'b': 30, 'r': 10, 't': 10},
                              xaxis_title_text=title, yaxis_title_text='Count')
            
            '''
            fig.update_layout(
                title_text='Sampled Results', # title of plot
                xaxis_title_text='Value', # xaxis label
                yaxis_title_text='Count', # yaxis label
                bargap=0.2, # gap between bars of adjacent location coordinates
                bargroupgap=0.1 # gap between bars of the same location coordinates
            )
            '''

            return fig


        @callback(
            Output('histogram-main', 'figure'),
            Output('crossfilter-feature-hist', 'value'),
            Output('crossfilter-agg-type-hist', 'value'),
            Output('crossfilter-uniprot-hist', 'value'),
            Output('crossfilter-chain-hist', 'value'),
            Output('crossfilter-chain-hist', 'options'),
            Output('crossfilter-resid-hist', 'value'),
            Output('crossfilter-resid-hist', 'options'),
            Input('crossfilter-feature-hist', 'value'),
            Input('crossfilter-agg-type-hist', 'value'),
            Input('crossfilter-uniprot-hist', 'value'),
            Input('crossfilter-chain-hist', 'value'),
            Input('crossfilter-chain-split-type', 'value'),
            Input('crossfilter-resid-hist', 'value'),
            Input('crossfilter-resid-split-type', 'value'))
        def update_main_hist(feature, agg_type, uniprot,
                             chain, chain_split, resid, resid_split):
            if feature is None:
                feature = self.feature_labels_reverse[self.features[0]]
            if agg_type is None:
                agg_type = 'None'
            if uniprot is None:
                uniprot = 'All'
            if chain is None:
                chain = 'All' 
            if resid is None:
                resid = 'All'
            
            pot_chains = ['All'] + list(self.df_measures['Chain'].unique())
            pot_resids = ['All'] + list(self.df_measures['Resid'].unique())

            if agg_type == 'None':
                feat = self.feature_labels[feature]
                cols_remain = [a for a in self.non_feat_cols if a in self.df_scalar_measures.columns] + [feat]
                dff = self.df_scalar_measures[cols_remain]
                if uniprot != 'All':
                    dff = dff[dff['Uniprot_Entry'] == uniprot]
                    pot_chains = list(dff['Chain'].unique())
                    pot_resids = list(dff['Resid'].unique())
                    if chain != 'All' and chain in list(dff['Chain']):
                        dff = dff[dff['Chain'] == chain]
                    if resid != 'All' and resid in list(dff['Resid']):
                        dff = dff[dff['Resid'] == resid]
                return create_main_hist(dff, feat, feature), feature, agg_type, uniprot, chain, pot_chains, resid, pot_resids

            else:
                feat = self.feature_labels[feature]
                agg = self.aggregation_types[agg_type]
                feat_col = f'{feat}_{agg}'
                cols_remain = [a for a in self.non_feat_cols if a in self.df_agg.columns] + [feat_col]
                dff = self.df_agg[cols_remain]
                if uniprot != 'All':
                    dff = dff[dff['Uniprot_Entry'] == uniprot]
                    pot_chains = list(dff['Chain'].unique())
                    pot_resids = list(dff['Resid'].unique())
                    if chain != 'All' and chain in list(dff['Chain']):
                        dff = dff[dff['Chain'] == chain]
                    if resid != 'All' and resid in list(dff['Resid']):
                        dff = dff[dff['Resid'] == resid]
                return create_main_hist(dff, feat_col, feature), feature, agg_type, uniprot, chain, pot_chains, resid, pot_resids



if __name__ == '__main__':
    outdir = 'demo'
    df_measures = f'{outdir}{os.sep}measures.csv'
    V = Viewer(df_measures=df_measures,
                outdir=outdir)
    V.launch_viewer()

