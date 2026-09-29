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
        self.outdir = outdir

        self.non_feat_cols = ['Uniprot_Entry', 'PDB_Code', 'Chain', 'Resid',
                         'Method', 'Resolution', 'PLDDT', 'class']
        self.features = [a for a in self.df_measures.columns if a not in self.non_feat_cols]

        # create the aggregated table for investigations
        self.A = Aggregation(df_measurements=df_measures,
                             aggregation_method='all',
                             features_to_include='all',
                             aev_red_method='pca')
        self.df_agg = self.A.aggregate_data()

        self.df_meas_stack = self._stack_data(self.df_measures, ['sasa', 0, 100], ['das', 0, 100])
        print(self.df_meas_stack)
        self.df_agg_stack = self._stack_data(self.df_agg, ['sasa', 0, 100], ['das', 0, 100])

        self.data_viewer_app = Dash(__name__)
        self._setup_html()


    def _stack_data(self, df, feat_a=['propka', 7.0, 11.5], feat_b=['das', 2, 40]):
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
        df_stack = selected_df.stack().reset_index()

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
            dcc.Tabs(id="tabs_scalar_analysis", value='tabs_scalar_analysis', children=[
                dcc.Tab(label='Scalar Measurements Analysis', value='tab_scalar_analysis'),
                dcc.Tab(label='Scalar Aggregation Analysis', value='tab_aggregation_analysis'),
            ]),
            html.Div(id='tabs_content_scalar_analysis')
        ])

        @callback(Output('tabs_content_scalar_analysis', 'children'),
                  Output('x_axis_feat_select', 'value'),
                  Output('scalar_feats_graph', 'figure'),
                  Input('tabs_scalar_analysis', 'value'),
                  Input('scalar_feats_graph', 'figure'),
                  Input('scalar_feats_graph', 'relayoutData'))
        def render_content(tab):
            if tab == 'tab_scalar_analysis':
                # investigate 2D and 3D plots of aggregated analysis
                return html.Div([
                    html.Div([
                        html.H3('Scalar Measurements Analysis'),
                        html.Div([
                            dcc.Dropdown(
                                self.features,
                                'depth',
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
                            dcc.Dropdown(
                                self.features,
                                'propka',
                                id='crossfilter-yaxis-column'
                            ),
                            dcc.RadioItems(
                                ['Linear', 'Log'],
                                'Linear',
                                id='crossfilter-yaxis-type',
                                labelStyle={'display': 'inline-block', 'marginTop': '5px'}
                            )
                        ],
                        style={'width': '32%', 'float': 'right', 'display': 'inline-block'}),

                        html.Div([
                            dcc.Dropdown(
                                self.features,
                                'average',
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
                            hoverData={'points': [{'customdata': 'Japan'}]}
                        )
                    ], style={'width': '49%', 'display': 'inline-block', 'padding': '0 20'}),
                    html.Div([
                        dcc.Graph(id='x-feat-hist'),
                        dcc.Graph(id='y-feat-hist'),
                    ], style={'display': 'inline-block', 'width': '49%'}),

                    html.Div(dcc.Slider(
                        self.df_measures['sasa'].min(),
                        self.df_measures['sasa'].max(),
                        step=None,
                        id='crossfilter-sasa--slider',
                        value=self.df_measures['sasa'].max(),
                        marks={str(sasa): str(sasa) for sasa in self.df_measures['sasa'].unique()}
                    ), style={'width': '49%', 'padding': '0px 20px 20px 20px'})
                ])
            elif tab == 'tab_aggregation_analysis':
                # investigate feature specific histogram and violins of different aggregation types
                return html.Div([
                    html.H3('Scalar Aggregation Analysis', style={'text-align': 'center'}),

                    dcc.Dropdown(id='x_axis_feat_select',
                                 options=self.features,
                                 optionHeight=25,
                                 multi=False,
                                 placeholder='Enter feature for x-axis...',
                                 clearable=True,
                                 value=['all'],
                                 style={'width': '100%'}
                                 ),

                    dcc.Graph(
                        id='scalar_feats_graph',
                        figure=self._scalar_feats_graph()
                    ),
                    dcc.Slider(
                        id='param_a_graph2_slider'
                    ),
                    html.Div([
                        dcc.Graph(id='tmp_graph')
                    ])
                ])

        @callback(
            Output('crossfilter-indicator-scatter', 'figure'),
            Input('crossfilter-xaxis-column', 'value'),
            Input('crossfilter-yaxis-column', 'value'),
            Input('crossfilter-xaxis-type', 'value'),
            Input('crossfilter-yaxis-type', 'value'),
            Input('crossfilter-agg-col', 'value'),
            Input('crossfilter-agg-type', 'value'),
            Input('crossfilter-sasa--slider', 'value'))
        def update_graph(xaxis_column_name, yaxis_column_name,
                        xaxis_type, yaxis_type,
                        sasa_value):
            dff = self.df_agg[self.df_agg['sasa_avg'] == sasa_value]

            fig = px.scatter(x=dff[dff['Feature'] == xaxis_column_name]['Value'],
                    y=dff[dff['Feature'] == yaxis_column_name]['Value'],
                    hover_name=dff[dff['Feature'] == yaxis_column_name][['Uniprot_Entry', 'Chain', 'Resid']]
                    )

            fig.update_traces(customdata=dff[dff['Feature'] == yaxis_column_name]['Uniprot_Entry'])

            fig.update_xaxes(title=xaxis_column_name, type='linear' if xaxis_type == 'Linear' else 'log')

            fig.update_yaxes(title=yaxis_column_name, type='linear' if yaxis_type == 'Linear' else 'log')

            fig.update_layout(margin={'l': 40, 'b': 40, 't': 10, 'r': 0}, hovermode='closest')

            return fig


        def create_feature_hist(dff, feature, title):

            fig = px.histogram(dff, x=feature)

            fig.update_traces(mode='lines+markers')

            fig.update_xaxes(showgrid=False)

            fig.update_yaxes(type='linear' if axis_type == 'Linear' else 'log')

            fig.add_annotation(x=0, y=0.85, xanchor='left', yanchor='bottom',
                            xref='paper', yref='paper', showarrow=False, align='left',
                            text=title)

            fig.update_layout(height=225, margin={'l': 20, 'b': 30, 'r': 10, 't': 10})

            return fig


        @callback(
            Output('x-feat-hist', 'figure'),
            Input('crossfilter-indicator-scatter', 'hoverData'),
            Input('crossfilter-xaxis-column', 'value'),
            Input('crossfilter-xaxis-type', 'value'))
        def update_x_hist(hoverData, xaxis_column_name, axis_type):
            country_name = hoverData['points'][0]['customdata']
            dff = self.df_measures[self.df_measures['Country Name'] == country_name]
            dff = dff[dff['Feature'] == xaxis_column_name]
            title = '<b>{}</b><br>{}'.format(country_name, xaxis_column_name)
            return create_feature_hist(dff, axis_type, title)


        @callback(
            Output('y-feat-hist', 'figure'),
            Input('crossfilter-indicator-scatter', 'hoverData'),
            Input('crossfilter-yaxis-column', 'value'),
            Input('crossfilter-yaxis-type', 'value'))
        def update_y_hist(hoverData, yaxis_column_name, axis_type):
            dff = self.df_measures[self.df_measures['Country Name'] == hoverData['points'][0]['customdata']]
            dff = dff[dff['Feature'] == yaxis_column_name]
            return create_feature_hist(dff, axis_type, yaxis_column_name)


if __name__ == '__main__':
    outdir = 'demo'
    df_measures = f'{outdir}{os.sep}measures.csv'
    V = Viewer(df_measures=df_measures,
                outdir=outdir)
    V.launch_viewer()


# https://dash.plotly.com/interactive-graphing
