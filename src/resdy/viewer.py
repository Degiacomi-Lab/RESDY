import os
import io
import numpy as np
import pandas as pd
import plotly.express as px
from dash import Dash, dcc, html, Input, Output, ctx, callback, no_update
from contextlib import redirect_stdout
from .aggregation import Aggregation
from .helper import METADATA_COLUMNS
from .analysis import Analysis


class Viewer(object):
    '''
    Class to investigate the measurements data produced from measures using the plotly
    dash functionality
    '''

    def __init__(self, df_measures, outdir = 'result', render_mode='auto',
                 suppress_unimportant_output=True):
        '''
        Initialise the viewer module for the RESDY package.

        :param df_measures: Dataframe of measurements created through measure.py to
            view in the viewer window.
        :type df_measures: Pandas Dataframe
        :param outdir: Location in which the measurements are stored and the place you
            would like any graphs saved to go to.
        :type outdir: string
        :param render_mode: Rendering mode of the graphs within the viewer window. Default
            is auto which allows plotly to choose the best option between webgl and svg depending
            on the size of the dataset in the plot. Can either be set to 'svg' or 'webgl' if user
            wants more custom control.
        :type str
        :param suppress_unimportant_output: Toggleable option to 
        '''
        if isinstance(df_measures, str):
            self.df_measures = pd.read_csv(df_measures)
        else:
            self.df_measures = df_measures

        self.df_measures.columns.name = 'Feature'
        self.outdir = outdir
        self.render_mode = render_mode
        if self.render_mode not in ['auto', 'svg', 'webgl']:
            print(f'>> Given render mode is not one accepted, setting render_mode to auto')
            self.render_mode = 'auto'
        self.suppress_unimportant_output = suppress_unimportant_output

        self.non_feat_cols = list(METADATA_COLUMNS)
        self.features = [a for a in self.df_measures.columns if a not in self.non_feat_cols]

        self.vector_features = ['aev', 'evolution']
        self.df_scalar_measures = self.df_measures.copy()
        self.df_scalar_measures['Chain'] = self.df_scalar_measures['Chain'].astype(str).str.strip()
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

        self.axis_labels = {'depth': 'Depth (Å)',
                        'sasa': 'Solvent Accessible Surface Area (Å\u00b2)',
                        'propka': 'pKa (PROPKA3)',
                        'pkaani': 'pKa (pKaANI)',
                        'flexibility': 'Flexibility (B-factor) (Å\u00b2)',
                        'curvature': 'Curvature',
                        'arc_length': 'Arc Length',
                        'das': 'Dynamically Accessible Surface Area',
                        'legolas': '15N NMR Backbone Shift (ppm) (LEGOLAS)',
                        'phi': 'Phi',
                        'psi': 'Psi',
                        'seqcharge': 'Sequence Charge',
                        'torsion': 'Torsion Angle',
                        'writhing': 'Writhing',
                        'frustration': 'Frustration',
                        'rmsf': 'Root Mean Square Fluctuation',
                        'density': 'Density',
                        'aev': 'Atomic Environment Vector'}

        self.feature_labels_reverse = {v: k for k, v in self.feature_labels.items()}
        self.feature_labels_list = [self.feature_labels_reverse[a] for a in self.features]

        self.colour_cols = [a for a in ['PLDDT', 'Method', 'Source', 'class'] if a in self.df_scalar_measures.columns]

        with redirect_stdout(io.StringIO()) as f:
            self.A = Aggregation(df_measurements=self.df_scalar_measures,
                                aggregation_method='all',
                                features_to_include='all',
                                aev_red_method='pca')
            self.df_agg = self.A.aggregate_data()
            self.df_agg.columns.name = 'Feature'
            self.res_key = [a for a in ['Uniprot_Entry', 'Chain', 'Resid']
                            if a in self.df_agg.columns]
            
            self.Analysis = Analysis(df=self.df_scalar_measures,
                                    outdir=self.outdir,
                                    features_to_analyse=self.features)
            self.Analysis.GO_get_data()
            self.GO_codes = {k: f'{k}: {v}' for k,v in self.Analysis.code_to_name.items()}
            self.GO_codes_reverse = {v:k for k,v in self.GO_codes.items()}
            self.GO_dict = self.Analysis.GO_dict
            self.GO_list = list(self.GO_codes.keys())
            self.GO_display_names = list(self.GO_codes.values())

        if not self.suppress_unimportant_output:
            for line in f:
                print(line)

        self.data_viewer_app = Dash(__name__, suppress_callback_exceptions=True)
        self._setup_html()


    def launch_viewer(self):
        '''
        Launch a plotly web browser window which can then be used for analysis of the
        measurements data. Will launch on a local address which can be opened in a browser
        or preview in many code editors.
        '''
        self.data_viewer_app.run(debug=True)


    def _setup_html(self):
        '''
        Houses all the html layout for the viewer dash app. Operates on a tab based system.
        '''
        self.data_viewer_app.layout = html.Div([
            html.H1('RESDY Measurements Analysis', style={'text-align': 'center'}),
            dcc.Tabs(id="tabs_analysis", value='tab_2d_scalar_analysis', children=[
                dcc.Tab(label='2D Scalar Measurements Analysis', value='tab_2d_scalar_analysis'),
                dcc.Tab(label='3D Scalar Measurements Analysis', value='tab_3d_scalar_analysis'),
                dcc.Tab(label='Feature Histogram Analysis', value='tab_histogram_analysis')
            ]),
            html.Div(id='tabs_content_analysis')
        ])

        @callback(Output('tabs_content_analysis', 'children'),
                  Input('tabs_analysis', 'value'))
        def render_content(tab):
            if tab == 'tab_2d_scalar_analysis':
                return html.Div([
                    html.Div([
                        html.Div([
                            html.H4('X-Axis Feature', style={'text-align': 'center'}),

                            dcc.Dropdown(
                                self.feature_labels_list,
                                self.feature_labels_list[0],
                                id='2d-axis-column',
                            ),

                            dcc.RadioItems(
                                ['Linear', 'Log'],
                                'Linear',
                                id='2d-xaxis-type',
                                labelStyle={'display': 'inline-block', 'marginTop': '5px'}
                            ),

                            html.Div([
                                html.H4('pKa (PROPKA3) Range',
                                        style={'text-align': 'center'},
                                        id='xaxis-feature-slider-header'),

                                dcc.RangeSlider(
                                min=round(self.df_scalar_measures['propka'].min(), 1),
                                max=round(self.df_scalar_measures['propka'].max(), 1),
                                step=0.1,
                                id='2d-xaxis-slider',
                                value=[round(self.df_scalar_measures['propka'].min(), 1),
                                       round(self.df_scalar_measures['propka'].max(), 1)],
                                marks=(int(((self.df_scalar_measures['propka'].max() -
                                             self.df_scalar_measures['propka'].min())/ 7)) or 1)),
                            ], style={'width': '95%', 'horizontal-align': 'center'})
                        ],
                        style={'width': '32%', 'display': 'inline-block'}),

                        html.Div([
                            html.H4('Y-Axis Feature', style={'text-align': 'center'}),

                            dcc.Dropdown(
                                self.feature_labels_list,
                                self.feature_labels_list[1],
                                id='2d-yaxis-column'
                            ),

                            dcc.RadioItems(
                                ['Linear', 'Log'],
                                'Linear',
                                id='2d-yaxis-type',
                                labelStyle={'display': 'inline-block', 'marginTop': '5px'}
                            ),

                            html.Div([
                                html.H4('Solvent Accessible Surface Area Range',
                                        style={'text-align': 'center'},
                                        id='yaxis-feature-slider-header'),

                                dcc.RangeSlider(
                                min=round(self.df_scalar_measures['sasa'].min(), 1),
                                max=round(self.df_scalar_measures['sasa'].max(), 1),
                                step=0.1,
                                id='2d-yaxis-slider',
                                value=[round(self.df_scalar_measures['sasa'].min(), 1),
                                       round(self.df_scalar_measures['sasa'].max(), 1)],
                                marks=(int(((self.df_scalar_measures['sasa'].max() -
                                             self.df_scalar_measures['sasa'].min())/ 7)) or 1)),
                            ], style={'width': '95%', 'horizontal-align': 'center'})
                        ],
                        style={'width': '32%', 'horizontal-align': 'center',
                               'display': 'inline-block', 'justify-content': 'center',
                               'align-items': 'center'}),

                        html.Div([
                            html.H4('Measurements Aggregation Type',
                                    style={'text-align': 'center'}),

                            dcc.Dropdown(
                                list(self.aggregation_types.keys()),
                                'Average',
                                id='2d-agg-col'
                            ),

                            html.Div(id='aggon-radio-container',
                                children=[
                                dcc.RadioItems(
                                    ['On', 'Off'],
                                    'On',
                                    id='2d-agg-type',
                                    labelStyle={'display': 'inline-block', 'marginTop': '5px'}
                                ),
                            ], style={'width': '48%', 'float': 'left', 'display': 'inline-block'}),

                            html.Div(id='sd-radio-container',
                                children=[
                                dcc.RadioItems(
                                    ['+ st dev', '- st dev'],
                                    '- st dev',
                                    id='2d-stdev-toggle',
                                    labelStyle={'display': 'inline-block', 'marginTop': '5px'}
                                ),
                            ], style={'width': '48%', 'float': 'right', 'display': 'inline-block'}),

                            html.Div([
                                html.H4('Scatter Colour Section', style={'text-align': 'center'}),

                                dcc.Dropdown(
                                    ['None'] + self.colour_cols,
                                    'None',
                                    id='2d-colour-col'
                                ),])
                        ], style={'width': '32%', 'float': 'right', 'display': 'inline-block'})
                    ],),

                    html.Div([
                        dcc.Graph(
                            id='2d-indicator-scatter',
                            hoverData={'points': [{'customdata': self.df_measures[self.res_key].iloc[0]}]}
                        )
                    ], style={'width': '49%', 'display': 'inline-block', 'padding': '0 20'}),

                    html.Div([
                        dcc.Graph(id='x-feat-hist'),

                        dcc.Graph(id='y-feat-hist'),

                        html.Div([
                            html.H5('Subset By:', style={'text-align': 'center'}),

                            dcc.Dropdown(
                                ['Uniprot Code', 'Uniprot Code + Resid', 'Uniprot Code + Chain + Resid'],
                                'Uniprot Code',
                                id='2d-little_hist_subset_dropdown'
                            ),
                        ], style={'display': 'inline-block', 'width': '48%'}),

                        html.Div([
                            html.H5('Colour By:', style={'text-align': 'center'}),

                            dcc.Dropdown(
                                ['Uniform', 'Chain', 'Resid', 'Chain + Resid'] + self.colour_cols,
                                'Uniform',
                                id='2d-little_hist_colour_dropdown'
                            ),
                        ], style={'display': 'inline-block', 'width': '48%'})
                    ], style={'display': 'inline-block', 'width': '49%', 'float': 'right'}),

                    html.Div([
                        html.H3('GO Term Subset', style={'text-align': 'center'}),

                        dcc.Dropdown(
                            ['All'] + self.GO_display_names,
                            'All',
                            id='2d-GOterm-dropdown',
                        ),
                    ],
                    style={'width': '95%', 'display': 'inline-block'}),
                ])


            elif tab == 'tab_histogram_analysis':
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
                                id='2d-agg-type-hist'
                            ),
                        ],
                        style={'width': '49%', 'horizontal-align': 'right',
                               'display': 'inline-block'}),
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
                        style={'width': '32%', 'horizontal-align': 'center',
                               'display': 'inline-block'}),

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

                    html.Div(id='mainhist-warningtext',
                             children=[
                        html.H4('WARNING: Plotting data without an aggregation type may be biased due to the number of available structures for each protein.')
                    ], style={'width': '100%', 'display': 'inline-block', 'float': 'center',
                              'text-align': 'center', 'padding': '0 20'}),

                    html.Div([
                        dcc.Graph(
                            id='histogram-main'
                        )
                    ], style={'width': '100%', 'display': 'inline-block', 'padding': '0 20'}),

                    html.Div([
                        html.H4('Number of Bins (Fits to Nearest Neat Splitting Bin Size)', style={'text-align': 'center'}),
                        dcc.Slider(
                        min=0,
                        max=100,
                        step=1,
                        id='main-hist-bins-slider',
                        marks=10),
                    ], style={'width': '95%', 'horizontal-align': 'center'}),

                    html.Div([
                        html.H4('pKa (PROPKA3) Range', style={'text-align': 'center'},
                                id='main-hist-feature-slider-header'),
                        dcc.RangeSlider(
                        min=round(self.df_scalar_measures['propka'].min(), 1),
                        max=round(self.df_scalar_measures['propka'].max(), 1),
                        step=0.1,
                        id='main-hist-feature-slider',
                        value=[round(self.df_scalar_measures['propka'].min(), 1), round(self.df_scalar_measures['propka'].max(), 1)],
                        marks=int(((self.df_scalar_measures['propka'].max() - self.df_scalar_measures['propka'].min())/ 5)) or 1),
                    ], style={'width': '95%', 'horizontal-align': 'center'}),

                    html.Div([
                        html.H3('GO Term Subset', style={'text-align': 'center'}),
                        dcc.Dropdown(
                            ['All'] + self.GO_display_names,
                            'All',
                            id='mainhist-GOterm-dropdown',
                        ),
                    ], style={'width': '95%', 'display': 'inline-block'}),

                    html.Div([
                        html.H4('Scatter Colour Section', style={'text-align': 'center'}),
                        dcc.Dropdown(
                            ['None'] + self.colour_cols,
                            'None',
                            id='mainhist-colour-col'
                        ),
                    ], style={'width': '95%', 'display': 'inline-block'})
                ])

            elif tab == 'tab_3d_scalar_analysis':
                return html.Div([
                    html.Div([
                        html.Div([
                            html.H4('X-Axis Feature', style={'text-align': 'center'}),

                            dcc.Dropdown(
                                self.feature_labels_list,
                                self.feature_labels_list[0],
                                id='3d-xaxis-column',
                            ),

                            dcc.RadioItems(
                                ['Linear', 'Log'],
                                'Linear',
                                id='3d-xaxis-type',
                                labelStyle={'display': 'inline-block', 'marginTop': '5px'}
                            ),

                            html.Div([
                                html.H4('pKa (PROPKA3) Range',
                                        style={'text-align': 'center'},
                                        id='3d-xaxis-feature-slider-header'),

                                dcc.RangeSlider(
                                min=round(self.df_scalar_measures['propka'].min(), 1),
                                max=round(self.df_scalar_measures['propka'].max(), 1),
                                step=0.1,
                                id='3d-xaxis-slider',
                                value=[round(self.df_scalar_measures['propka'].min(), 1), round(self.df_scalar_measures['propka'].max(), 1)],
                                marks=(int(((self.df_scalar_measures['propka'].max() - self.df_scalar_measures['propka'].min())/ 7)) or 1)),
                            ], style={'width': '95%', 'horizontal-align': 'center'})
                        ],
                        style={'width': '32%', 'display': 'inline-block'}),

                        html.Div([
                            html.H4('Y-Axis Feature', style={'text-align': 'center'}),

                            dcc.Dropdown(
                                self.feature_labels_list,
                                self.feature_labels_list[1],
                                id='3d-yaxis-column'
                            ),

                            dcc.RadioItems(
                                ['Linear', 'Log'],
                                'Linear',
                                id='3d-yaxis-type',
                                labelStyle={'display': 'inline-block', 'marginTop': '5px'}
                            ),

                            html.Div([
                                html.H4('Solvent Accessible Surface Area Range',
                                        style={'text-align': 'center'},
                                        id='3d-yaxis-feature-slider-header'),

                                dcc.RangeSlider(
                                min=round(self.df_scalar_measures['sasa'].min(), 1),
                                max=round(self.df_scalar_measures['sasa'].max(), 1),
                                step=0.1,
                                id='3d-yaxis-slider',
                                value=[round(self.df_scalar_measures['sasa'].min(), 1), round(self.df_scalar_measures['sasa'].max(), 1)],
                                marks=(int(((self.df_scalar_measures['sasa'].max() - self.df_scalar_measures['sasa'].min())/ 7)) or 1)),
                            ], style={'width': '95%', 'horizontal-align': 'center'})
                        ],
                        style={'width': '32%', 'horizontal-align': 'center', 'display': 'inline-block', 'justify-content': 'center', 'align-items': 'center'}),

                        html.Div([
                            html.H4('Z-Axis Feature', style={'text-align': 'center'}),

                            dcc.Dropdown(
                                self.feature_labels_list,
                                self.feature_labels_list[2],
                                id='3d-zaxis-column',
                            ),

                            dcc.RadioItems(
                                ['Linear', 'Log'],
                                'Linear',
                                id='3d-zaxis-type',
                                labelStyle={'display': 'inline-block', 'marginTop': '5px'}
                            ),

                            html.Div([
                                html.H4('Depth', style={'text-align': 'center'},
                                        id='3d-zaxis-feature-slider-header'),

                                dcc.RangeSlider(
                                min=round(self.df_scalar_measures['depth'].min(), 1),
                                max=round(self.df_scalar_measures['depth'].max(), 1),
                                step=0.1,
                                id='3d-zaxis-slider',
                                value=[round(self.df_scalar_measures['depth'].min(), 1), round(self.df_scalar_measures['depth'].max(), 1)],
                                marks=(int(((self.df_scalar_measures['depth'].max() - self.df_scalar_measures['depth'].min())/ 7)) or 1)),
                            ], style={'width': '95%', 'horizontal-align': 'center'})
                        ],
                        style={'width': '32%', 'display': 'inline-block', 'horizontal-align': 'right'}),
                    ],),

                    html.Div([
                        html.H4('Measurements Aggregation Type', style={'text-align': 'center'}),

                        dcc.Dropdown(
                            ['None'] + list(self.aggregation_types.keys()),
                            'Average',
                            id='3d-agg-col'
                        ),
                    ], style={'width': '50%', 'float': 'left', 'display': 'inline-block'}),

                    html.Div([
                        html.Div([
                            html.H4('Scatter Colour Section',
                                    style={'text-align': 'center'}),

                            dcc.Dropdown(
                                ['None'] + self.colour_cols,
                                'None',
                                id='3d-colour-col'
                            ),])
                    ], style={'width': '50%', 'float': 'right', 'display': 'inline-block'}),

                    html.Div([
                        dcc.Graph(
                            id='3d-indicator-scatter',
                            hoverData={'points': [{'customdata': self.df_measures[self.res_key].iloc[0]}]}
                        )
                    ], style={'width': '98%', 'display': 'inline-block', 'padding': '0 20'}),

                    html.Div([
                        dcc.Graph(id='3d-x-feat-hist'),
                    ], style={'display': 'inline-block', 'width': '32%'}),

                    html.Div([
                        dcc.Graph(id='3d-y-feat-hist'),
                    ], style={'display': 'inline-block', 'width': '32%'}),

                    html.Div([
                        dcc.Graph(id='3d-z-feat-hist'),
                    ], style={'display': 'inline-block', 'width': '32%'}),


                    html.Div([
                        html.H5('Individual Histograms Subset By:',
                                style={'text-align': 'center'}),

                        dcc.Dropdown(
                            ['Uniprot Code', 'Uniprot Code + Resid', 'Uniprot Code + Chain + Resid'],
                            'Uniprot Code',
                            id='3d-little_hist_subset_dropdown'
                        ),
                    ], style={'display': 'inline-block', 'width': '48%'}),

                    html.Div([
                        html.H5('Individual Histograms Colour By:',
                                style={'text-align': 'center'}),

                        dcc.Dropdown(
                            ['Uniform', 'Chain', 'Resid', 'Chain + Resid'] + self.colour_cols,
                            'Uniform',
                            id='3d-little_hist_colour_dropdown'
                        ),
                    ], style={'display': 'inline-block', 'width': '48%'}),

                    html.Div([
                        html.H3('GO Term Overall Subset',
                                style={'text-align': 'center'}),

                        dcc.Dropdown(
                            ['All'] + self.GO_display_names,
                            'All',
                            id='3d-GOterm-dropdown',
                        ),
                    ],
                    style={'width': '95%', 'display': 'inline-block'}),
                ])


        @callback(
            Output('2d-indicator-scatter', 'figure'),
            Output('2d-axis-column', 'value'),
            Output('2d-yaxis-column', 'value'),
            Output('2d-agg-col', 'value'),
            Output('xaxis-feature-slider-header', 'children'),
            Output('yaxis-feature-slider-header', 'children'),
            Output('2d-xaxis-slider', 'min'),
            Output('2d-xaxis-slider', 'max'),
            Output('2d-yaxis-slider', 'min'),
            Output('2d-yaxis-slider', 'max'),
            Output('2d-xaxis-slider', 'value'),
            Output('2d-yaxis-slider', 'value'),
            Output('sd-radio-container', 'style'),
            Output('aggon-radio-container', 'style'),
            Input('2d-axis-column', 'value'),
            Input('2d-yaxis-column', 'value'),
            Input('2d-xaxis-type', 'value'),
            Input('2d-yaxis-type', 'value'),
            Input('2d-agg-col', 'value'),
            Input('2d-agg-type', 'value'),
            Input('2d-xaxis-slider', 'value'),
            Input('2d-yaxis-slider', 'value'),
            Input('2d-colour-col', 'value'),
            Input('2d-GOterm-dropdown', 'value'),
            Input('2d-stdev-toggle', 'value')
            )
        def update_graph_2d(xaxis_column_name, yaxis_column_name, xaxis_type, yaxis_type,
                         agg_type, agg_on, xaxis_range, yaxis_range, colour_col, go_term, st_dev_tog):
            trig_id = ctx.triggered_id

            x_low, x_high = xaxis_range
            y_low, y_high = yaxis_range

            if xaxis_column_name is None:
                xaxis_column_name = self.feature_labels_list[0]
            if yaxis_column_name is None:
                yaxis_column_name = self.feature_labels_list[1]
            if agg_type is None:
                agg_type = 'Average'
            if colour_col is None:
                colour_col = 'None'
            if go_term is None:
                go_term = 'All'

            if agg_on == 'On':
                dff = self.df_agg

                xaxis = f'{self.feature_labels[xaxis_column_name]}_{self.aggregation_types[agg_type]}'
                yaxis = f'{self.feature_labels[yaxis_column_name]}_{self.aggregation_types[agg_type]}'
            else:
                dff = self.df_scalar_measures

                xaxis = self.feature_labels[xaxis_column_name]
                yaxis = self.feature_labels[yaxis_column_name]

            if go_term != 'All':
                associated_uniprots = self.GO_dict[self.GO_codes_reverse[go_term]]
                dff = dff[dff['Uniprot_Entry'].isin(associated_uniprots)]

            xaxis_min = round(dff[xaxis].min(), 1)
            xaxis_max = round(dff[xaxis].max(), 1)
            yaxis_min = round(dff[yaxis].min(), 1)
            yaxis_max = round(dff[yaxis].max(), 1)
            xaxis_range_header = f'{xaxis_column_name} Range'
            yaxis_range_header = f'{yaxis_column_name} Range'

            if trig_id == '2d-axis-column':
                xaxis_range = [xaxis_min, xaxis_max]
            else:
                xaxis_range = no_update

            if trig_id == '2d-yaxis-column':
                yaxis_range = [yaxis_min, yaxis_max]
            else:
                yaxis_range = no_update

            dff = dff[(dff[xaxis] >= x_low) & (dff[xaxis] <= x_high)]
            dff = dff[(dff[yaxis] >= y_low) & (dff[yaxis] <= y_high)]

            res_key = self.res_key.copy()
            if 'Chain' in dff.columns:
                res_key.append('Chain')

            if colour_col != 'None' and (agg_on == 'Off' or agg_type == 'None'):
                aggon_radio = {'width': '98%', 'float': 'left', 'display': 'inline-block'}
                sd_radio = {'display': 'None'}
                fig = px.scatter(dff,
                                x=xaxis,
                                y=yaxis,
                                color=colour_col,
                                hover_name='Uniprot_Entry',
                                hover_data=res_key)
            elif agg_type == 'Average':
                if agg_on == 'On':
                    aggon_radio = {'width': '48%', 'float': 'left', 'display': 'inline-block'}
                    sd_radio = {'width': '48%', 'float': 'right', 'display': 'inline-block'}
                    if st_dev_tog == '+ st dev':
                        xaxis_error = f'{self.feature_labels[xaxis_column_name]}_sd'
                        yaxis_error = f'{self.feature_labels[yaxis_column_name]}_sd'
                        fig = px.scatter(dff,
                            x=xaxis,
                            y=yaxis,
                            color_discrete_sequence=['#682860'],
                            hover_name='Uniprot_Entry',
                            hover_data=res_key,
                            error_x=xaxis_error,
                            error_y=yaxis_error)
                    else:
                        fig = px.scatter(dff,
                            x=xaxis,
                            y=yaxis,
                            color_discrete_sequence=['#682860'],
                            hover_name='Uniprot_Entry',
                            hover_data=res_key)
                else:
                    aggon_radio = {'width': '98%', 'float': 'left', 'display': 'inline-block'}
                    sd_radio = {'display': 'None'}
                    fig = px.scatter(dff,
                        x=xaxis,
                        y=yaxis,
                        color_discrete_sequence=['#682860'],
                        hover_name='Uniprot_Entry',
                        hover_data=res_key)

            else:
                aggon_radio = {'width': '98%', 'float': 'left', 'display': 'inline-block'}
                sd_radio = {'display': 'None'}
                fig = px.scatter(dff,
                    x=xaxis,
                    y=yaxis,
                    color_discrete_sequence=['#682860'],
                    hover_name='Uniprot_Entry',
                    hover_data=res_key)

            fig.update_traces(mode='markers', customdata=dff[res_key])

            fig.update_layout(
                hoverlabel=dict(
                    bgcolor="white",
                    font_size=16
                )
            )

            fig.update_xaxes(title=self.axis_labels[self.feature_labels[xaxis_column_name]],
                             type='linear' if xaxis_type == 'Linear' else 'log')

            fig.update_yaxes(title=self.axis_labels[self.feature_labels[yaxis_column_name]],
                             type='linear' if yaxis_type == 'Linear' else 'log')

            fig.update_layout(margin={'l': 40, 'b': 40, 't': 10, 'r': 0}, hovermode='closest')

            return fig, xaxis_column_name, yaxis_column_name, agg_type, xaxis_range_header, yaxis_range_header, xaxis_min, xaxis_max, yaxis_min, yaxis_max, xaxis_range, yaxis_range, sd_radio, aggon_radio


        def create_feature_hist(dff, feature, title, axis_type, subset, colour):

            if colour == 'Uniform':
                fig = px.histogram(dff, x=feature, color_discrete_sequence=['#682860'])
            elif colour == 'Chain + Resid':
                if 'Chain' in dff.columns:
                    dfff = dff.copy()
                    dfff['CombColour'] = dfff[['Chain', 'Resid']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
                    fig = px.histogram(dfff, x=feature, color='CombColour')
                else:
                    fig = px.histogram(dff, x=feature, color_discrete_sequence=['#682860'])
            elif colour == 'Chain':
                if 'Chain' in dff.columns:
                    fig = px.histogram(dff, x=feature, color=colour)
                else:
                    fig = px.histogram(dff, x=feature, color_discrete_sequence=['#682860'])
            else:
                if colour in dff.columns:
                    if colour == 'PLDDT':
                        dfff = dff.copy()
                        plddt_bin_gap = 5
                        bin_edges = [float(a) for a in list(range(0, 101, plddt_bin_gap))]
                        bin_labels = [f'{int(a)}-{int(a)+5}' for a in bin_edges[:-1]]
                        dfff['PLDDT'] = pd.cut(dff['PLDDT'],
                                                    bins=bin_edges,
                                                    labels=bin_labels,
                                                    include_lowest=True,
                                                    ordered=True).astype(str)
                        fig = px.histogram(dfff, x=feature, color='PLDDT',
                                            category_orders={'PLDDT': bin_labels})
                    else:
                        fig = px.histogram(dff, x=feature, color=colour)
                else:
                    fig = px.histogram(dff, x=feature, color_discrete_sequence=['#682860'])

            fig.update_xaxes(showgrid=False)

            fig.update_yaxes(type='linear' if axis_type == 'Linear' else 'log')

            fig.add_annotation(x=0, y=0.85, xanchor='left', yanchor='bottom',
                            xref='paper', yref='paper', showarrow=False, align='left',
                            text=title)

            fig.update_layout(height=225, margin={'l': 20, 'b': 30, 'r': 10, 't': 10},
                              xaxis_title_text=self.axis_labels[feature],
                              yaxis_title_text='Count')

            return fig


        def update_little_hist(hoverData, col_name, subset, colour):
            if col_name is None:
                col_name = self.feature_labels_list[0]
            if subset is None:
                subset = 'Uniprot Code'
            if colour is None:
                colour = 'Uniform'

            axis_name = self.feature_labels[col_name]
            custom_hover_data = hoverData['points'][0]['customdata']
            if len(custom_hover_data) == 3:
                uniprot_entry = custom_hover_data[0]
                chain = custom_hover_data['customdata'][1]
                resid = custom_hover_data['customdata'][2]
                if subset == 'Uniprot Code + Chain + Resid':
                    dff = self.df_scalar_measures[(self.df_scalar_measures['Uniprot_Entry'] == uniprot_entry) &
                                                  (self.df_scalar_measures['Chain'] == chain) &
                                                  (self.df_scalar_measures['Resid'] == resid)]
                    title = '<b>{} {} {}</b><br>{}'.format(uniprot_entry, chain, resid, col_name)
                elif subset == 'Uniprot Code + Resid':
                    dff = self.df_scalar_measures[(self.df_scalar_measures['Uniprot_Entry'] == uniprot_entry) &
                                                    (self.df_scalar_measures['Resid'] == resid)]
                    title = '<b>{} {}</b><br>{}'.format(uniprot_entry, resid, col_name)
                else:
                    dff = self.df_scalar_measures[(self.df_scalar_measures['Uniprot_Entry'] == uniprot_entry)]
                    title = '<b>{}</b><br>{}'.format(uniprot_entry, col_name)

            else:
                uniprot_entry = custom_hover_data[0]
                resid = custom_hover_data[1]
                if subset == 'Uniprot Code + Resid':
                    dff = self.df_scalar_measures[(self.df_scalar_measures['Uniprot_Entry'] == uniprot_entry) &
                                                  (self.df_scalar_measures['Resid'] == resid)]
                    title = '<b>{} {}</b><br>{}'.format(uniprot_entry, resid, col_name)
                else:
                    dff = self.df_scalar_measures[(self.df_scalar_measures['Uniprot_Entry'] == uniprot_entry)]
                    title = '<b>{}</b><br>{}'.format(uniprot_entry, col_name)

            return dff, axis_name, title, subset, colour


        @callback(
            Output('x-feat-hist', 'figure'),
            Output('2d-little_hist_subset_dropdown', 'value'),
            Output('2d-little_hist_colour_dropdown', 'value'),
            Input('2d-indicator-scatter', 'hoverData'),
            Input('2d-axis-column', 'value'),
            Input('2d-xaxis-type', 'value'),
            Input('2d-little_hist_subset_dropdown', 'value'),
            Input('2d-little_hist_colour_dropdown', 'value'))
        def update_x_hist(hoverData, xaxis_column_name, axis_type, subset, colour):
            dff, xaxis, title, subset, colour = update_little_hist(hoverData, xaxis_column_name, subset, colour)
            return create_feature_hist(dff, xaxis, title, axis_type, subset, colour), subset, colour


        @callback(
            Output('y-feat-hist', 'figure'),
            Input('2d-indicator-scatter', 'hoverData'),
            Input('2d-yaxis-column', 'value'),
            Input('2d-yaxis-type', 'value'),
            Input('2d-little_hist_subset_dropdown', 'value'),
            Input('2d-little_hist_colour_dropdown', 'value'))
        def update_y_hist(hoverData, yaxis_column_name, axis_type, subset, colour):
            dff, yaxis, title, subset, colour = update_little_hist(hoverData, yaxis_column_name, subset, colour)
            return create_feature_hist(dff, yaxis, title, axis_type, subset, colour)


        # tab 2 functions
        def create_main_hist(dff, feature, title, nbins, chain_split, resid_split, colour_col):

            if colour_col != 'None':
                if colour_col in dff.columns:
                    if colour_col == 'PLDDT':
                        dfff = dff.copy()
                        plddt_bin_gap = 5
                        bin_edges = [float(a) for a in list(range(0, 101, plddt_bin_gap))]
                        bin_labels = [f'{int(a)}-{int(a)+5}' for a in bin_edges[:-1]]
                        dfff['PLDDT'] = pd.cut(dff['PLDDT'],
                                                   bins=bin_edges,
                                                   labels=bin_labels,
                                                   include_lowest=True,
                                                   ordered=True).astype(str)
                        fig = px.histogram(dfff, x=feature, marginal='rug',
                                            nbins=nbins, color='PLDDT',
                                            category_orders={'PLDDT': bin_labels})
                    else:
                        fig = px.histogram(dff, x=feature, marginal='rug',
                                            nbins=nbins, color=colour_col)
                else:
                    fig = px.histogram(dff, x=feature, marginal='rug',
                                        nbins=nbins, color_discrete_sequence=['#682860'])
            elif chain_split == 'Together' and resid_split == 'Together':
                fig = px.histogram(dff, x=feature, marginal='rug', color_discrete_sequence=['#682860'], nbins=nbins)
            elif chain_split == 'Together' and resid_split == 'Seperate':
                fig = px.histogram(dff, x=feature, marginal='rug',
                                   nbins=nbins, color='Resid')
            elif chain_split == 'Seperate' and resid_split == 'Together':
                if 'Chain' in dff.columns:
                    fig = px.histogram(dff, x=feature, marginal='rug',
                                        nbins=nbins, color='Chain')
                else:
                    fig = px.histogram(dff, x=feature, marginal='rug',
                                       color_discrete_sequence=['#682860'], nbins=nbins)
            else:
                dfff = dff.copy()
                dfff['CombColour'] = dfff[['Chain', 'Resid']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
                fig = px.histogram(dfff, x=feature, marginal='rug',
                                    nbins=nbins, color='CombColour')

            fig.update_xaxes(showgrid=False)
            fig.update_yaxes(showgrid=False)

            fig.update_layout(height=225, margin={'l': 20, 'b': 30, 'r': 10, 't': 10},
                              xaxis_title_text=self.axis_labels[self.feature_labels[title]], yaxis_title_text='Count')

            return fig


        @callback(
            Output('histogram-main', 'figure'),
            Output('crossfilter-feature-hist', 'value'),
            Output('2d-agg-type-hist', 'value'),
            Output('crossfilter-uniprot-hist', 'value'),
            Output('crossfilter-chain-hist', 'value'),
            Output('crossfilter-chain-hist', 'options'),
            Output('crossfilter-resid-hist', 'value'),
            Output('crossfilter-resid-hist', 'options'),
            Output('main-hist-feature-slider-header', 'children'),
            Output('main-hist-feature-slider', 'value'),
            Output('main-hist-feature-slider', 'min'),
            Output('main-hist-feature-slider', 'max'),
            Output('mainhist-warningtext', 'style'),
            Input('crossfilter-feature-hist', 'value'),
            Input('2d-agg-type-hist', 'value'),
            Input('crossfilter-uniprot-hist', 'value'),
            Input('crossfilter-chain-hist', 'value'),
            Input('crossfilter-chain-split-type', 'value'),
            Input('crossfilter-resid-hist', 'value'),
            Input('crossfilter-resid-split-type', 'value'),
            Input('main-hist-bins-slider', 'value'),
            Input('main-hist-feature-slider', 'value'),
            Input('mainhist-GOterm-dropdown', 'value'),
            Input('mainhist-colour-col', 'value'))
        def update_main_hist(feature, agg_type, uniprot, chain, chain_split,
                             resid, resid_split, nbins, feat_range, go_term, colour_col):

            trig_id = ctx.triggered_id

            if feature is None:
                feature = self.feature_labels_reverse[self.features[0]]
            if agg_type is None:
                agg_type = 'None'
            if uniprot is None:
                uniprot = 'All'
                chain = 'All'
                resid = 'All'
            if chain is None:
                chain = 'All'
            if resid is None:
                resid = 'All'
            if go_term is None:
                go_term = 'All'
            if colour_col is None:
                colour_col = 'None'

            pot_chains = ['All'] + list(self.df_measures['Chain'].unique())
            pot_resids = ['All'] + list(self.df_measures['Resid'].unique())

            feat_low, feat_high = feat_range
            header = f'{feature} Range'

            if agg_type == 'None':
                warning_label = {'width': '100%', 'display': 'inline-block',
                                 'float': 'center', 'text-align': 'center',
                                 'padding': '0 20'}
                
                feat = self.feature_labels[feature]

                cols_remain = [a for a in self.non_feat_cols if a in self.df_scalar_measures.columns] + [feat]
                dff = self.df_scalar_measures[cols_remain]
                dff = dff[(dff[feat] > feat_low) & (dff[feat] < feat_high)]
                if uniprot != 'All':
                    dff = dff[dff['Uniprot_Entry'] == uniprot]
                    pot_chains = list(dff['Chain'].unique())
                    pot_resids = list(dff['Resid'].unique())
                    if chain != 'All' and chain in list(dff['Chain']):
                        dff = dff[dff['Chain'] == chain]
                    if resid != 'All' and resid in list(dff['Resid']):
                        dff = dff[dff['Resid'] == resid]

                if go_term != 'All':
                    associated_uniprots = self.GO_dict[self.GO_codes_reverse[go_term]]
                    dff = dff[dff['Uniprot_Entry'].isin(associated_uniprots)]

                feat_min = round(self.df_scalar_measures[feat].min(), 1)
                feat_max = round(self.df_scalar_measures[feat].max(), 1)

                if trig_id == 'crossfilter-feature-hist':
                    feat_range = [feat_min, feat_max]
                else:
                    feat_range = no_update

                return create_main_hist(dff, feat, feature, nbins, chain_split, resid_split, colour_col), feature, agg_type, uniprot, chain, pot_chains, resid, pot_resids, header, feat_range, feat_min, feat_max, warning_label

            else:
                warning_label = {'display': 'None'}
                feat = self.feature_labels[feature]
                agg = self.aggregation_types[agg_type]
                feat_col = f'{feat}_{agg}'
                cols_remain = [a for a in self.non_feat_cols if a in self.df_agg.columns] + [feat_col]
                dff = self.df_agg[cols_remain]
                dff = dff[(dff[feat_col] > feat_low) & (dff[feat_col] < feat_high)]
                if uniprot != 'All':
                    dff = dff[dff['Uniprot_Entry'] == uniprot]
                    pot_chains = list(dff['Chain'].unique())
                    pot_resids = list(dff['Resid'].unique())
                    if chain != 'All' and chain in list(dff['Chain']):
                        dff = dff[dff['Chain'] == chain]
                    if resid != 'All' and resid in list(dff['Resid']):
                        dff = dff[dff['Resid'] == resid]

                if go_term != 'All':
                    associated_uniprots = self.GO_dict[self.GO_codes_reverse[go_term]]
                    dff = dff[dff['Uniprot_Entry'].isin(associated_uniprots)]

                feat_min = round(self.df_agg[feat_col].min(), 1)
                feat_max = round(self.df_agg[feat_col].max(), 1)

                if trig_id == 'crossfilter-feature-hist':
                    feat_range = [feat_min, feat_max]
                else:
                    feat_range = no_update

                return create_main_hist(dff, feat_col, feature, nbins, chain_split, resid_split, colour_col), feature, agg_type, uniprot, chain, pot_chains, resid, pot_resids, header, feat_range, feat_min, feat_max, warning_label



        @callback(
            Output('3d-indicator-scatter', 'figure'),
            Output('3d-xaxis-column', 'value'),
            Output('3d-yaxis-column', 'value'),
            Output('3d-zaxis-column', 'value'),
            Output('3d-agg-col', 'value'),
            Output('3d-xaxis-feature-slider-header', 'children'),
            Output('3d-yaxis-feature-slider-header', 'children'),
            Output('3d-zaxis-feature-slider-header', 'children'),
            Output('3d-xaxis-slider', 'min'),
            Output('3d-xaxis-slider', 'max'),
            Output('3d-yaxis-slider', 'min'),
            Output('3d-yaxis-slider', 'max'),
            Output('3d-zaxis-slider', 'min'),
            Output('3d-zaxis-slider', 'max'),
            Output('3d-xaxis-slider', 'value'),
            Output('3d-yaxis-slider', 'value'),
            Output('3d-zaxis-slider', 'value'),
            Input('3d-xaxis-column', 'value'),
            Input('3d-yaxis-column', 'value'),
            Input('3d-zaxis-column', 'value'),
            Input('3d-xaxis-type', 'value'),
            Input('3d-yaxis-type', 'value'),
            Input('3d-zaxis-type', 'value'),
            Input('3d-agg-col', 'value'),
            Input('3d-xaxis-slider', 'value'),
            Input('3d-yaxis-slider', 'value'),
            Input('3d-zaxis-slider', 'value'),
            Input('3d-colour-col', 'value'),
            Input('3d-GOterm-dropdown', 'value')
            )
        def update_graph_3d(xaxis_column_name, yaxis_column_name, zaxis_column_name, xaxis_type,
                            yaxis_type, zaxis_type, agg_type, xaxis_range, yaxis_range,
                            zaxis_range, colour_col, go_term):
            trig_id = ctx.triggered_id

            x_low, x_high = xaxis_range
            y_low, y_high = yaxis_range
            z_low, z_high = zaxis_range

            if xaxis_column_name is None:
                xaxis_column_name = self.feature_labels_list[0]
            if yaxis_column_name is None:
                yaxis_column_name = self.feature_labels_list[1]
            if zaxis_column_name is None:
                zaxis_column_name = self.feature_labels_list[2]
            if agg_type is None:
                agg_type = 'Average'
            if colour_col is None:
                colour_col = 'None'
            if go_term is None:
                go_term = 'All'

            if agg_type != 'None':
                dff = self.df_agg

                xaxis = f'{self.feature_labels[xaxis_column_name]}_{self.aggregation_types[agg_type]}'
                yaxis = f'{self.feature_labels[yaxis_column_name]}_{self.aggregation_types[agg_type]}'
                zaxis = f'{self.feature_labels[zaxis_column_name]}_{self.aggregation_types[agg_type]}'
            else:
                dff = self.df_scalar_measures

                xaxis = self.feature_labels[xaxis_column_name]
                yaxis = self.feature_labels[yaxis_column_name]
                zaxis = self.feature_labels[zaxis_column_name]

            xaxis_min = round(dff[xaxis].min(), 1)
            xaxis_max = round(dff[xaxis].max(), 1)
            yaxis_min = round(dff[yaxis].min(), 1)
            yaxis_max = round(dff[yaxis].max(), 1)
            zaxis_min = round(dff[zaxis].min(), 1)
            zaxis_max = round(dff[zaxis].max(), 1)
            xaxis_range_header = f'{xaxis_column_name} Range'
            yaxis_range_header = f'{yaxis_column_name} Range'
            zaxis_range_header = f'{zaxis_column_name} Range'

            if trig_id == '3d-xaxis-column':
                xaxis_range = [xaxis_min, xaxis_max]
            else:
                xaxis_range = no_update

            if trig_id == '3d-yaxis-column':
                yaxis_range = [yaxis_min, yaxis_max]
            else:
                yaxis_range = no_update

            if trig_id == '3d-zaxis-column':
                zaxis_range = [zaxis_min, zaxis_max]
            else:
                zaxis_range = no_update

            if go_term != 'All':
                associated_uniprots = self.GO_dict[self.GO_codes_reverse[go_term]]
                dff = dff[dff['Uniprot_Entry'].isin(associated_uniprots)]

            dff = dff[(dff[xaxis] >= x_low) & (dff[xaxis] <= x_high)]
            dff = dff[(dff[yaxis] >= y_low) & (dff[yaxis] <= y_high)]
            dff = dff[(dff[zaxis] >= z_low) & (dff[zaxis] <= z_high)]

            res_key = self.res_key.copy()
            if 'Chain' in dff.columns:
                res_key.append('Chain')

            if colour_col != 'None' and agg_type == 'None':
                fig = px.scatter_3d(dff,
                                x=xaxis,
                                y=yaxis,
                                z=zaxis,
                                color=colour_col,
                                hover_name='Uniprot_Entry',
                                hover_data=res_key)
            else:
                fig = px.scatter_3d(dff,
                    x=xaxis,
                    y=yaxis,
                    z=zaxis,
                    color_discrete_sequence=['#682860'],
                    hover_name='Uniprot_Entry',
                    hover_data=res_key)

            fig.update_traces(mode='markers', customdata=dff[res_key])

            fig.update_layout(
                hoverlabel=dict(
                    bgcolor="white",
                    font_size=16
                ),
                autosize=True,
                height=700,
                scene=dict(
                    aspectmode='cube',
                    xaxis=dict(title=self.axis_labels[self.feature_labels[xaxis_column_name]],
                               type='linear' if xaxis_type == 'Linear' else 'log'),
                    yaxis=dict(title=self.axis_labels[self.feature_labels[yaxis_column_name]],
                               type='linear' if yaxis_type == 'Linear' else 'log'),
                    zaxis=dict(title=self.axis_labels[self.feature_labels[zaxis_column_name]],
                               type='linear' if zaxis_type == 'Linear' else 'log'),
                )
            )

            fig.update_layout(margin={'l': 40, 'b': 40, 't': 10, 'r': 0}, hovermode='closest')

            return fig, xaxis_column_name, yaxis_column_name, zaxis_column_name, agg_type, xaxis_range_header, yaxis_range_header, zaxis_range_header, xaxis_min, xaxis_max, yaxis_min, yaxis_max, zaxis_min, zaxis_max, xaxis_range, yaxis_range, zaxis_range


        @callback(
            Output('3d-x-feat-hist', 'figure'),
            Output('3d-little_hist_subset_dropdown', 'value'),
            Output('3d-little_hist_colour_dropdown', 'value'),
            Input('3d-indicator-scatter', 'hoverData'),
            Input('3d-xaxis-column', 'value'),
            Input('3d-xaxis-type', 'value'),
            Input('3d-little_hist_subset_dropdown', 'value'),
            Input('3d-little_hist_colour_dropdown', 'value'))
        def update_3d_x_hist(hoverData, xaxis_column_name, axis_type, subset, colour):
            dff, xaxis, title, subset, colour = update_little_hist(hoverData, xaxis_column_name, subset, colour)
            return create_feature_hist(dff, xaxis, title, axis_type, subset, colour), subset, colour

        @callback(
            Output('3d-y-feat-hist', 'figure'),
            Input('3d-indicator-scatter', 'hoverData'),
            Input('3d-yaxis-column', 'value'),
            Input('3d-yaxis-type', 'value'),
            Input('3d-little_hist_subset_dropdown', 'value'),
            Input('3d-little_hist_colour_dropdown', 'value'))
        def update_3d_y_hist(hoverData, yaxis_column_name, axis_type, subset, colour):
            dff, yaxis, title, subset, colour = update_little_hist(hoverData, yaxis_column_name, subset, colour)
            return create_feature_hist(dff, yaxis, title, axis_type, subset, colour)

        @callback(
            Output('3d-z-feat-hist', 'figure'),
            Input('3d-indicator-scatter', 'hoverData'),
            Input('3d-zaxis-column', 'value'),
            Input('3d-zaxis-type', 'value'),
            Input('3d-little_hist_subset_dropdown', 'value'),
            Input('3d-little_hist_colour_dropdown', 'value'))
        def update_3d_z_hist(hoverData, zaxis_column_name, axis_type, subset, colour):
            dff, zaxis, title, subset, colour = update_little_hist(hoverData, zaxis_column_name, subset, colour)
            return create_feature_hist(dff, zaxis, title, axis_type, subset, colour)



if __name__ == '__main__':
    outdir = 'demo'
    df_measures = f'{outdir}{os.sep}measures.csv'
    V = Viewer(df_measures=df_measures,
                outdir=outdir)
    V.launch_viewer()
