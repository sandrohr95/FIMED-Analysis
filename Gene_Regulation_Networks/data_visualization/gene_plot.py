#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 09:48:09 2019

@author: Sandro-Hurtado
"""
from grn_algorithms.grn_arbotero import run_arboreto
import networkx as nx
from scipy.spatial.distance import pdist, squareform
import pandas as pd
import seaborn as sns
import plotly.graph_objs as go
import plotly.plotly
import plotly.figure_factory as ff
import plotly
from plotly.offline import plot
from plotly.graph_objs import Figure, Data, Layout, Scatter3d, XAxis, YAxis, ZAxis, Scene, Margin, Annotations, Font, \
    Annotation, Line, Marker
from math import pi
from bokeh.io import show
from bokeh.embed import components
from bokeh.plotting import figure
from bokeh.models import (
    ColumnDataSource,
    HoverTool,
    LinearColorMapper,
    BasicTicker,
    PrintfTickFormatter,
    ColorBar,
)

plotly.tools.set_credentials_file(username='sandrohr', api_key='geLPhCC0443kFfOFbeU5')


class GenePlots:

    def __init__(self, expression_data):
        """

        :param expression_data: gene expression data
        """
        self.df = expression_data

    def plot_heatmap(self):

        self.df.index.name = 'Genes'
        self.df.columns.name = 'Muestras'

        dfej = pd.DataFrame(self.df.stack(), columns=['percentage']).reset_index()

        genes = list(self.df.index)
        samples = list(self.df.columns)

        colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]
        mapper = LinearColorMapper(palette=colors, low=dfej.percentage.min(), high=dfej.percentage.max())

        source = ColumnDataSource(dfej)  # Contiene los datos que le hemos pasado en forma de columnas

        tools = "hover,save,pan,box_zoom,reset,wheel_zoom"

        p = figure(title="Gene expression Heatmap",
                   x_range=genes, y_range=samples,
                   x_axis_location="above", plot_width=900, plot_height=400,
                   tools=tools, toolbar_location='below')

        p.grid.grid_line_color = None
        p.axis.axis_line_color = None
        p.axis.major_tick_line_color = None
        p.axis.major_label_text_font_size = "7pt"
        p.axis.major_label_standoff = 0
        p.xaxis.major_label_orientation = pi / 3

        # Formamos un rectangulo
        p.rect(x='Genes', y='Muestras', width=1, height=1,
               source=source,
               fill_color={'field': 'percentage', 'transform': mapper},
               line_color=None)

        color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="5pt",
                             ticker=BasicTicker(desired_num_ticks=len(colors)),
                             formatter=PrintfTickFormatter(format="%d"),
                             label_standoff=6, border_line_color=None, location=(0, 0))
        p.add_layout(color_bar, 'right')

        # rate: percentage
        p.select_one(HoverTool).tooltips = [
            ('value', '@percentage'),
            ('Gen', '@Genes'),
            ('Sample', '@Muestras')
        ]

        # show(p)      # show the plot
        # Me devuelve el script que contiene la gráfica y el div asociado
        script, div = components(p)
        print(script, div)

        return div

    def cluster_heat_map(self):

        # Compute the correlation matrix
        genes = list(self.df.index)
        # Initialize figure by creating upper dendrogram
        figure = ff.create_dendrogram(self.df, orientation='bottom', labels=genes)

        for i in range(len(figure['data'])):
            figure['data'][i]['yaxis'] = 'y2'

        # Create Side Dendrogram
        dendro_side = ff.create_dendrogram(self.df, orientation='right')
        for i in range(len(dendro_side['data'])):
            dendro_side['data'][i]['xaxis'] = 'x2'

        # Add Side Dendrogram Data to Figure
        for data in dendro_side['data']:
            figure.add_trace(data)

        # Create Heatmap
        dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
        dendro_leaves = list(map(int, dendro_leaves))
        # Distancias por parejas entre observaciones en el espacio n-dimensional.
        data_dist = pdist(self.df)
        # Convierte un vector de distancia de forma vectorial en una matriz de distancia de forma cuadrada.
        heat_data = squareform(data_dist)
        heat_data = heat_data[dendro_leaves, :]
        heat_data = heat_data[:, dendro_leaves]

        heat_map = [
            go.Heatmap(
                x=dendro_leaves,
                y=dendro_leaves,
                z=heat_data,
                colorscale='Blues'
            )
        ]

        heat_map[0]['x'] = figure['layout']['xaxis']['tickvals']
        heat_map[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

        # Add Heatmap Data to Figure
        for data in heat_map:
            figure.add_trace(data)

        # Edit Layout
        figure['layout'].update({'width': 800, 'height': 800,
                                 'showlegend': False, 'hovermode': 'closest',
                                 'title': 'Clustering Gene Expression Data'
                                 })
        # Edit xaxis
        figure['layout']['xaxis'].update({'domain': [.15, 1],
                                          'mirror': False,
                                          'showgrid': False,
                                          'showline': False,
                                          'zeroline': False,
                                          'ticks': ""})
        # Edit xaxis2
        figure['layout'].update({'xaxis2': {'domain': [0, .15],
                                            'mirror': False,
                                            'showgrid': False,
                                            'showline': False,
                                            'zeroline': False,
                                            'showticklabels': False,
                                            'ticks': ""}})

        # Edit yaxis
        figure['layout']['yaxis'].update({'domain': [0, .85],
                                          'mirror': False,
                                          'showgrid': False,
                                          'showline': False,
                                          'zeroline': False,
                                          'showticklabels': False,
                                          'ticks': ""})
        # Edit yaxis2
        figure['layout'].update({'yaxis2': {'domain': [.825, .975],
                                            'mirror': False,
                                            'showgrid': False,
                                            'showline': False,
                                            'zeroline': False,
                                            'showticklabels': False,
                                            'ticks': ""}})

        # plotly.offline.plot(figure, filename='ClusterMap.html', auto_open=True)
        script = plot(figure, output_type='div', include_plotlyjs=False, show_link=True)
        print(script)

        return script

    def sns_clusterMap(self):
        self.df.index.name = 'Genes'
        self.df.columns.name = 'Muestras'
        sns.set(color_codes=True)
        sns.clustermap(self.df)
