import networkx as nx
import plotly
from plotly.offline import plot
from plotly.graph_objs import Figure, Data, Layout, Scatter3d, XAxis, YAxis, ZAxis, Scene, Margin, Annotations, Font, \
    Annotation, Line, Marker
from collections import Counter
import pandas as pd
import codecs

plotly.tools.set_credentials_file(username='sandrohr', api_key='geLPhCC0443kFfOFbeU5')

"""  Gene Regulatory Network Plot  """


def gene_regulatory_network_plot(network, net_threshold, config, algorithm_name):
    """

    :param network:
    :param net_threshold: double
    :param config: int between 1 or 2 to select layout
    :param algorithm_name: 1, 2, 3, 4
    :return:
    """
    # We put a threshold to obtain a clear graph with the most representatives genes
    limit = network.index.size * net_threshold

    graph = nx.from_pandas_edgelist(network.head(int(limit)), 'TF', 'target', ['importance'],
                                    create_using=nx.Graph(directed=True))

    n_nodes = len(list(graph.node()))  # number of genes n_nodes
    list_nodes = list(graph.node())  # list of genes n_nodes

    layout = {
        1: nx.fruchterman_reingold_layout(graph, dim=3),
        2: nx.circular_layout(graph, dim=3)
    }.get(config, nx.circular_layout(graph, dim=3))

    layout_network = list(layout.values())

    edges = list(graph.edges())

    xn = [layout_network[k][0] for k in range(n_nodes)]  # x-coordinates of n_nodes
    yn = [layout_network[k][1] for k in range(n_nodes)]  # y-coordinates
    zn = [layout_network[k][2] for k in range(n_nodes)]  # z-coordinates
    xe = []
    ye = []
    ze = []
    for e in edges:
        xe += [layout[e[0]][0], layout[e[1]][0], None]  # x-coordinates of edge ends
        ye += [layout[e[0]][1], layout[e[1]][1], None]
        ze += [layout[e[0]][2], layout[e[1]][2], None]

    trace1 = Scatter3d(x=xe,
                       y=ye,
                       z=ze,
                       mode='lines',
                       line=Line(color='rgb(125,125,125)', width=2)
                       # hoverinfo='text'
                       )

    trace2 = Scatter3d(x=xn,
                       y=yn,
                       z=zn,
                       mode='markers+text',
                       textposition='top center',
                       name='genes',
                       marker=Marker(symbol='circle',
                                     size=6,
                                     color='#6959CD',
                                     colorscale='Viridis',
                                     line=Line(color='rgb(50,50,50)', width=1)
                                     ),
                       text=list_nodes,
                       hoverinfo='text'
                       )

    axis = dict(showbackground=False,
                showline=False,
                zeroline=False,
                showgrid=False,
                showticklabels=False,
                title=''
                )

    fig = Figure(data=Data([trace1, trace2]),
                 layout=Layout(
                     title="Gene Regulation Network (" + algorithm_name + ")",
                     width=1000,
                     height=1000,
                     showlegend=False,
                     scene=Scene(
                         xaxis=XAxis(axis),
                         yaxis=YAxis(axis),
                         zaxis=ZAxis(axis),
                     ),
                     margin=Margin(
                         t=100
                     ),
                     hovermode='closest',
                     annotations=Annotations([
                         Annotation(
                             showarrow=False,
                             text="Khaos Research Group",
                             xref='paper',
                             yref='paper',
                             x=0,
                             y=0.1,
                             xanchor='left',
                             yanchor='bottom',
                             font=Font(
                                 size=20
                             )
                         )
                     ]),
                 ))

    plotly.offline.plot(fig, filename=algorithm_name + '.html', auto_open=False)
    script = plot(fig, output_type='div', include_plotlyjs=False, show_link=True)
    print(script)
    return script


def new_gene_regulatory():
    f = codecs.open("/home/antonio/Anaconda_Projects/Python_Projects/FIMED-Analysis/VisualizationTools/bokeh_plots/Gene Regulatory Network.html", 'r')
    print(f.read())
    return f.read()
