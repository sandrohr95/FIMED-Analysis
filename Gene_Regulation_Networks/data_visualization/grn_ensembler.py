import networkx as nx
import plotly
from plotly.offline import plot
from plotly.graph_objs import Figure, Data, Layout, Scatter3d, XAxis, YAxis, ZAxis, Scene, Margin, Annotations, Font, \
    Annotation, Line, Marker
from collections import Counter
import pandas as pd

plotly.tools.set_credentials_file(username='sandrohr', api_key='geLPhCC0443kFfOFbeU5')

"""  Gene Regulatory Network Plot  """


def gene_regulatory_network_plot_emsembler(network1, network2, network3, network4, net_threshold, config,
                                           algorithm_name):
    """

    :param network:
    :param net_threshold: double
    :param config: int between 1 or 2 to select layout
    :param algorithm_name: 1, 2, 3, 4
    :return:
    """
    network = get_most_important_arcs(network1, network2, network3, network4, net_threshold)
    # We put a threshold to obtain a clear graph with the most representatives genes

    graph = nx.from_pandas_edgelist(network, 'TF', 'target', ['frequency'],
                                    create_using=nx.Graph())

    layout = {
        1: nx.fruchterman_reingold_layout(graph, dim=3),
        2: nx.circular_layout(graph, dim=3)
    }.get(config, nx.circular_layout(graph, dim=3))

    freq = network['frequency']
    tf_mapping = network['TF']
    target_mapping = network['target']

    # Almaceno la lista de arcos en tres listas distintas dependiendo de la frecuencia con la que aparecen los arcos
    list_edges = []
    edges1 = []
    edges2 = []
    edges3 = []
    edges4 = []
    for i, f in enumerate(freq):
        if f == 1:
            edges1.append((tf_mapping[i], target_mapping[i], f))
        elif f == 2:
            edges2.append((tf_mapping[i], target_mapping[i], f))
        elif f == 3:
            edges3.append((tf_mapping[i], target_mapping[i], f))
        elif f == 4:
            edges4.append((tf_mapping[i], target_mapping[i], f))

    list_edges.append(edges1)
    list_edges.append(edges2)
    list_edges.append(edges3)
    list_edges.append(edges4)

    layout_network = list(layout.values())

    n_nodes = len(list(graph.node()))  # number of genes n_nodes
    list_nodes = list(graph.node())  # list of genes n_nodes

    Xn = [layout_network[k][0] for k in range(n_nodes)]  # x-coordinates of nodes
    Yn = [layout_network[k][1] for k in range(n_nodes)]  # y-coordinates
    Zn = [layout_network[k][2] for k in range(n_nodes)]  # z-coordinates

    trace = []
    for i, Edges in enumerate(list_edges):
        Xe = []
        Ye = []
        Ze = []
        for e in Edges:
            Xe += [layout[e[0]][0], layout[e[1]][0], None]  # x-coordinates of edge ends
            Ye += [layout[e[0]][1], layout[e[1]][1], None]
            Ze += [layout[e[0]][2], layout[e[1]][2], None]
        #            width = e[2]/5 #Cogemos la frecuencia con la que aparecen
        if i == 0:  # freq = 1
            color = 'rgb(125,125,125)'
            width = 2
        elif i == 1:  # freq = 2
            color = 'blue'
            width = 4
        elif i == 2:  # freq = 3
            color = 'green'
            width = 5
        else:  # freq = 4
            color = 'red'
            width = 6

        trace1 = Scatter3d(x=Xe,
                           y=Ye,
                           z=Ze,
                           mode='lines',
                           line=Line(color=color,
                                     width=width,
                                     ),
                           hoverinfo='text',
                           text="frequency:" + str(i + 1),
                           showlegend=True
                           )
        trace.append(trace1)

    trace2 = Scatter3d(x=Xn,
                       y=Yn,
                       z=Zn,
                       mode='markers+text',
                       textposition='top center',
                       name='genes',
                       marker=Marker(symbol='circle',
                                     size=12,
                                     color='#6959CD',
                                     colorscale='Viridis',
                                     line=Line(color='rgb(50,50,50)', width=2)
                                     ),
                       text=list_nodes,
                       hoverinfo='text'
                       )
    trace.append(trace2)

    axis = dict(showbackground=False,
                showline=False,
                zeroline=False,
                showgrid=False,
                showticklabels=False,
                title=''
                )

    fig = Figure(data=Data(trace),
                 layout=Layout(
                     title="Gene Regulatory Network Emsembler",
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
                     ]), ))

    plotly.offline.plot(fig, filename="ensembler" + '.html', auto_open=False)
    script = plot(fig, output_type='div', include_plotlyjs=False, show_link=True)
    print(script)
    return script


def get_most_important_arcs(network1, network2, network3, network4, net_threshold):
    limit1 = network1.index.size * net_threshold
    network1 = network1.head(int(limit1))

    limit2 = network2.index.size * net_threshold
    network2 = network2.head(int(limit2))

    limit3 = network3.index.size * net_threshold
    network3 = network3.head(int(limit3))

    limit4 = network4.index.size * net_threshold
    network4 = network4.head(int(limit4))

    # TAKE MOST IMPORTANT ARCS FROM NETWORK
    tf_mapping_n1 = network1['TF']
    target_mapping_n1 = network1['target']
    arc_n1 = list(zip(tf_mapping_n1, target_mapping_n1))

    tf_mapping_n2 = network2['TF']
    target_mapping_n2 = network2['target']
    arc_n2 = list(zip(tf_mapping_n2, target_mapping_n2))

    tf_mapping_n3 = network3['TF']
    target_mapping_n3 = network3['target']
    arc_n3 = list(zip(tf_mapping_n3, target_mapping_n3))

    tf_mapping_n4 = network4['TF']
    target_mapping_n4 = network4['target']
    arc_n4 = list(zip(tf_mapping_n4, target_mapping_n4))

    arc_n4.extend(arc_n1)
    arc_n4.extend(arc_n2)
    arc_n4.extend(arc_n3)

    common_arcs = Counter(arc_n4).most_common(40)

    # network TF, TG, Freq
    tf_list = []
    target_gen = []
    frequency = []

    for arcs in common_arcs:
        tf = arcs[0][0]
        tg = arcs[0][1]
        tf_list.append(tf)
        target_gen.append(tg)
        frequency.append(arcs[1])

    data_tuples = list(zip(tf_list, target_gen, frequency))
    network_arcs = pd.DataFrame(data_tuples, columns=['TF', 'target', 'frequency'])
    # print(network_arcs)
    return network_arcs
