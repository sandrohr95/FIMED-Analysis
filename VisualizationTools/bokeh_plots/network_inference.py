import pandas as pd
import networkx
import matplotlib.pyplot as plt
import numpy as np
from bokeh.io import output_notebook, show, save
from bokeh.io import output_notebook, show, save
from bokeh.models import Range1d, Circle, ColumnDataSource, MultiLine, EdgesAndLinkedNodes, NodesAndLinkedEdges
from bokeh.plotting import figure
from bokeh.plotting import from_networkx
from bokeh.palettes import Blues8, Reds8, Purples8, Oranges8, Viridis8, Spectral8
from bokeh.transform import linear_cmap
from networkx.algorithms import community
from bokeh.models import EdgesAndLinkedNodes, NodesAndLinkedEdges
from bokeh.embed import components


network_file = "/home/antonio/Anaconda_Projects/Python_Projects/FIMED_ANALYSIS/Gene_Regulation_Networks/data_dir/Results00afb767-2295-11ec-abbf-fa163e4623e5/grn_panda_network.csv"
df = pd.read_csv(network_file, header=None, sep=" ")
df.columns = ["TF", "target", "importance"]
# We put a threshold to obtain a clear graph with the most representatives genes
limit = df.index.size * 0.2
# Create Graph
G = networkx.from_pandas_edgelist(df.head(int(limit)), 'TF', 'target', 'importance')
degrees = dict(networkx.degree(G))
networkx.set_node_attributes(G, name='degree', values=degrees)
# Slightly adjust degree so that the nodes with very small degrees are still visible
number_to_adjust_by = 4
adjusted_node_size = dict([(node, degree + number_to_adjust_by) for node, degree in networkx.degree(G)])
networkx.set_node_attributes(G, name='adjusted_node_size', values=adjusted_node_size)


# Calculate communities
communities = community.greedy_modularity_communities(G)

# Add modularity class and color as attributes to network graph
# Create empty dictionaries
modularity_class = {}
modularity_color = {}
# Loop through each community in the network
for community_number, community in enumerate(communities):
    # For each member of the community, add their community number and a distinct color
    for name in community:
        modularity_class[name] = community_number
        modularity_color[name] = Spectral8[community_number]

print(modularity_class)
# Add modularity class and color as attributes from the network above
networkx.set_node_attributes(G, modularity_class, 'modularity_class')
networkx.set_node_attributes(G, modularity_color, 'modularity_color')

# Choose colors for node and edge highlighting
node_highlight_color = 'white'
edge_highlight_color = 'black'

# Choose attributes from G network to size and color by —
# setting manual size (e.g. 10) or color (e.g. 'skyblue') also allowed
size_by_this_attribute = 'adjusted_node_size'
color_by_this_attribute = 'modularity_color'

# Pick a color palette — Blues8, Reds8, Purples8, Oranges8, Viridis8
color_palette = Blues8

# Choose a title!
title = 'Gene Regulatory Network'

# Establish which categories will appear when hovering over each node
HOVER_TOOLTIPS = [
    ("Gen", "@index"),
    ("Degree", "@degree"),
    ("Modularity Class", "@modularity_class"),
    ("Modularity Color", "$color[swatch]:modularity_color")
]

# Create a plot — set dimensions, toolbar, and title
plot = figure(tooltips=HOVER_TOOLTIPS,
              tools="pan,wheel_zoom,save,reset", active_scroll='wheel_zoom',
              x_range=Range1d(-10.1, 10.1), y_range=Range1d(-10.1, 10.1), title=title)

# Create a network graph object
# https://networkx.github.io/documentation/networkx-1.9/reference/generated/networkx.drawing.layout.spring_layout.html
network_graph = from_networkx(G, networkx.spring_layout, scale=10, center=(0, 0))

# Set node sizes and colors according to node degree (color as category from attribute)
network_graph.node_renderer.glyph = Circle(size=size_by_this_attribute, fill_color=color_by_this_attribute)
# Set node highlight colors
network_graph.node_renderer.hover_glyph = Circle(size=size_by_this_attribute, fill_color=node_highlight_color,
                                                 line_width=2)
network_graph.node_renderer.selection_glyph = Circle(size=size_by_this_attribute, fill_color=node_highlight_color,
                                                     line_width=2)

# Set edge opacity and width
network_graph.edge_renderer.glyph = MultiLine(line_alpha=0.5, line_width=1)
# Set edge highlight colors
network_graph.edge_renderer.selection_glyph = MultiLine(line_color=edge_highlight_color, line_width=2)
network_graph.edge_renderer.hover_glyph = MultiLine(line_color=edge_highlight_color, line_width=2)

# Highlight nodes and edges
network_graph.selection_policy = NodesAndLinkedEdges()
network_graph.inspection_policy = NodesAndLinkedEdges()

plot.renderers.append(network_graph)

show(plot)
script, div = components(plot)
print(div,script)

# save(plot, filename=f"{title}.html")
