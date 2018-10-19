"""

Example illustrating how to construct and
plot a random geometric graph.

"""

import rgg_ensemble_analysis.RGG as RGG
import numpy as np

np.random.seed(1)

N=100
d=2
r = 0.13
positions = np.random.uniform(0.0,1.0,(N,d) ) 
rgg_graph = RGG.RGG(N,r,d,Positions=positions)


network_plot = RGG.graph_plot(rgg_graph.graph,positions=positions)
network_plot.save_plot("graph_1",form="png")

#We can also adjust parameters such as the node shape, size and transparency and colour the nodes in terms of properties such as their degree:
network_plot.node_color = rgg_graph.Degrees()
network_plot.node_transparency = 0.6
network_plot.node_shape = '<'
network_plot.node_size = 200
network_plot.edge_style = 'solid'
network_plot.save_plot("graph_2",form="png")
