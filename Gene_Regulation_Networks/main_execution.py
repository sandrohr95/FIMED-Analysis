import pre_processing.gene_expression as ge
import run_grn_algorithms as run
import data_visualization.gene_plot as plots
import run_go_enrichment as go_enrich

genEx = ge.GeneExpressions(["5dcd5365f534772c8292a6bd", "5dcd5381f534772c8292a6c0", "5dcd5396f534772c8292a6c3","5ddd26f3f534772c8292a6d0","5e3151efc4aa3b3ad03a3d80"], ['sample1','sample2','sample3','s4','s5'], 0.05)
data = genEx.load_gene_expressions()

# Execute Go enrichment
# go_enrich.run_go_enrich(data)



# Execute GRN with Panda Lioness (4)
run.run_grn_algorithms(data, 0.1, 1, 4)
# plot = plots.GenePlots(data)
# plot.plot_heatmap()
# plot.cluster_heat_map()



