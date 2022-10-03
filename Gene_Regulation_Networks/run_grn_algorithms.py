import grn_algorithms.grn_panda as grn_panda
import grn_algorithms.grn_arbotero as grn_arboreto
import data_visualization.grn_plot as grn_plot
import data_visualization.grn_ensembler as ensembler

import pandas as pd


def run_grn_algorithms(expression_data: pd.DataFrame, net_threshold: float, config: int, algorithm_code: int):
    """
        RUN GRN ALGORITHMS

        :param algorithm_code: 1:(genie3), 2:(grnboost2), 3:(panda), 4:(lioness)
        :param config: 1:(force-layout), 2:(circular-layout)
        :param net_threshold: By default --> 0.05 (only 5% of interactions. The strongest interactions.
        :type expression_data: expression matrix dataframe
    """

    result = ""
    if algorithm_code == 1:
        algorithm_name = "GENIE3"
        network = grn_arboreto.run_arboreto(expression_data, net_threshold, config, algorithm_name)
        result = grn_plot.gene_regulatory_network_plot(network, net_threshold, config, algorithm_name)
    elif algorithm_code == 2:
        algorithm_name = "GRNBOOST2"
        network = grn_arboreto.run_arboreto(expression_data, net_threshold, config, algorithm_name)
        result = grn_plot.gene_regulatory_network_plot(network, net_threshold, config, algorithm_name)

    elif algorithm_code == 3:
        algorithm_name = "Panda"
        network = grn_panda.run_grn_panda(expression_data, net_threshold, config, algorithm_name)
        result = grn_plot.gene_regulatory_network_plot(network, net_threshold, config, algorithm_name)

    elif algorithm_code == 4:
        algorithm_name = "Lioness"
        network = grn_panda.run_grn_panda(expression_data, net_threshold, config, algorithm_name)
        result = grn_plot.gene_regulatory_network_plot(network, net_threshold, config, algorithm_name)

    elif algorithm_code == 5:
        genie3_net = grn_arboreto.run_arboreto(expression_data, net_threshold, config, "GENIE3")
        grnboost2_net = grn_arboreto.run_arboreto(expression_data, net_threshold, config, "GRNBOOST2")
        panda_net = grn_panda.run_grn_panda(expression_data, net_threshold, config, "Panda")
        lioness_net = grn_panda.run_grn_panda(expression_data, net_threshold, config, "Lioness")
        result = ensembler.gene_regulatory_network_plot_emsembler(genie3_net, grnboost2_net, panda_net, lioness_net,
                                                                  net_threshold, config, "Ensembler")
    else:
        print("error")

    return result
