import pandas as pd
from arboreto.algo import grnboost2, genie3
from arboreto.utils import load_tf_names
from distributed import Client, LocalCluster
import data_visualization.grn_plot as grn_plot

""" Run arboreto algorithms and return GRN graph """


def run_arboreto(expression_matrix: pd.DataFrame, net_threshold: float = 0.05, config: int = 1, algorithm_name: str = "GRNBOOST2"):
    """
    :param config: configuration layout (1: force-layout, 2: circular-layout)
    :param net_threshold: filter genes for better visualization
    :param expression_matrix: gene expression data
    :param algorithm_name:
    :return:
    """
    # Transpose the dataframe to get correct format to create the network
    df = expression_matrix.transpose()

    # Get all the TF Gene names
    tf_names = list(df)

    # Create a Dask Client, just in case we want parellalize the algorithm
    client = Client(processes=False)

    network = None
    # create network with columns --> TF, target Gene, Importance
    if algorithm_name == "GRNBOOST2":
        network = grnboost2(expression_data=df, tf_names=tf_names, client_or_address=client)
    elif algorithm_name == "GENIE3":
        network = genie3(expression_data=df, tf_names=tf_names, client_or_address=client)
    else:
        print("Wrong algorithm name. Should either be GENIE3 or GRNBoost2.")

    script_html = grn_plot.gene_regulatory_network_plot(network, net_threshold, config, algorithm_name)
    return script_html
