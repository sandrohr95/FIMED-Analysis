import grn_algorithms.grn_panda as grn_panda
import grn_algorithms.grn_arbotero as grn_arboreto
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
        result = grn_arboreto.run_arboreto(expression_data, net_threshold, config, algorithm_name)
    elif algorithm_code == 2:
        algorithm_name = "GRNBOOST2"
        result = grn_arboreto.run_arboreto(expression_data, net_threshold, config, algorithm_name)
    elif algorithm_code == 3:
        algorithm_name = "Panda"
        result = grn_panda.run_grn_panda(expression_data, net_threshold, config, algorithm_name)
    elif algorithm_code == 4:
        algorithm_name = "Lioness"
        result = grn_panda.run_grn_panda(expression_data, net_threshold, config, algorithm_name)
    else:
        print("Wrong algorithm name.")

    return result
