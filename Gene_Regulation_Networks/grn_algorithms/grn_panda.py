import subprocess
import data_visualization.grn_plot as grn_plot
import data_visualization.grn_ensembler as grn_plot_improve

# import data_visualization.grn_bokeh as grn_bokeh
import pandas as pd
import uuid
import os
import shutil
from config import Config


def run_grn_panda(expression_data: pd.DataFrame, net_threshold: float = 0.05, config: int = 1,
                  algorithm_name: str = "Panda", motif_data=None, ppi_data=None):
    """ Run grn PANDAS. Example: python run_panda.py ToyExpressionData.txt ToyMotifData.txt ToyPPIData.txt Panda"""

    # Generate folders with unique id for each execution
    expression_data_dir = Config.config("GENE_EXPRESSION_DATA") + str(uuid.uuid1())
    network_dir = Config.config("GRN_PANDA_DIR") + str(uuid.uuid1())

    try:
        os.mkdir(expression_data_dir)
        os.mkdir(network_dir)
    except OSError:
        print("Creation of the directory %s failed")

    # Save gene expression data file as csv and send directory to read it in others grn scripts
    expression_data_file = expression_data_dir + "/gene_expression_data.csv"
    expression_data.to_csv(expression_data_file, header=False)

    # Generate grn_panda.csv random folder for each execution and send it to others grn scripts to save the network here
    network_file = network_dir + "/grn_panda_network.csv"

    # Execute subprocess to run grn algorithm
    command = Config.config("GRN_PANDA_ENV") + ' ' + Config.config("GRN_PANDA_SCRIPT") + ' ' + expression_data_file + ' ' + \
              algorithm_name + ' ' + network_file

    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True, executable='/bin/bash')

    # Wait until subprocess has finished
    p.communicate()

    # read_csv of temporary network
    network = pd.read_csv(network_file, header=None, sep=" ")
    network.columns = ["TF", "target", "importance"]

    # Run Script Plotly
    # script_html = grn_plot.gene_regulatory_network_plot(network, net_threshold, config, algorithm_name)


    # Delete Folders
    try:
        shutil.rmtree(expression_data_dir)
        shutil.rmtree(network_dir)
    except OSError as e:
        print("Error: %s : %s" % (network_dir, e.strerror))
        print("Error: %s : %s" % (expression_data_dir, e.strerror))

    # return script_html
    return network
