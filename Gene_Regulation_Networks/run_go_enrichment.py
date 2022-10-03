import go_enrichment_tools.go_enrichment as go_en
import pandas as pd


def run_go_enrich(expression_data: pd.DataFrame, tax_id: int = 9096):
    """

    :param expression_data:
    :param tax_id: gene id --> by default gene_id = 9606 for humans
    """
    # Cogemos la lista de genes resultado del pre-procesamiento y le aplicamos el enriquecimiento de genes
    gene_list = expression_data.index.values
    print(len(gene_list))
    # tax_id for human = 9096
    # we can load more than one specie taxids=[9606, 7227]
    go_en.execute_go_enrichment_analysis(gene_list, tax_id)
