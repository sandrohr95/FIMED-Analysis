from __future__ import print_function
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
import collections as cx
from goatools.godag_plot import plot_results
# For Humans:
import go_enrichment_tools.gene_protein_coding.genes_ncbi_9096_proteincoding as GeneID2nt_human


def execute_go_enrichment_analysis(gene_list: list, tax_id: int = 9606):
    # Download Ontologies
    download_go_basic_obo("go_enrichment_tools/ontology/go-basic.obo")
    # Download Associations
    fin_gene2go = download_ncbi_associations("go_enrichment_tools/ontology/gene2go")
    # Load Ontologies
    obodag = GODag("go_enrichment_tools/ontology/go-basic.obo")

    # Load Associations
    # Read NCBI's gene2go. Store annotations in a list of namedtuples
    obj_anno = Gene2GoReader(fin_gene2go, taxids=[9606])

    # Get namespace2association where:
    #    namespace is:
    #        BP: biological_process
    #        MF: molecular_function
    #        CC: cellular_component
    #    association is a dict:
    #        key: NCBI GeneID
    #        value: A set of GO IDs associated with that gene
    ns2assoc = obj_anno.get_ns2assc()

    for nspc, id2gos in ns2assoc.items():
        print("{NS} {N:,} annotated humans genes".format(NS=nspc, N=len(id2gos)))

    # 3. Initialize a GOEA object

    # The GO-EA object holds the Ontologies, Associations, and background.
    # Numerous studies can then be run without needing to re-load the above items.
    # In this case, we only run one GO-EA.

    list_human_genes = GeneID2nt_human.GENEID2NT

    goeaobj = GOEnrichmentStudyNS(
        # GeneID2nt_mus.keys(),
        list_human_genes.keys(),  # List of human protein-coding genes
        ns2assoc,  # gene_id/GO associations
        obodag,  # Ontology
        propagate_counts=False,
        alpha=0.05,  # default significance cut-off
        methods=['fdr_bh'])  # default multiple test correction method

    # 4. Read study genes
    # g[0] = tad_id
    # g[2] = GENE_ID
    # g[5] = Symbol
    gene_ids_study = []
    for g_human in list_human_genes.values():
        for g_name in gene_list:
            if g_name == g_human[5]:
                gene_ids_study.append(g_human[2])

    # 5. Run Gene Ontology Enrichment Analysis (GOEA)
    # You may choose to keep all results or just the significant results.
    # In this example, we choose to keep only the significant results.
    # 'p_' means "p-value". 'fdr_bh' is the multiple test method we are currently using.
    go_ea_results_all = goeaobj.run_study(gene_ids_study)
    go_ea_results_sig = [r for r in go_ea_results_all if r.p_fdr_bh < 0.05]

    # Example 1: Significant v All GOEA results
    print('{N} of {M:,} results were significant'.format(N=len(go_ea_results_sig), M=len(go_ea_results_all)))

    # Example 2: Enriched v Purified GOEA results
    print('Significant results: {E} enriched, {P} purified'.format(
        E=sum(1 for r in go_ea_results_sig if r.enrichment == 'e'),
        P=sum(1 for r in go_ea_results_sig if r.enrichment == 'p')))

    # Example 3: Significant GOEA results by namespace
    ctr = cx.Counter([r.NS for r in go_ea_results_sig])
    print('Significant results[{TOTAL}] = {BP} BP + {MF} MF + {CC} CC'.format(
        TOTAL=len(go_ea_results_sig),
        BP=ctr['BP'],  # biological_process
        MF=ctr['MF'],  # molecular_function
        CC=ctr['CC']))  # cellular_component

    # 6. Write results to an Excel file and to a text file
    goeaobj.wr_xlsx("go_enrichment_tools/go_results/result.xlsx", go_ea_results_sig)
    goeaobj.wr_txt("go_enrichment_tools/go_results/result.txt", go_ea_results_sig)

    # 7. Plot all significant GO terms
    # The "{NS}" in "nbt3102_{NS}.png" indicates that you will see three plots,
    # one for "biological_process"(BP), "molecular_function"(MF), and "cellular_component"(CC)
    plot_results("go_enrichment_tools/go_results/result_{NS}.png", go_ea_results_sig)
