import networkx as nx
import pandas as pd
from typing import Optional, Set, List, Union
import random
import numpy as np
import goenrichclone as goenrich

random.seed(1)
np.random.seed(1)


def printa_go():
    goenrich.obo.printa_goenrich()


def check_namespace(func):
    def wrapper(*args, **kwargs):
        valid_namespace = ['p', 'f', 'c', 'all', None]
        if not kwargs['namespace'] in valid_namespace:
            raise ValueError(f'namespace must be one of {valid_namespace}.')
        return func(*args, **kwargs)
    return wrapper


@check_namespace
def build_ontology(obo_file: str,
                   namespace: Optional[str] = None) -> nx.DiGraph:
    '''Returns gene ontology as a networkx DiGraph (using goenrich coveragetools).
    Can be filtered by 'namespace', which can be 'p' for biological process,
    'f' for molecular function, 'c' for cellular component, or 'all'.
    '''

    if not namespace:
        namespace = 'all'

    ont = goenrich.obo.ontology(obo_file)

    # Filter ontology by namespace
    ns_dict = {'p': 'biological_process',
               'f': 'molecular_function',
               'c': 'cellular_component'}
    if namespace in ['p', 'f', 'c']:
        ns_nodes = []
        for node, ns in ont.nodes.data('namespace'):
            if ns == ns_dict[namespace]:
                ns_nodes.append(node)
        ont = ont.subgraph(ns_nodes)

    return ont


@check_namespace
def read_gaf(gaf_file: str,
             namespace: Optional[str] = None,
             qualifiers: Optional[bool] = None) -> pd.DataFrame:
    '''Returns pandas dataframe containing the info from gaf file, for each
    or all of the three GO branches. Branch is set with 'namespace', which can
    be 'p' for biological process, 'f' for molecular function, 'c' for cellular
    component, or 'all'. Annotations with qualifiers are removed by default,
    but can be left by setting 'qualifiers=True'.
    '''

    column_names = ('db', 'db_object_id', 'db_object_symbol',
                    'qualifier', 'go_id', 'db_reference',
                    'evidence_code', 'with_from', 'aspect',
                    'db_object_name', 'db_object_synonym',
                    'db_object_type', 'taxon', 'date', 'assigned_by',
                    'annotation_extension', 'gene_product_form_id')

    # read gaf file and rename columns
    gaf_df = pd.read_csv(gaf_file,
                         sep='\t',
                         header=None,
                         comment='!',
                         names=column_names)

    # keep only rows without qualifiers
    if not qualifiers:
        gaf_df = gaf_df[gaf_df['qualifier'].isnull()]

    namespace = namespace.capitalize()

    # keep only terms associated with namespace
    if namespace != 'All':
        gaf_df = gaf_df[gaf_df['aspect'] == namespace]

    cols = ['db_object_id', 'db_object_symbol', 'qualifier', 'go_id',
            'evidence_code', 'aspect']

    gaf_df = gaf_df[cols]

    return gaf_df


def build_bg_ontology(bg_genes: Set[str], annots: pd.DataFrame,
                      ont: nx.DiGraph) -> nx.DiGraph:
    '''Get sub-ontology (subgraph) containing only terms (nodes) to which
    background genes (given by their uniprot ids) have been annotated.
    '''

    # Propagate gene annotations through hierarchy
    bg_ont = propagate(bg_genes, annots, ont)

    # Remove terms without annotations to keep only background terms
    no_annot = []
    for node, annots in bg_ont.nodes.data('annots'):
        if annots:
            no_annot.append(node)
    bg_ont = bg_ont.subgraph(no_annot)

    return bg_ont


def propagate(annots: pd.DataFrame,
              ont: nx.DiGraph,
              uniprot_ids: Optional[Set[str]] = None) -> nx.DiGraph:
    '''Adds and propagates annotations to terms in ontology for the given
    genes (identified by their uniprot ids).
    '''

    # Get subset of human annotations that contain the background genes
    if uniprot_ids:
        annots = annots[annots['db_object_id'].isin(uniprot_ids)]

    # Get dict with GO terms as keys and set of annotated genes as values
    values = {k: set(v) for k, v in
              annots.groupby('go_id')['db_object_symbol']}

    # Propagate annotations through hierarchy
    goenrich.enrich.propagate(ont, values, 'annots')

    return ont


def build_qry_ontology(query_genes: Set[str], bg_genes: Set[str],
                       annots: pd.DataFrame,
                       bg_ont:  nx.DiGraph) -> nx.DiGraph:

    # Remove query genes that are not present in background genes
    if diff := query_genes.difference(bg_genes):
        print(f'Warning! The following {len(diff)} genes are query genes '
              f'but are not contained in background. \n{diff}')

    query_genes = bg_genes.intersection(query_genes)

    # Create new graph to add query annotations and delete all previous
    # annotations
    qry_ont = nx.DiGraph(bg_ont)
    for node in qry_ont:
        del qry_ont.node[node]['annots']

    # Propagate gene annotations through hierarchy
    propagate(query_genes, annots, qry_ont)

    # Get fraction of query vs background genes that are annotated to a term
    for node in qry_ont:
        len_bg = len(bg_ont.nodes[node]['annots'])
        len_query = len(qry_ont.nodes[node]['annots'])
        qry_ont.nodes[node]['ratio'] = len_query/len_bg
        qry_ont.nodes[node]['frac'] = f'{len_query}/{len_bg}'
        qry_ont.nodes[node]['term_size'] = len_bg
        qry_ont.nodes[node]['term_id'] = node

        # Get annotations that have been included and excluded in screens
        qry_ont.nodes[node]['inc_annots'] = list(bg_ont.nodes[node]['annots']
                                                 .intersection(qry_ont.nodes
                                                               [node]['annots']))
        qry_ont.nodes[node]['exc_annots'] = list(bg_ont.nodes[node]['annots']
                                                 .difference(qry_ont.nodes
                                                             [node]['annots']))

    return qry_ont


def filter_by_depth(ont: nx.DiGraph, depths: List[int]) -> nx.DiGraph:

    # Filter terms by depth
    nodes = []
    for node, depth in ont.nodes.data('depth'):
        if depth in depths:
            nodes.append(node)
    d_ont = ont.subgraph(nodes)

    return d_ont


def ont_to_df(ont: nx.DiGraph) -> pd.DataFrame:
    '''Write ontology annotation data as dataframe, ordering terms by their
    size and their screen coverage ratio.
    '''

    ont_df = pd.DataFrame()
    ont_df['term_id'] = [node for node in ont]

    cols = ['name', 'namespace', 'depth', 'term_size',
            'frac', 'ratio', 'inc_annots', 'exc_annots']
    for col in cols:
        ont_df[col] = [ont.nodes[node][col] for node in ont]

    cols_lists = ['inc_annots', 'exc_annots']
    for col in cols_lists:
        ont_df[col] = [','.join(ont.nodes[node][col]) for node in ont]

    ont_df = ont_df.sort_values(['ratio', 'term_size'],
                                ascending=[True, False])

    ont_df = ont_df.rename(columns={'frac': 'fraction_covered',
                                    'ratio': 'ratio_covered',
                                    'inc_annots': 'covered_genes',
                                    'exc_annots': 'not_covered_genes'})

    return ont_df


def get_gene_stats(ont_df: pd.DataFrame) -> pd.DataFrame:
    '''Get stats on non-covered genes, on how many terms they are annotated,
    the average rank of those terms, and a score based on term # / term rank.
    '''

    def stats_per_gene(group):
        gene_stats = pd.Series()
        gene_stats['mean_term_rank'] = group.term_rank.mean()
        gene_stats['term_count'] = int(len(group))
        gene_stats['score'] = gene_stats.term_count/gene_stats.mean_term_rank
        gene_stats['term_id_list'] = ','.join(group.term_id.to_list())
        gene_stats['term_name_list'] = ','.join(group['name'].to_list())

        return gene_stats

    genes = ont_df.copy(deep=True)
    genes.not_covered_genes = genes.not_covered_genes.str.split(',')
    genes['term_rank'] = genes.reset_index().index + 1
    genes = genes.explode(column='not_covered_genes')
    genes = genes.query('not_covered_genes != ""')

    gene_stats = genes.groupby('not_covered_genes').apply(stats_per_gene)

    gene_stats = gene_stats.sort_values(['score', 'term_count',
                                         'mean_term_rank'],
                                        ascending=[False, False, True])

    return gene_stats
