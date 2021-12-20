import pandas as pd
import umap
from sklearn.cluster import DBSCAN
from reactome2py import analysis

# Data matrix with genes and screens
gene_data = [['CD274', 5, -4, 3, 8, -1], ['IFNGR1', 6, -4, 4, 7, -2],
             ['STAT1', 5, -3, 3, 8, -2], ['JAK1', 6, -3, 4, 7, -1],
             ['ATP5F1A', -3, 2, -5, 5, 3], ['ATP5F1E', -4, 3, -4, 4, 3],
             ['ATP5PF', -3, 3, -5, 4, 4]]
gene_data = pd.DataFrame(gene_data, columns=['gene', 'screenA', 'screenB',
                                             'screenC', 'screenD', 'screenE']
                         ).set_index('gene')

# Get UMAP using pearson correlation as similarity measure between genes
reducer = umap.UMAP(metric='correlation')
embedding = reducer.fit_transform(gene_data.values)

# Cluster the result from UMAP using DBSCAN (parameters depend on how the data
# looks like. eps should be similar to the distance between points in the
# clusters. min_samples is the minimum number of genes a cluster should
# contain)
db = DBSCAN(eps=1, min_samples=2).fit(embedding)
labels = db.labels_

# Organize data into dataframe
umap_data = pd.DataFrame({'genes': gene_data.index,
                          'x_coord': [x[0] for x in embedding],
                          'y_coord': [x[1] for x in embedding]})
# Add id of clusters found by DBSCAN
umap_data['cluster_id'] = db.labels_

# Group genes within the same cluster together to get a dataframe where
# every row corresponds to a cluster
clusters = umap_data.groupby('cluster_id')['genes'].apply(list).reset_index()


# Function for Reactome enrichment analysis per cluster
def cluster_enrichment(cluster):

    genes_cluster = cluster['genes']
    react = analysis.identifiers(ids=','.join(genes_cluster),
                                 species='Homo Sapiens', include_disease=False,
                                 page_size=1)['pathways']
    react_df = (pd.DataFrame
                .from_dict({'id': [x['stId'] for x in react],
                            'name': [x['name'] for x in react],
                            'p': [x['entities']['pValue'] for x in react]
                            }).sort_values('p'))

    cluster_data = pd.Series({'reactome_term': react_df['id'].values[0],
                              'reactome_name': react_df['name'].values[0],
                              'reactome_pval': react_df['p'].values[0]})

    return cluster_data


# Add enrichment data to our cluster dataframe
cols = ['reactome_term', 'reactome_name', 'reactome_pval']
clusters[cols] = clusters.apply(cluster_enrichment, axis=1)
