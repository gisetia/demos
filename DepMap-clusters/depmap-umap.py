#%%
import umap
import pandas as pd
from sklearn.cluster import DBSCAN
from reactome2py import analysis

#%%
# Load gene effect data
crispr = pd.read_parquet('gs://gisetia-ccle/processed_data/cell_cols/'
                         'DepMap23Q4/CRISPRGeneEffect.pq').fillna(100)

#%%
# Get UMAP using pearson correlation as similarity measure between genes
reducer = umap.UMAP(metric='correlation', n_neighbors=2,
                    min_dist=0, random_state=0,)
embedding = reducer.fit_transform(crispr.values)

reducer_plot = umap.UMAP(metric='correlation',
                         n_neighbors=2, min_dist=.6, random_state=0,)
embedding_plot = reducer.fit_transform(crispr.values)

#%%
# Cluster the result from UMAP using DBSCAN (parameters depend on how the data
# looks like. eps should be similar to the distance between points in the
# clusters. min_samples is the minimum number of genes a cluster should
# contain)
db = DBSCAN(eps=0.1, min_samples=2).fit(embedding)
labels = db.labels_

# Organize data into dataframe
umap_data = pd.DataFrame({'genes': crispr.index,
                          'x_coord': [x[0] for x in embedding],
                          'y_coord': [x[1] for x in embedding],
                          'x_coord_plot': [x[0] for x in embedding_plot],
                          'y_coord_plot': [x[1] for x in embedding_plot], })
# Add id of clusters found by DBSCAN
umap_data['cluster_id'] = db.labels_

umap_data.to_csv('depmap_umap_data.csv')

#%%
# Group genes within the same cluster together to get a dataframe where
# every row corresponds to a cluster
clusters = umap_data.groupby('cluster_id')['genes'].apply(list).reset_index()
clusters['size'] = clusters['genes'].apply(len)
clusters.sort_values('size', ascending=False)

def cluster_enrichment(cluster):
    # Function for Reactome enrichment analysis per cluster
    genes_cluster = cluster['genes']
    print('Processing: ', cluster['cluster_id'])
    react = analysis.identifiers(ids=','.join(genes_cluster),
                                 species='Homo Sapiens', include_disease=False,
                                 page_size=1)['pathways']
    
    try:
        react_df = (pd.DataFrame
                    .from_dict({'id': [x['stId'] for x in react],
                                'name': [x['name'] for x in react],
                                'p': [x['entities']['pValue'] for x in react]
                                }).sort_values('p'))

        cluster_data = pd.Series({'reactome_term': react_df['id'].values[0],
                                'reactome_name': react_df['name'].values[0],
                                'reactome_pval': react_df['p'].values[0]})
    except IndexError:

        cluster_data = pd.Series({'reactome_term': 'None',
                                  'reactome_name': 'None',
                                  'reactome_pval': 1})

    return cluster_data


# Add enrichment data to our cluster dataframe
cols = ['reactome_term', 'reactome_name', 'reactome_pval']
clusters[cols] = clusters.apply(cluster_enrichment, axis=1)

clusters.to_csv('depmap_umap_clusters.csv')

# %%
