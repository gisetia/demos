# %%
import json
import pandas as pd
import numpy as np
import multiprocessing as mp
from itertools import combinations
from functools import partial

with open('processed_data/gene_terms.json') as json_file:
    gene_terms = json.load(json_file)

genes = list(gene_terms.keys())
gene_pairs = combinations(genes, 2)

ont = pd.read_csv('processed_data/ontology.csv.gz')


def semantic_similarity(pair, ont, gene_terms, count=[0]):
    # as defined in https://www.nature.com/articles/s41467-019-13058-9#Sec9

    if not count[0] % 50000:
        print(count)
    count[0] += 1

    common_terms = (set(gene_terms.get(pair[0]))
                    .intersection(gene_terms.get(pair[1])))

    if common_terms == {}:
        sem_sim = 0
        min_term = 'None'
        min_num = np.inf
        min_name = 'None'
    else:
        common_ont = (ont.query('term_id in @common_terms')
                      .sort_values('size').iloc[0])

        min_term = common_ont['term_id']
        min_num = common_ont['size']
        sem_sim = 2/min_num
        min_name = common_ont['name']

    return *pair, sem_sim, min_term, min_num, min_name

#%%
partial_func = partial(semantic_similarity, ont=ont, gene_terms=gene_terms)

cores = 30
with mp.Pool(cores) as p:
    result = p.imap_unordered(partial_func, gene_pairs, chunksize=10000)
    result = list(result)

data = pd.DataFrame(result, columns=['gene1', 'gene2', 'sem_sim', 'term',
                                     'term_size', 'term_name'])
# data.to_parquet('processed_data/sem-sim.pq')



# %%

# data = pd.read_parquet('processed_data/sem-sim.pq')
# data

# # %%
# sample = data.sample(100, random_state=0)
# sample.to_parquet(
#         'processed_data/sample_sem-sim.pq')
# sample_grouped = sample.groupby(['sem_sim', 'term', 'term_size', 'term_name']).apply(
#     lambda x: str(list(zip(x['gene1'], x['gene2'])))).to_frame().rename(
#         columns={0: 'gene_pairs'})

# sample_grouped.to_parquet('processed_data/sample-grouped_sem-sim.pq')

# sample_grouped

# %%

from fastparquet import write
write('processed_data/chunks_sem-sim.pq', data,
      row_group_offsets=range(0, len(data), 5000000),
      compression='SNAPPY', file_scheme='hive')

# %%
print(data)
