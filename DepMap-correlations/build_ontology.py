# %%
import json
import numpy as np
import pandas as pd
from gotools import go

save_dir = 'processed_data'

# % Build ontology and read human annotations
namespace = 'p'
obo_file = 'raw_data/go-basic.obo'
gaf_file = 'raw_data/goa_human.gaf.zip'

# Construct ontology
ont = go.build_ontology(obo_file, namespace=namespace)

# Read annotations and remove some evidence codes to annotations
annots = go.read_gaf(gaf_file, namespace=namespace)
exclude_ev_codes = ['IEA', 'ND']
annots = annots.query('evidence_code not in @exclude_ev_codes')

# Propagate annotations
go.propagate(annots, ont)

# Look for root term of ontology branch and get its size
root = [n for n, d in ont.out_degree() if d == 0]
max_annots = len(ont.nodes[root[0]]['annots'])

# Get information content (IC) for each ontology term
for node in ont:
    len_annots = len(ont.nodes[node]['annots'])
    ont.nodes[node]['size'] = len_annots
    ont.nodes[node]['IC'] = len_annots/max_annots

# Convert ontology to dataframe and save
ont_df = pd.DataFrame()
ont_df['term_id'] = [node for node in ont]

cols = ['name', 'namespace', 'depth', 'annots',
        'size', 'IC']
for col in cols:
    ont_df[col] = [ont.nodes[node][col] for node in ont]

ont_df['annots'] = ont_df['annots'].apply(lambda x: list(x))

# Keep only terms that have annotations
ont_df = ont_df.query('size > 0')
ont_df.to_csv(f'{save_dir}/ontology.csv.gz', compression='gzip')

# Get genes and the terms to which they are annotated
gene_terms = ont_df.explode('annots').groupby('annots')['term_id'].agg(list)
with open(f'{save_dir}/gene_terms.json', 'w') as outfile:
    json.dump(gene_terms.to_dict(), outfile)
# gene_terms.to_csv(f'{save_dir}/gene_terms.csv.gz', compression='gzip')

# %%
